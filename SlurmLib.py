'''
Created on Jun 25, 2015

Slurm Library for submitting Jobs from Python

@author: lgoff
'''
import os, re
import subprocess
import time
import sys

from misc import pp

#Constants
slurm_mem = 32
slurm_default_queue = "shared" # normal_parallel  since it has less users 

#######################
#Error Handling
#######################
class SlurmError(Exception):
	"""Base class for exceptions in this module."""
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

#################
#Base Class
#################
class SlurmJob(object):
	'''
	Slurm Job
	'''


	def __init__(self,cmd_str,job_name=None,job_group=None,blocking=False,outfilename=None,errfilename=None,queue_name=None,job_mem=None,job_cores=1,notify=None):
		'''
		Creates instance of SlurmJob
		
		'''
		self.cmd_str = cmd_str
		
		global slurm_default_queue
		if queue_name == None:
			self.queue = slurm_default_queue
		else:
			self.queue = queue_name

		if outfilename == None:
			self.outfile = tmp_name("sbatch_out_")
		else:
			self.outfile = outfilename
		if errfilename == None:
			self.errfile = tmp_name("sbatch_err_")
		else:
			self.errfile = errfilename
		
		self.job_name = job_name
		self.group = job_group
		self.job_mem = job_mem
		self.submit_flag = False
		self.complete = False
		self.status = 'NOT SUBMITTED'
		self.jobID= -999
				
		self.submit_time = ""
		self.exec_host = ""
		self.submit_host = ""
		
		if blocking != False:
			sbatch_str = ["srun"]
		else:
			sbatch_str = ["sbatch"]
		
		if notify:
			sbatch_str.extend(["--mail-type=END"])
		
		sbatch_str.extend(["-p", self.queue])
		
		if self.job_name != None:
			sbatch_str.extend(["--job-name=%s" % self.job_name])
		
		#TODO: Do I really need this?
		#if self.group != None:
		#	sbatch_str.extend(['-g', self.group])
		
		global slurm_mem
		if job_mem != None and slurm_mem != None:
			self.job_mem = min(self.job_mem, slurm_mem)
			sbatch_str.extend(["--mem=%d" % self.job_mem])
		
		#sbatch_str.extend(["-R span[hosts=1]"])
		
		sbatch_str.extend(["-o", self.outfile])
		sbatch_str.extend(["-e", self.errfile])
		sbatch_str.extend(['--wrap="%s"' % self.cmd_str])
		
		self.sbatch_str = sbatch_str
		
		#Handle if queue == "local"
		if self.queue == "local":
			local_str = [""]
			local_str.extend([">", self.outfile])
			local_str.extend(["2>", self.errfile])
			
			#TODO: Add self.cmd_str to sbatch_str so command actually gets run.
			self.sbatch_str = local_str

	def __repr__(self):
		return "Instance of class SlurmJob:\n\t%s\n\tSubmitted: %s\n\t Complete: %s\n" % (self.cmd_str,self.submit_flag,self.complete) + str(pp(self.__dict__))
	
	def __str__(self):
		return " ".join(self.sbatch_str)

	def submit(self): # wait pend
		if self.submit_flag == True:
			print >>sys.stderr, "Job already submitted"
			return 0# what do you return here?
		
		self.submit_proc = subprocess.Popen(self.sbatch_str,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		
		#Handle local jobs
		if self.queue == "local":
			self.submit_flag = True
			self.status = 'RUNNING'
			self.submit
			self.jobID = self.submit_proc.pid
			print >>sys.stderr, "Job running locally with pid %d" % self.jobID
			return 0
		
		#Handle queued jobs
		if self.submit_proc.wait() != 0:
			raise SlurmError("Could not submit to Slurm. Error %d" % self.submit_proc.poll())
		else:
			self.submit_flag = True
			self.status = 'SUBMITTED'
			self.submit_status = self.submit_proc.stdout.read().rstrip() 
			self.getJobId()
			#Wait until job switched from submitted to pend/run
			while self.status in ['SUBMITTED'] :
				try:
					self.poll()
				except Exception , e:
					print >> sys.stderr,'Exception poll error: %s\n' %e
					
		print >>sys.stderr, self.submit_status
		return self.submit_proc.wait()
	
	def poll(self):
		"""This will poll using scontrol view job for the specific jobID for a given instance of SlurmJob"""
		if not self.submit_flag:
			return "Job not yet submitted"
		elif self.complete:
			return "Job completed"
		else:
			#Handle local jobs
			if self.queue == "local":
				if self.submit_proc.poll() == 0:
					self.complete = True
					self.status = 'DONE'
					return self.status
				if self.submit_proc.poll() == None:
					return self.status
				if self.submit_proc.poll() > 0 or self.submit_proc.poll() < 0:
					raise SlurmError("Problem with local job %d. Error %d" % (self.jobID,int(self.submit_proc.poll())))
				return
			tmp = subprocess.Popen('scontrol show job %d' % self.jobID,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			tmp_err = tmp.stderr.read().rstrip()
			notfoundpat = re.compile("Invalid job id specified")
			failedpat = "Failed in a Slurm library call" #TODO:  THis is wrong
				
			#wait until the scontrol show job query returns  (not until the job itself is finished)
			while tmp.wait() > 0:
				if tmp_err.count(failedpat) > 0:
					print >>sys.stderr, tmp_err
					time.sleep(20)
					tmp = subprocess.Popen('scontrol show job %d' % self.jobID,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
					tmp_err = tmp.stderr.read().rstrip()
					continue
				else:
					raise SlurmError(tmp_err) # never caught
			if tmp.wait() == 0:
				#If job is not found and was previously pending/running assume it is finished
				if notfoundpat.match(tmp_err):
					if self.status in ['RUNNING','PENDING']:
						self.status = 'DONE'
						self.complete = True
						return self.status
					else: # was never run
						print >>sys.stderr, "waited, job did not run " + tmp_err
						return tmp_err
				#else: job still exists, update its status	
				tmp_lines = [x.rstrip() for x in tmp.stdout.readlines()]
				kv_pairs = [x.split(" ") for x in tmp_lines]
				keys,values = [x.split("=") for x in kv_pairs]
				tmpDict = dict(zip(keys,values))
				#pp(tmpDict)
				self.status = tmpDict['JobState']
				self.submit_time = tmpDict['SubmitTime']
				self.exec_host = tmpDict['ExcNodeList']
				self.submit_host = tmpDict['AllocNode:Sid']
				return self.status
			else:
				#Should not reach this line... CONSIDER erasing and doing while tmp.wait!=0 
				raise SlurmError("Problem with scontrol show job polling. Error %s" % tmp_err)
	
	def getJobId(self):
		if self.submit_flag:
			jobID_search = re.search("[0-9]+",self.submit_status)
			#self.jobID = int(jobID_search.group().strip("><"))
			self.jobID = int(jobID_search.group().rstrip())
			return
		else:
			print "Job not yet submitted."
			return
	
	def kill(self):
		#Added this to fix cases were kill fails because there is no job id
		if self.status in ['NOT SUBMITTED'] or self.jobID== -999 :
			self.status = 'NOT SUBMITTED'
			return
		tmp = subprocess.Popen('scancel %d' % self.jobID,shell=True) #This is not working : we kill but the tmp wait doesn't return 0
		while tmp.wait() > 0:
			time.sleep(3)
			if tmp.wait()< 0: #if were not able to kill , try again
				tmp = subprocess.Popen('scancel %d' % self.jobID,shell=True) #This is not working : we kill but the tmp wait doesn't return 0
		if tmp.wait() == 0:
			self.status = 'KILLED'
			self.submit_flag = False
			self.complete = False
			self.status = 'NOT SUBMITTED'
		return
	
	def wait(self):
		self.poll()
		if not self.submit_flag:
			print "Job not yet submitted"
			return
		while self.status in['SUBMITTED','PENDING','RUNNING','SUSPENDED']:
			time.sleep(30)
			self.poll()
			if self.status in ['SUSPENDED']:
				print >> sys.stderr,'SUSPENDED : %d \n' % self.jobID 
		self.status = 'DONE'
		self.complete = True
		return
			

##############
#Helper functions
##############
def tmp_name(prefix):
	tmp_root = "tmp/"
	if os.path.exists(tmp_root):
		pass
	else:
		os.mkdir(tmp_root)
	return tmp_root + prefix + os.tmpnam().split('/')[-1]
