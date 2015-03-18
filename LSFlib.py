'''
Created on Jun 29, 2011

@author: lgoff
'''
import os, re
import subprocess
import time
import sys

from misc import pp

#Constants
lsf_mem = 32
lsf_default_queue = "normal_parallel" # normal_parallel  since it has less users 

#######################
#Error Handling
#######################
class LSFError(Exception):
	"""Base class for exceptions in this module."""
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

#################
#Base Class
#################
class LSFJob(object):
	'''
	LSF Job
	'''


	def __init__(self,cmd_str,job_name=None,job_group=None,blocking=False,outfilename=None,errfilename=None,queue_name=None,job_mem=None,job_cores=1,notify=None):
		'''
		Creates instance of LSFJob
		#Don't use blocking because this is a limiting resource on Odyssey LSF
		'''
		self.cmd_str = cmd_str
		
		global lsf_default_queue
		if queue_name == None:
			self.queue = lsf_default_queue
		else:
			self.queue = queue_name

		if outfilename == None:
			self.outfile = tmp_name("bsub_out_")
		else:
			self.outfile = outfilename
		if errfilename == None:
			self.errfile = tmp_name("bsub_err_")
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
		
		bsub_str = ["bsub"]
		
		if notify:
			bsub_str.extend(["-N"])
		
		bsub_str.extend(["-q", self.queue])
		
		if self.job_name != None:
			bsub_str.extend(["-J", self.job_name])
		
		if self.group != None:
			bsub_str.extend(['-g', self.group])
		
		if blocking != False:
			bsub_str.extend(["-K"])
		
		global lsf_mem
		if job_mem != None and lsf_mem != None:
			self.job_mem = min(self.job_mem, lsf_mem)
			bsub_str.extend(["-R rusage[mem=%d]" % self.job_mem])
		
		bsub_str.extend(["-R span[hosts=1]"])
		
		bsub_str.extend(["-oo", self.outfile])
		bsub_str.extend(["-eo", self.errfile])
		bsub_str.extend(["%s" % self.cmd_str])
		
		self.bsub_str = bsub_str
		
		#Handle if queue == "local"
		if self.queue == "local":
			local_str = [""]
			local_str.extend([">", self.outfile])
			local_str.extend(["2>", self.errfile])
			
			#TODO: Add self.cmd_str to bsub_str so command actually gets run.
			self.bsub_str = local_str

	def __repr__(self):
		return "Instance of class LSF Job:\n\t%s\n\tSubmitted: %s\n\t Complete: %s\n" % (self.cmd_str,self.submit_flag,self.complete) + str(pp(self.__dict__))
	
	def __str__(self):
		return " ".join(self.bsub_str)

	def submit(self): # wait pend
		if self.submit_flag == True:
			print >>sys.stderr, "Job already submitted"
			return 0# what do you return here?
		
		self.submit_proc = subprocess.Popen(self.bsub_str,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		
		#Handle local jobs
		if self.queue == "local":
			self.submit_flag = True
			self.status = 'RUN'
			self.submit
			self.jobID = self.submit_proc.pid
			print >>sys.stderr, "Job running locally with pid %d" % self.jobID
			return 0
		
		#Handle queued jobs
		if self.submit_proc.wait() != 0:
			raise LSFError("Could not submit to LSF. Error %d" % self.submit_proc.poll())
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
		"""This will poll using bjobs for the specific jobID for a given instance of LSFJob"""
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
					raise LSFError("Problem with local job %d. Error %d" % (self.jobID,int(self.submit_proc.poll())))
				return
			tmp = subprocess.Popen('bjobs -a -w %d' % self.jobID,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			tmp_err = tmp.stderr.read().rstrip()
			notfoundpat = re.compile("Job \<[0-9]+\> is not found")
			failedpat = "Failed in an LSF library call"
				
			#wait until the bjobs query returns  (not until the job itself is finished)
			while tmp.wait() > 0:
				if tmp_err.count(failedpat) > 0:
					print >>sys.stderr, tmp_err
					time.sleep(20)
					tmp = subprocess.Popen('bjobs -w %d' % self.jobID,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
					tmp_err = tmp.stderr.read().rstrip()
					continue
				else:
					raise LSFError(tmp_err) # never caught
			if tmp.wait() == 0:
				#If job is not found and was previously pending/running assume it is finished
				if notfoundpat.match(tmp_err):
					if self.status in ['RUN','PEND']:
						self.status = 'DONE'
						self.complete = True
						return self.status
					else: # was never run
						print >>sys.stderr, "waited, job did not run " + tmp_err
						return tmp_err
				#else: job still exists, update its status	
				tmp_lines = [x.rstrip() for x in tmp.stdout.readlines()]
				keys,values = [x.split() for x in tmp_lines]
				tmpDict = dict(zip(keys,values))
				#pp(tmpDict)
				self.status = tmpDict['STAT']
				self.submit_time = tmpDict['SUBMIT_TIME']
				self.exec_host = tmpDict['EXEC_HOST']
				self.submit_host = tmpDict['FROM_HOST']
				return self.status
			else:
				#Should not reach this line... CONSIDER erasing and doing while tmp.wait!=0 
				raise LSFError("Problem with bjobs polling. Error %s" % tmp_err)
	
	def getJobId(self):
		if self.submit_flag:
			jobID_search = re.search("\<[0-9]+\>",self.submit_status)
			self.jobID = int(jobID_search.group().strip("><"))
			return
		else:
			print "Job not yet submitted."
			return
	
	def kill(self):
		#Added this to fix cases were kill fails because there is no job id
		if self.status in ['NOT SUBMITTED'] or self.jobID== -999 :
			self.status = 'NOT SUBMITTED'
			return
		tmp = subprocess.Popen('bkill %d' % self.jobID,shell=True) #This is not working : we kill but the tmp wait doesn't return 0
		while tmp.wait() > 0:
			time.sleep(3)
			if tmp.wait()< 0: #if were not able to kill , try again
				tmp = subprocess.Popen('bkill %d' % self.jobID,shell=True) #This is not working : we kill but the tmp wait doesn't return 0
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
		while self.status in['SUBMITTED','PEND','RUN','SUSP']:
			time.sleep(30)
			self.poll()
			if self.status in ['SUSP']:
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
