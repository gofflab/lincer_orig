#!/usr/bin/env python
# encoding: utf-8
"""
build-transcript-db.py

Created by Cole Trapnell on 2010-4-11.
Copyright (c) 2010 Cole Trapnell. All rights reserved.

last Revision: Feb 25th 2011, Moran Cabili
"""

import sys
try:
    import psyco
    psyco.full()
except ImportError:
    pass

#import sys
import getopt
import subprocess
import errno
import os
import tempfile
import warnings
import shutil
import copy
from datetime import datetime, date, time
import ConfigParser
import math
import time
from types import *

# Needed to handle LSF submissions and polling
import LSFlib # Must be replaced/avoided in MARCC slurm environment

# Needed for linc naming conventions
import GTF2Genbank
import newGTF
import intervallib
import dbConn
import bisect
import sys,getopt
from misc import rstrips

use_message = '''
 build-transcript-db is a script to build a single "meta-assembly" from many Cufflinks assemblies

 Usage:
     build-transcript-db [options] <primary_asm_list.txt>

 Options:
     -o/--output-dir                <string>    [ default: ./ ]
     -p/--threads-per-node          <int>       [ default: 1  ]
     -n/--lsf-nodes                 <int>       [ default: 10  ]
     --system                       <string> ("Broad" or "Valor" or "MARCC") [ default: MARCC]
     --sample-reads                 <filename>
     --ref-sequence                 <dir>
     -g/--genome                    <string>    [ default: hg19 ]
     --lsf-queue                    <string>
     --cufflinks-queue              <string>
     --lsf-mem                      <int>
     --min-transcript-cov           <int>
     --min-linc-length              <int>
     --min-linc-exons               <int>
     --no-pfam
     --ref-gtf                      <string>
     --mask-gtf                     <string>
     --pseudogene                   <string>
     --external-linc-gtf            <string>
     --csf-score                    <int>
     --out-prefix                   <String>
     --run-non-assembly-mode
     --run-merge-catalog-mode
     --run-get-intergenic-mode
     --run-meta-assembly
     --keep-tmp
     --version
     --no-notify                                [ Turns off email Notification for LSF jobs ]
     --help

'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

output_dir = "./"
logging_dir = output_dir + "logs/"
run_log = None
run_cmd = None
run_meta_assembly = 0
non_assembly_mode = 0
merge_catalog_mode = 0
get_intergenic_mode = 0
csf_only_mode = 0
csf_score_threshold= 100
tmp_dir = output_dir + "tmp/"
bin_dir = sys.path[0] + "/"

lsf_mem = 64
#lsf_queue = "compbio-week"
#cufflinks_queue="compbio-week"

#Valor-specific global settings
lsf_queue = "normal_parallel"   # Changed from short_parallel , no need for parallel
cufflinks_queue="normal_parallel"
pfam_db = "/n/rinn_data1/indexes/Pfam/"
aa_converter =  "perl /n/rinn_data1/user-supported/bin/seqToAminoAcid.pl"
gffread = "gffread"
pfammer = "perl /n/rinn_data1/indexes/Pfam/Tools/PfamScan/pfam_scan.pl"
csf_bin = "/n/rinn_data1/users/lgoff/sw/PhyloCSF/batchPhyloCSF_Loyal.ml.exe"

#ok_str = "\t\t\t\t[OK]\n"
fail_str = "\t[FAILED]\n"

def whereis(program):
   for path in os.environ.get('PATH', '').split(':'):
       if os.path.exists(os.path.join(path, program)) and \
          not os.path.isdir(os.path.join(path, program)):
           return os.path.join(path, program)
   return None

awk_bin = whereis("awk")


class TestParams:

    class SystemParams:
        def __init__(self,
                     threads,
                     keep_tmp):
            self.threads = threads
            self.keep_tmp = keep_tmp
            self.lsf_nodes = 10
            self.system = "valor"
            self.notify = True


        def parse_options(self, opts):
            for option, value in opts:
                if option in ("-p", "--num-threads"):
                    self.threads = int(value)
                if option in ("--keep-tmp"):
                    self.keep_tmp = True
                if option in ("-n", "--lsf-nodes"):
                    self.lsf_nodes = int(value)
                if option in ("--system"):
                    self.system = value
                if option in ("--no-notify"):
                    self.notify = False
        def check(self):
            pass

    def __init__(self):
        self.system_params = self.SystemParams(1,               # threads
                                               False)           # keep_tmp
        self.ref_gtf = None
        self.external_linc_gtf = None
        self.external_true_lincs_gtf = None
        self.sample_reads_list = None
        self.mask_gtf = None
        self.fasta = None
        self.genome = 'hg19'
        self.min_linc_length = 250
        self.min_linc_exons = 2
        self.min_transcript_cov = 2
        self.no_csf = False
        self.no_pfam = False
        self.out_prefix=""
        self.pseudogene = None

    def check(self):
        self.system_params.check()

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:],
                                       "hvp:o:G:M:s:q:g:",
                                       ["version",
                                        "help",
                                        "ref-sequence=",
                                        "min-linc-length=",
                                        "min-linc-exons=",
                                        "sample-reads=",
                                        "min-transcript-cov=",
                                        "ref-gtf=",
                                        "mask-gtf=",
                                        "pseudogene=",
                                        "external-linc-gtf=",
                                        "external-true-lincs-gtf=",
                                        "output-dir=",
                                        "num-threads=",
                                        "lsf-nodes=",
                                        "no-notify",
                                        "system=",
                                        "lsf-queue=",
                                        "cufflinks-queue=",
                                        "lsf-mem=",
                                        "csf-score-threshold=",
                                        "out-prefix=",
                                        "run-meta-assembly",
                                        "run-non-assembly-mode",
                                        "run-merge-catalog-mode",
                                        "run-get-intergenic-mode",
                                        "keep-tmp",
                                        "genome",
                                        "no-pfam",
                                        "no-csf",
                                        "csf-only"])
        except getopt.error, msg:
            raise Usage(msg)

        self.system_params.parse_options(opts)

        global lsf_queue
        global output_dir
        global logging_dir
        global tmp_dir
        global lsf_mem
        #Moran added q for cufflinks runs
        global cufflinks_queue
        global run_meta_assembly
        global non_assembly_mode
        global merge_catalog_mode
        global get_intergenic_mode
        global csf_only_mode

        # option processing
        for option, value in opts:
            if option in ("-v", "--version"):
                print "linc builder v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(use_message)
            if option in ("-G", "--ref-gtf"):
                self.ref_gtf = value
            if option in ("-M", "--mask-gtf"):
                self.mask_gtf = value
            if option == "--external-linc-gtf":
                self.external_linc_gtf = value
            if option == "--external-true-lincs-gtf":
                self.external_true_lincs_gtf = value
            if option in ("-s", "--ref-sequence"):
                self.fasta = value
            if option in ("-g", "--genome"):
                self.genome = value
            if option in "--sample-reads":
                self.sample_reads_list = value
            if option in "--min-transcript-cov":
                self.min_transcript_cov = float(value)
            if option == "--min-linc-length":
                self.min_linc_length = int(value)
            if option == "--min-linc-exons":
                self.min_linc_exons = int(value)
            if option == "--out-prefix":
                self.out_prefix = value
            if option in ("-q", "--lsf-queue"):
                lsf_queue = value
            if option == "--cufflinks-queue":
                cufflinks_queue = value
            if option == "--lsf-mem":
                lsf_mem = int(value)
            if option == "--csf-score-threshold":
                csf_score_threshold = int(value)
            if option == "--run-meta-assembly":
                run_meta_assembly = 1
            if option == "--run-non-assembly-mode":
                non_assembly_mode = 1
            if option == "--run-merge-catalog-mode":
                merge_catalog_mode = 1
            if option == "--run-get-intergenic-mode":
                get_intergenic_mode = 1
            if option == "--pseudogene" :
                self.pseudogene = value
            if option in ("-o", "--output-dir"):
                output_dir = value + "/"
                logging_dir = output_dir + "logs/"
                tmp_dir = output_dir + "tmp/"
            if option == '--no-pfam':
                self.no_pfam = True
            if option == '--no-csf':
                self.no_csf = True
            if option == '--csf-only':
                csf_only_mode = 1

        return args


def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def prepare_output_dir():

    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)

    if os.path.exists(logging_dir):
        pass
    else:
        os.mkdir(logging_dir)

    if os.path.exists(tmp_dir):
        pass
    else:
        os.mkdir(tmp_dir)

def formatTD(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)

def tmp_name(prefix):
    tmp_root = output_dir + "tmp/"
    if os.path.exists(tmp_root):
        pass
    else:
        os.mkdir(tmp_root)
    return tmp_root + prefix + os.tmpnam().split('/')[-1]

def get_version():
   return "0.0.1"

# Job group can be specified with a string like "/cufflinks_grp"
# Calling with blocking=True will embed a -K in the bsub command,
# which will cause the shell to block until the submitted command finishes
# Running such a command string with sys.executable and then waiting
# until all shell commands are done will block the pipeline until the whole
# job group returns
def bsub_cmd(cmd_str, job_group=None, blocking=False, outfilename=None, errfilename=None, queue_name=None, job_mem=None, job_cores=1, notify=None):
    #return cmd_str
    global lsf_queue

    if lsf_queue == "local" or queue_name == "local":
        return "%s 1>%s 2>%s" % (cmd_str, outfilename, errfilename)
        
    if queue_name == None:
        queue_name = lsf_queue

    if outfilename == None:
        out_file = tmp_name("bsub_out_")
    else:
        out_file = outfilename
    if errfilename == None:
        err_file = tmp_name("err_")
    else:
        err_file = errfilename
    
    if notify == None:
        notify = params.system_params.notify
    
    bsub_str = ["bsub"]
    if notify:
        bsub_str.extend(["-N"])
    bsub_str.extend(["-q", queue_name])
    if job_group != None:
        bsub_str.extend(["-g", job_group])
    if blocking != False:
        bsub_str.extend(["-K"])

    if job_mem != None and lsf_mem != None:
        global lsf_mem
        job_mem = lsf_mem
        bsub_str.extend(["-R rusage[mem=%d]" % job_mem])

    bsub_str.extend(["-R span[hosts=1]"])
    
    #bsub_str.extend(["-oo", out_file, "-f", '\"%s <\"' % out_file])
    bsub_str.extend(["-oo", out_file])
    bsub_str.extend(["-eo", err_file])
    bsub_str.extend(["'%s'" % cmd_str])
    return " ".join(bsub_str)

def sbatch_cmd(cmd_str, job_group=None, blocking=False, outfilename=None, errfilename=None, queue_name=None, job_mem=None, job_cores=1, notify=None):
    #return cmd_str
    global lsf_queue

    if lsf_queue == "local" or queue_name == "local":
        return "%s 1>%s 2>%s" % (cmd_str, outfilename, errfilename)
        
    if queue_name == None:
        queue_name = lsf_queue

    if outfilename == None:
        out_file = tmp_name("bsub_out_")
    else:
        out_file = outfilename
    if errfilename == None:
        err_file = tmp_name("err_")
    else:
        err_file = errfilename
    
    if notify == None:
        notify = params.system_params.notify
    
    if blocking == True:
        bsub_str = ["srun"]
    else:
        bsub_str = ["sbatch"]

    if notify:
        bsub_str.extend(["--mail=type=END"])
    
    bsub_str.extend(["--partition=%s" % queue_name])
    
    if job_group != None:
        #bsub_str.extend(["-g", job_group])
        pass
    
    if job_mem != None and lsf_mem != None:
        global lsf_mem
        job_mem = lsf_mem
        bsub_str.extend(["--mem=%d" % job_mem])
    
    #bsub_str.extend(["-oo", out_file, "-f", '\"%s <\"' % out_file])
    bsub_str.extend(["-o", out_file])
    bsub_str.extend(["-e", err_file])
    bsub_str.extend(["'%s'" % cmd_str])
    return " ".join(bsub_str)

def cufflinks(params,
              out_dir,
              sam_file,
              gtf_file=None,
              extra_opts=["-q", "--num-bootstrap-samples", "0", "--max-bundle-frags", "200000", "--overhang-tolerance", "200", "--library-type=transfrags",  "-A","0.0", "-F", "0.0", "-j 0.0"],
              lsf=False,
              curr_queue=None):
    if gtf_file != None:
        print >> sys.stderr, "[%s] Quantitating transcripts" % (right_now())
    else:
        print >> sys.stderr, "[%s] Assembling transcripts" % (right_now())

    global cufflinks_queue
    if curr_queue == None:
        curr_queue = cufflinks_queue

    cmd = ["cufflinks"]

    if out_dir != None and out_dir != "":
        cmd.extend(["-o", out_dir])

    if gtf_file != None:
        cmd.extend(["-G", gtf_file])

    if extra_opts != None:
        cmd.extend(extra_opts)

    # Run Cufflinks with more than one thread?
    cmd.extend(["-p", str(params.system_params.threads)])

    cmd.append(sam_file)

    try:

        if lsf:
            global lsf_mem
            cmd = bsub_cmd(" ".join(cmd), "/cov_cufflinks", blocking=True, job_cores=params.system_params.threads, queue_name=curr_queue, job_mem=lsf_mem)
            print >> run_log, cmd
            p = subprocess.Popen(cmd, shell=True)
            return p
        else:
            print >> run_log, " ".join(cmd)
            ret = subprocess.call(cmd)
            if ret != 0:
                print >> sys.stderr, fail_str, "Error: could not execute cufflinks"
                exit(1)
    # cufflinks not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cufflinks not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def cuffdiff(params, out_dir, gtf_file, sam_files):

    print >> sys.stderr, "[%s] Testing for differences with cuffdiff" % (right_now())

    cmd = ["cuffdiff"]

    if out_dir != None and out_dir != "":
        cmd.extend(["-o", out_dir])

    # Run Cuffdiff with more than one thread?
    cmd.extend(["-p", str(params.system_params.threads)])

    cmd.append(gtf_file)
    cmd.extend(sam_files)

    try:
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffdiff"
            exit(1)
    # cuffdiff not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffdiff not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def cuffcompare(params, prefix, ref_gtf, fasta, cuff_gtf):

    print >> sys.stderr, "[%s] Comparing reference %s to assembly %s" % (right_now(), ref_gtf, cuff_gtf)
    cmd = ["cuffcompare"]

    if  prefix != None:
        cmd.extend(["-o", prefix])
    if  ref_gtf != None:
        cmd.extend(["-r", ref_gtf])
    if  fasta != None:
        cmd.extend(["-s", fasta])
    if type(cuff_gtf) == ListType:
        for g in cuff_gtf:
            cmd.extend([g])
    else:
        cmd.extend([cuff_gtf])

    try:
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def gtf_to_sam(gtf_filename):

    sam_out = tmp_name("gtf2sam_")

    cmd = ["gtf_to_sam"]
    cmd.append(gtf_filename)
    cmd.append(sam_out)
    try:
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute gtf_to_sam"
            exit(1)
    # gtf_to_sam not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: gtf_to_sam not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
    return sam_out

def test_input_files(filename_list):
    """This function takes a file that contains a list of GTF files,
       tests accessibility of each, and returns the list of filenames"""

    OK = True
    input_files = []
    for line in filename_list:
        line = line.strip()

        # Skip comment line
        if len(line) == 0 or line[0] == "#":
            continue
        try:
            g = open(line,"r")
            input_files.append(line)

        except OSError, o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                print >> sys.stderr, fail_str, "Error: could not open %s" % line
            OK = False
    if OK == False:
        sys.exit(1)
    return input_files

def convert_gtf_to_sam(gtf_filename_list):
    """This function takes a list of GTF files, converts them all to
       temporary SAM files, and returns the list of temporary file names."""
    print >> sys.stderr, "[%s] Converting GTF files to SAM" % (right_now())
    OK = True
    sam_input_filenames = []
    for line in gtf_filename_list:
        try:
            sam_out = gtf_to_sam(line)
            sam_input_filenames.append(sam_out)
        except OSError, o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                print >> sys.stderr, fail_str, "Error: could not open %s" % line
            OK = False
    if OK == False:
        sys.exit(1)
    return sam_input_filenames

def merge_sam_inputs(sam_input_list):
    sorted_map_name = tmp_name( "mergeSam_")

    sorted_map = open(sorted_map_name, "w")

    sort_cmd =["sort",
               "-k",
               "3,3",
               "-k",
               "4,4n",
               "--temporary-directory="+tmp_dir]
    sort_cmd.extend(sam_input_list)

    print >> run_log, " ".join(sort_cmd), ">", sorted_map_name
    subprocess.call(sort_cmd,
                   stdout=open(sorted_map_name, "w"))
    return sorted_map_name

def compare_to_reference(meta_asm_gtf, ref_gtf, fasta):
    print >> sys.stderr, "[%s] Comparing against reference file %s" % (right_now(), ref_gtf)
    if fasta != None:
        comp_cmd = '''cuffcompare -o meta_asm_vs_GENCODE -r %s -s %s %s;''' % (ref_gtf, fasta, meta_asm_gtf)
    else:
        comp_cmd = '''cuffcompare -o meta_asm_vs_GENCODE -r %s %s;''' % (ref_gtf, meta_asm_gtf)

    #cmd = bsub_cmd(comp_cmd, "/gencode_cmp", True, job_mem=8)
    cmd = comp_cmd

    try:
        print >> run_log, cmd
        ret = subprocess.call(cmd,shell=True)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
        tmap_out = "meta_asm_vs_GENCODE." + meta_asm_gtf.split("/")[-1] + ".tmap"
        shutil.copy2(tmap_out, "./meta_asm_vs_GENCODE.tmap")
        return ("./meta_asm_vs_GENCODE.combined.gtf", tmap_out)
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def subtract_mask_transcripts(meta_asm_gtf, mask_gtf):
    #This is super dangerous and needs to be reeled in...
    print >> sys.stderr, "[%s] Comparing against mask file %s" % (right_now(), mask_gtf)
    comp_cmd = '''cuffcompare -o meta_asm_vs_mask -r %s %s;''' % (mask_gtf, meta_asm_gtf)
    
    #This should be a much better output format to avoid overwriting existing files.
    tmap_out = "./" + "/".join(meta_asm_gtf.split("/")[:-1]) + "/meta_asm_vs_mask." + meta_asm_gtf.split("/")[-1] + ".tmap" 
    print >> sys.stderr, meta_asm_gtf, tmap_out
    selected_ids = "%s.selected.ids" % tmap_out
    selected_gtf = "%s.selected.gtf" % tmap_out

    #cmd = bsub_cmd(comp_cmd, "/mask_gtf", True, job_mem=8)
    cmd = comp_cmd

    try:
        print >> run_log, cmd
        ret = subprocess.call(cmd,shell=True)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not mask meta assembly"
            exit(1)

    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

    #This will select out 'u' and 'x' class codes from output of cuffcompare.  This is where filtering actually happens.
    try:
        comp_cmd = "awk 'NR>1 {if (($3 ~\"u|x\")) print $5}' %s > %s" % (tmap_out, selected_ids)
        print >> run_log, comp_cmd
        ret = subprocess.call(comp_cmd, shell=True)
    # awk or not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: awk not found on this system.  Did you forget to include it in your PATH?"
        print o
        exit(1)

    #Select records from meta_asm_gtf that have an identifier in selected_ids and write to selected_gtf. 
    try:
        comp_cmd = whereis("select_gtf.py") + " %s %s > %s;" % (selected_ids, meta_asm_gtf, selected_gtf)
        print >> run_log, comp_cmd
        ret = subprocess.call(comp_cmd, shell=True)
        return selected_gtf
    # awk or not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: select_gtf.py not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

#get a specific marked class from a cuffcompare output : use u or X for lincs or lancs respectivaly
def get_class_transcripts(meta_asm_gtf,potential_ids,tmap_out,class_code,out_file_name):

    class_code_name = class_code
    if (class_code.find(".") ==1) :
        class_code_name="point"
    class_ids_file = "selected.%s.%s.ids" % ( tmap_out,class_code_name)

    try:
        comp_cmd = "awk 'NR>1 {if (($3 == \"%s\")) print $5}' %s > %s" % (class_code,tmap_out, class_ids_file)
        print >> run_log, comp_cmd
        ret = subprocess.call(comp_cmd, shell=True)
    # awk or not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: awk not found on this system.  Did you forget to include it in your PATH?"
        print o
        exit(1)

    class_ids=set ([])
    f = open(class_ids_file)
    for line in f:
        line = line.strip()
        class_ids.add(line)

    selected_ids=set([])
    selected_ids=potential_ids & class_ids
    select_gtf ( meta_asm_gtf,selected_ids, out_file_name)

    os.remove(class_ids_file)


# def subtract_short_transcripts(meta_asm_gtf, min_linc_length):
#     print >> sys.stderr, "[%s] Excluding short transfrags %s" % (right_now(), mask_gtf)
#
#     tmap_out = meta_asm_gtf.split("/")[-1] + ".tmap"
#     selected_ids = "%s.selected.ids" % tmap_out
#     selected_gtf = "%s.selected.gtf" % tmap_out
#
#     try:
#         comp_cmd = "awk 'NR>1 {if (($3 ~\"u|x\")) print $5}' %s > %s" % (tmap_out, selected_ids)
#         print >> run_log, comp_cmd
#         ret = subprocess.call(comp_cmd, shell=True)
#     # awk or not found
#     except OSError, o:
#         if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
#             print >> sys.stderr, fail_str, "Error: awk not found on this system.  Did you forget to include it in your PATH?"
#         print o
#         exit(1)
#
#     try:
#         comp_cmd = whereis("select_gtf.py") + " %s %s > %s;" % (selected_ids, meta_asm_gtf, selected_gtf)
#         print >> run_log, comp_cmd
#         ret = subprocess.call(comp_cmd, shell=True)
#         return selected_gtf
#     # awk or not found
#     except OSError, o:
#         if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
#             print >> sys.stderr, fail_str, "Error: select_gtf.py not found on this system.  Did you forget to include it in your PATH?"
#         exit(1)


def get_gtf_ids(gtf_filename):
    f_gtf = open(gtf_filename)
    ids = set([])
    for line in f_gtf:
        line = line.strip()
        transcript = get_gtf_line_id(line)
        if transcript != "":
            ids.add(transcript)
            
        #print line
        #cols = line.split('\t')
        #if len(cols) < 9:
         #   continue
        #if cols[2] != "exon":
         #   continue
        #attrs = cols[8]
        #attr_cols = attrs.split(';')

        #for col in attr_cols:
         #   if col.find("transcript_id") != -1:
          #      first_quote = col.find('"')
           #     last_quote = col.find('"', first_quote + 1)
            #    transcript = col[first_quote + 1:last_quote]
                #print transcript, t_len
                #ids.add(transcript)
    return ids

def get_gtf_line_id(line, transcript_id=True):
    cols = line.split('\t')
    transcript = ""
    id_ret = None
    if len(cols) >= 9 and cols[2] == "exon":
        attrs = cols[8]
        attr_cols = attrs.split(';')
        for col in attr_cols:
            if transcript_id == True:
                if col.find("transcript_id") != -1:
                    first_quote = col.find('"')
                    last_quote = col.find('"', first_quote + 1)
                    transcript = col[first_quote + 1:last_quote]
                    id_ret = transcript
            else:
                if col.find("gene_id") != -1:
                    first_quote = col.find('"')
                    last_quote = col.find('"', first_quote + 1)
                    gene = col[first_quote + 1:last_quote]
                    id_ret = gene
    return id_ret
    
def expand_ids_by_gene_grouping(gtf_filename, transcript_ids):
    #Returns expanded list of transcript_ids (all t_ids for a gene that includes and t_id in transcript_ids).
    f_gtf = open(gtf_filename)
    t_id_to_g_id = {}
    g_id_to_t_ids = {}

    expanded_ids = set(transcript_ids)
    for line in f_gtf:
        t_id = get_gtf_line_id(line)
        g_id = get_gtf_line_id(line, False)
        t_id_to_g_id[t_id] = g_id
        t_ids_for_g_id = g_id_to_t_ids.setdefault(g_id, [])
        t_ids_for_g_id.append(t_id)
        #g_id_to_t_ids[g_id] = t_ids_for_g_id

    for t_id in transcript_ids:
        g_id = t_id_to_g_id.get(t_id)
        if g_id != None:
            t_ids = g_id_to_t_ids[g_id] 
            expanded_ids |= set(t_ids)
        else:
            print >> sys.stderr, "Error: can't find %s" % g_id

    print >> sys.stderr, "Expanded list with %d ids to one with %d ids" % (len(transcript_ids), len(expanded_ids))
    return expanded_ids


def extend_gtf_line (line, new_keys, new_vals):
    str="";
    for i in range(len(new_keys)):
        k = new_keys[i]
        v = new_vals[i]
        v2=v
        if type(v) is int or type(v) is float:
            v2 = "%d"%v
        str += " " + k + " \"" + v2 + "\";"
    line += str
    return line    

def get_long_transcripts(gtf_filename, min_length):
    #Returns a set of ids for transcripts in gtf_filename that are >= min_length
    
    print >> sys.stderr, "[%s] Selecting long transcripts from meta-assembly" % (right_now())
    f_gtf = open(gtf_filename)
    transcript_lengths = {}
    for line in f_gtf:
        line = line.strip()
        #print line
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        if cols[2] != "exon":
            continue
        attrs = cols[8]
        attr_cols = attrs.split(';')

        for col in attr_cols:
            if col.find("transcript_id") != -1:
                first_quote = col.find('"')
                last_quote = col.find('"', first_quote + 1)
                transcript = col[first_quote + 1:last_quote]
                t_len = transcript_lengths.setdefault(transcript,0)
                t_len += int(cols[4]) - int(cols[3])
                #print transcript, t_len
                transcript_lengths[transcript] = t_len

    ids = []
    for (t_id, t_len) in transcript_lengths.iteritems():
        if t_len >= min_length:
            ids.append(t_id)
    return set(ids)

def get_multiexonic_transcripts(gtf_filename, min_exons):
    #Returns a set of transcript ids with >= min_exons number of exons.
    print >> sys.stderr, "[%s] Selecting multiexoninc transcripts from meta-assembly" % (right_now())
    f_gtf = open(gtf_filename)
    transcript_lengths = {}
    for line in f_gtf:
        line = line.strip()
        #print line
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        if cols[2] != "exon":
            continue
        attrs = cols[8]
        attr_cols = attrs.split(';')

        for col in attr_cols:
            if col.find("transcript_id") != -1:
                first_quote = col.find('"')
                last_quote = col.find('"', first_quote + 1)
                transcript = col[first_quote + 1:last_quote]
                t_len = transcript_lengths.setdefault(transcript,0)
                t_len += 1
                #print transcript, t_len
                transcript_lengths[transcript] = t_len

    ids = []
    for (t_id, t_len) in transcript_lengths.iteritems():
        if t_len >= min_exons:
            ids.append(t_id)
    return set(ids)



def gtf_to_bed(gtf_filename):
    f_gtf = open(gtf_filename)
    f_bed_filename = tmp_name("gtf2bed_")

    transcript_lengths = {}
    for line in f_gtf:
        line = line.strip()
        #print line
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        if cols[2] != "exon":
            continue
        attrs = cols[8]
        attr_cols = attrs.split(';')

        for col in attr_cols:
            if col.find("transcript_id") != -1:
                first_quote = col.find('"')
                last_quote = col.find('"', first_quote + 1)
                transcript = col[first_quote + 1:last_quote]
                t_lines = transcript_lengths.setdefault(transcript,[])
                t_lines.append(cols)

    f_bed = open(f_bed_filename, "w")

    for (t_id, t_lines) in transcript_lengths.iteritems():
        t_exons = []
        chrom = t_lines[0][0]
        t_start = int(t_lines[0][3]) - 1
        t_strand = t_lines[0][6]
        t_lines.sort(lambda x,y: int(x[3]) - int(y[3]))
        for cols in t_lines:
            exon_start = int(cols[3]) - 1
            exon_end = int(cols[4]) - 1
            exon_len = exon_end - exon_start + 1
            exon_rel_start = exon_start - t_start
            t_exons.append((exon_rel_start, exon_len))
        t_end = t_start + t_exons[-1][0] + t_exons[-1][1]
        t_block_starts = ",".join([str(x[0]) for x in t_exons])
        t_block_sizes = ",".join([str(x[1]) for x in t_exons])
        #print t_exons
        print >> f_bed, "%s\t%s\t%s\t%s\t%d\t%c\t%d\t%d\t0\t%d\t%s\t%s" % (chrom, t_start, t_end, t_id, 1000, t_strand, t_start, t_end, len(t_exons), t_block_sizes, t_block_starts)

    return f_bed_filename


def split_bed(input_bed_filename, lines_per_sub_bed, max_files):
    if lines_per_sub_bed == 0:
        return [input_bed_filename]
    input_bed = open(input_bed_filename, "r")

    i = 1
    curr_sub_bed_filename = tmp_name( "split_bed_")
    sub_beds = [curr_sub_bed_filename]
    curr_sub_bed = open(curr_sub_bed_filename, "w")

    for line in input_bed:
        line = line.strip()
        cols = line.split()
        if len(cols) != 12:
            continue
        print >> curr_sub_bed, line
        if i % lines_per_sub_bed == 0 and len(sub_beds) != max_files:
            curr_sub_bed_filename = tmp_name("split_bed_")
            sub_beds.append(curr_sub_bed_filename)
            curr_sub_bed = open(curr_sub_bed_filename, "w")
        i = i + 1
    return sub_beds

def count_fasta_recs(fasta_filename):
    f_fasta = open(fasta_filename)
    i = 0
    for line in f_fasta:
        line = line.strip()
        if line.startswith(">"):
            i = i + 1
    return i

def split_fasta(fasta_filename, records_per_fasta_file, max_files):
    f_fasta = open(fasta_filename)

    i = 0
    sub_fastas = []
    curr_sub_fasta = None

    for line in f_fasta:
        line = line.strip()
        if line.startswith(">"):
            if i % records_per_fasta_file == 0 and len(sub_fastas) < (max_files):
                curr_sub_fasta_filename = tmp_name("split_fasta_")
                #print curr_sub_fasta_filename
                sub_fastas.append(curr_sub_fasta_filename)
                curr_sub_fasta = open(curr_sub_fasta_filename, "w")
            i = i + 1
            print >> curr_sub_fasta, line
        else:
            print >> curr_sub_fasta, line
    print len(sub_fastas), max_files
    return sub_fastas

def calculate_csf(bed_filename, num_lsf_nodes, genome):
    print >> sys.stderr, "[%s] Calculating CSF scores for possible lincs" % (right_now())

    #csf_bin = "~mlin/bin/batchPhyloCSF"
    #csf_bin = "/rinn_data1/user-supported/PhyloCSF/batchPhyloCSF_Loyal.ml.exe"
    min_csf_orf = 30

    num_bed_lines = len(open(bed_filename).readlines())

    lines_per_chunk = (num_bed_lines / num_lsf_nodes)

    sub_beds = split_bed(bed_filename, lines_per_chunk, num_lsf_nodes)
    print "processing %d bed files with ~%d transcripts per file" % (len(sub_beds), lines_per_chunk)
    jobs = []
    outfiles = {}
    global lsf_queue
    for bed in sub_beds:
        csf_score_filename = tmp_name("csf_")

        ## If using LSF for the CSF calculation, uncomment these lines
        if params.system_params.system=="Broad":
            csf_cmd = csf_bin + ' -sp %s -strategy GuessLik -bestORF %d  -best3  -bed %s' % (genome, min_csf_orf, bed)
        else:
            #edited for updated batchPhyloCSF on Valor LAG 3-14-11
            csf_cmd = csf_bin + ' -s %s -p"--minCodons %d -f3" -f BED --exe=/n/rinn_data1/user-supported/PhyloCSF/PhyloCSF %s' % (genome, min_csf_orf, bed)
        
        #Cole's original Method
        #csf_bsub_cmd = bsub_cmd(csf_cmd, "/csf_calc", blocking=True, outfilename=csf_score_filename, job_cores=1, job_mem=2)
        #p = subprocess.Popen(csf_bsub_cmd, shell=True)
        #
        #while p.poll() != None:
        #    print >> sys.stderr,"Resubmitting Job"
        #    p = subprocess.Popen(cmd,shell = True)
        #    time.sleep(5)
        
        #Using LSFlib
        csf_bsub_job = LSFlib.LSFJob(csf_cmd,job_group = "/csf_calc",blocking=False, notify = False, queue_name=lsf_queue, outfilename = csf_score_filename,job_cores=1,job_mem=4)

        ret_val=-1
        while ret_val != 0:
            try:
                ret_val = csf_bsub_job.submit()
            except Exception , e:
                print >> sys.stderr, 'Submit Exception %s\n' %e
                print >>sys.stderr, "Resubmitting job %s" % str(csf_bsub_job)
                time.sleep(10)
                
        # Otherwise, just use the generic version
        # csf_cmd = csf_bin + ' -sp hg19 -strategy GuessLik -bestORF %d  -best3  -bed %s' % (min_csf_orf, bed)
        #  csf_bsub_cmd = csf_cmd
        #  csf_tmp_filename = tmp_name()
        #  #print csf_bsub_cmd, ">", csf_score_filename, "2>", csf_tmp_filename
        #  p = subprocess.Popen(csf_bsub_cmd, stdout=open(csf_score_filename,"w"), stderr=open(csf_tmp_filename, "w"), shell=True)
        #
        outfiles[csf_score_filename] = csf_bsub_job
        print >> run_log, str(csf_bsub_job)

        jobs.append(csf_bsub_job)

    #for j in jobs:
    #    ret = j.wait()
        #assert ret == 0
    check_lsf_output(outfiles)
    
    #time.sleep(60)
    csf_score_filename = tmp_name("csf_merged_")
    csf_score_file = open(csf_score_filename, "w")
    print >> run_log, "#cat *.csf > %s" % csf_score_filename
    for outfilename in outfiles.keys():
        try:
            outfile = open(outfilename)
        except IOError:
            print >>run_log, "re-running job: %s" % outfiles[outfilename]
            ret_val = outfiles[outfilename].submit()
            ret = outfiles[outfilename].wait()
            time.sleep(60)
            outfile = open(outfilename)
        
        for line in outfile:
            print >> csf_score_file, line,
            
    return csf_score_filename

#Return : negative_hits,positive_hits
def get_negative_csf_transfrags(gtf_filename, num_lsf_nodes, annotated_gtf_f_name, genome):
    global csf_score_threshold
    print >> sys.stderr, "[%s] Screening lincRNAs for codon-level conservation using CSF threshold = %d" % (right_now(),csf_score_threshold)
    bed_filename = gtf_to_bed(gtf_filename)
    csf_scorefilename = calculate_csf(bed_filename, num_lsf_nodes, genome)
    csf_scorefile = open(csf_scorefilename)
    csf_hash = dict ([]);
    
    positive_csf_ids = set([])
    for line in csf_scorefile:
        line = line.strip()
        cols = line.split('\t')
        #print len(cols)
        if len(cols) < 15:
            continue
        csf_score = float(cols[12])
        csf_start = float(cols[13])
        csf_end = float(cols[14])
        csf_hash[cols[3]]= [csf_score,csf_start,csf_end]
        # print cols[3], cols[12]
        if csf_score >= csf_score_threshold:
            positive_csf_ids.add(cols[3])
            
    positive_csf_ids = expand_ids_by_gene_grouping(gtf_filename, positive_csf_ids)
    
    all_hits = get_gtf_ids(gtf_filename)
    negative_hits = all_hits - positive_csf_ids
    positive_hits = all_hits & positive_csf_ids
    
    #Generate a new version of the input GTF that has CSF annotations
    #Step 1 - generate a dictinary of maximal csf per transcript
    #csf_gene_hash = dict ([])
    #for t_id in csf_hash.keys():
    #    csf_vals = csf_hash[t_id]
    #    g_id = t_id
    #    dot = t_id.rfind(".")
    #    if dot != -1:
    #        g_id = t_id[:dot] 
    #    if csf_gene_hash.has_key(g_id):
    #        curr_csf_vals = csf_gene_hash[g_id]
    #        if (curr_csf_vals.get(0) > csf_vals.get(0)):
    #            csf_vals = curr_csf_vals
    #    csf_gene_hash[g_id] = csf_vals
    
    #Step 2 - go through the gtf and extend
    gtf_file = open(gtf_filename)
    new_gtf_file = open(annotated_gtf_f_name, 'w')
    new_keys = ['csf','csf_start','csf_end']
    for line in gtf_file:
        new_vals = [0,0,0]
        line = line.strip()
        t_id = get_gtf_line_id(line)
        if csf_hash.has_key(t_id):
            new_vals = csf_hash[t_id]
        new_line = extend_gtf_line(line, new_keys, new_vals)
        print >> new_gtf_file, new_line 
    
    return negative_hits,positive_hits


def get_pfam_hit_transcripts(gtf_filename, fasta, num_lsf_nodes, annotated_gtf_f_name):
    print >> sys.stderr, "[%s] Screening possible lincs for hits to Pfam" % (right_now())

    dna_fastafilename = tmp_name("pfam_dna_")

    # Get fasta sequences for possible_lincs
    cmd = gffread + " %s -g %s -w %s" % (gtf_filename, fasta, dna_fastafilename)
    print >> run_log, cmd
    ret = subprocess.call(cmd, shell=True)

    num_fasta_recs = count_fasta_recs(dna_fastafilename)
    recsPerNode= (num_fasta_recs / num_lsf_nodes);
    print >> sys.stderr, "%d fasta recs per node " % recsPerNode

    if  recsPerNode > 0:
        sub_fastas = split_fasta(dna_fastafilename,  recsPerNode, num_lsf_nodes)
    else:
        sub_fastas = [dna_fastafilename]

    pfam_jobs = []
    hmmer_out_filenames = {}

    #print len(sub_fastas)
    global lsf_queue
    for fasta in sub_fastas:
        aa_fastafilename = tmp_name("pfam_aa_")

        # Convert these to amino acid sequences in all three frames
        #cmd = aa_converter + " %s %s" % (dna_fastafilename, aa_fastafilename) % this was an error  (big one:)
        cmd = aa_converter + " %s %s" % (fasta, aa_fastafilename)
        print >> run_log, cmd
        ret = subprocess.call(cmd, shell=True)
        
        hmmer_out_filename = tmp_name("pfam_hmm_") #Key for dictionary of LSFJobs
        
        cmd = pfammer + " -fasta %s -dir %s -outfile %s  -e_seq 0.1 -pfamB -cpu %d" % (aa_fastafilename, pfam_db, hmmer_out_filename, 1)
        ##############BROAD###############
        
        #Cole's method
        #
        #cmd = bsub_cmd(cmd, "/pfam_search", blocking=True, job_cores=1, job_mem=2)
        #
        #hmmer_out_filenames[hmmer_out_filename] = cmd
        #
        
        #Using LSFlib for Job sumission
        #pfam_job = LSFlib.LSFJob(cmd,job_group="/pfam_search",outfilename = hmmer_out_filename, notify = False, queue_name=lsf_queue, blocking=False,job_cores=1,job_mem=4)
        pfam_job = LSFlib.LSFJob(cmd,job_group="/pfam_search", notify = False, queue_name=lsf_queue, blocking=False,job_cores=1,job_mem=4)
        #Contains file names as keys and jobs as values
        hmmer_out_filenames[hmmer_out_filename] = pfam_job
        
        print >> run_log, "[%s] Submitting Job: %s" % (right_now(),str(pfam_job))
        
        ret_val=-1
        while ret_val != 0:
            try:
                ret_val = pfam_job.submit()
            except Exception , e:
                print >> sys.stderr, 'Submit Exception %s\n' %e
                print >>sys.stderr, "Resubmitting job %s" % str(pfam_job)
                time.sleep(10)
        
        #Cole's method (prior to LSFLib)
        #Wait and then poll to make sure job was submitted
        #while p.poll() != None:
        #    print >> sys.stderr,"Resubmitting Job"
        #    p = subprocess.Popen(cmd,shell = True)
        #    time.sleep(5)
        
        pfam_jobs.append(pfam_job)
    
    #This is my 'wait' for all jobs to complete
    check_lsf_output(hmmer_out_filenames)
    
    hmmer_hits = set([])

    #Store hits info
    pfam_hit_info = dict()
    for hmmer_out_filename in hmmer_out_filenames.keys():
        #try:
        hmmer_out = open(hmmer_out_filename)
        #except IOError:
        #    ret_val = hmmer_out_filenames[hmmer_out_filename].submit()
        #    ret=p.wait()
        #    time.sleep(60)
        #    hmmer_out = open(hmmer_out_filename)
        for line in hmmer_out:
            line = line.strip()
            cols = line.split()
            if len(cols) < 14:
                continue
            sig = cols[13]
            domain = cols[6]
            if sig == "1" or sig == "NA" :
                uscore = cols[0].rfind("_")
                if uscore != -1:
                    t_id = cols[0][:uscore]
                    hmmer_hits.add(t_id)
                    oldDat=""
                    if pfam_hit_info.has_key(t_id):
                        oldDat = pfam_hit_info[t_id]
                    newDat = oldDat + domain + "_" +sig +","
                    pfam_hit_info[t_id] = newDat
                else:
                    print >> sys.stderr, "Warning: could not parse meta-assembly transcript id %s" % cols[0]

    hmmer_hits = expand_ids_by_gene_grouping(gtf_filename, hmmer_hits)
    
    # for hit in hmmer_hits:
    #     print hit
    all_hits = get_gtf_ids(gtf_filename)
    negative_hits = all_hits - hmmer_hits
    positive_hits = all_hits & hmmer_hits
    
    
    #Generate a new version of the input GTF that has Domain annotations
    #Step 1 - generate a dictinary of maximal csf per transcript
    #pfam_gene_hash = dict()
    #for t_id in pfam_hit_info.keys():
    #    pfam_vals = pfam_hit_info[t_id]
    #    g_id = t_id
    #    dot = t_id.rfind(".")
    #    if dot != -1:
    #        g_id = t_id[:dot]
    #    if pfam_gene_hash.has_key(g_id):
    #        curr_vals = pfam_gene_hash[g_id]
    #        pfam_vals = curr_vals + pfam_vals
    #    pfam_gene_hash[g_id] = pfam_vals
    
    #Step 2 - go through the gtf and extend
    gtf_file = open(gtf_filename)
    new_gtf_file = open(annotated_gtf_f_name, 'w')
    new_keys = ['PFAM_Domains']
    for line in gtf_file:
        new_vals = ["None"]
        line = line.strip()
        t_id = get_gtf_line_id(line)
        if pfam_hit_info.has_key(t_id):
            new_vals = [pfam_hit_info[t_id]]
        new_line = extend_gtf_line (line, new_keys, new_vals)
        print >> new_gtf_file, new_line 
        
    
    return (negative_hits,positive_hits) #Returns gene-level expansion of transcript_ids for transcripts with a Pfam domain (in positive). Negative is 'none of the transcripts from a given gene have PFam hits'.

#THIS WAS TRAPPED IN AN INFINITE LOOP
def check_lsf_output(fileDict):
    """Checks a dictionary of k=ouput_file_names v=LSFlib.LSFJob Instances"""
    print >> sys.stderr, "Moran debug: check LSF output"
    reruns = {}
    for fname in fileDict.keys():
        try:
            fileDict[fname].wait()
        except:
            time.sleep(30)
            fileDict[fname].wait()
        if not os.path.exists(fname):
            print >>sys.stderr, "Waiting 60 sec..."
            time.sleep(60)
        if not os.path.exists(fname):
            print >>sys.stderr, "Resubmitting job %s" % str(fileDict[fname])
            reruns[fname] = fileDict[fname]
            reruns[fname].kill()
            reruns[fname].submit()
    if len(reruns.keys()) == 0: #Base case
        return True 
    else:
        check_lsf_output(reruns) #recursive call

def check_lsf_output2(jobList):
    #base case
    newList = []
    for job in jobList:
        if job.poll() != 'DONE':
            newList.append(job)
    if len(newList) == 0:
        return
    else:
        time.sleep(60)
        check_lsf_output2(newList)


def annot_overlap_transcripts (gtf_filename,overlap_filename, gtf_field, annotated_gtf_f_name):
    all_ids = get_gtf_ids(gtf_filename)
    masked_gtf = subtract_mask_transcripts(gtf_filename, overlap_filename)
    non_overlap_ids = get_gtf_ids(masked_gtf)
    overlap_ids = all_ids - non_overlap_ids 
    overlap_ids = expand_ids_by_gene_grouping(gtf_filename, overlap_ids)
    non_overlap_ids = all_ids - overlap_ids
    
    gtf_file = open(gtf_filename)
    new_gtf_file = open(annotated_gtf_f_name, 'w')
    new_keys = [gtf_field]
    for line in gtf_file:
       new_vals = [0]
       line = line.strip()
       t_id = get_gtf_line_id(line)
       if t_id in overlap_ids:
           new_vals = [1]
       new_line = extend_gtf_line (line, new_keys, new_vals)
       print >> new_gtf_file, new_line 
    
    return (overlap_ids , non_overlap_ids) 
        

def select_gtf(gtf_in_filename, ids, gtf_out_filename):
    f_gtf = open(gtf_in_filename)
    #print >> sys.stderr, "Select GTF: Ids are: "
    #print >> sys.stderr, ids
    #print >> sys.stderr, "reference gtf file name:"
    #print >> sys.stderr, gtf_in_filename
    out_gtf = open(gtf_out_filename, "w")
    for line in f_gtf:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        attrs = cols[8]
        attr_cols = attrs.split(';')
        for col in attr_cols:
            if col.find("transcript_id") != -1:
                first_quote = col.find('"')
                last_quote = col.find('"', first_quote + 1)
                transcript = col[first_quote + 1:last_quote]
                #print >> sys.stderr, transcript
                if transcript in ids:
                    print >> out_gtf, line

# Originally from Loyal Goff's naming script
overlapThreshold = 0.20
extensionLength = 500 #grow 5'end of lincRNA by this many bases to test for Bidirectional promoter
strandLookup = {'+':'-','-':'+'}

def test5PrimeOverlap(lincInt,geneInt):
    """May need to validate this.  I'm not sure this works when a lincRNA completely covers a PC gene on the opposite strand"""
    assert lincInt.overlaps(geneInt)
    if lincInt.strand == "+":
        if lincInt.start <= geneInt.end and lincInt.end > geneInt.end:
            return True
        else:
            return False
    elif lincInt.strand == "-":
        if geneInt.start <= lincInt.end and geneInt.end > lincInt.end:
            return True
        else:
            return False
    else:
        raise ValueError("Could not determine")

def bpOverlap(lincInt,geneInt):
    assert lincInt.overlaps(geneInt)
    bounds = [lincInt.start,lincInt.end,geneInt.start,geneInt.end]
    bounds.sort()
    #range = bounds[3]-bounds[0]
    overlap = bounds[2]-bounds[1]
    return overlap

def printLincs(handle,lincs):
    for linc in lincs:
        print >>handle, linc.getGTF(),

def name_lincs(gtfFile, genome, verbose=False):
    print >> sys.stderr, "[%s] Attaching names to lincRNAs" % (right_now())

    #Parse GTF File for lincs
    lincIter = newGTF.GTFGeneIterator(gtfFile,verbose=verbose)

    #Retrieve and index RefSeq genes
    refSeqs = dbConn.fetchRefSeqIntervalsIndexed(genome=genome,proteinCodingOnly=True,verbose=verbose)

    #Results container
    res = set([])

    #Container for gene:linc assoc.
    geneLincs = {}

    #Loop through lincRNAs
    for linc in lincIter:
        flag = False
        bdFlag = False #True if linc is bidirectional
        asFlag = False #True if linc is antisense
        #Convert to Interval
        interval = linc.toInterval()

        #Test for weird chromosome (ie. not in refSeqs.keys() )
        if not interval.chr in refSeqs.keys():
            res.add(linc)
            continue

        #Bug tracking only
        if verbose:
            sys.stderr.write(str(interval)+"\n")

        #Get list of gene positions that are relevant
        senseGeneStarts = [x.start for x in refSeqs[interval.chr][interval.strand]]
        senseGeneEnds = [x.end for x in refSeqs[interval.chr][interval.strand]]

        #Get opposite strand to test
        testStrand = strandLookup[interval.strand]

        #Test overlap with genes on opposite strand
        for gene in refSeqs[interval.chr][testStrand]:
            extendedInterval = copy.copy(interval)
            extendedInterval.grow5_prime(extensionLength)
            
            if extendedInterval.overlaps(gene):
                #If 5' end of linc overlaps the 5' of a coding gene on the opposite strand,
                #by more than 0bp but less than min(BP_THRESH * length(L), BP_THRESH * length(coding gene))
                #THEN name linc "linc-[HUGO_GENE_NAME]-BP"
                overlap = bpOverlap(extendedInterval,gene)
                fivePrime = test5PrimeOverlap(extendedInterval,gene)
                cutoff = min(len(extendedInterval)*overlapThreshold,gene.intervalLen()*overlapThreshold)
                if fivePrime and overlap <= cutoff:
                    linc.propogateLincName("linc-%s-BP" % gene.name)
                    linc.addAttribute("bidirectional_prom",gene.name)
                    res.add(linc)
                    flag = True
                    bdFlag = True
                    #break
                    continue
                
                #TODO FIX this so that ANY overlap that is not a BP becomes and -AS
                if not bdFlag:
                    linc.propogateLincName("linc-%s-AS" % gene.name)
                linc.addAttribute("antisense",gene.name)
                res.add(linc)
                flag = True
                asFlag = True
                break

        #ELSE find the closest coding gene on the same strand as the L, starting from the 3' end of the linc.
        #Suppose its HUGO name is NCG1.Add L to a list of lincs to be named after NCG1.
        if not flag:
            if interval.strand == "+":
                nearestGeneIdx = bisect.bisect(senseGeneStarts,interval.end) #choose most adjacent gene 3' to lincRNA
            elif interval.strand == "-":
                nearestGeneIdx = bisect.bisect(senseGeneEnds,interval.start)-1
            try:
                nearestGene = refSeqs[interval.chr][interval.strand][nearestGeneIdx]
            except IndexError:
                #If I cannot find the nearestGene (e.g. end of chromosome or something, just push linc to results
                #and deal with them later. (for now)

                #print nearestGeneIdx
                #print interval.toBed()
                res.add(linc)
                continue
            geneLincs.setdefault(nearestGene.name,[]).append(linc)

    #Evaluate container for linc:gene assocs
    """
    FOREACH coding gene G in the table above:
    IF there's only one linc to be named after G THEN
     name that linc "linc-G"
    ELSE
     sort the list of lincs by proximity to G, with the closest linc at the front of the list
     FOR i = 1 to #number of lincs named after G
         name linc i "linc-G-i"
    """
    for k,v in geneLincs.iteritems():
        if len(v) == 1:
            v[0].propogateLincName("linc-%s" % (k))
            res.add(v[0])
        elif len(v) >1:
            if v[0].strand == "+":
                v.sort(reverse=True)
            elif v[0].strand == "-":
                v.sort()
            for i in xrange(len(v)):
                v[i].propogateLincName("linc-%s-%d" % (k,i+1))
                res.add(v[i])
    return res

def merge_gtfs(gtf_filenames, merged_gtf, ref_gtf=None):
    print >> sys.stderr, "[%s] Merging linc gtf files with cuffcompare" % (right_now())
    cmd = ["cuffcompare"]

    cmd.extend(["-o", merged_gtf])
    if ref_gtf != None:
        cmd.extend(["-r", ref_gtf])

    cmd.extend(gtf_filenames)
    cmd = " ".join(cmd)
    #cmd = bsub_cmd(cmd, "/merge_gtf", True, job_mem=8)

    try:
        print >> run_log, cmd
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
        return merged_gtf + ".combined.gtf"
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def gtf_to_genbank(gtf_filename, genbank_out, genome):
    print >> sys.stderr, "[%s] Building Genbank records for linc catalog" % (right_now())

    cmd = ["GTF2Genbank.py", "-g", genome, "-o", genbank_out, gtf_filename]
    if params.system_params.system =="Broad":
        cmd = ["GTF2Genbank.py", "-g", genome, "-s", "Broad", "-o", genbank_out, gtf_filename]
    #print cmd
    GTF2Genbank.main(cmd)
    
    # try:
    #     print >> run_log, " ".join(cmd)
    #     ret = subprocess.call(cmd)
    #     if ret != 0:
    #         print >> sys.stderr, fail_str, "Error: could not execute GTF2Genbank.py"
    #         exit(1)
    #     print >> sys.stderr
    #     return genbank_out
    # # GTF2Genbank.py not found
    # except OSError, o:
    #     if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
    #         print >> sys.stderr, fail_str, "Error: GTF2Genbank.py not found on this system.  Did you forget to include it in your PATH?"
    #     exit(1)

def cuffcompare_all_assemblies(gtf_input_files):
    print >> sys.stderr, "Cuffcmpare all assemblies GTFs"
    cmd = ["cuffcompare"]
    gtf_lst=" ".join(gtf_input_files)
    cuffcompare(None, "cuffcmpr1" , None,None, gtf_input_files)
    print >> sys.stderr, "Cuffcmpare all assemblies GTFs : run second cuffcompare"
    cuffcompare(None, "cuffcmpr2" , "cuffcmpr1.combined.gtf",None, "cuffcmpr1.combined.gtf")
    tmap="cuffcmpr2.cuffcmpr1.combined.gtf.tmap" #!!!!!! CUFFCMPR CHANGES ITS OUT PUT STRING !!!#
    print >> sys.stderr, "Cuffcmpare all assemblies GTFs : filter c"
    selected_ids= set([])
    f_tmap = open(tmap)
    out = open("cuffcmpr1_selectedIds.txt", "w")
    for line in f_tmap:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 5:
            continue
        class_code = cols[2]
        name =cols[4]
        if class_code != "c" :
            print >> out, name
        selected_ids.add(name)

    global output_dir
    asm_dir = output_dir + "merged_asm/"
    print >> sys.stderr, asm_dir
    if os.path.exists(asm_dir):
        pass
    else:
        os.mkdir(asm_dir)
    current_asm_gtf = "cuffcmpr1.combined.gtf"
    select_gtf(current_asm_gtf, selected_ids, "merged_asm/transcripts.gtf")


def run_non_assembly_mode(current_asm_gtf,params):

    current_asm_ids=get_gtf_ids(current_asm_gtf)
    unmasked_transcript_ids = current_asm_ids

    ##Testing
    #tmap_out = current_asm_gtf.split("/")[-1] + ".tmap"
    #lincRNA_ids=get_gtf_ids("possible_lincs_final.gtf")
    #print >> sys.stderr, "Extracting intergenic transcripts"
    #get_class_transcripts(current_asm_gtf,lincRNA_ids,tmap_out,"u","lincs_final.gtf")
    #return
    ##testing

    # If a mask GTF file is provided, cuffcompare against it to import reference metadata
    if params.mask_gtf != None:
       unmasked_gtf = subtract_mask_transcripts(current_asm_gtf, params.mask_gtf)
       unmasked_transcript_ids = get_gtf_ids(unmasked_gtf)
       tmap_out = "meta_asm_vs_mask." + current_asm_gtf.split("/")[-1] + ".tmap"

    long_transcript_ids = set([])
    # Exclude possible noncoding transcripts shorter than 250bp or which are unspliced.
    # If a mask GTF file is provided, cuffcompare against it to import reference metadata
    if params.min_linc_length > 0:
        long_transcript_ids = get_long_transcripts(current_asm_gtf, params.min_linc_length)
    min_exon_transcript_ids = set([])
    # Exclude possible noncoding transcripts shorter than 250bp or which are unspliced.
    # If a mask GTF file is provided, cuffcompare against it to import reference metadata
    if params.min_linc_exons > 1:
        min_exon_transcript_ids = get_multiexonic_transcripts(current_asm_gtf, params.min_linc_exons)

    lincRNA_ids = current_asm_ids

    lincRNA_ids &= long_transcript_ids
    print >> sys.stderr, "Current lincIDS after masking SHORT %d" % len(lincRNA_ids)
    lincRNA_ids &= unmasked_transcript_ids
    print >> sys.stderr, "Current lincIDS after masking CODING %d"  % len(lincRNA_ids)
    lincRNA_ids &= min_exon_transcript_ids
    print >> sys.stderr, "Current lincIDS after masking single EXON %d" % len(lincRNA_ids)

    #this file was empty in last test run!
    possible_lincs_filename = "possible_lincs.gtf"
    select_gtf(current_asm_gtf, lincRNA_ids, possible_lincs_filename)
    no_pfam_hit_ids = set([])

    ##** Need to 
    if params.fasta != None and 1==1:
        (no_pfam_hit_ids , positive_pfam_hit_ids) = get_pfam_hit_transcripts(possible_lincs_filename, params.fasta, params.system_params.lsf_nodes,"pfam_annot.gtf")
        lincRNA_ids &= no_pfam_hit_ids
        possible_lincs_filename2 = "possible_lincs_noPfam.gtf"
        select_gtf(possible_lincs_filename, lincRNA_ids, possible_lincs_filename2)
        # Calculate CSF score for possible long noncoding transcripts
        possible_lincs_filename2_csf = "possible_lincs_noPfam_extendCsf.gtf"
        (negative_csf_ids,positive_csf_ids) = get_negative_csf_transfrags(possible_lincs_filename2, params.system_params.lsf_nodes,possible_lincs_filename2_csf,params.genome)
        lincRNA_ids &= negative_csf_ids

    print >> sys.stderr, "Found %d possible lincRNAs" % len(lincRNA_ids)
    # Extract the final catalog of lincRNAs
    possible_lincs_final_filename="possible_lincs_final.gtf"
    select_gtf(current_asm_gtf, lincRNA_ids, possible_lincs_final_filename)
    #retrieve just lincs (not lancs)
    print >> sys.stderr, "Extracting intergenic transcripts"
    get_class_transcripts(current_asm_gtf,lincRNA_ids,tmap_out,"u","lincs_final.gtf")

    print >> sys.stderr,"-----------------------------------------------"
    print >> sys.stderr, "Trancriptome catalog build complete "

def run_merge_catalogs( merge_gtf_list,out_prefix, genome ):
    merged_lincs = merge_gtfs(merge_gtf_list, "merge_catalog")
    print >> sys.stderr,"name of merged file %s", merged_lincs
    # Run Loyal's linc naming script to annotate the lincRNA catalog with sensible names
    lincs = name_lincs(merged_lincs, genome=genome, verbose=False)
    fgtf="%s_linc_catalog.gtf" % out_prefix
    fgbk="%s_linc_catalog.gbk" % out_prefix
    printLincs(open(fgtf, "w"), lincs)
    gtf_to_genbank(fgtf, fgbk, genome)

def run_get_intergenic(current_asm_gtf,params,out_name):
    unmasked_gtf = subtract_mask_transcripts(current_asm_gtf, params.mask_gtf)
    unmasked_transcript_ids = get_gtf_ids(unmasked_gtf)
    tmap_out = "meta_asm_vs_mask." +  current_asm_gtf.split("/")[-1] + ".tmap"
    get_class_transcripts(current_asm_gtf,unmasked_transcript_ids,tmap_out,"u",out_name)
    #get_class_transcripts(current_asm_gtf,unmasked_transcript_ids,tmap_out,".",out_name+"2")
    #cmd="cat %s >> %s" % (out_name+"2",out_name)
    #os.system (cmd)
    #os.remove(out_name+"2")

    if 0:
     os.remove(tmap_out)
     pre = current_asm_gtf.split("/")[-1]
     if os.path.exists(pre+".refmap"):
            os.remove(pre+".refmap")
     if os.path.exists(pre+".tmap.selected.gtf"):
            os.remove(pre+".tmap.selected.gtf")
     if os.path.exists(pre+".tmap.selected.ids"):
            os.remove(pre+".tmap.selected.ids")
     os.system("rm meta_asm_vs_mask*")

    bed_tmp_name=gtf_to_bed(out_name)
    bed_name=out_name.replace("gtf","bed")
    str="cp %s %s" % (bed_tmp_name,bed_name)
    os.system(str)
    os.remove(bed_tmp_name)


def meta_assemble_gtfs(params, gtf_filename_manifest, ref_gtf=None):
    print >> sys.stderr, "[%s] Merging linc gtf files with cuffmerge" % (right_now())
    cmd = ["cuffmerge"]

    
    cmd.extend(["-p", str(params.system_params.threads)])
    if ref_gtf != None:
        cmd.extend(["-g", ref_gtf])

    cmd.extend(["-F", "0.20"])
        
    cmd.append(gtf_filename_manifest)
    #print cmd    
    cmd = " ".join(cmd)
    #cmd = bsub_cmd(cmd, "/merge_gtf", True, job_mem=8)

    try:
        print >> run_log, cmd
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffmerge"
            exit(1)
        return "merged_asm/merged.gtf"
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)

def cuffmerge(params,gtf_filename_manifest,ref_gtf=None,outDir="./tmp"):
    print >> sys.stderr, "[%s] Merging linc gtf files with cuffmerge" % (right_now())
    cmd = ["cuffmerge"]
    
    cmd.extend(["-p", str(params.system_params.threads)])
    if ref_gtf != None:
        cmd.extend(["-g", ref_gtf])
    
    cmd.extend(["-o", outDir])
    
    cmd.extend(["-F", "0.20"])
        
    cmd.append(gtf_filename_manifest)
    #print cmd    
    cmd = " ".join(cmd)
    #cmd = bsub_cmd(cmd, "/merge_gtf", True, job_mem=8)

    try:
        print >> run_log, cmd
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffmerge"
            exit(1)
        return outDir+"/merged.gtf"
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)


def main(argv=None):
    warnings.filterwarnings("ignore", "tmpnam is a potential security risk")
    
    # Initialize default parameter values
    global params
    params = TestParams()
    
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()
        
        if params.system_params.system == "Broad": 
            global lsf_queue
            global cufflinks_queue
            global pfam_db
            global aa_converter
            global gffread
            global pfammer
            global csf_bin
            
            csf_bin = "~mlin/bin/batchPhyloCSF"
            #lsf_queue = "regevlab"
            #cufflinks_queue="regevlab"
            pfam_db = "/ahg/regev/users/nmcabili/Programs/Hmmer/pfam"
            aa_converter =  "perl /ahg/regev/users/nmcabili/Src/scripts/pfamAnalysis/seqToAminoAcid.pl"
            gffread = "gffread"
            pfammer = "perl /ahg/regev/users/nmcabili/Programs/Hmmer/PfamScan/pfam_scan.pl"
        
        
        global run_log
        global run_cmd
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning transcriptome build (suite v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------"
        print >> sys.stderr
        print >> sys.stderr, "  Configuration:"
        print >> sys.stderr, "     system: %s" % params.system_params.system
        print >> sys.stderr, "     genome: %s" % params.genome
        print >> sys.stderr, "     output_dir: %s" % output_dir
        print >> sys.stderr, "     logging_dir: %s" % logging_dir
        print >> sys.stderr, "     run log: %s" % run_log
        #print >> sys.stderr, "     Meta assembly: %s" % ("yes" if run_meta_assembly else "no")
        #print >> sys.stderr, "     Merged catalog: %s" % ("yes" if merge_catalog_mode else "no")
        print >> sys.stderr, "     CSF threshold: %f" % csf_score_threshold
        print >> sys.stderr, "     LSF nodes: %d" % params.system_params.lsf_nodes
        if lsf_mem != None:
            print >> sys.stderr, "     LSF memory per node: %d" % lsf_mem
        else:
            print >> sys.stderr, "     LSF memory per node: Not specified" 
        print >> sys.stderr
        
        start_time = datetime.now()
        prepare_output_dir()
        #now=datetime.datetime.now()
        #curr_time=now.strftime("%Y%m%d_%H%M")

        
        run_log = open(logging_dir + "run.log", "w", 0)
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd

        if len(args) < 1:
            raise Usage(use_message)

        transfrag_list_filename = args[0]
        transfrag_list_file = open(args[0], "r")

        if params.ref_gtf != None:
            test_input_files([params.ref_gtf])
        else:
            print >> sys.stderr, "Warning: no reference GTF provided!"
        if params.fasta != None:
            if os.path.isdir(params.fasta) == False:
                test_input_files([params.fasta])

        if params.mask_gtf != None:
            test_input_files([params.mask_gtf])

        if params.external_linc_gtf != None:
            test_input_files([params.external_linc_gtf])

        #############33
        
        
        global non_assembly_mode
        if non_assembly_mode==1:
            known_transcript_gtf=args[0]
            print >> sys.stderr, "Running non assembly mode on %s \n" % (known_transcript_gtf)
            run_non_assembly_mode(known_transcript_gtf,params)
            return 0
        global merge_catalog_mode
        if merge_catalog_mode==1:
            gtf_lst_file_name=args[0]
            gtf_lst_file=open(gtf_lst_file_name)
            gtf_lst = test_input_files(gtf_lst_file)
            print >> sys.stderr, "Running  merge catalog mode on %s \n" % (gtf_lst_file_name)
            run_merge_catalogs(gtf_lst,params.out_prefix )
            return 0
        global get_intergenic_mode
        if get_intergenic_mode==1:
            transcript_gtf=args[0]
            print >> sys.stderr, "Running get intergenic mode on %s \n" % (transcript_gtf)
            out_name="%sIntergenic_%s" % (params.out_prefix,transcript_gtf)
            run_get_intergenic(transcript_gtf,params,out_name)
            return 0


        
        #### Part 1: Cufflinks Runs
        runCuff=1


        if runCuff == 1:
            print >> sys.stderr, "Running Cufflinks section \n"

            # Check that all the primary assemblies are accessible before starting the time consuming stuff
            gtf_input_files = test_input_files(transfrag_list_file)

            if params.sample_reads_list != None:
                sample_reads_file_list = open(params.sample_reads_list)
                sample_reads = test_input_files(sample_reads_file_list)
            else:
                sample_reads = None

            #Meta assembly option:
            global run_meta_assembly
            if run_meta_assembly == 1:
                print >> sys.stderr, "Run Meta-Assembly \n"
                meta_assemble_gtfs(params, transfrag_list_filename, params.external_linc_gtf)
            #Meta Cuffcompare option:
            else:
                cuffcompare_all_assemblies(gtf_input_files);


            current_asm_gtf = "merged_asm/transcripts.gtf"  #This should probably be merged.gtf instead of transcripts.gtf
            current_asm_ids = get_gtf_ids(current_asm_gtf)

            #Subtract mask file (before running cuffdiff)
            unmasked_transcript_ids = current_asm_ids
            # # If a mask GTF file is provided, cuffcompare against it to import reference metadata
            if params.mask_gtf != None:
                 unmasked_gtf = subtract_mask_transcripts(current_asm_gtf, params.mask_gtf)
                 unmasked_transcript_ids = get_gtf_ids(unmasked_gtf)
                 current_asm_gtf = unmasked_gtf

            # Coverage-based filtering of meta-assembly
            if sample_reads != None  :
                deeply_covered = set([])
                out_dirs = []
                procs = []

                # For each BAM file we want to use to calculate coverage on transcripts in the meta-assembly
                for reads in sample_reads:
                    out_dir_name = tmp_name("cufflinks_")
                    out_dirs.append(out_dir_name)

                    # Spin in this loop until we have a free node at our disposal
                    while len(procs) == params.system_params.lsf_nodes:
                        not_done = []
                        for p in procs:
                            ret = p.poll()
                            if ret != None:
                                if ret != 0:
                                    print >> sys.stderr, "Error: Cufflinks return error on coverage calculation"
                                    sys.exit(1)
                                else:
                                    not_done.append(p)
                        if len(not_done) < params.system_params.lsf_nodes:
                            break
                        else:
                            time.sleep(10)

                    p = cufflinks(params, out_dir_name, reads, current_asm_gtf,extra_opts=["-F", "0.05", "--max-bundle-frags", "300000"], lsf=True, curr_queue=cufflinks_queue)
                    procs.append(p)

                # Wait until all the jobs we spun off have finished.
                for p in procs:
                    ret = p.wait()
                    if ret != None and ret != 0:
                        print >> sys.stderr, "Error: Cufflinks return error on coverage calculation"
                        sys.exit(1)

                # Extract IDs of sufficiently covered transcripts from the Cufflinks expression files
                for out_dir_name in out_dirs:
                    expr = open(out_dir_name+"/isoforms.fpkm_tracking")
                    expr.readline()
                    for line in expr:
                        line = line.strip()
                        cols = line.split('\t')
                        if len(cols) < 11:
                            continue
                        cov = float(cols[8])
                        fpkm_status = cols[12].rstrip()
                        if cov > params.min_transcript_cov or fpkm_status == "HIDATA":
                            deeply_covered.add(cols[0])
                print "full-length transcripts:", len(deeply_covered)
                deeply_covered_gtf_filename = "probably_full_length.gtf"
                select_gtf(current_asm_gtf, deeply_covered, deeply_covered_gtf_filename)
                current_asm_gtf = deeply_covered_gtf_filename

            current_asm_ids = get_gtf_ids(current_asm_gtf)

            #gtf_to_bed(current_asm_gtf)

            # If a reference GTF file is provided, cuffcompare against it to import reference metadata
            # if params.ref_gtf != None:
            #    (current_asm_gtf, current_asm_ids) = compare_to_reference(current_asm_gtf, params.ref_gtf, params.fasta)

            #unmasked_transcript_ids = current_asm_ids
            # If a mask GTF file is provided, cuffcompare against it to import reference metadata
            #if params.mask_gtf != None:
            #    unmasked_gtf = subtract_mask_transcripts(current_asm_gtf, params.mask_gtf)
            #    unmasked_transcript_ids = get_gtf_ids(unmasked_gtf)

            # Extract (nonrepetive) intergenic transcripts, antisense transcripts, annotated noncoding transcripts,
            # and any "orphan ORF" transcripts that are either very small peptide coding genes or lincRNAs.  Let these
            # be "possible noncoding" transcripts
            # TODO

            long_transcript_ids = set([])
            # Exclude possible noncoding transcripts shorter than 250bp or which are unspliced.
            # If a mask GTF file is provided, cuffcompare against it to import reference metadata
            if params.min_linc_length > 0:
                long_transcript_ids = get_long_transcripts(current_asm_gtf, params.min_linc_length)
            #print unmasked_transcript_ids

            min_exon_transcript_ids = set([])

            # Exclude possible noncoding transcripts shorter than 250bp or which are unspliced.
            # If a mask GTF file is provided, cuffcompare against it to import reference metadata
            if params.min_linc_exons > 1:
                min_exon_transcript_ids = get_multiexonic_transcripts(current_asm_gtf, params.min_linc_exons)


            lincRNA_ids = current_asm_ids

            print >> sys.stderr, "Current asmIds %d" % len(lincRNA_ids)
            print >> sys.stderr, "Current long transcripts  %d" % len(long_transcript_ids )
            print >> sys.stderr, "Current unmasked %d" % len(unmasked_transcript_ids)
            print >> sys.stderr, "Current min exon  %d" % len(min_exon_transcript_ids)

            lincRNA_ids &= long_transcript_ids
            print >> sys.stderr, "Current lincIDS after masking SHORT %d" % len(lincRNA_ids)
            #lincRNA_ids &= unmasked_transcript_ids
            #print >> sys.stderr, "Current lincIDS after masking CODING %d"  % len(lincRNA_ids)
            lincRNA_ids &= min_exon_transcript_ids
            print >> sys.stderr, "Current lincIDS after masking single EXON %d" % len(lincRNA_ids)

           
            possible_lincs_filename = "possible_lincs.gtf"
            select_gtf(current_asm_gtf, lincRNA_ids, possible_lincs_filename)


        #########################
        ##This is code that is used in case the mother job had crushed and one wants
        ##to restart from the pfam,csf step..
        if (runCuff == 0):
            print >> sys.stderr,"Skipped Cufflinks section updating parameters"
            possible_lincs_filename = "possible_lincs.gtf"
            deeply_covered_gtf_filename = "probably_full_length.gtf"
            current_asm_gtf = deeply_covered_gtf_filename
            current_asm_ids = get_gtf_ids(current_asm_gtf)
            lincRNA_ids = get_gtf_ids(possible_lincs_filename)
            
            ##Subtract mask file (before running cuffdiff)
            unmasked_transcript_ids = current_asm_ids
            # # If a mask GTF file is provided, cuffcompare against it to import reference metadata
            if params.mask_gtf != None:
                print >>sys.stderr, "[%s] Masking transcripts" % right_now()
                unmasked_gtf = subtract_mask_transcripts(current_asm_gtf, params.mask_gtf)
                unmasked_transcript_ids = get_gtf_ids(unmasked_gtf)
                current_asm_gtf = unmasked_gtf
                current_asm_ids = unmasked_transcript_ids
                lincRNA_ids = get_gtf_ids(current_asm_gtf)

        """  Old pipeline filtered CSF and pfam prior to naming lincRNAs
        no_pfam_hit_ids = set([])
        possible_lincs_filename2=possible_lincs_filename
        #We should run this on all possible lincs rather than the current_asm_gtf
        if params.fasta != None :
            #no_pfam_hit_ids = get_pfam_hit_transcripts(current_asm_gtf, params.fasta, params.system_params.lsf_nodes)
            possible_lincs_filename_exPfam = "possible_lincs_extendPFAM.gtf"
            (no_pfam_hit_ids , positive_pfam_hit_ids) = get_pfam_hit_transcripts(possible_lincs_filename, params.fasta, params.system_params.lsf_nodes,possible_lincs_filename_exPfam)
            lincRNA_ids &= no_pfam_hit_ids
            possible_lincs_filename2 = "possible_lincs_noPfam.gtf"
            select_gtf(possible_lincs_filename_exPfam, lincRNA_ids, possible_lincs_filename2)

        # Calculate CSF score for possible long noncoding transcripts
        #negative_csf_ids = get_negative_csf_transfrags(possible_lincs_filename, params.system_params.lsf_nodes)
        possible_lincs_filename_exPfamCsf = "possible_lincs_noPfam_extendPFAMCsf.gtf"
        (negative_csf_ids, positive_csf_ids) = get_negative_csf_transfrags(possible_lincs_filename_exPfam, params.system_params.lsf_nodes,possible_lincs_filename_exPfamCsf)
        lincRNA_ids &= negative_csf_ids
        """
        
        print >> sys.stderr, "Found %d possible lincRNAs" % len(lincRNA_ids)
                 
        # Extract the final catalog of lincRNAs
        possible_lincs_final_filename="possible_lincs_final.gtf"
        select_gtf(current_asm_gtf, lincRNA_ids, possible_lincs_final_filename)
        

        merge_gtf_list = []
        
        if params.external_linc_gtf != None and run_meta_assembly == 0:
            merge_gtf_list = [possible_lincs_final_filename, params.external_linc_gtf]
        else:
            merge_gtf_list = [possible_lincs_final_filename, possible_lincs_final_filename]

        if params.ref_gtf != None:
            merged_lincs = merge_gtfs(merge_gtf_list, "merged", params.ref_gtf)
        else:
            merged_lincs = merge_gtfs(merge_gtf_list, "merged")

        # Collect some statistics about how the discovered lincRNA set compares to the reference GTF.
        # compare_to_reference

        # Run Loyal's linc naming script to annotate the lincRNA catalog with sensible names
        lincs = name_lincs(merged_lincs, genome=params.genome, verbose=False)
        printLincs(open("linc_catalog.gtf", "w"), lincs)

        gtf_to_genbank("linc_catalog.gtf", "linc_catalog.gbk",params.genome)


        # *** #
        #Generate the level 1 (no CSF no PFAM) and level 2 sets 
        current_asm_gtf = "linc_catalog.gtf"
        current_asm_ids = get_gtf_ids(current_asm_gtf)
        lincRNA_ids = current_asm_ids
        
        
        no_pfam_hit_ids = set([])
        possible_lincs_filename = current_asm_gtf
        
        #We should run this on all possible lincs rather than the current_asm_gtf
        neg_pfam_hit_ids = lincRNA_ids
        if params.fasta != None :
            #no_pfam_hit_ids = get_pfam_hit_transcripts(current_asm_gtf, params.fasta, params.system_params.lsf_nodes)
            possible_lincs_filename_exPfam = "linc_catalog_extendPFAM.gtf"
            (neg_pfam_hit_ids , pos_pfam_hit_ids) = get_pfam_hit_transcripts(possible_lincs_filename, params.fasta, params.system_params.lsf_nodes,possible_lincs_filename_exPfam)
        
        # Calculate CSF score for possible long noncoding transcripts
        possible_lincs_filename_exPfamCsf = "linc_catalog_extendPFAMCsf.gtf"
        (neg_csf_ids, pos_csf_ids) = get_negative_csf_transfrags(possible_lincs_filename_exPfam, params.system_params.lsf_nodes,possible_lincs_filename_exPfamCsf, params.genome)
        
        curr_asm_gtf = possible_lincs_filename_exPfamCsf
        
        #Cross with pseudogenes
        neg_pseudo_ids = set([])
        possible_lincs_filename_exPfamCsfPseudo = "linc_catalog_extendPFAMCsfPseudo.gtf"
        if params.pseudogene != None:
            (pos_pseudo_ids,neg_pseudo_ids) = annot_overlap_transcripts (possible_lincs_filename_exPfamCsf,params.pseudogene, "pseudogene" ,possible_lincs_filename_exPfamCsfPseudo)
            curr_asm_gtf = possible_lincs_filename_exPfamCsfPseudo
        else:
            neg_pseudo_ids = lincRNA_ids
        
        #Cross with merged annotations
        external_transcript_ids = set([])
        if params.external_linc_gtf != None:
            non_external_gtf = subtract_mask_transcripts(current_asm_gtf, params.external_linc_gtf )
            non_external_transcript_ids = get_gtf_ids(non_external_gtf)
            external_transcript_ids = lincRNA_ids - non_external_transcript_ids
        
        # Level 1 : csf_neg, pfam_neg, pseudo_neg and merge_file_neg
        lincRNA_ids1 = set([])
        lincRNA_ids1 = lincRNA_ids1 | lincRNA_ids
        lincRNA_ids1 &= neg_pfam_hit_ids
        lincRNA_ids1 &= neg_csf_ids
        lincRNA_ids1 &= neg_pseudo_ids 
        lincRNA_ids1 = lincRNA_ids1 | external_transcript_ids
        linc_catalog_level1 = "linc_catalog.level1.gtf"
        select_gtf(curr_asm_gtf, lincRNA_ids1, linc_catalog_level1)
        
        # Level 2 : others
        linc_catalog_level2 = "linc_catalog.level2.gtf"
        lincRNA_ids2 = lincRNA_ids - lincRNA_ids1
        #print >> sys.stderr, "length of lincRNA ids : %d" % len(lincRNA_ids)
        #print >> sys.stderr, "length of lincRNA ids1 : %d" % len(lincRNA_ids1)
        #print >> sys.stderr, lincRNA_ids
        #print >> sys.stderr, lincRNA_ids1
        #print >> sys.stderr, lincRNA_ids2
        
        select_gtf(curr_asm_gtf, lincRNA_ids2, linc_catalog_level2) 
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Trancriptome catalog build complete [%s elapsed]" %  formatTD(duration)

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        #print >> sys.stderr, "    for detailed help see http://spats.cbcb.umd.edu/manual.html"
        return 2


def main2(argv=None):
    warnings.filterwarnings("ignore", "tmpnam is a potential security risk") #
    
    # Initialize default parameter values
    global params #
    params = TestParams()#
    
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()
        
        if params.system_params.system == "Broad": 
            global lsf_queue
            global cufflinks_queue
            global pfam_db
            global aa_converter
            global gffread
            global pfammer
            global csf_bin
            
            csf_bin = "~mlin/bin/batchPhyloCSF"
            #lsf_queue = "regevlab"
            #cufflinks_queue="regevlab"
            pfam_db = "/ahg/regev/users/nmcabili/Programs/Hmmer/pfam"
            aa_converter =  "perl /ahg/regev/users/nmcabili/Src/scripts/pfamAnalysis/seqToAminoAcid.pl"
            gffread = "gffread"
            pfammer = "perl /ahg/regev/users/nmcabili/Programs/Hmmer/PfamScan/pfam_scan.pl"
        
        
        global run_log
        global run_cmd
        
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning transcriptome build (suite v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------"
        print >> sys.stderr
        print >> sys.stderr, "  Configuration:"
        print >> sys.stderr, "     system: %s" % params.system_params.system
        print >> sys.stderr, "     genome: %s" % params.genome
        print >> sys.stderr, "     output_dir: %s" % output_dir
        print >> sys.stderr, "     logging_dir: %s" % logging_dir
        print >> sys.stderr, "     run log: %s" % run_log
        #print >> sys.stderr, "     Meta assembly: %s" % ("yes" if run_meta_assembly else "no")
        #print >> sys.stderr, "     Merged catalog: %s" % ("yes" if merge_catalog_mode else "no")
        print >> sys.stderr, "     CSF threshold: %f" % csf_score_threshold
        print >> sys.stderr, "     LSF nodes: %d" % params.system_params.lsf_nodes
        if lsf_mem != None:
            print >> sys.stderr, "     LSF memory per node: %d" % lsf_mem
        else:
            print >> sys.stderr, "     LSF memory per node: Not specified" 
        print >> sys.stderr
        
        start_time = datetime.now()
        prepare_output_dir()
        #now=datetime.datetime.now()
        #curr_time=now.strftime("%Y%m%d_%H%M")
        
        run_log = open(logging_dir + "run.log", "w", 0)
        run_cmd = " ".join(argv)

        print >> run_log, run_cmd

        if len(args) < 1:
            raise Usage(use_message)

        transfrag_list_filename = args[0]
        transfrag_list_file = open(args[0], "r")

        if params.ref_gtf != None:
            test_input_files([params.ref_gtf])
        else:
            print >> sys.stderr, "Warning: no reference GTF provided!"
        if params.fasta != None:
            if os.path.isdir(params.fasta) == False:
                test_input_files([params.fasta])

        if params.mask_gtf != None:
            test_input_files([params.mask_gtf])

        if params.external_linc_gtf != None:
            test_input_files([params.external_linc_gtf])

        ###############
        #Non-assembly mode:
        # There is an external GTF file of annotations that you want to process with PFam, CSF, and the rest of the lincer characterization pipeline. It excludes abundance estimations.
        global non_assembly_mode
        if non_assembly_mode==1:
            known_transcript_gtf=args[0]
            print >> sys.stderr, "Running non assembly mode on %s \n" % (known_transcript_gtf)
            run_non_assembly_mode(known_transcript_gtf,params)
            return 0
        
        #Merge Catalog Mode
        # You have several lincRNA catalogs, and you want to make them a unique set of lincRNAs. Goes through ID assignment, creates Genbank.  Doesn't do PFam and CSF?
        global merge_catalog_mode
        if merge_catalog_mode==1:
            gtf_lst_file_name=args[0]
            gtf_lst_file=open(gtf_lst_file_name)
            gtf_lst = test_input_files(gtf_lst_file)
            print >> sys.stderr, "Running  merge catalog mode on %s \n" % (gtf_lst_file_name)
            run_merge_catalogs(gtf_lst,params.out_prefix )
            return 0
        
        #Get Intergenic Mode
        # A run that extracts intergenic lincRNAs from a complete lincRNA catalog output (must be used in the directory with accesory files from lincer pipeline.
        global get_intergenic_mode
        if get_intergenic_mode==1:
            transcript_gtf=args[0]
            print >> sys.stderr, "Running get intergenic mode on %s \n" % (transcript_gtf)
            out_name="%sIntergenic_%s" % (params.out_prefix,transcript_gtf)
            run_get_intergenic(transcript_gtf,params,out_name)
            return 0

        #CSF only mode:
        # There is an external GTF file of annotations that you want to process with CSF.
        global csf_only_mode
        if csf_only_mode==1:
            known_transcript_gtf=args[0]
            print >> sys.stderr, "Running CSF only mode on %s\n" % known_transcript_gtf
            known_transcript_bed = gtf_to_bed(known_transcript_gtf)
            csf_scorefilename = calculate_csf(known_transcript_bed, params.system_params.lsf_nodes, params.genome)
            return 0

        
        #### Part 1: Cufflinks Runs
        runCuff=1

        if runCuff == 1:
            print >> sys.stderr, "Running Cufflinks section \n"

            # Check that all the primary assemblies are accessible before starting the time consuming stuff
            gtf_input_files = test_input_files(transfrag_list_file)

            if params.sample_reads_list != None:
                sample_reads_file_list = open(params.sample_reads_list)
                sample_reads = test_input_files(sample_reads_file_list)
            else:
                sample_reads = None

            #Meta assembly option:
            global run_meta_assembly
            if run_meta_assembly == 1:
                print >> sys.stderr, "Run Meta-Assembly \n"
                meta_assemble_gtfs(params, transfrag_list_filename, params.external_linc_gtf) #This returns a pointer string to 'merged_asm/merged.gtf'
                current_asm_gtf = "merged_asm/merged.gtf"
            #Meta Cuffcompare option:
            else:
                cuffcompare_all_assemblies(gtf_input_files); #This doesn not return a value but writes results to 'merged_asm/transcripts.gtf'
                current_asm_gtf = "merged_asm/transcripts.gtf"

            current_asm_ids = get_gtf_ids(current_asm_gtf)
            
            #TODO: Write an ouputfile of the correct current_asm_gtf with standardized name to preserve current_asm (post_merge_transcripts.gtf as e.g.)
            
            #Now we have a correct pointer for current_asm_gtf and a CORRECT list of ids from current_asm_gtf (LAG + NMC)

            #Subtract mask file (before running cuffdiff)
            unmasked_transcript_ids = current_asm_ids
            # # If a mask GTF file is provided, cuffcompare against it to import reference metadata
            if params.mask_gtf != None:
                 unmasked_gtf = subtract_mask_transcripts(current_asm_gtf, params.mask_gtf) #Returns string pointer for unmasked gtf file. (ie. transcripts overlapping masked transcripts are removed)
                 unmasked_transcript_ids = get_gtf_ids(unmasked_gtf) #Gets ids from unmasked_gtf
                 current_asm_gtf = unmasked_gtf #This moves unmasked gtf file to main current_asm_gtf

            # Coverage-based filtering of meta-assembly
            if sample_reads != None  :
                deeply_covered = set([])
                out_dirs = []
                procs = []

                # For each BAM file we want to use to calculate coverage on transcripts in the meta-assembly
                for reads in sample_reads:
                    out_dir_name = tmp_name("cufflinks_")
                    out_dirs.append(out_dir_name)

                    # Spin in this loop until we have a free node at our disposal
                    while len(procs) == params.system_params.lsf_nodes:
                        not_done = []
                        for p in procs:
                            ret = p.poll()
                            if ret != None:
                                if ret != 0:
                                    print >> sys.stderr, "Error: Cufflinks return error on coverage calculation"
                                    sys.exit(1)
                                else:
                                    not_done.append(p)
                        if len(not_done) < params.system_params.lsf_nodes:
                            break
                        else:
                            time.sleep(10)

                    p = cufflinks(params, out_dir_name, reads, current_asm_gtf,extra_opts=["-F", "0.05", "--max-bundle-frags", "300000"], lsf=True, curr_queue=cufflinks_queue)
                    procs.append(p)

                # Wait until all the jobs we spun off have finished.
                for p in procs:
                    ret = p.wait()
                    if ret != None and ret != 0:
                        print >> sys.stderr, "Error: Cufflinks return error on coverage calculation"
                        sys.exit(1)

                # Extract IDs of sufficiently covered transcripts from the Cufflinks expression files
                for out_dir_name in out_dirs:
                    expr = open(out_dir_name+"/isoforms.fpkm_tracking")
                    expr.readline()
                    for line in expr:
                        line = line.strip()
                        cols = line.split('\t')
                        if len(cols) < 11:
                            continue
                        cov = float(cols[8])
                        fpkm_status = cols[12].rstrip()
                        if cov > params.min_transcript_cov or fpkm_status == "HIDATA":
                            deeply_covered.add(cols[0])
                print "full-length transcripts:", len(deeply_covered)
                deeply_covered_gtf_filename = "probably_full_length.gtf"
                select_gtf(current_asm_gtf, deeply_covered, deeply_covered_gtf_filename) #Writing ouput of deeply covered  and masked transcripts to 'probably_full_length.gtf'
                current_asm_gtf = deeply_covered_gtf_filename #Replace current_asm_gtf with deeply_covered_gtf ('probably_full_length.gtf')

            current_asm_ids = get_gtf_ids(current_asm_gtf) #Update current_asm_ids so that they match current_asm_gtf
            
            ###########
            #Filter for lincRNA size
            ##########
            
            long_transcript_ids = set([])
            # Exclude possible noncoding transcripts shorter than 250bp or which are unspliced.
            if params.min_linc_length > 0:
                long_transcript_ids = get_long_transcripts(current_asm_gtf, params.min_linc_length)
            
            ###########
            #Filter for minimum number of exons
            ###########
            
            min_exon_transcript_ids = set([])

            # Exclude possible noncoding transcripts shorter than 250bp or which are unspliced.
            if params.min_linc_exons > 1:
                min_exon_transcript_ids = get_multiexonic_transcripts(current_asm_gtf, params.min_linc_exons)

            #####
            #Begin Set intersections for size and min_exon filtering
            #####
            
            lincRNA_ids = current_asm_ids

            print >> sys.stderr, "Current asmIds %d" % len(lincRNA_ids)
            print >> sys.stderr, "Current long transcripts  %d" % len(long_transcript_ids )
            print >> sys.stderr, "Current unmasked %d" % len(unmasked_transcript_ids)
            print >> sys.stderr, "Current min exon  %d" % len(min_exon_transcript_ids)

            lincRNA_ids &= long_transcript_ids
            print >> sys.stderr, "Current lincIDS after masking SHORT %d" % len(lincRNA_ids)
          
            lincRNA_ids &= min_exon_transcript_ids
            print >> sys.stderr, "Current lincIDS after masking single EXON %d" % len(lincRNA_ids)

            #Possible lincs are those that have minimum coverage, masked via mask_gtf, are >=min_length and >=min_exons. These are in 'possible_lincs.gtf' 
            possible_lincs_filename = "possible_lincs.gtf"
            select_gtf(current_asm_gtf, lincRNA_ids, possible_lincs_filename)


        #########################
        ##This is code that is used in case the mother job had crashed and one wants
        ##to restart from the pfam,csf step..
        if (runCuff == 0):
            print >> sys.stderr,"Skipped Cufflinks section updating parameters"
            possible_lincs_filename = "possible_lincs.gtf" # This is deeply covered, masked, filtered for length and n exons
            #deeply_covered_gtf_filename = "probably_full_length.gtf" # This is just deeply covered and masked
            #current_asm_gtf = deeply_covered_gtf_filename
            #current_asm_ids = get_gtf_ids(current_asm_gtf)
            lincRNA_ids = get_gtf_ids(possible_lincs_filename)
            
        print >> sys.stderr, "Found %d possible lincRNAs" % len(lincRNA_ids)

        ##########
        #Bring in External linc file (gtf)
        ##########
        #Must use cuffmerge with -g for external_linc_gtf
        #Output will be a gtf that guarantees inclusion of external_linc_gtfs
        
        #Write file_manifest for cuffmerge (should only include possible_lincs.gtf)
        if params.external_linc_gtf != None:
            merge_manifest_filename = "tmp/merge_manifest.txt"
            merge_manifest = open(merge_manifest_filename,'w')
            print >> merge_manifest, possible_lincs_filename
            merge_manifest.close()
        
            external_linc_merge_dir = "tmp/external_linc_merge"
            os.mkdir(external_linc_merge_dir)
        
            merged_lincs = cuffmerge(params, merge_manifest_filename, params.external_linc_gtf, outDir=external_linc_merge_dir)
        else:
            merged_lincs = possible_lincs_filename

        ##########
        #name lincRNAs
        ##########
        # Run Loyal's linc naming script to annotate the lincRNA catalog with sensible names
        lincs = name_lincs(merged_lincs, genome=params.genome, verbose=False)
        printLincs(open("linc_catalog.gtf", "w"), lincs)
        
        ##########
        #Create and write genbank file
        ##########
        #gtf_to_genbank("linc_catalog.gtf", "linc_catalog.gbk", params.genome)
        
        ##########
        #Collect initial info
        ##########
        current_asm_gtf = "linc_catalog.gtf"
        current_asm_ids = get_gtf_ids(current_asm_gtf)
        lincRNA_ids = current_asm_ids
        
        ##########
        #PFAM analysis
        ##########
        neg_pfam_hit_ids = lincRNA_ids
        if not params.no_pfam and params.fasta != None:
            #no_pfam_hit_ids = get_pfam_hit_transcripts(current_asm_gtf, params.fasta, params.system_params.lsf_nodes)
            possible_lincs_filename_exPfam = "linc_catalog_extendPFAM.gtf"
            (neg_pfam_hit_ids , pos_pfam_hit_ids) = get_pfam_hit_transcripts(current_asm_gtf, params.fasta, params.system_params.lsf_nodes, possible_lincs_filename_exPfam)
            current_asm_gtf = possible_lincs_filename_exPfam

        ###########
        #CSF analysis
        ###########
        # Calculate CSF score for possible long noncoding transcripts
        neg_csf_ids = lincRNA_ids
        if not params.no_csf:
            possible_lincs_filename_exPfamCsf = "linc_catalog_extendPFAMCsf.gtf"
            (neg_csf_ids, pos_csf_ids) = get_negative_csf_transfrags(possible_lincs_filename_exPfam, params.system_params.lsf_nodes, possible_lincs_filename_exPfamCsf, params.genome)
            current_asm_gtf = possible_lincs_filename_exPfamCsf
        
        ###########
        #Cross with external pseudogene list
        #############
        neg_pseudo_ids = lincRNA_ids
        if params.pseudogene != None:
            possible_lincs_filename_exPfamCsfPseudo = "linc_catalog_extendPFAMCsfPseudo.gtf"
            (pos_pseudo_ids,neg_pseudo_ids) = annot_overlap_transcripts(current_asm_gtf, params.pseudogene, "pseudogene", possible_lincs_filename_exPfamCsfPseudo)
            curr_asm_gtf = possible_lincs_filename_exPfamCsfPseudo
        
        ############
        #Cross with merged annotations
        ############
        #TODO: Add a parameter for external_true_lincs to be included by force at this point (even if they have PFam hits, high CSF, etc).
        
        external_true_transcript_ids = set([])
        if params.external_true_lincs_gtf != None:
            non_external_gtf = subtract_mask_transcripts(current_asm_gtf, params.external_true_lincs_gtf)
            non_external_transcript_ids = get_gtf_ids(non_external_gtf)
            external_true_transcript_ids = lincRNA_ids - non_external_transcript_ids
        
        ###########
        #Level 1 / 2 classification (set operations)
        ###########
        lincRNA_ids1 = set([])
        lincRNA_ids1 = lincRNA_ids1 | lincRNA_ids
        lincRNA_ids1 &= neg_pfam_hit_ids
        lincRNA_ids1 &= neg_csf_ids
        lincRNA_ids1 &= neg_pseudo_ids 
        lincRNA_ids1 = lincRNA_ids1 | external_true_transcript_ids
        
        linc_catalog_level1 = "linc_catalog.level1.gtf"
        select_gtf(curr_asm_gtf, lincRNA_ids1, linc_catalog_level1)
        
        linc_catalog_level2 = "linc_catalog.level2.gtf"
        lincRNA_ids2 = lincRNA_ids - lincRNA_ids1
        select_gtf(curr_asm_gtf, lincRNA_ids2, linc_catalog_level2) 
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Trancriptome catalog build complete [%s elapsed]" %  formatTD(duration)

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        #print >> sys.stderr, "    for detailed help see http://spats.cbcb.umd.edu/manual.html"
        return 2


if __name__ == "__main__":
    sys.exit(main2())
