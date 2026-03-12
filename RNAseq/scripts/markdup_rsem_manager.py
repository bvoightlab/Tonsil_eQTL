import sys
import numpy as np
import os
import pdb
import gzip
import argparse
import subprocess
import os.path
import string
from time import sleep

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="list of filepaths",type=str)
parser.add_argument("-n", help="group name",default='',type=str)
parser.add_argument("-j", help="number of jobs to run at the same time",default='5',type=int)
parser.add_argument("-c", help="# seconds to wait before checking",default='300',type=int)
args = parser.parse_args()

Infile = args.i
Name = args.n
Njobs = 1+args.j
Wait = args.c


###################################################
##definitions
###################################################

def runstuff(cmd):
	subprocess.call(cmd, shell=True)
	print("submitted "+cmd)
	
def get_out(cmd):
	stuff = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	print("returned "+cmd)
	stdout = []
	while True:
		l = stuff.stdout.readline()
		line = l.rstrip()
		if line == b'' and stuff.poll() != None:
			break
		else:
			stdout.append(line)
	return stdout

#checks how many jobs are running and waits until less than Njobs are running. 
def wait5(cmd):
	while True:
		stuff = get_out("bjobs -g /lorenzk/"+cmd+" | wc -l")
		count = stuff[0]
		if "No unfinished job found" in str(count):
			print("returned: No unfinished job found")
			break
		elif (int(count) < Njobs):
			print("count less than njobs")
			break
		sleep(Wait)

###################################################
##read in inputs & make tuples of all outpaths and inpaths
###################################################
Romp = set()

for line in open(Infile,'r'):
	#remove white space at end of line
	sample = line.rstrip()
	#save as set to prevent duplication
	Romp.add(sample)

###################################################
##step through list and submit Njobs at a time
###################################################


for sample in Romp:
	wait5(Name)
	
	#assemble mark duplicates command
	markdup = "python -u /project/voight_ML/lorenzk/Romberg_eQTL/scripts/run_MarkDuplicates.py "+sample+".Aligned.sortedByCoord.out.bam "+sample
	
	#submit job
	runstuff("bsub -g /lorenzk/"+Name+" -o "+sample+"markdup_lsf_log.out "+ markdup)
	
	#assemble rsem command
	rsem = "python /project/voight_ML/lorenzk/Romberg_eQTL/scripts/run_RSEM.py --max_frag_len 1000 --estimate_rspd true --is_stranded true --threads 2 /project/voight_ML/lorenzk/Romberg_eQTL/reference_files "+sample+".Aligned.toTranscriptome.out.bam "+sample
	
	#submit job
	runstuff("bsub -g /lorenzk/"+Name+" -o "+sample+"rsem_lsf_log.out "+ rsem)

	
	sleep(10) #pause before submitting more.
