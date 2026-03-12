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
	
	#remove unneeded bam file. 
	remove = "rm "+sample+".Aligned.sortedByCoord.out.bam*"
	runstuff(remove)
	
	#assemble bam sample rename command
	rename = "\'samtools view -H "+sample+".Aligned.sortedByCoord.out.md.bam  | sed \"s/SM:[^\t]*/SM:"+sample+"/g\" | samtools reheader - "+sample+".Aligned.sortedByCoord.out.md.bam > "+sample+".Aligned.sortedByCoord.out.id.bam\'"
	
	#submit job
	runstuff("bsub -g /lorenzk/"+Name+" -o "+sample+"bam_lsf_log.out "+ rename)
	
	
	sleep(10) #pause before submitting more.
