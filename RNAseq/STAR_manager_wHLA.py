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
Romp = []

for line in open(Infile,'r'):
	#remove white space at end of line
	line2 = line.rstrip()
	#split input lines by tab
	outpath,inpath = line2.split("\t")
	#save as dictionaries to prevent duplication
	Romp.append((outpath,inpath))

###################################################
##step through list and submit Njobs at a time
###################################################

#assemble STAR command

#control script and ref file location
pt1 = "python /project/voight_ML/lorenzk/Romberg_eQTL/scripts/run_STAR.py /project/voight_ML/lorenzk/Romberg_eQTL/reference_files/hg38_human_gencodev34_HLA "

#pt2 is variable by input/output filepaths
#TC257NaiveT_S146_L001_R1_001.fastq.gz,TC257NaiveT_S146_L002_R1_001.fastq.gz,TC257NaiveT_S146_L003_R1_001.fastq.gz,TC257NaiveT_S146_L004_R1_001.fastq.gz \
#TC257NaiveT_S146_L001_R2_001.fastq.gz,TC257NaiveT_S146_L002_R2_001.fastq.gz,TC257NaiveT_S146_L003_R2_001.fastq.gz,TC257NaiveT_S146_L004_R2_001.fastq.gz \
#TC257_NaiveT \

#a ton of settings
pt3 =  " --output_dir star_out --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 --outSAMstrandField intronMotif --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMattrRGline ID:rg1 SM:sm1 --outSAMattributes NH HI AS nM NM ch --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType Junctions WithinBAM SoftClip --chimMainSegmentMultNmax 1 --threads 8"


for (outpath, inpath) in Romp:
	wait5(Name)
	
	#assemble pt2
	pt2 = inpath+"_L001_R1_001.fastq.gz,"+inpath+"_L002_R1_001.fastq.gz,"+inpath+"_L003_R1_001.fastq.gz,"+inpath+"_L004_R1_001.fastq.gz "+inpath+"_L001_R2_001.fastq.gz,"+inpath+"_L002_R2_001.fastq.gz,"+inpath+"_L003_R2_001.fastq.gz,"+inpath+"_L004_R2_001.fastq.gz "+outpath
	
	#submit job
	runstuff("bsub -m vengeance -g /lorenzk/"+Name+" -o "+outpath+"_lsf_log.out "+ pt1 + pt2 + pt3)
	sleep(10) #pause before submitting more.
