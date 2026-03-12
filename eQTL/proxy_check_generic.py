import argparse
import os
import numpy as np
import subprocess
import gzip
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="LD base filename",type=str)
parser.add_argument("-l", help="list of rsIDs filename",type=str)
parser.add_argument("-o", help="outfile name",type=str)
args = parser.parse_args()

LDbase = args.i
snps = args.l
Out = args.o

###################################################
##step through each LD file and save proxies
###################################################

proxies = {}

loop = list(range(1,23))
loop.append("X")

for i in loop:
	LDfile = LDbase+str(i)+".ld"
	with open(LDfile,'r') as t:
		t.readline() #remove header
		for line in t:
			stuff= line.rstrip().split()
			
			#set up info
			snp1pos = "chr"+str(stuff[0])+":"+str(stuff[1])
			snp1rs = stuff[2]
			snp2pos = "chr"+str(stuff[3])+":"+str(stuff[4])
			snp2rs = stuff[5]
			LD = float(stuff[6])
			
			#if rsid is compound (ie contains two rsids), skip it. 
			if ";" in snp1rs or ";" in snp2rs:
				continue
			if "ss" in snp1rs or "ss" in snp2rs:
				continue
			if "esv" in snp1rs or "esv" in snp2rs:
				continue
			#test = snp1rs+" "+snp2rs+"\n"
			#print(test)
			#if proxy isn't in list, add it & create set, then add lead
			if snp2rs not in proxies:
				proxies[snp2rs] = set()
				proxies[snp2rs].add(snp1rs)
				
			else:
				proxies[snp2rs].add(snp1rs)

###################################################
##step through each DICE file and check for proxies
###################################################

present = {}


with open(snps,'rt') as t:
	t.readline() #remove header
	for line in t:
		rsID= line.rstrip()
		
		#print("test rsid: "+rsID)
		if rsID in proxies:
			leads = proxies[rsID]
			
			#add leads to present dict
			
			for i in leads:
				if i in present:
					present[i].add(rsID)
				else:
					present[i]=set()
					present[i].add(rsID)
			
		

###################################################
##step through present list and output to file
###################################################

outfile = open(Out, "w")
outfile.write("tonsil_lead\tproxy_hits\n")

for lead in present:
	other = present[lead]
	outfile.write(lead+"\t"+", ".join(x for x in other)+"\n")


