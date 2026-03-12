import argparse
import os
import numpy as np
import subprocess
import gzip
import re

parser = argparse.ArgumentParser()
parser.add_argument("-e", help="egenes file",type=str)
parser.add_argument("-i", help="LD base filename",type=str)
parser.add_argument("-o", help="outfile name",type=str)
args = parser.parse_args()

egenes_file = args.e
LDbase = args.i
Out = args.o

###################################################
##read in the egenes file
###################################################
print("reading in egenes file")

egenes = {}
#structure here is going to be
#egenes[ENSG]["rsID"]=set(SNPs)
#egenes[ENSG]["info"]=outlist

genebareENSG = {}
#dict of ENSG ids without a .## pointing to ones with a .##

snpdict = {}
#lead SNP pointing to gene name

with open(egenes_file,'rt') as f:
	header = f.readline().rstrip().split("\t")

	for line in f:
	#remove white space at end of line, split by tab
		stuff = line.rstrip().split("\t")
		
		#initialize gene in dictionary
		ENSGID = stuff[0]
		egenes[ENSGID]={}
		egenes[ENSGID]["rsid"]=set()
		egenes[ENSGID]["info"]= stuff
		#default answer for is SNP sig in other data is no
		egenes[ENSGID]["DICE"]="no"
		egenes[ENSGID]["eQTLgen"]="no"
		
		genebare, toss = ENSGID.split('.')
		genebareENSG[genebare]=ENSGID
		
		rsid = stuff[17]
		if rsid in snpdict:
			snpdict[rsid].add(ENSGID)
		else:
			snpdict[rsid] = set()
			snpdict[rsid].add(ENSGID)



###################################################
##step through each LD file and save proxies
###################################################
print("getting LD info")

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
			#if lead is in snpdict, add proxy to correct gene SNPlists
			if snp1rs  in snpdict:
				genelist = snpdict[snp1rs]
				for ENSGID in genelist:
					egenes[ENSGID]["rsid"].add(snp2rs)

###################################################
##step through each DICE file and check for sig proxies
###################################################
print("checking DICE for sig proxies")

DICE = ["/project/voight_datasets_01/eQTL/DICE/B_CELL_NAIVE_filtered.txt", "/project/voight_datasets_01/eQTL/DICE/TFH_filtered.txt", "/project/voight_datasets_01/eQTL/DICE/CD4_NAIVE_filtered.txt"]

for i in DICE:
	with open(i,'r') as t:
		t.readline() #remove header
		for line in t:
			stuff= line.rstrip().split()
			
			rsID = stuff[4]
			genebare = stuff[0]
			
			#only care if gene is an egene in our tissue
			if genebare in genebareENSG:
				ENSGID = genebareENSG[genebare]
				
				#then if rsID is in proxy list, change flag to yes
				if rsID in egenes[ENSGID]["rsid"]:
					egenes[ENSGID]["DICE"]="yes"

###################################################
##step through eQTLgen file and check for sig proxies
###################################################
print("checking eQTLgen for sig proxies")

eQTLgen = ["/project/voight_datasets/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"]

for i in eQTLgen:
	with open(i,'r') as t:
		t.readline() #remove header
		for line in t:
			stuff= line.rstrip().split()
			
			rsID = stuff[1]
			genebare = stuff[7]
			
			#only care if gene is an egene in our tissue
			if genebare in genebareENSG:
				ENSGID = genebareENSG[genebare]
				
				#then if rsID is in proxy list, change flag to yes
				if rsID in egenes[ENSGID]["rsid"]:
					egenes[ENSGID]["eQTLgen"]="yes"
			


###################################################
##step through present list and output to file
###################################################
print("writing output file")

outfile = open(Out, "w")

outfile.write("\t".join(str(x) for x in header)+"\tlead_sig_DICE\tlead_sig_eQTLgen\n")

for ENSGID in egenes:
	stuff = egenes[ENSGID]["info"]
	DICE_out = egenes[ENSGID]["DICE"]
	eQTLgen_out = egenes[ENSGID]["eQTLgen"]

	outfile.write("\t".join(str(x) for x in stuff)+"\t"+DICE_out+"\t"+eQTLgen_out+"\n")


