import sys
import numpy as np
import os
import pdb
import gzip
import argparse
import subprocess
import os.path
import string
import random
from collections import defaultdict
import re

parser = argparse.ArgumentParser()
parser.add_argument("-x", help="config file",type=str)
parser.add_argument("-o", help="out file",type=str)

args = parser.parse_args()

trait_config = args.x
Out = args.o
tab = "\t"


###################################################
##script goal
###################################################

#script designed to take in
# config file which defines a slew of input files to consolidate by egene
# see /project/voight_viz/lorenzk/Romberg_eQTL/GCB_egenes_overlaps_config.txt as example

#and outputs 
# based on the egenes file, lots of info about the egene. 
#call with
#python /project/voight_ML/lorenzk/scripts/eQTL_miniscoreboard_egenes.py \
#-x /project/voight_viz/lorenzk/Romberg_eQTL/GCB_egenes_overlaps_config.txt \
#-o /project/voight_viz/lorenzk/Romberg_eQTL/GCB_egenes_scoreboard.txt  

###################################################
##read in trait config file, save all the files
###################################################


with open(trait_config,'r') as f:
	# dice tested leads
	DICE = f.readline().rstrip().split(": ")
	
	# dice egenes
	DICE_egenes = f.readline().rstrip().split(": ")
	
	# eQTLgen tested leads
	eQTLgen = f.readline().rstrip().split(": ")
	
	# eQTLgen egenes
	eQTLgen_egenes = f.readline().rstrip().split(": ")

	# cAID leads
	cAID_snps = f.readline().rstrip().split(": ")
	
	# # basepairs in promoter
	prom_length = f.readline().rstrip().split(": ")
	
	# egenes file
	egenes_file = f.readline().rstrip().split(": ")
	
	# significant pairs file
	sigpairs = f.readline().rstrip().split(": ")
	
	# IUIS genes
	IUIS = f.readline().rstrip().split(": ")
	
	# top haploinsufficient genes
	haplo_top = f.readline().rstrip().split(": ")
	
	# all haploinsufficient genes
	haplo_all = f.readline().rstrip().split(": ")
	
	# ATAC regions
	ATAC = f.readline().rstrip().split(": ")
	
	# Capture C 1frag
	onefrag = f.readline().rstrip().split(": ")
	
	# Capture C 4frag
	fourfrag = f.readline().rstrip().split(": ")
	
	# colocalization trait list, space separated
	traitlist = f.readline().rstrip().split(": ")
	traits = traitlist[1].split(" ")
	
	# coloc summary file base
	colocsummarybase = f.readline().rstrip().split(": ")
	
	# tonsil focus trait
	tonsil_focus = f.readline().rstrip().split(": ")
	
	# DICE focus trait
	DICE_focus = f.readline().rstrip().split(": ")
	

###################################################
##read in the egenes file
###################################################
print("reading in egenes file")

egenes = {}
#structure here is going to be
#egenes[ENSG]["rsIDs"]=set(SNPs)
#egenes[ENSG]["info"]=outlist

genenames = {}
#dict of gene names pointing to ENSG ids

genebareENSG = {}
#dict of ENSG ids without a .## pointing to ones with a .##

with open(egenes_file[1],'rt') as f:
	header = f.readline().rstrip().split("\t")
	
	header_out = header[0:7]
	header_out.pop(5)
	
	for line in f:
	#remove white space at end of line, split by tab
		stuff = line.rstrip().split("\t")
		
		keep = stuff[0:7]
		keep.pop(5)
		
		eQTLgen_sig = stuff[-1]
		DICE_sig = stuff[-2]
		
		#initialize gene in dictionary
		#stuff[0] is ENSG ID
		egenes[stuff[0]]={}
		egenes[stuff[0]]["rsid"]=set()
		egenes[stuff[0]]["info"]= keep
		#default answer for is SNP in promoter is no
		egenes[stuff[0]]["check"]="no"
		#save whether has sig proxy in DICE or eQTLgen
		egenes[stuff[0]]["DICE_sig"]= DICE_sig
		egenes[stuff[0]]["eQTLgen_sig"]= eQTLgen_sig
		#add to genenames dict
		#stuff[1] is gene name
		genenames[stuff[1]] = stuff[0]
		
		genebare, toss = stuff[0].split('.')
		genebareENSG[genebare]=stuff[0]

#then step through sig pairs file and add rsIDs & check promotor
print("adding sig SNPs")

with open(sigpairs[1],'rt') as f:
	header = f.readline().rstrip().split("\t")
	for line in f:
		stuff = line.rstrip().split("\t")
		
		rsid = stuff[12]
		otter = stuff[1].split(":")
		pos = int(otter[1])
		
		#add rsid
		egenes[stuff[0]]["rsid"].add(rsid)
		
		#is rsid in promoter?
		#define promoter
		info = egenes[stuff[0]]["info"]
		#if start < end, promoter is start-1000 to start
		if int(info[3]) < int(info[4]):
			promstart = int(info[3])-1000
			promend = int(info[3])
		else:
			promstart = int(info[4])
			promend = int(info[4])+1000
		#if snp pos is in promoter, flip check to yes
		if promstart <= pos and pos <= promend:
			egenes[stuff[0]]["check"]="yes"
			#print(str(pos)+" is between "+str(promstart)+" and "+str(promend))
		
#add # significant SNPs & y/n on promoter to out list
header_out.append("sig_var")
header_out.append("promoter")
for ENSG in egenes:
	value = len(egenes[ENSG]["rsid"])
	egenes[ENSG]["info"].append(value)
	check = egenes[ENSG]["check"]
	egenes[ENSG]["info"].append(check)
	
	#reset check value to no
	egenes[ENSG]["check"] = "no"
	
###################################################
##check if lead SNP/proxy tested in DICE or eQTLgen
###################################################
print("checking DICE and eQTLgen SNPs & egenes")

#make set of DICE and eQTLgen SNPs
DICE_snps = set()
with open(DICE[1],'rt') as f:
	for line in f:
		SNP = line.rstrip()
		
		DICE_snps.add(SNP)

DICE_setegenes = set()
#for GCB, no DICE testing done
if DICE_egenes[1] == "NA":
	print("no file for DICE egenes")
else:
	with open(DICE_egenes[1],'rt') as f:
		header = f.readline()
		for line in f:
			stuff = line.rstrip().split("\t")
			DICE_setegenes.add(stuff[0])

eQTLgen_snps = set()
with open(eQTLgen[1],'rt') as f:
	for line in f:
		SNP = line.rstrip()
		
		eQTLgen_snps.add(SNP)

eQTLgen_setegenes = set()
with open(eQTLgen_egenes[1],'rt') as f:
	header = f.readline()
	for line in f:
		stuff = line.rstrip().split("\t")
		
		eQTLgen_setegenes.add(stuff[7])

#step through egenes and add summary
header_out.append("DICE")
header_out.append("DICE_egene")
header_out.append("lead_sig_DICE")
header_out.append("eQTLgen")
header_out.append("eQTLgen_egene")
header_out.append("lead_sig_eQTLgen")
for ENSG in egenes:
	sig_SNPs = egenes[ENSG]["rsid"]
	bareENSG, toss = ENSG.split('.')
	#if no overlap between significant SNPs and DICE SNPs, not tested
	if sig_SNPs.isdisjoint(DICE_snps):
		egenes[ENSG]["info"].append("not_tested")
	else:
		egenes[ENSG]["info"].append("tested")
	#check egenes in DICE
	if bareENSG in DICE_setegenes:
		egenes[ENSG]["info"].append("egene")
	else:
		egenes[ENSG]["info"].append("NA")
	#add sig lead in DICE
	DICE_sig = egenes[ENSG]["DICE_sig"]
	egenes[ENSG]["info"].append(DICE_sig)
	#if no overlap between significant SNPs and eQTLgen SNPs, not tested
	if sig_SNPs.isdisjoint(eQTLgen_snps):
		egenes[ENSG]["info"].append("not_tested")
	else:
		egenes[ENSG]["info"].append("tested")
	#check egenes in eQTLgen
	if bareENSG in eQTLgen_setegenes:
		egenes[ENSG]["info"].append("egene")
	else:
		egenes[ENSG]["info"].append("NA")
	#add sig lead in eQTLgen
	eQTLgen_sig = egenes[ENSG]["eQTLgen_sig"]
	egenes[ENSG]["info"].append(eQTLgen_sig)

###################################################
##get info for haploinsuffciency and IUIS data
###################################################
print("checking for IUIS and haploinsufficient egenes")

#make set of DICE and eQTLgen SNPs
IUIS_genes = set()
with open(IUIS[1],'rt') as f:
	for line in f:
		stuff = line.rstrip().split("\t")
		
		IUIS_genes.add(stuff[0])

haplo_top_genes = set()
with open(haplo_top[1],'rt') as f:
	for line in f:
		stuff = line.rstrip().split("\t")
		
		haplo_top_genes.add(stuff[0])

haplo_all_genes = set()
with open(haplo_all[1],'rt') as f:
	for line in f:
		stuff = line.rstrip().split("\t")
		
		haplo_all_genes.add(stuff[0])

#step through egenes and add summary
header_out.append("IUIS")
header_out.append("top_haploinsufficient")
header_out.append("all_haploinsufficient")
for ENSG in egenes:
	#check for ENSG in IUIS set
	if ENSG in IUIS_genes:
		egenes[ENSG]["info"].append("yes")
	else:
		egenes[ENSG]["info"].append("no")
	#check for ENSG in top haplo set
	if ENSG in haplo_top_genes:
		egenes[ENSG]["info"].append("yes")
	else:
		egenes[ENSG]["info"].append("no")
	#check for ENSG in all haplo set
	if ENSG in haplo_all_genes:
		egenes[ENSG]["info"].append("yes")
	else:
		egenes[ENSG]["info"].append("no")

###################################################
##cAID lookups
###################################################
print("checking for cAID proxies")

#make set of cAID SNPs
cAID_set_snps = set()
with open(cAID_snps[1],'rt') as f:
	header = f.readline()
	for line in f:
		stuff = line.rstrip().split("\t")
		
		cAID_set_snps.add(stuff[0])
header_out.append("cAID_snp")
for ENSG in egenes:
	sig_SNPs = egenes[ENSG]["rsid"]
	#if no overlap between significant SNPs and cAID SNPs, no
	if sig_SNPs.isdisjoint(cAID_set_snps):
		egenes[ENSG]["info"].append("no")
	else:
		egenes[ENSG]["info"].append("yes")

###################################################
##parse ATAC overlaps
###################################################
print("checking for ATAC overlaps")

ATAC_genes = set()
with open(ATAC[1],'rt') as f:
	for line in f:
		stuff = line.rstrip().split("\t")
		
		ENSG = stuff[7]
		ATAC_genes.add(ENSG)
		
header_out.append("ATAC")
for ENSG in egenes:
	#if gene in set, yes
	if ENSG in ATAC_genes:
		egenes[ENSG]["info"].append("yes")
	else:
		egenes[ENSG]["info"].append("no")


###################################################
##parse CaptureC overlaps
###################################################
print("checking for CaptureC links")

CaptureC_genes = set()
with open(onefrag[1],'rt') as f:
	for line in f:
		stuff = line.rstrip().split("\t")
		#have to parse to see if our egene is the link 
		ENSG = stuff[14]
		link = stuff[9].split("|")
		
		for l in link:
			transcript, gene = l.split("+")
			
			#get ENSG id from gene names if there. 
			testENSG = genenames.get(gene, "NA")
			
			#if not in genenames, skip
			if testENSG == "NA":
				continue
			else:
				#otherwise test if genes match. if they do, add to set
				if testENSG == ENSG:
					CaptureC_genes.add(ENSG)
		
with open(fourfrag[1],'rt') as f:
	for line in f:
		stuff = line.rstrip().split("\t")
		#have to parse to see if our egene is the link 
		ENSG = stuff[14]
		link = stuff[9].split("|")
		
		for l in link:
			transcript, gene = l.split("+")
			
			#get ENSG id from gene names if there. 
			testENSG = genenames.get(gene, "NA")
			
			#if not in genenames, skip
			if testENSG == "NA":
				continue
			else:
				#otherwise test if genes match. if they do, add to set
				if testENSG == ENSG:
					CaptureC_genes.add(ENSG)

header_out.append("CaptureC")
#in this cycle, also initialize the DICE and tonsil trait sets
for ENSG in egenes:
	#if gene in set, yes
	if ENSG in CaptureC_genes:
		egenes[ENSG]["info"].append("yes")
	else:
		egenes[ENSG]["info"].append("no")
	egenes[ENSG]["DICE"] = set()
	egenes[ENSG]["tonsil"] = set()
	
###################################################
##parse COLOCs
###################################################
print("checking for COLOCs")

for trait in traits:
	#make summary filename
	sumfile = colocsummarybase[1].replace("TRAIT", trait)
	with open(sumfile,'rt') as f:
		#skip header
		header = f.readline()
		for line in f:
			stuff = line.rstrip().split("\t")
			
			otter = stuff[2].split("_")
			ENSG = otter.pop(0)
			tissue = "_".join(str(x) for x in otter)
			PP4 = stuff[8]
			
			#if colocalization
			if float(PP4) > 0.8:
				#if tonsil focus add to tonsil set
				if tissue == tonsil_focus[1]:
					egenes[ENSG]["tonsil"].add(trait)
				#if DICE focus
				if tissue == DICE_focus[1]:
					if ENSG in genebareENSG:
						fullENSG = genebareENSG[ENSG]
						egenes[fullENSG]["DICE"].add(trait)
			
header_out.append(tonsil_focus[1]+"_colocs")
header_out.append("DICE_"+tonsil_focus[1]+"_colocs")
for ENSG in egenes:
	#make set of tonsil traits a list, add as column
	value = len(egenes[ENSG]["tonsil"])
	if value == 0:
		egenes[ENSG]["tonsil"].add("NA")
	setlist = ", ".join(str(x) for x in sorted(egenes[ENSG]["tonsil"]))
	egenes[ENSG]["info"].append(setlist)
	
	#make set of DICE traits a list, add as column
	value = len(egenes[ENSG]["DICE"])
	if value == 0:
		egenes[ENSG]["DICE"].add("NA")
	setlist = ", ".join(str(x) for x in sorted(egenes[ENSG]["DICE"]))
	egenes[ENSG]["info"].append(setlist)

###################################################
##output everything
###################################################

outfile = open(Out, "w")

#make header
outfile.write(tab.join(str(x) for x in header_out)+"\n")


for ENSG in egenes:
	outlist = egenes[ENSG]["info"]
	#output the full line
	outfile.write(tab.join(str(x) for x in outlist)+"\n")


