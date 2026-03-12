# Tonsil_eQTL

This repository contains code relating to the following paper  
Tonsillar expression quantitative trait loci verify and expand genetic contributors to childhood atopic diseases  
Kim Lorenz, Samuel Yoon, Carole Le Coz, Karen Zur, Andrew Wells, Neil Romberg, Benjamin F. Voight  


##Genotyping 

explains modifications made to raw array calls, imputation, and final vcf QCs


##RNAseq

explains how we got from fastq files to .gct files

##differential_expression

takes as input the .gct file outputs of RNAseq and a covariates file  
outputs log fold change differences across naive vs active, naive vs naive, and active vs active cell type comparisons, as well as corrsponding upset plots

##eQTL

takes as input the .gct file outputs of RNAseq, a covariates file, the genotypes vcf from genotyping, and the collapse gene .gtf file from GTEx  
outputs eQTL summary statistics and Upset plots  
with additional annotation files, code is provided to recreate annotations in supplementary table 1  

##colocalization
uses GWAS trait summary statistics and the [ColocQuial pipeline](https://github.com/bvoightlab/ColocQuiaL)  
outputs list of colocalization probabilities, locuscompare plots, and stacked locusplots

