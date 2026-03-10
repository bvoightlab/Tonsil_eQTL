#programs required:
rnaseqc/2.3.6 bcftools/1.9 picard/2.23.3 htslib/1.9 
samtools/1.9 rsem/1.3 python/2.7.10 multiqc/1.4
R/4.2

#this pipeline starts from a list of .fastq files and results in .gct files, both tpm and readcount.
#it is run separately by cell type, for clarity scripts are only provided for GCB cell type
#we started from https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md 
#in Jan 2022 so listing all code here as there have been changes to their pipeline since. 

--------------------------------------------------------
####make reference & input files
#GCB_path_list_first.txt lines look like the following:
#/scratch/lorenzk/TC202GCB	TC202GCB_ds.d2792d10621c4edf8946eaf884bd6eca/TC202GCB_S187
#tab separated, first path is where I want output, second path is location of fastq files & base filename

grep "GCB" file_list.txt > file_list_GCB.txt
awk -F'/' '{ split($7,a,"_"); print "/scratch/lorenzk/" a[1] "\t" "/" $2 "/" $3 "/" $4 "/" $5 "/" $6 "/" a[1] "_" a[2] }' file_list_GCB.txt | awk '!x[$0]++' > GCB_path_list.txt
head -n 52 GCB_path_list.txt  > GCB_path_list_first.txt
tail -n 53 GCB_path_list.txt  > GCB_path_list_second.txt
awk -F'[/\t]' '{print $4}' GCB_path_list_first.txt > GCB_path_list_first_markdup.txt
awk -F'[/\t]' '{print $4}' GCB_path_list_second.txt > GCB_path_list_second_markdup.txt

#sample participant files look like
#sample_id	participant_id
#TC202GCB	TC202_TC202
awk '{print $1 "GCB\t" $2}' sample_participant_lookup.txt > sample_participant_lookup_GCB.txt
sed -i '1 i\sample_id\tparticipant_id' sample_participant_lookup_GCB.txt

--------------------------------------------------------

#fastQC
bsub -o lsf_log.out python scripts/fastqc_manager.py -i GCB_path_list_first.txt -n fastqc -j 20

bsub -o lsf_log.out python scripts/fastqc_manager.py -i GCB_path_list_second.txt -n fastqc -j 20  

#once fastQC finishes, run multiqc from fastqc directory
cd /pathto/fastqc
module unload python/2
module load python/2.7.10
module load multiqc/1.4
multiqc .

--------------------------------------------------------

#RNA alignment

#run alignment with STAR, confirm successes
bsub -o lsf_log.out -m vengeance python scripts/STAR_manager_wHLA.py -i GCB_path_list_first.txt -n eQTL -j 10
ls TC*GCB_lsf_log.out | wc -l
grep "Successfully" TC*GCB_lsf_log.out | wc -l
mv * /pathto/STAR_alignments

bsub -o lsf_log.out -m vengeance python scripts/STAR_manager_wHLA.py -i GCB_path_list_second.txt -n eQTL -j 10
ls TC*GCB_lsf_log.out | wc -l
grep "Successfully" TC*GCB_lsf_log.out | wc -l
mv * /pathto/STAR_alignments

cd /pathto/STAR_alignments


#run mark duplicates and run rsem
bsub -o lsf_log.out python scripts/markdup_rsem_manager.py -i GCB_path_list_first_markdup.txt -n markdup -j 20
ls *GCBmarkdup_lsf_log.out | wc -l
grep "Successfully" *GCBmarkdup_lsf_log.out | wc -l
ls *GCBrsem_lsf_log.out | wc -l
grep "Successfully" *GCBrsem_lsf_log.out | wc -l
bsub -o lsf_log.out python scripts/markdup_rsem_manager.py -i GCB_path_list_second_markdup.txt -n markdup -j 20
ls *GCBmarkdup_lsf_log.out | wc -l
grep "Successfully" *GCBmarkdup_lsf_log.out | wc -l
ls *GCBrsem_lsf_log.out | wc -l
grep "Successfully" *GCBrsem_lsf_log.out | wc -l

#rename bam files
bsub -o lsf_log.out python scripts/bam_manager.py -i GCB_path_list_first_markdup.txt -n bam -j 10
ls *GCBbam_lsf_log.out | wc -l
grep "Successfully" *GCBbam_lsf_log.out | wc -l
bsub -o lsf_log.out python scripts/bam_manager.py -i GCB_path_list_second_markdup.txt -n bam -j 10
ls *GCBbam_lsf_log.out | wc -l
grep "Successfully" *GCBbam_lsf_log.out | wc -l

#run rnaseqc to make individual gct files
bsub -o lsf_log.out python scripts/rnaseqc_manager.py -i GCB_path_list_first_markdup.txt -n markdup -j 10
ls *GCBrnaseqc_lsf_log.out | wc -l
grep "Successfully" *GCBrnaseqc_lsf_log.out | wc -l
bsub -o lsf_log.out python scripts/rnaseqc_manager.py -i GCB_path_list_second_markdup.txt -n markdup -j 10
ls *GCBrnaseqc_lsf_log.out | wc -l
grep "Successfully" *GCBrnaseqc_lsf_log.out | wc -l

--------------------------------------------------------

run multiqc on aligned reads 

cd /pathto/STAR_alignments
module unload python/2
module load python/2.7.10
module load multiqc/1.4
multiqc .

----------------------------------
#see genotyping readme for details on comparing RNA sequencing genotypes with array genotypes
#in this process, we determined some sample names were swapped. the following changes were made:

#cell type swaps were further confirmed using gene expressions of genes encoding the proteins we sorted cells upon (shown below with a *) AND lineage defining transcripts (shown below with a ^) 
#Naïve B cells: CD19*, IGD*
#GC B cells: CD19*, CD38*, AIDCA^, MKI67^, UNG^
#Naïve T cells: CD3*, CD4*, CCR7^
#Tfh cells: CD3*, CD4*, PDCD1*, CXCR5*, BCL6^, ICOS^, CD40LG^, SH2D1A^, IL21^

#remove TC335 TFH sample from collection list
ls TC335TFH*
rm TC335TFH*

#swap NaiveB 274 and 273; 220 and 221
mv TC274NaiveB.gene_reads.gct.gz placeholder.gene_reads.gct.gz
mv TC273NaiveB.gene_reads.gct.gz TC274NaiveB.gene_reads.gct.gz
mv placeholder.gene_reads.gct.gz TC273NaiveB.gene_reads.gct.gz

mv TC274NaiveB.gene_tpm.gct.gz placeholder.gene_tpm.gct.gz
mv TC273NaiveB.gene_tpm.gct.gz TC274NaiveB.gene_tpm.gct.gz
mv placeholder.gene_tpm.gct.gz TC273NaiveB.gene_tpm.gct.gz 

mv TC220NaiveB.gene_reads.gct.gz placeholder.gene_reads.gct.gz
mv TC221NaiveB.gene_reads.gct.gz TC220NaiveB.gene_reads.gct.gz
mv placeholder.gene_reads.gct.gz TC221NaiveB.gene_reads.gct.gz

mv TC220NaiveB.gene_tpm.gct.gz placeholder.gene_tpm.gct.gz
mv TC221NaiveB.gene_tpm.gct.gz TC220NaiveB.gene_tpm.gct.gz
mv placeholder.gene_tpm.gct.gz TC221NaiveB.gene_tpm.gct.gz

#swap GCB 220 and 221
mv TC220GCB.gene_reads.gct.gz placeholder.gene_reads.gct.gz
mv TC221GCB.gene_reads.gct.gz TC220GCB.gene_reads.gct.gz
mv placeholder.gene_reads.gct.gz TC221GCB.gene_reads.gct.gz

mv TC220GCB.gene_tpm.gct.gz placeholder.gene_tpm.gct.gz
mv TC221GCB.gene_tpm.gct.gz TC220GCB.gene_tpm.gct.gz
mv placeholder.gene_tpm.gct.gz TC221GCB.gene_tpm.gct.gz

#relabel sample currently labeled TC304 TFH as TC303 NaiveB 
#and sample currently labeled TC303 NaiveB as TC304 TFH

mv TC304TFH.gene_reads.gct.gz placeholder.gene_reads.gct.gz
mv TC303NaiveB.gene_reads.gct.gz TC304TFH.gene_reads.gct.gz
mv placeholder.gene_reads.gct.gz TC303NaiveB.gene_reads.gct.gz

mv TC304TFH.gene_tpm.gct.gz placeholder.gene_tpm.gct.gz
mv TC303NaiveB.gene_tpm.gct.gz TC304TFH.gene_tpm.gct.gz
mv placeholder.gene_tpm.gct.gz TC303NaiveB.gene_tpm.gct.gz

----------------------------------

#after all samples were appropriately labeled, they were combined into .gct file by cell type

ls -1 $PWD/*GCB.gene_tpm.gct.gz > GCB_tpm_gct_list.txt
python scripts/combine_GCTs.py GCB_tpm_gct_list.txt GCB_tpm

ls -1 $PWD//*GCB.gene_reads.gct.gz > GCB_readcount_gct_list.txt
python scripts/combine_GCTs.py GCB_readcount_gct_list.txt GCB_readcount
