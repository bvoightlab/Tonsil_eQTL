This colocalization pipeline uses Colocquial https://github.com/bvoightlab/ColocQuiaL
As input it request eQTL allpairs files and gwas summary statistics
It outputs pp4 and conditional pp4, as well as locuscompare, locuszoom and looping plots
plots require pygenometracks and looping data
--------------------------------------------------------

#GWAS summary statistics were downloaded from sources described in supplementary table 4

#We then harmonized data by lifting over to hg38 and assigning rsIDs based on dbSNP155. 
#an example .yml config file is provided. 

python process_gwas_summary_stats.py -c Allergic.yml

#then reformat for locus ID script
for i in Asthma pAsthma Atopic_Asthma Child_Asthma ADE Allergic;
do bsub -o reformat.out Rscript scripts/reformat_allergic.R "$i".quant.processed.bed;
done

#Find intervals not in 
for i in Asthma pAsthma Atopic_Asthma Child_Asthma ADE Allergic;
do bsub -o locus_def_out.txt python3 scripts/define_sig_loci_for_trait_extra_padding.py \
-g "$i".temp.txt \
-o "$i"_intervals.txt \
-t 50000 \
-c grch38_chr_bounds_mhc_only.bed \
-p 5e-8 \
-s 5e-5 \
-n TRAIT.ANC.txt;
done

#combine any intervals that overlap (not an issue for this data)

#add chr to make bed compatible
sed -i 's/^/chr/' *intervals.bed

#for each study, pulled list of leads from supplementary material, lifting over to hg38 if needed
#then overlapped with above, only used intervals with a published lead
#this was done because some summary statistics were data subsets (excluding 23&me data) and some lead lists included non genome-wide significant leads. 
#Note: for Atopic Asthma and Child Asthma studies, summary statistics were in hg19, but published leads were already in hg38
#for ADE study, many duplicate SNP positions were present in the summary statistics. To facilitate analysis, these were removed. 

for i in Asthma pAsthma Atopic_Asthma Child_Asthma ADE Allergic;
do bedtools intersect -a "$i"/"$i"_intervals.bed -b "$i"/"$i"_leads_hg38.bed -loj > "$i"/"$i"_leads_hg38.signal_check.bed;
echo "# intervals that don't overlap";
bedtools intersect -a "$i"/"$i"_intervals.bed -b "$i"/"$i"_leads_hg38.bed -loj -v | wc -l ;
echo "leads";
bedtools intersect -a "$i"/"$i"_leads_hg38.bed -b "$i"/"$i"_intervals.bed | wc -l ;
bedtools intersect -a "$i"/"$i"_leads_hg38.bed -b "$i"/"$i"_intervals.bed -v;
done

for i in Asthma pAsthma Atopic_Asthma Child_Asthma ADE Allergic;
do grep -v '\.' "$i"/"$i"_leads_hg38.signal_check.bed | awk '!x[$1, $2, $3]++' > "$i"/"$i"_overlapping_intervals.bed;
wc -l "$i"/"$i"_overlapping_intervals.bed
done

#pull out lead SNP (minimum p value)
for i in Asthma pAsthma Atopic_Asthma Child_Asthma ADE Allergic;
do bsub -o leads.out Rscript scripts/minp_regions_auto.R "$i"/"$i"_overlapping_intervals.bed "$i".quant.processed.bed "$i"/;
done

#the lead SNPs and intervals tested are provided in supplementary table 11

--------------------------------------------------------

#we reformatted eQTL data as described by Colocquial
#then ran colocquial per GWAS in folders using
bash colocquial_wrapper.sh

#after completion, results were collected and processed as follows

for i in Asthma pAsthma Atopic_Asthma Child_Asthma ADE Allergic;
do cd "$i";
bash summarize_results.bsub;
cd ..;
done

head -n 1 ADE/ADE_eqtl_coloc_results_all_summaryspread_condPP4.txt > all_sig_colocs.txt
for i in Asthma pAsthma Atopic_Asthma Child_Asthma ADE Allergic;
do sed -e '1d' "$i"/"$i"_eqtl_coloc_results_all_summaryspread_condPP4.txt | awk '$9 > 0.8' >> all_sig_colocs.txt;
done;

--------------------------------------------------------

#plots

#locuscompare plots were made using R and plink

#the eQTL x trait version requires as input the colocquial input file, and the genotyping plink files (for getting pairwise LD)
#Figure 3A
lsf_locuscompare_colocquial_insample_paper_TRAF3.R

#Figure 4A
lsf_locuscompare_colocquial_insample_paper_ZBTB10.R

#Figure 5A
lsf_locuscompare_colocquial_insample_paper_JAZF1.R

#the eQTL x eQTL version greps gene info from the cell type all pairs files, and again uses the genotyping plink files for pairwise LD
#Figure 5B
eQTL_eQTL_locuscompare_paper_JAZF1_NB_NT.R

#Figure 5C
eQTL_eQTL_locuscompare_paper_JAZF1_NB_GCB.R

#Figure 5D
eQTL_eQTL_locuscompare_paper_JAZF1_NT_TFH.R


#stacked locus plots were made using pygenometracks and bedtools2 in python 3.8

#bigwig ATAC and bedpe looping files were provided by their paper authors
#genes are plotted according to gencode.v34.annotation.gtf
#also needs the credible set bounds output by eQTL_ABF.R in eQTL pipeline

#other input files for these plots were created using
bsub -o pygenometracks_prep_out.txt Rscript scripts/prep_locuszoom_for_pygenometracks.R NaiveT GCB ENSG00000131323.16 TRAF3 pAsthma 102777476 102911500
bsub -o pygenometracks_prep_out.txt Rscript scripts/prep_locuszoom_for_pygenometracks.R TFH GCB ENSG00000205189.12 ZBTB10 Asthma 80485619 80526265    
#these require the GWAS trait summary statistics, the all pairs files for relevant cell types, and use the genotyping plink file for pairwise LD

#Figure 3B:
pyGenomeTracks --dpi 300 --fontSize 20 --plotWidth 41 --tracks TRAF3_NaiveT_GCB_pAsthma_lzloops3.ini -o TRAF3_NaiveT_GCB_pAsthma_lzloops_paper2.png --region chr14:102764707-102915467

#Figure 4B:
pyGenomeTracks --dpi 300 --fontSize 20 --plotWidth 41 --tracks ZBTB10_TFH_pAsthma_paper.ini -o ZBTB10_TFH_pAsthma_paper.png --region chr8:80370736-80530468

the .ini files for the above are provided 
