We used Genome Studio 2.0.551 to filter SNPs from the raw array calls if they had the following characteristics: AB T mean > 0.15, AB R mean < 0.3, call frequency < 0.5, cluster separation < 0.2748, het excess < -0.6 or > 0.4. 
Additionally, we removed all X chromosome SNPs with call frequency < 0.95. 
We filtered manually if the call frequency was between 0.5 and 0.75, het excess was between -0.6 and -0.35, or 0.2 and 0.4; removing SNPs in these zones without clearly clustered calls. 
We modified the recommended heterozygote excess filter due to our mixed ancestry population. 
Our final call set included array called genotypes for 566,408 SNPs.

This file is genotypes_22_02_23_array.vcf and is archived on dbGAP

----------------------------------------------------

IBD calculations

#using plink 1.90Beta6.18, was run on array level genotypes (no imputation)

plink --noweb --allow-no-sex --bfile genotypes_22_02_23 --indep-pairwise 10000kb 1 0.05 --out genotypes_22_02_23_prune

plink --noweb --allow-no-sex --bfile genotypes_22_02_23 --extract genotypes_22_02_23_prune.prune.in --genome --out genotypes_22_02_23_pruned

#then compared pihat column of the output genotypes_22_02_23_pruned.genome file

----------------------------------------------------

#using bcftools/1.9, vcftools and htslib/1.9 I split files into individual chromosomes for submission to Topmed

for i in {1..23};
do vcftools --vcf genotypes_22_02_23_all.vcf --recode --chr "$i" --out genotypes_22_02_23_chr"$i";
done

#need to rename chr 23 X. 

echo "23 X" >> chr_name_conv.txt

module load bcftools/1.9

bcftools annotate --rename-chrs chr_name_conv.txt genotypes_22_02_23_chr23.recode.vcf | bgzip > genotypes_22_02_23_chr23.rename.vcf.gz

for i in {1..22};
do bgzip genotypes_22_02_23_chr"$i".recode.vcf;
done

#submitted all 23 chromosomes to TopMed at https://imputation.biodatacatalyst.nhlbi.nih.gov/#! on 02/24/2022

#submitted without allele frequency check, but with r2 filter of .3

QC returned:
23 valid VCF file(s) found.

Samples: 105
Chromosomes: 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 X 3 4 5 6 7 8 9
SNPs: 734073
Chunks: 309
Datatype: unphased
Build: hg19
Reference Panel: apps@topmed-r2@1.0.0 (hg38)
Population: mixed
Phasing: eagle
Mode: imputation
Rsq filter: 0.3
Quality Control

Uploaded data is hg19 and reference is hg38.

Lift Over

Skip allele frequency check.

Calculating QC Statistics

Statistics:
Alternative allele frequency > 0.5 sites: 0
Reference Overlap: 96.26 %
Match: 457,943
Allele switch: 83,776
Strand flip: 1
Strand flip and allele switch: 0
A/T, C/G genotypes: 831
Filtered sites:
Filter flag set: 0
Invalid alleles: 167,704
Multiallelic sites: 0
Duplicated sites: 1,590
NonSNP sites: 0
Monomorphic sites: 0
Allele mismatch: 470
SNPs call rate < 90%: 1,267

Excluded sites in total: 171,032
Remaining sites in total: 541,283
See snps-excluded.txt for details
Typed only sites: 21,084
See typed-only.txt for details

----------------------------------------------------

confirm array genotyping and expression sequencing match - must be done AFTER RNAseq pipeline
(checked all cell type expression data, only showing one)

#pull list of variants from genotyping array
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' genotypes_22_02_23_all.vcf.gz > als.tsv

#add chr so they match bamfiles
sed 's/^/chr/' als.tsv > als_chr.tsv

#sample list in same format as in RNAseq pipeline

#make vcf from bam files
ls -1 $PWD/*Aligned.sortedByCoord.out.id.bam > GCB_bamlist_for_variantcalling.txt
bsub -o bcftools_variant_calling.out 'bcftools mpileup -Ou -d 8000 -f reference_files/Homo_sapiens_assembly38_noALT_noDecoy.fasta -T als_chr.tsv -b GCB_bamlist_for_variantcalling.txt | bcftools call -mv -T als_chr.tsv -Oz -o GCB_variant_calls.vcf.gz'

#split vcf by sample & run gtcheck to compare sample to array genotypes
while IFS= read -r line;
do bcftools view -Oz -s "$line" /pathto/STAR_alignments/GCB_variant_calls.vcf.gz > /pathto/STAR_alignments/"$line".vcf.gz;
tabix -p vcf /pathto/STAR_alignments/"$line".vcf.gz;
bsub -o lsfout_GCB_gtcheck.txt "bcftools gtcheck -g /pathto/ImputationII_withX/all_chr.maf01.Rsq5.vcf.gz /pathto/STAR_alignments/"$line".vcf.gz > "$line"_gtcheck.txt";
done < GCB_path_list_first_markdup.txt

#check outputs
while IFS= read -r line;
do echo -n $line;
tail -n 1 "$line"_gtcheck.txt | awk -F '[\t_]' '{print " " $5}';
done < GCB_path_list_first_markdup.txt > GCB_gtcheck_output.txt

#gtcheck was run before (to identify sample swaps) and after (to confirm new labels were correct, using final vcf)

----------------------------------------------------

#imputation file was filtered using bcftools

filter -i 'INFO/R2 > 0.5 && INFO/MAF > 0.01' all_chr.dose.vcf.gz

#after imputation we discovered some sample swaps. they were corrected as follows using bcftools reheader

swap the following samples:
TC331 for TC332
TC279 for TC285
TC280 for TC284
TC282 for TC283

#using bcftools/1.9
bcftools reheader -s samples_file.txt -o all_chr.maf01.Rsq5.vcf.gz all_chr.maf01.Rsq5.vcf.gz

#during initial IBD analysis we identified a trio of siblings. 2 were removed using bcftools

bcftools view -s^TC277_TC277,TC278_TC278 all_chr.maf01.Rsq5.vcf.gz -Oz -o all_chr.final.vcf.gz 

#final imputation vcf is all_chr.final.vcf.gz and is also archived at dbGAP

----------------------------------------------------

Genotyping PCs 

#were obtained as follows using plink 1.90Beta4.5
 
plink --vcf all_chr.maf01.Rsq5.vcf.gz --indep-pairwise 10000kb 1 0.05 --out ldp1

plink --vcf all_chr.maf01.Rsq5.vcf.gz --extract ldp1.prune.in --pca --out all_chr_LDpruned

#then used plot_PCs.R

