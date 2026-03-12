library(tidyverse)
library(DESeq2)
 

#use args to get the input files
args = commandArgs(TRUE)
#testing purposes:
#args <- c("NaiveB,/project/voight_viz/lorenzk/Romberg_eQTL/fixed_samplenames/NaiveB_readcount.gct.gz,/project/voight_viz/lorenzk/Romberg_eQTL/tensorqtl/NaiveB/NaiveB_peer12.combined_covariates_header.txt","GCB,/project/voight_viz/lorenzk/Romberg_eQTL/fixed_samplenames/GCB_readcount.gct.gz,/project/voight_viz/lorenzk/Romberg_eQTL/tensorqtl/GCB/GCB_peer11.combined_covariates_header.txt")
#args <- c("NaiveT,/project/voight_viz/lorenzk/Romberg_eQTL/fixed_samplenames/NaiveT_readcount.gct.gz,/project/voight_viz/lorenzk/Romberg_eQTL/tensorqtl/NaiveT/NaiveT_peer13.combined_covariates_header.txt","NaiveB,/project/voight_viz/lorenzk/Romberg_eQTL/fixed_samplenames/NaiveB_readcount.gct.gz,/project/voight_viz/lorenzk/Romberg_eQTL/tensorqtl/NaiveB/NaiveB_peer12.combined_covariates_header.txt")
trait1stuff <- strsplit(args[1], ",")[[1]]
trait2stuff <- strsplit(args[2], ",")[[1]]
#assumes each arguement is traitname,path/to/gctfile,path/to/covfile

trait1 <- read.delim(trait1stuff[2], header = TRUE, skip=2)
trait2 <- read.delim(trait2stuff[2], header = TRUE, skip=2)

trait1_meta <- read.delim(trait1stuff[3], header = TRUE)
trait2_meta <- read.delim(trait2stuff[3], header = TRUE)

##merge counts files & print
merged <- inner_join(trait1, select(trait2, !Description), "Name")
rownames(merged) <- merged$Name
merged[,2] <- NULL
merged[,1] <- NULL

##make metadata file
meta <- data.frame(cbind(colnames(merged)))
colnames(meta) <- "header"
meta$header2 <- meta$header
meta <- separate(meta, header2, c("sample", "celltype"), 5)

#run DESeq2

#import data
DE_data <- DESeqDataSetFromMatrix(countData = merged, 
                                  colData = meta,
                                  design = formula(paste("~ sample + celltype")))

#prefilter
keep <- rowSums(counts(DE_data)) >= 10
DE_data <- DE_data[keep,]

#run DESeq
#DE_out <- DESeq(DE_data)
DE_data <- estimateSizeFactors(DE_data)
DE_data <- estimateDispersions(DE_data)
DE_data <- nbinomWaldTest(DE_data, maxit=1000)
res <- results(DE_data)

resultsNames(DE_data)
#resLFC <- lfcShrink(DE_data, coef="celltype_NaiveB_vs_GCB", type="apeglm")
#plotMA(resLFC, ylim=c(-2,2))

#check degrees of freedom
ncol(model.matrix(design(DE_data), colData(DE_data)))


resultsfile <- str_c("/project/voight_viz/lorenzk/Romberg_eQTL/DESeq2/", trait1stuff[1], "_", trait2stuff[1], "_DESEQ_results_paired.txt")
write.table(res, file= resultsfile, sep="\t", row.names=TRUE, quote=FALSE)

resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.0025)

resultsfile <- str_c("/project/voight_viz/lorenzk/Romberg_eQTL/DESeq2/", trait1stuff[1], "_", trait2stuff[1], "_DESEQ_results_paired_sig.txt")
write.table(resSig, file= resultsfile, sep="\t", row.names=TRUE, quote=FALSE)

