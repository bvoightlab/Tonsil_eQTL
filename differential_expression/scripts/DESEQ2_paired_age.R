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

#match covariate files to meta file rownames
t1t <- data.frame(t(trait1_meta))
colnames(t1t) <- trait1_meta$cov
t1t <- t1t[-(1:1),]
t1t$celltype <- trait1stuff[1]
t1t <- rownames_to_column(t1t, "sample")
t1t <- separate(t1t, sample, c("sample", "header"), 5)
t1t$header <- str_c(t1t$sample, t1t$celltype)

t2t <- data.frame(t(trait2_meta))
colnames(t2t) <- trait2_meta$cov
t2t <- t2t[-(1:1),]
t2t$celltype <- trait2stuff[1]
t2t <- rownames_to_column(t2t, "sample")
t2t <- separate(t2t, sample, c("sample", "header"), 5)
t2t$header <- str_c(t2t$sample, t2t$celltype)

both_cov <- bind_rows(t1t, t2t)
both_cov <- select(both_cov, -c(sample, celltype))
full_meta <- inner_join(meta, both_cov , "header")

#change explicit cov to factors
makefac <- c("Age", "Sex", "Infectious", "Obstructive")
full_meta[makefac] <- sapply(full_meta[makefac],as.numeric)


#make M and F splits
meta4 <- filter(full_meta, Age <= 4)
meta5 <- filter(full_meta, Age >= 5)

merged4  <- merged[,meta4$header]
merged5  <- merged[,meta5$header]

#run DESeq2

#import data
DE_data <- DESeqDataSetFromMatrix(countData = merged4, 
                                  colData = meta4,
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


resultsfile <- str_c("/project/voight_viz/lorenzk/Romberg_eQTL/DESeq2/", trait1stuff[1], "_", trait2stuff[1], "_DESEQ_results_4paired.txt")
write.table(res, file= resultsfile, sep="\t", row.names=TRUE, quote=FALSE)

resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.0025)

resultsfile <- str_c("/project/voight_viz/lorenzk/Romberg_eQTL/DESeq2/", trait1stuff[1], "_", trait2stuff[1], "_DESEQ_results_4paired_sig.txt")
write.table(resSig, file= resultsfile, sep="\t", row.names=TRUE, quote=FALSE)

#import data
DE_data <- DESeqDataSetFromMatrix(countData = merged5, 
                                  colData = meta5,
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


resultsfile <- str_c("/project/voight_viz/lorenzk/Romberg_eQTL/DESeq2/", trait1stuff[1], "_", trait2stuff[1], "_DESEQ_results_5paired.txt")
write.table(res, file= resultsfile, sep="\t", row.names=TRUE, quote=FALSE)

resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.0025)

resultsfile <- str_c("/project/voight_viz/lorenzk/Romberg_eQTL/DESeq2/", trait1stuff[1], "_", trait2stuff[1], "_DESEQ_results_5paired_sig.txt")
write.table(resSig, file= resultsfile, sep="\t", row.names=TRUE, quote=FALSE)
