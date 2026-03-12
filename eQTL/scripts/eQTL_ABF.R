library(tidyverse)


#Read in 
#args <- c("TFH", "ENSG00000205189", "ZBTB10")
#args <- c("NaiveT", "ENSG00000101246", "ARFRP1", "pAsthma", "63698642", "63708025")
#args <- c("NaiveB", "ENSG00000107485", "GATA3", "Asthma", "8053604", "8075198")
#args <- c("GCB", "ENSG00000228463.10", "AP006222.1")

args <- commandArgs(trailingOnly = TRUE)

cell1 <- args[1]
ENSG <- args[2]
gene <- args[3]

file_1 <- paste0("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/gene_allpairs/", ENSG, "_", cell1, "_allpairs.txt")

#check if files exist, run grep if not
if (file.exists(file_1)){
  df <- read.table(file_1, header=FALSE, fill=TRUE)

} else {
  system(paste0("zgrep ", ENSG, " /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/", cell1, "_all.cis_qtl_pairs.txt.gz > /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/gene_allpairs/", ENSG, "_", cell1, "_allpairs.txt"))
  Sys.sleep(10)
  df <- read.table(file_1, header=FALSE, fill=TRUE)
  
}


colnames(df) <- c("phenotype_id", "variant_id", "tss_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")


#separate variant ID to get pos and chr
df <- separate(df, variant_id, c("chr", "pos", "a1", "a2"), remove = FALSE)
df$pos <- as.numeric(df$pos)
#chr <- df$chr[1]

#remove rows with NA
df = df[complete.cases(df), ]


#calculate Bayes Factor & PIP
#variance prior set to 0.15 as eQTL are quantitative traits
df$W2 <- (0.15)^2
df$SNP_BF <- (sqrt(df$slope_se^2 / (df$slope_se^2 + df$W2)) * exp((df$W2 / (df$slope_se^2 + df$W2)) * ((df$slope / df$slope_se)^2 / 2)))
df$PIP <- df$SNP_BF/sum(df$SNP_BF)

ABF <- df %>% mutate(rank = rank(1-PIP)) %>% arrange(rank)
ABF$cum_PIP <- cumsum(ABF$PIP)

#pull SNPs in credible set
#find index of 95% cred set
CS_index <- max(which(ABF$cum_PIP < .95)) +1
if (CS_index == -Inf){
  CS_index <- 1
}
credset <- ABF[0:CS_index,]
credset <- select(credset, chr, pos, a1, a2, pval_nominal, slope, slope_se, PIP)

#add gene and tissue info for eventual consolidating
credset <- credset %>% add_column(cell_type = cell1, .before = "chr")
credset <- credset %>% add_column(gene_name = gene, .before = "chr")
credset <- credset %>% add_column(ENSG_ID = ENSG, .before = "chr")


#output abf credible set
#filename
out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/ABF/", cell1,  "/", gene, "_", ENSG, "_", cell1, "_abf_credset.txt", sep='')
#write table
write.table(credset, file= out, sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)

#output line for bedfile with first and last SNP in credset, by position
CS_min <- min(credset$pos)
CS_max <- max(credset$pos)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/ABF/", cell1,  "/", gene, "_", cell1, "_abf_credset_bounds.bed", sep='')
bounds <- data.frame(chr, CS_min, CS_max)
write.table(bounds, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

