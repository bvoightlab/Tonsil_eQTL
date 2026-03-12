library(EnsDb.Hsapiens.v86)
library(RColorBrewer)
library(metafor)
library(ggforestplot)
library(tidyverse)
library(ggforce)

args <- c("Mpaired", "Fpaired", "NaiveB_GCB")
args <- c("4paired", "5paired", "NaiveB_GCB")
args <- c("Mpaired", "Fpaired", "NaiveT_TFH")
args <- c("4paired", "5paired", "NaiveT_TFH")
args <- c("Mpaired", "Fpaired", "NaiveT_NaiveB")
args <- c("4paired", "5paired", "NaiveT_NaiveB")
args <- c("Mpaired", "Fpaired", "GCB_TFH")
args <- c("4paired", "5paired", "GCB_TFH")

analysis1 <- args[1]
analysis2 <- args[2]
pair <- args[3]


#fancy the analysis name for printing
analysisp <- gsub("paired", " only", analysis1)
analysisp <- gsub("4", "<4", analysisp)
analysisp <- gsub("M", "Male", analysisp)


analysis2p <- gsub("paired", " only", analysis2)
analysis2p <- gsub("5", ">5", analysis2p)
analysis2p <- gsub("F", "Female", analysis2p)

#make the printed plotname pretty
pairp <- str_replace_all(pair, "_", " vs ")
pairp <- str_replace_all(pairp, "NaiveT", "Naive T")
pairp <- str_replace_all(pairp, "NaiveB", "Naive B")
pairp <- str_replace_all(pairp, "TFH", "Tfh")
pairp <- str_replace_all(pairp, "GCB", "GC B")


afile <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair, "_DESEQ_results_", analysis1, ".txt", sep='')
a2file <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair, "_DESEQ_results_", analysis2, ".txt", sep='')

a1 <- read.delim(afile, header = TRUE)
a1 <- rownames_to_column(a1, "ENSG")
base <- read.delim(a2file, header = TRUE)
base <- rownames_to_column(base, "ENSG")

#combine into one table:

merged <- inner_join(base, select(a1, ENSG, log2FoldChange, lfcSE, padj), c("ENSG" = "ENSG"))
#.x is a2 and .y is a1

#filter so only have genes sig in one subset or the other
sig <- filter(merged, (padj.x < 0.0025 & (log2FoldChange.x >= 1.5 | log2FoldChange.x <= -1.5)) | (padj.y < 0.0025 & (log2FoldChange.y >= 1.5 | log2FoldChange.y <= -1.5)))

#add gene names
sig <- separate(sig, col = "ENSG", into = c("ENSG", "transcript?"), sep = "[.]")
genes <- as.data.frame(genes(EnsDb.Hsapiens.v86))
sig <- left_join(sig, select(genes, gene_id, gene_name, gene_biotype, seqnames, start, end), c("ENSG" = "gene_id") )

#rearrange things to be sensible
sig <- unite(sig, ENSG, c("ENSG", "transcript?"), sep = ".")
sig <- relocate(sig, c("gene_name", "gene_biotype", "seqnames", "start", "end"), .after = "ENSG")
sig <- rename(sig, "chr" = "seqnames")

sig$Q <- NA
sig$Qp <- NA
sig$I2 <- NA

#step through dataset and run meta on each gene result
for (i in 1:nrow(sig)){
  #i <- 1
  
  #put data from this row into a table
  formeta <- data.frame(dataset = c("all_data", "subset"),
                        lfc = c(sig[[i, "log2FoldChange.x"]], sig[[i, "log2FoldChange.y"]]),
                        SE = c(sig[[i, "lfcSE.x"]], sig[[i, "lfcSE.y"]]))
  
  #run meta
  results <- rma(yi=lfc, sei=SE, method="FE", data=formeta)
  
  #add to table
  sig[[i, "Q"]] <- results$QE
  sig[[i, "Qp"]] <- results$QEp
  sig[[i, "I2"]] <- results$I2
  
  
}
#filter for signficiant heterogenetity
SE_interest <- filter(sig, Qp < (0.05/nrow(sig)))
#calculate difference in lfc and remove X
SE_interest$lfc_delta <- abs(SE_interest$log2FoldChange.y - SE_interest$log2FoldChange.x)
SE_interest <- filter(SE_interest, chr != "X")
SE_interest <- filter(SE_interest, lfc_delta >= .5)
#tag with orientation & split
SE_interest$orientation <- ifelse(SE_interest$log2FoldChange.x < SE_interest$log2FoldChange.y, paste(analysis2p, "<", analysisp,  sep = " "), paste(analysis2p, ">", analysisp, sep = " "))
dfx <- select(SE_interest, gene_name, chr, log2FoldChange.x, lfcSE.x, orientation)
dfy <- select(SE_interest, gene_name, chr, log2FoldChange.y, lfcSE.y, orientation)
#rename
analysisq <- gsub("paired", "_only", analysis2)
dfx$Subset <- analysis2p
dfx <-rename(dfx, log2FoldChange = log2FoldChange.x)
dfx <-rename(dfx, lfcSE = lfcSE.x)

analysist <- gsub("paired", "_only", analysis1)
dfy$Subset <- analysisp
dfy <-rename(dfy, log2FoldChange = log2FoldChange.y)
dfy <-rename(dfy, lfcSE = lfcSE.y)
#merge & reorder columns
df_combo <- rbind(dfx, dfy)
df_combo <- df_combo %>% 
  right_join(select(SE_interest, gene_name, orientation), by = "gene_name") %>% 
  mutate(df_combo, group = factor(.data$orientation, levels = unique(.data$orientation)))
df_combo <- arrange(df_combo, log2FoldChange)

png(paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair, "_DESEQ_results_", analysis1, "_vs_", analysis2, "_lfc_Bfforest2.png", sep=''), width=5, height=10, unit="in", res=300);
ggforestplot::forestplot(
  df = df_combo,
  name = gene_name,
  estimate = log2FoldChange,
  se = lfcSE,
  colour = Subset,
  xlab = "Log2FoldChange",
  title = paste(pairp, " DE genes", sep = ""), subtitle = paste("Subsets: ", analysisp, " vs ",  analysis2p, sep = "")) +
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
dev.off()



#output file
sig$interest <- ifelse(sig$Qp < (0.05/nrow(sig)), TRUE, FALSE)
sig$lfc_delta <- abs(sig$log2FoldChange.y - sig$log2FoldChange.x)

#rename lfc columns so distinguishable
analysist <- gsub("paired", "_only", analysis1)
newname <- paste(analysist, "_lfc", sep = "")
sig <-rename(sig, (!!newname) := log2FoldChange.y)
analysisq <- gsub("paired", "_only", analysis2)
newname <- paste(analysisq, "_lfc", sep = "")
sig <-rename(sig, (!!newname) := log2FoldChange.x)

newname <- paste(analysist, "_lfcSE", sep = "")
sig <-rename(sig, (!!newname) := lfcSE.y)
newname <- paste(analysisq, "_lfcSE", sep = "")
sig <-rename(sig, (!!newname) := lfcSE.x)

newname <- paste(analysist, "_padj", sep = "")
sig <-rename(sig, (!!newname) := padj.y)
newname <- paste(analysisq, "_padj", sep = "")
sig <-rename(sig, (!!newname) := padj.x)

sig <- select(sig, !c(baseMean, pvalue, stat))

signals_out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair, "_DESEQ_results_", analysis1, "_vs_", analysis2, "_lfc_het.txt", sep='')
write.table(sig, file= signals_out, sep="\t", row.names=FALSE, quote=FALSE)


sig <- filter(sig, interest == "TRUE")

signals_out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair, "_DESEQ_results_", analysis1, "_vs_", analysis2, "_lfc_het_sig.txt", sep='')
write.table(sig, file= signals_out, sep="\t", row.names=FALSE, quote=FALSE)

sig <- filter(sig, chr != "X")
sig <- filter(sig, lfc_delta >= .5)

signals_out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair, "_DESEQ_results_", analysis1, "_vs_", analysis2, "_lfc_het_sig_noX_lfc5.txt", sep='')
write.table(sig, file= signals_out, sep="\t", row.names=FALSE, quote=FALSE)

  