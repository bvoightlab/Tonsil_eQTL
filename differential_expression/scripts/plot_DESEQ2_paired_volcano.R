library(EnsDb.Hsapiens.v86)
library(tidyverse)


args = commandArgs(TRUE)
#args <- c("NaiveT_NaiveB","paired")
celltype <- args[1]
analysis <- args[2]


res <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", celltype, "_DESEQ_results_", analysis, "_sig.txt", sep='')

results <- read.delim(res, header = TRUE)

results <- rownames_to_column(results)

#add gene names?

DE <- separate(results, col = "rowname", into = c("ENSG", "transcript?"), sep = "[.]")

genes <- as.data.frame(genes(EnsDb.Hsapiens.v86))

DE_named <- left_join(DE, select(genes, gene_id, gene_name, gene_biotype, seqnames, start, end), c("ENSG" = "gene_id") )

#rearrange things to be pretty
DE_named <- unite(DE_named, ENSG, c("ENSG", "transcript?"), sep = ".")
DE_named <- relocate(DE_named, c("gene_name", "gene_biotype", "seqnames", "start", "end"), .after = "ENSG")
DE_named <- rename(DE_named, "chr" = "seqnames")

#remove genes with FC <|1.5|
DE_final <- DE_named %>% filter(log2FoldChange <= -1.5 | log2FoldChange >= 1.5)

#make the printed plotname pretty
printcelltype <- str_replace_all(celltype, "_", " vs ")
printcelltype <- str_replace_all(printcelltype, "NaiveT", "naive T")
printcelltype <- str_replace_all(printcelltype, "NaiveB", "naive B")
printcelltype <- str_replace_all(printcelltype, "TFH", "Tfh")
printcelltype <- str_replace_all(printcelltype, "GCB", "GC B")

#volcano plot might be better?
png(paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", celltype, "_DESEQ_results_", analysis, "_volcano_bonf.png", sep=''), width=7, height=7, unit="in", res=300);
ggplot(DE_final, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point() +
  ggtitle(paste(printcelltype, " Differential Expression", sep=''))
dev.off()


