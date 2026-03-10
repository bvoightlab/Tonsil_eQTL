library(tidyverse)
library(UpSetR)


pair1 <- "NaiveB_GCB"
pair2 <- "NaiveT_TFH"

args = commandArgs(TRUE)
#args <- c("paired")
analysis <- args[1]

#analysis <- "paired"
#

Bcellsfile <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair1, "_DESEQ_results_", analysis, "_sig_named.txt", sep='')
Tcellsfile <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/", pair2, "_DESEQ_results_", analysis, "_sig_named.txt", sep='')

Bcells1 <- read.delim(Bcellsfile, header = TRUE)
Tcells1 <- read.delim(Tcellsfile, header = TRUE)

#straight compare:

l <- list(
  Bcells = (Bcells1[[1]]),
  Tcells = (Tcells1[[1]])
)

png(paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/Upset_", pair1, "_vs_", pair2, "_DESEQ_", analysis, ".png", sep=''), width=10, height=5.5, unit="in", res=300);
upset(fromList(l), 
      order.by = "freq",
      nsets = 2, nintersects = NA,
      #text.scale = 1.85,
      text.scale = c(1.85, 1.85, 1.85, 1.85, 1.85, 1.85), 
      #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
      point.size = 3,
      line.size = 1.5,
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total DE genes",
      mainbar.y.label = "DE genes in common"
);
dev.off();

#direction of effect

#B cell data is inverted, so flip the sign here
bcell_pos1 <- filter(Bcells1, log2FoldChange < 0)
bcell_neg1 <- filter(Bcells1, log2FoldChange > 0)

tcell_neg1 <- filter(Tcells1, log2FoldChange < 0)
tcell_pos1 <- filter(Tcells1, log2FoldChange > 0)

#title so pos = increased expression in activated cell type
m <- list(
  "GC B down" = (bcell_neg1[[1]]),
  "GC B up" = (bcell_pos1[[1]]),
  "Tfh down" = (tcell_neg1[[1]]),
  "Tfh up" = (tcell_pos1[[1]])
)

png(paste("/project/voightlab_01/lorenzk/Romberg_eQTL/DESeq2/Upset_", pair1, "_vs_", pair2, "_DESEQ_", analysis, "_signed.png", sep=''), width=10, height=5.5, unit="in", res=300);
upset(fromList(m), 
      queries = list(
        list(query = intersects, params = list("GC B down", "Tfh down"), color = "orange", active = T),
        list(query = intersects, params = list("GC B up", "Tfh up"), color = "orange", active = T)
      ),
      order.by = "freq",
      nsets = 4,  nintersects = NA,
      #text.scale = 1.85,
      text.scale = c(1.85, 1.85, 1.85, 1.4, 1.85, 1.85), 
      #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
      point.size = 3,
      line.size = 1.5,
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total DE genes",
      mainbar.y.label = "DE genes in common"
);
dev.off();
