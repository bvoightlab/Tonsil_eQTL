library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("two arguments are required, scoreboardfile, celltype", call.=FALSE)
}
scoreboardfile <- args[1]
celltype <- args[2]



#read in files

#scoreboardfile <- "/project/voightlab_01/lorenzk/Romberg_eQTL/TFH_egenes_sig_scoreboard.txt"
#celltype <- "T"

base <- read_delim(scoreboardfile)


#WGCNA modules and lfcs
DE_out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/All_DESEQ_WGCNAjaccard_eQTL_results.txt", sep='')

DF <-read.delim(DE_out, header = TRUE)
#change module names to yes/no based on modules of interest (B, C, E, F, H, N, S)
DF <- mutate(DF, NB_module = if_else(NB_module %in% c("B", "C", "E", "F", "H", "N", "S"), "yes", "no"))
DF <- mutate(DF, NT_module = if_else(NT_module %in% c("B", "C", "E", "F", "H", "N", "S"), "yes", "no"))
DF <- mutate(DF, GCB_module = if_else(GCB_module %in% c("B", "C", "E", "F", "H", "N", "S"), "yes", "no"))
DF <- mutate(DF, TFH_module = if_else(TFH_module %in% c("B", "C", "E", "F", "H", "N", "S"), "yes", "no"))

#add overall (yes if yes in any celltype)
DF <- mutate(DF, all_module = if_else((NB_module =="yes" | NT_module =="yes" | GCB_module =="yes" | TFH_module =="yes"), "yes", "no"))


#panther gene lists
PanImm <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/reference_files/panther_gene_lists/panther_immunesystem_process.txt", sep='')
ImmSys <-read_delim(PanImm, col_names = FALSE, delim = ";")

#GSEA gene lists
GSEA_T1 <- "/project/voightlab_01/lorenzk/Romberg_eQTL/GSEA_gene_sets/GSEA_allgenes_TFH.txt"
GSEA_B1 <- "/project/voightlab_01/lorenzk/Romberg_eQTL/GSEA_gene_sets/GSEA_allgenes_GCB.txt"

GSEA_T <-read_delim(GSEA_T1, col_names = FALSE, delim = ";")
GSEA_B <-read_delim(GSEA_B1, col_names = FALSE, delim = ";")


#add WGCNA modules and lfcs
#right now set to T cell type, need to split later for both
if (celltype == "T"){
  scoreboard <- left_join(base, select(DF, c("ENSG", "NTT_lfc", "NT_module", "TFH_module", "all_module" )), c("phenotype_id" = "ENSG"))
}else {
  scoreboard <- left_join(base, select(DF, c("ENSG", "NBG_lfc", "NB_module", "GCB_module", "all_module" )), c("phenotype_id" = "ENSG"))
  
}

#add panther immune system designation
scoreboard <- mutate(scoreboard, ImmuneSys_Panther = if_else(gene_name %in% ImmSys$X2, "yes", "no"))

#add appropriate GSEA overlaps
if (celltype == "T"){
  scoreboard <- mutate(scoreboard, GSEA_TFH = if_else(gene_name %in% GSEA_T$X1, "yes", "no"))
}else {
  scoreboard <- mutate(scoreboard, GSEA_GCB = if_else(gene_name %in% GSEA_B$X1, "yes", "no"))
  
}

#output those

scoreboard_out <- gsub(".txt", "_DE.txt", scoreboardfile)
write.table(scoreboard, file= scoreboard_out, sep="\t", row.names=FALSE, quote=FALSE)

