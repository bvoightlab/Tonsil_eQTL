library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/project/voightlab_01/lorenzk/Romberg_eQTL/coloc/sumstats/Atopic_Asthma.quant.processed.bed")
afile <- args[1]


data <- read.table(afile, header = TRUE)

#remove chr from chrom column
data <- data %>% 
  mutate(across('chrom', str_replace, 'chr', ''))

#make ID column
data <- data %>% 
  unite("ID", c("chrom", "chromStart", "a1", "a2"), sep = ":", remove = FALSE)

#rename columns
data <- data %>% 
  rename(pos = chromStart ) %>% 
  rename(rsid = name ) 
  
#rearrange columns
data <- select(data, !chromEnd)  %>% 
  relocate(ID, rsid, chrom, pos )


genes_out <- str_replace(afile, ".quant.processed.bed", ".temp.txt")
write.table(data, file= genes_out, sep="\t", row.names=FALSE, quote=FALSE)

