library(tidyverse)
library(grid)

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/Atopic_Asthma/Atopic_Asthma_overlapping_intervals.bed", "/project/voightlab_01/lorenzk/Romberg_eQTL/coloc/sumstats/Atopic_Asthma.quant.processed.bed","/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/Atopic_Asthma/region_plots/")
data <- args[1]
trait_int <- args[2]
outdir <- args[3]

#leadsfile <- str_replace(data, "intervals.bed", "leads_hg38.bed")

df <- read.table(data, header = FALSE)
#leads <- read.table(leadsfile, header = FALSE)
trait_df <- read.table(trait_int, header = TRUE)

#make df to hold leads
output <- filter(trait_df, sample_size == 1 )
#make colocquial output file
colocquial <- select(df, V1, V2, V3)
colnames(colocquial) <- c("chr", "start", "stop")
colocquial$rsID <- "NA"
colocquial$pos <- "NA"
colocquial <- colocquial %>% 
  mutate(across('chr', str_replace, 'chr', '')) %>% 
  relocate(rsID, chr, pos, start, stop)



#region = 1
for (region in 1:nrow(df)) {
  #pull values we need
  start <- df[[region, "V2"]]
  end <- df[[region, "V3"]]
  chr <- df[[region, "V1"]]
  

  #add 50kb flanks
  startf = start - 50000
  endf = end + 50000
  
    #get size 
  size = end - start + 1
  
  #pull SNPs in interval
  fullint <- filter(trait_df, chrom == chr & chromStart >= startf & chromStart <= endf)
  
  #pull minp SNP and add to output dataframe
  labels <- fullint[which.min(fullint$pval),]
  output <- rbind(output, labels)
  
  #also add to df in lead info columns
  colocquial[region, "rsID"] = labels$name
  colocquial[region, "pos"] = labels$chromStart
  
  
}


out <- str_replace(data, "overlapping_intervals.bed", "colocquial_input.txt")
write.table(colocquial, file= out, sep="\t", row.names=FALSE, col.names = TRUE, quote=FALSE)

out <- str_replace(data, "overlapping_intervals.bed", "colocquial_leads_info.txt")
write.table(output, file= out, sep="\t", row.names=FALSE, col.names = TRUE, quote=FALSE)






