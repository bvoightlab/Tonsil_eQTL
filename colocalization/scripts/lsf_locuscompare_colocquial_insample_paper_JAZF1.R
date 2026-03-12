library(tidyverse)
library(ggtext)
library(ggrepel)
library(cowplot)

#get inputs
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) != 2) {
#  stop("two arguments are required, coloc input filepath, plink dir", call.=FALSE)
#}
#input <- args[1]
#plink <- args[2]

input <- "rs4722758/JAZF1_ENSG00000153814.13_NaiveB_pAsthma_coloc_input_data.txt"
plink <- "/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plink"
setwd("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/pAsthma")
df <- read.table(input, header = TRUE)

#parse input to pull out needed things
stuff <- str_split_1(input, "[/_]")
snp <- stuff[1]
gene <- stuff[2]
celltype <- "naive B"
trait <- stuff[length(stuff) -3]
chr <- paste0("chr", df$chrom[1])

#make chr:pos:a1:a2 column
df <- df %>% unite("chrpos", c(chrom_b38, chromEnd_b38, a1_trait, a2_trait), sep = ":" , remove = FALSE)
#label is lowest SNP in 
label <- df[which.min(df$pvalue_eQTL),]
chrpos <- label$chrpos
chrpos_print <- str_replace_all(chrpos, ":", "_")

#get LD info
LD_file <- paste(plink, "/",chrpos_print,".ld", sep='')
#run plink
#system(paste("plink --bfile /project/voightlab_01/lorenzk/Romberg_eQTL/ImputationII_withX/plink/Romberg --r2 --ld-window-r2 0 --ld-window 100000 --ld-window-kb 750 --ld-snp ", chrpos, " --out ", plink, "/", chrpos_print,  sep=''))
#print(paste("plink --bfile /project/voightlab_01/lorenzk/Romberg_eQTL/ImputationII_withX/plink/Romberg --r2 --ld-window-r2 0 --ld-window 100000 --ld-window-kb 750 --ld-snp ", chrpos, " --out ", plink, "/", chrpos_print,  sep=''))

#if LD file exists, add it
if (file.exists(LD_file)) {
  LD <- read.table(LD_file, header = TRUE)
  
  #merge that in
  df <- left_join(df, select(LD, c("SNP_B", "R2")), c("chrpos" = "SNP_B"))
  
  #add color
  df$color <- as.character(cut(df$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('0','0.2','0.4','0.6','0.8'), include.lowest=TRUE))
  df$color <- factor(df$color, levels = c(0.8,0.6,0.4,0.2,0))

  
} else {
  #otherwise add column of NAs
  df$color <- NA
}



##make locuscompare plot


plot1 <- ggplot(df, aes(x = -log10(pvalue_eQTL), y = -log10(pval))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  #geom_text_repel(data = label, label = snp) +
  xlab(paste(celltype, " ", gene, " eQTL\n-log10(P)", sep=''))+ ylab(paste(trait, " GWAS\n-log10(P)", sep='')) + 
  guides(fill = guide_legend(title=bquote(R^2), override.aes = list(size = 6, shape = 22)))+
  ggtitle(chrpos) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 1), legend.title = element_text(hjust = 1), legend.position=c(.9,.4), legend.background = element_rect(fill = "NA"))

lz1 <- ggplot(df, aes(x = chromStart, y = -log10(pvalue_eQTL))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  #geom_text_repel(data = label, label = snp) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}) +
  xlab("")+ ylab(paste(celltype, " ", gene, " eQTL\n-log10(P)", sep='')) + 
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.title.y=(element_text(size = rel(1))))
lz2 <- ggplot(df, aes(x = chromStart, y = -log10(pval))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2 ) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  #geom_text_repel(data = label, label = snp) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}) +
  xlab(paste(chr,  " (Mb)", sep=''))+ ylab(paste(trait, " GWAS\n-log10(P)", sep='')) + 
  theme_classic() +
  theme(axis.title.y=(element_text(size = rel(1))))

bothzooms <- cowplot::plot_grid(lz1, lz2, align = "v", nrow = 2, rel_heights = c(.9, 1))
allthree <- cowplot::plot_grid(plot1, bothzooms)

#outfile name
out_file <- paste("locuscompare/paper_",gene,"_", celltype, "_", trait, "_",snp, "_locuscompare.png", sep='')

save_plot(out_file, allthree )
