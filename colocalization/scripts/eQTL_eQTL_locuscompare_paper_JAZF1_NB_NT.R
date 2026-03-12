library(tidyverse)
library(ggtext)
library(ggrepel)
library(cowplot)


#Read in cell types to compare
args <- commandArgs(trailingOnly = TRUE)
#args <- c("GCB", "TFH", "ENSG00000165895", "ARHGAP42", "rs585183")
#args <- c("GCB", "NaiveT", "ENSG00000131323", "TRAF3", "chr14:102768675:G:A")

#args <- c("NaiveB", "NaiveT", "ENSG00000153814", "JAZF1")

cell1 <- args[1]
cell2 <- args[2]
ENSG <- args[3]
gene <- args[4]

cell1print <- "naive B"
cell2print <- "naive T"


file_1 <- paste0("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/colocs/", ENSG, "_", cell1, "_allpairs.txt")
file_2 <- paste0("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/colocs/", ENSG, "_", cell2, "_allpairs.txt")

#read in
eqtl_1 <- read.table(file_1, header=FALSE, fill=TRUE)
eqtl_2 <- read.table(file_2, header=FALSE, fill=TRUE)
  
colnames(eqtl_1) <- c("phenotype_id", "variant_id", "tss_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")
colnames(eqtl_2) <- c("phenotype_id", "variant_id", "tss_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")

#merge 
joint <- inner_join(eqtl_1, eqtl_2, join_by(phenotype_id, variant_id))

#separate variant ID to get pos and chr
joint <- separate(joint, variant_id, c("chr", "pos", "a1", "a2"), remove = FALSE)
joint$pos <- as.numeric(joint$pos)

#extract chr
chr <- joint$chr

#pull label as lowest pval in celltype1
label <- joint[which.min(joint$pval_nominal.x),]
chrpos <- label$variant_id
SNP_print <- str_replace_all(chrpos, ":", "_")

#get LD info
LD_file <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/plink/",SNP_print,".ld", sep='')
#run plink
#system(paste("plink --bfile /project/voightlab_01/lorenzk/Romberg_eQTL/ImputationII_withX/plink/Romberg --r2 --ld-window-r2 0 --ld-window 100000 --ld-window-kb 750 --ld-snp ", chrpos, " --out /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/plink/", SNP_print,  sep=''))
#print(paste("plink --bfile /project/voightlab_01/lorenzk/Romberg_eQTL/ImputationII_withX/plink/Romberg --r2 --ld-window-r2 0 --ld-window 100000 --ld-window-kb 750 --ld-snp ", chrpos, " --out /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/plink/", SNP_print,  sep=''))

#if LD file exists, add it
if (file.exists(LD_file)) {
  LD <- read.table(LD_file, header = TRUE)
  
  #merge that in
  joint <- left_join(joint, select(LD, c("SNP_B", "R2")), c("variant_id" = "SNP_B"))
  
  #add color
  joint$color <- as.character(cut(joint$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('0','0.2','0.4','0.6','0.8'), include.lowest=TRUE))
  joint$color <- factor(joint$color, levels = c(0.8,0.6,0.4,0.2,0))
  
  
} else {
  #otherwise add column of NAs
  joint$color <- NA
}

#reduce to window of eQTL x GWAS plot

joint <- filter(joint, pos > 28000137 & pos < 28298646)


##make locuscompare plot


plot1 <- ggplot(joint, aes(x = -log10(pval_nominal.x), y = -log10(pval_nominal.y))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  #geom_text_repel(data = label, label = snp) +
  xlab(paste(cell1print, " ", gene, " eQTL\n-log10(P)", sep=''))+ ylab(paste(cell2print, " ", gene, " eQTL\n-log10(P)", sep='')) + 
  guides(fill = guide_legend(title=bquote(R^2), override.aes = list(size = 6, shape = 22)))+
  ggtitle(chrpos) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 1), legend.title = element_text(hjust = 1), legend.position=c(.9,.4), legend.background = element_rect(fill = "NA"))

lz1 <- ggplot(joint, aes(x = pos, y = -log10(pval_nominal.x))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}) +
  xlab("")+ ylab(paste(cell1print, " ", gene, " eQTL\n-log10(P)", sep='')) + 
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.title.y=(element_text(size = rel(1))))
lz2 <- ggplot(joint, aes(x = pos, y = -log10(pval_nominal.y))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2 ) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}) +
  xlab(paste(chr,  " (Mb)", sep=''))+ ylab(paste(cell2print, " ", gene, " eQTL\n-log10(P)", sep='')) + 
  theme_classic() +
  theme(axis.title.y=(element_text(size = rel(1))))

bothzooms <- cowplot::plot_grid(lz1, lz2, align = "v", nrow = 2, rel_heights = c(0.9, 1))
allthree <- cowplot::plot_grid(plot1, bothzooms)

#outfile name
out_file <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Paper_Figures/paper_", gene, "_", cell1, "_", cell2, "_locuscompare.png", sep='')

save_plot(out_file, allthree )


