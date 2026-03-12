library(tidyverse)
library(ggtext)
library(ggrepel)
library(cowplot)


#Read in cell types to compare
#args <- c("TFH", "GCB", "ENSG00000205189.12", "ZBTB10", "pAsthma", "80485619", "80526265")

args <- commandArgs(trailingOnly = TRUE)




cell1 <- args[1]
cell2 <- args[2]
ENSG <- args[3]
gene <- args[4]
trait <- args[5]
gene_start <- as.numeric(args[6])
gene_stop <- as.numeric(args[7])


file_1 <- paste0("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/colocs/", ENSG, "_", cell1, "_allpairs.txt")
file_2 <- paste0("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/colocs/", ENSG, "_", cell2, "_allpairs.txt")

#check if files exist, run grep if not
if (file.exists(file_1)){
  eqtl_1 <- read.table(file_1, header=FALSE, fill=TRUE)

} else {
  system(paste0("zgrep ", ENSG, " /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/", cell1, "_all.cis_qtl_pairs.txt.gz > /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/colocs/", ENSG, "_", cell1, "_allpairs.txt"))
  Sys.sleep(10)
  eqtl_1 <- read.table(file_1, header=FALSE, fill=TRUE)
  
}

if (file.exists(file_2)){
  eqtl_2 <- read.table(file_2, header=FALSE, fill=TRUE)
  
} else {
  system(paste0("zgrep ", ENSG, " /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/", cell2, "_all.cis_qtl_pairs.txt.gz > /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/colocs/", ENSG, "_", cell2, "_allpairs.txt"))
  Sys.sleep(10)
  eqtl_2 <- read.table(file_2, header=FALSE, fill=TRUE)
  
}

colnames(eqtl_1) <- c("phenotype_id", "variant_id", "tss_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")
colnames(eqtl_2) <- c("phenotype_id", "variant_id", "tss_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")


#merge 
joint <- full_join(eqtl_1, eqtl_2, join_by(phenotype_id, variant_id))

#separate variant ID to get pos and chr
joint <- separate(joint, variant_id, c("chr", "pos", "a1", "a2"), remove = FALSE)
joint$pos <- as.numeric(joint$pos)

#extract chr
chr <- joint$chr[1]

#extract lead SNPs and gene positions to figure out interval
label <- joint[which.min(joint$pval_nominal.x),]
othertype <- joint[which.min(joint$pval_nominal.y),]
key_pos <- c(gene_start, gene_stop, label$pos, othertype$pos)
key_min <- min(key_pos)
key_min <- min(key_pos)
key_max <- max(key_pos)
#add buffer based on size of region
buf <- (key_max - key_min)/36
left <- as.integer(key_min - buf)
right <- as.integer(key_max + buf)

#get trait sumstats too

traitfile <- paste0("/project/voightlab_01/lorenzk/Romberg_eQTL/coloc/sumstats/", trait, ".quant.processed.bed.gz")
traitdf <- read.table(traitfile, header=TRUE)
#reduce to chr and window asked for
traitdf$chromStart <- as.numeric(traitdf$chromStart)
traitdf2 <- filter(traitdf, chrom == chr & chromStart > left & chromStart < right)

#Print small table for trait
#out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", trait, ".bedgraph", sep='')
#traitbedgraph <- select(traitdf2, c(chrom, chromStart, chromEnd, pval))
#traitbedgraph$pval <- -log10(traitbedgraph$pval)
#write.table(traitbedgraph, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)


#merge into big table
joint <- left_join(joint, traitdf2, join_by(pos == chromStart, chr == chrom))

#pull label as lowest pval in celltype1
label <- joint[which.min(joint$pval_nominal.x),]
chrpos <- label$variant_id
SNP_print <- str_replace_all(chrpos, ":", "_")

#get LD info
LD_file <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/plink/",SNP_print,".ld", sep='')
#run plink
#commending this out as running locally and this call is already done. 
system(paste("plink --bfile /project/voightlab_01/lorenzk/Romberg_eQTL/ImputationII_withX/plink/Romberg --r2 --ld-window-r2 0 --ld-window 100000 --ld-window-kb 750 --ld-snp ", chrpos, " --out /project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/plink/", SNP_print,  sep=''))
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

#reduce entire matrix to window requested
joint <- filter(joint, chr == chr & pos > left & pos < right)

##step through label and R2 levels and output bedgraph for each, for trait and both eQTLs
############################################################
#bedgraphs for label
#filename
out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", trait, "_label.bedgraph", sep='')
#pull needed columns
LDlabel <- select(label, c(chr, pos, pval))
#make end position column
LDlabel <- add_column(LDlabel,end = LDlabel$pos+1, .after = "pos")
#negative log10 of the pval
LDlabel$pval <- -log10(LDlabel$pval)
#write table
write.table(LDlabel, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_label.bedgraph", sep='')
LDlabel <- select(label, c(chr, pos, pval_nominal.x))
LDlabel <- add_column(LDlabel,end = LDlabel$pos+1, .after = "pos")
LDlabel$pval_nominal.x <- -log10(LDlabel$pval_nominal.x)
write.table(LDlabel, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell2, "_label.bedgraph", sep='')
LDlabel <- select(label, c(chr, pos, pval_nominal.y))
LDlabel <- add_column(LDlabel,end = LDlabel$pos+1, .after = "pos")
LDlabel$pval_nominal.y <- -log10(LDlabel$pval_nominal.y)
write.table(LDlabel, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#for dark blue (R2 0-0.2)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", trait, "_LD0.bedgraph", sep='')
#pull color subset needed
LD0 <- filter(joint, color == 0 )
#select columns of interest
LD0 <- select(LD0, c(chr, pos, pval))
#make end position column
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
#remove rows with no pvalue
LD0 <- drop_na(LD0, pval)
#take negative log10 of the pvalue
LD0$pval <- -log10(LD0$pval)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_LD0.bedgraph", sep='')
LD0 <- filter(joint, color == 0 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.x))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.x)
LD0$pval_nominal.x <- -log10(LD0$pval_nominal.x)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell2, "_LD0.bedgraph", sep='')
LD0 <- filter(joint, color == 0 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.y))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.y)
LD0$pval_nominal.y <- -log10(LD0$pval_nominal.y)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#for light blue (R2 0.2-0.4)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", trait, "_LD.2.bedgraph", sep='')
LD0 <- filter(joint, color == 0.2 )
LD0 <- select(LD0, c(chr, pos, pval))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval)
LD0$pval <- -log10(LD0$pval)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_LD.2.bedgraph", sep='')
LD0 <- filter(joint, color == 0.2 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.x))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.x)
LD0$pval_nominal.x <- -log10(LD0$pval_nominal.x)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell2, "_LD.2.bedgraph", sep='')
LD0 <- filter(joint, color == 0.2 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.y))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.y)
LD0$pval_nominal.y <- -log10(LD0$pval_nominal.y)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#for green (R2 0.4-0.6)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", trait, "_LD.4.bedgraph", sep='')
LD0 <- filter(joint, color == 0.4 )
LD0 <- select(LD0, c(chr, pos, pval))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval)
LD0$pval <- -log10(LD0$pval)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_LD.4.bedgraph", sep='')
LD0 <- filter(joint, color == 0.4 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.x))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.x)
LD0$pval_nominal.x <- -log10(LD0$pval_nominal.x)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell2, "_LD.4.bedgraph", sep='')
LD0 <- filter(joint, color == 0.4 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.y))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.y)
LD0$pval_nominal.y <- -log10(LD0$pval_nominal.y)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#for yellow (R2 0.6-0.8)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", trait, "_LD.6.bedgraph", sep='')
LD0 <- filter(joint, color == 0.6 )
LD0 <- select(LD0, c(chr, pos, pval))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval)
LD0$pval <- -log10(LD0$pval)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_LD.6.bedgraph", sep='')
LD0 <- filter(joint, color == 0.6 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.x))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.x)
LD0$pval_nominal.x <- -log10(LD0$pval_nominal.x)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell2, "_LD.6.bedgraph", sep='')
LD0 <- filter(joint, color == 0.6 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.y))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.y)
LD0$pval_nominal.y <- -log10(LD0$pval_nominal.y)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#for red (R2 0.8-1)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", trait, "_LD.8.bedgraph", sep='')
LD0 <- filter(joint, color == 0.8 )
LD0 <- select(LD0, c(chr, pos, pval))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval)
LD0$pval <- -log10(LD0$pval)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_LD.8.bedgraph", sep='')
LD0 <- filter(joint, color == 0.8 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.x))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.x)
LD0$pval_nominal.x <- -log10(LD0$pval_nominal.x)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

out <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell2, "_LD.8.bedgraph", sep='')
LD0 <- filter(joint, color == 0.8 )
LD0 <- select(LD0, c(chr, pos, pval_nominal.y))
LD0 <- add_column(LD0,end = LD0$pos+1, .after = "pos")
LD0 <- drop_na(LD0, pval_nominal.y)
LD0$pval_nominal.y <- -log10(LD0$pval_nominal.y)
write.table(LD0, file= out, sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)


##make stacked locuszoom plot for reference

lz1 <- ggplot(joint, aes(x = pos, y = -log10(pval_nominal.x))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}) +
  xlab("")+ ylab(paste(cell1, " ", gene, " -log10(P)", sep='')) + 
  theme_classic() + 
  guides(fill = guide_legend(title = "R2", override.aes = list(size = 6, shape = 22)))+
  theme(plot.title = element_text(hjust = 1), legend.title = element_text(hjust = 1), legend.justification = "top", legend.background = element_rect(fill = "NA"), axis.text.x=element_blank())
lz2 <- ggplot(joint, aes(x = pos, y = -log10(pval_nominal.y))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}) +
  xlab("")+ ylab(paste(cell2, " ", gene, " -log10(P)", sep='')) + 
  theme_classic() + 
  theme(axis.text.x=element_blank())
lz3 <- ggplot(joint, aes(x = pos, y = -log10(pval))) +
  geom_point(aes(fill = color), alpha = 0.8, shape = 21, size = 2 ) +
  scale_fill_manual(values = c('red','orange','darkgreen','skyblue','blue4'), guide = "none") +
  geom_point(data = label, fill = "purple", shape = 23, size = 3) +
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}) +
  xlab(paste(chr,  " (Mb)", sep=''))+ ylab(paste(trait, " -log10(P)", sep='')) + 
  theme_classic()

zoomstack <- cowplot::plot_grid(lz1, lz2, lz3, align = "v", nrow = 3, rel_heights = c(0.95, 0.95, 1))
#outfile name
out_file <- paste("/project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_", cell2, "_", trait, "_stackedlocuszoom.png", sep='')

save_plot(out_file, zoomstack , base_height = 6, base_width = 6)

##################
#modify .ini template and call pygenometracks

#replace 'celltype1' in .ini, creating named .ini file
system(paste0("sed \"s/celltype1/", cell1, "/\" /project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/template_lzloops.ini > /project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini"))
#replace 'celltype2' 
system(paste0("sed -i \"s/celltype2/",  cell2, "/\" /project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini"))
#replace 'trait' in filenames
system(paste0("sed -i \"s/trait/", trait, "/\" /project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini"))
#replace 'gene' in filenames
system(paste0("sed -i \"s/GENE/", gene, "/\" /project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini"))

#calculate max -log10(pval) for eQTL and gwas
temp <- drop_na(joint, pval)
traitmax <- round(max(-log10(temp$pval)) +1 )
cell1max <- round(max(-log10(temp$pval_nominal.x)) +1 )
cell2max <- round(max(-log10(temp$pval_nominal.y)) +1 )
eQTLmax <- max(cell1max, cell2max)

#replace 'eQTL_value' in filenames
system(paste0("sed -i \"s/eQTL_value/", eQTLmax, "/\" /project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini"))
#replace 'trait_value' in filenames
system(paste0("sed -i \"s/TRAIT_value/", traitmax, "/\" /project/voightlab_01/lorenzk/Romberg_eQTL/Allergy_coloc/plots/", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini"))

print(paste0("pyGenomeTracks --dpi 300 --fontSize 16 --tracks ", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini -o ", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.png --region ", chr, ":", left, "-", right))

system(paste0("pyGenomeTracks --dpi 300 --fontSize 16 --tracks ", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.ini -o ", gene, "_", cell1, "_", cell2, "_", trait, "_lzloops.png --region ", chr, ":", left, "-", right))


