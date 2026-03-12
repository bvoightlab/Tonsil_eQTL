library(UpSetR)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(cowplot)

NaiveB1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveB.sig_egenes.txt", header = TRUE)
GCB1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/GCB.sig_egenes.txt", header = TRUE)
NaiveT1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveT.sig_egenes.txt", header = TRUE)
TFH1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/TFH.sig_egenes.txt", header = TRUE)


NaiveB1_m <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveB_sex1.sig_egenes.txt", header = TRUE)
GCB1_m <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/GCB_sex1.sig_egenes.txt", header = TRUE)
NaiveT1_m <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveT_sex1.sig_egenes.txt", header = TRUE)
TFH1_m <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/TFH_sex1.sig_egenes.txt", header = TRUE)

NaiveB1_f <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveB_sex2.sig_egenes.txt", header = TRUE)
GCB1_f <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/GCB_sex2.sig_egenes.txt", header = TRUE)
NaiveT1_f <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveT_sex2.sig_egenes.txt", header = TRUE)
TFH1_f <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/TFH_sex2.sig_egenes.txt", header = TRUE)

NaiveB1_4y <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveB_lt4.sig_egenes.txt", header = TRUE)
GCB1_4y <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/GCB_lt4.sig_egenes.txt", header = TRUE)
NaiveT1_4y <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveT_lt4.sig_egenes.txt", header = TRUE)
TFH1_4y <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/TFH_lt4.sig_egenes.txt", header = TRUE)

NaiveB1_5o <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveB_gt4.sig_egenes.txt", header = TRUE)
GCB1_5o <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/GCB_gt4.sig_egenes.txt", header = TRUE)
NaiveT1_5o <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveT_gt4.sig_egenes.txt", header = TRUE)
TFH1_5o <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/TFH_gt4.sig_egenes.txt", header = TRUE)

#age comparisons
l <- list(
  "Naive B" = (NaiveB1[[1]]),
  "Naive B <4" = (NaiveB1_4y[[1]]),
  "Naive B >5" = (NaiveB1_5o[[1]])
)
m <- list(
  "GC B" = (GCB1[[1]]),
  "GC B <4" = (GCB1_4y[[1]]),
  "GC B >5" = (GCB1_5o[[1]])
)
n <- list(
  "Naive T" = (NaiveT1[[1]]),
  "Naive T <4" = (NaiveT1_4y[[1]]),
  "Naive T >5" = (NaiveT1_5o[[1]])
)
o <- list(
  "Tfh" = (TFH1[[1]]),
  "Tfh <4" = (TFH1_4y[[1]]),
  "Tfh >5" = (TFH1_5o[[1]])
)

NB <- upset(fromList(l), 
      queries = list(
        list(query = intersects, params = list("Naive B <4"), color = "orange", active = T),
        list(query = intersects, params = list("Naive B >5"), color = "orange", active = T)
      ),
      order.by = c("freq"),
      nsets = 3, nintersects = NA,
      text.scale = 1.95,
      point.size = 3,
      line.size = 1.5,
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
);
GC <- upset(fromList(m), 
      queries = list(
        list(query = intersects, params = list("GC B <4"), color = "orange", active = T),
        list(query = intersects, params = list("GC B >5"), color = "orange", active = T)
      ),
      order.by = c("freq"),
      nsets = 3, nintersects = NA,
      text.scale = 1.95,
      point.size = 3,
      line.size = 1.5,
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
);
NT <- upset(fromList(n), 
      queries = list(
        list(query = intersects, params = list("Naive T <4"), color = "orange", active = T),
        list(query = intersects, params = list("Naive T >5"), color = "orange", active = T)
      ),
      order.by = c("freq"),
      nsets = 3, nintersects = NA,
      text.scale = 1.95,
      point.size = 3,
      line.size = 1.5,
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
);
TF <- upset(fromList(o), 
      queries = list(
        list(query = intersects, params = list("Tfh <4"), color = "orange", active = T),
        list(query = intersects, params = list("Tfh >5"), color = "orange", active = T),
        list(query = intersects, params = list("Tfh <4", "Tfh >5"), color = "orange", active = T)
      ),
      order.by = c("freq"),
      nsets = 3, nintersects = NA,
      text.scale = 1.95,
      point.size = 3,
      line.size = 1.5,
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
);
NB_c <- cowplot::plot_grid(NULL, NB$Main_bar, NB$Sizes, NB$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
GC_c <- cowplot::plot_grid(NULL, GC$Main_bar, GC$Sizes, GC$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
NT_c <- cowplot::plot_grid(NULL, NT$Main_bar, NT$Sizes, NT$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
TF_c <- cowplot::plot_grid(NULL, TF$Main_bar, TF$Sizes, TF$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
png("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/Upset_age.png", width=15, height=15, unit="in", res=300);
plot_grid(NB_c, GC_c, NT_c, TF_c, labels = c("A", "B", "C", "D"), ncol = 2, label_size = 25)

dev.off()

#print gene overlaps
cycle <- list(list(l, "NaiveB"), list(m,"GCB"), list(n,"NaiveT"), list(o,"TFH"))

for (i in cycle){
  upsetlist = i[[1]]
  celltype = i[[2]]
  df2 <- data.frame(gene=unique(unlist(upsetlist)))
  df1 <- lapply(upsetlist,function(x){
    data.frame(gene = x)
  }) %>% 
    bind_rows(.id = "celltype")
  
  df_int <- lapply(df2$gene,function(x){
    # pull the name of the intersections
    intersection <- df1 %>% 
      dplyr::filter(gene==x) %>% 
      arrange(celltype) %>% 
      pull("celltype") %>% 
      paste0(collapse = "|")
    
    # build the dataframe
    data.frame(gene = x,int = intersection)
  }) %>% 
    bind_rows()
  
  df_summary <- df_int %>% 
    group_by(int) %>% 
    summarise(n=n()) %>% 
    arrange(desc(n))
  
  #add gene names
  genes <- as.data.frame(genes(EnsDb.Hsapiens.v86))
  
  df_int <- separate(df_int, col = "gene", into = c("ENSG", "transcript?"), sep = "[.]")
  
  df_int_named <- left_join(df_int, select(genes, gene_id, gene_name, gene_biotype, seqnames, start, end), c("ENSG" = "gene_id") )
  
  write.table(df_summary, file= paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/",celltype,"_age_intersections_summary.txt",sep=''), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(df_int_named, file= paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/",celltype,"_age_intersections_named.txt",sep=''), sep="\t", row.names=FALSE, quote=FALSE)
}

#sex comparisons
l <- list(
  "Naive B" = (NaiveB1[[1]]),
  "Naive B female" = (NaiveB1_f[[1]]),
  "Naive B male" = (NaiveB1_m[[1]])
)
m <- list(
  "GC B" = (GCB1[[1]]),
  "GC B female" = (GCB1_f[[1]]),
  "GC B male" = (GCB1_m[[1]])
)
n <- list(
  "Naive T" = (NaiveT1[[1]]),
  "Naive T female" = (NaiveT1_f[[1]]),
  "Naive T male" = (NaiveT1_m[[1]])
)
o <- list(
  "Tfh" = (TFH1[[1]]),
  "Tfh female" = (TFH1_f[[1]]),
  "Tfh male" = (TFH1_m[[1]])
)

NB <- upset(fromList(l), 
            queries = list(
              list(query = intersects, params = list("Naive B female"), color = "orange", active = T),
              list(query = intersects, params = list("Naive B male"), color = "orange", active = T),
              list(query = intersects, params = list("Naive B female","Naive B male"), color = "orange", active = T)
            ),
            order.by = c("freq"),
            nsets = 3, nintersects = NA,
            text.scale = 1.95,
            point.size = 3,
            line.size = 1.5,
            main.bar.color = "gray23",
            matrix.color = "gray23",
            sets.bar.color = "gray23",
            sets.x.label = "total eGenes",
            mainbar.y.label = "eGenes in common"
);
GC <- upset(fromList(m), 
            queries = list(
              list(query = intersects, params = list("GC B female"), color = "orange", active = T),
              list(query = intersects, params = list("GC B male"), color = "orange", active = T),
              list(query = intersects, params = list("GC B female","GC B male"), color = "orange", active = T)
            ),
            order.by = c("freq"),
            nsets = 3, nintersects = NA,
            text.scale = 1.95,
            point.size = 3,
            line.size = 1.5,
            main.bar.color = "gray23",
            matrix.color = "gray23",
            sets.bar.color = "gray23",
            sets.x.label = "total eGenes",
            mainbar.y.label = "eGenes in common"
);
NT <- upset(fromList(n), 
            queries = list(
              list(query = intersects, params = list("Naive T female"), color = "orange", active = T),
              list(query = intersects, params = list("Naive T male"), color = "orange", active = T)
            ),
            order.by = c("freq"),
            nsets = 3, nintersects = NA,
            text.scale = 1.95,
            point.size = 3,
            line.size = 1.5,
            main.bar.color = "gray23",
            matrix.color = "gray23",
            sets.bar.color = "gray23",
            sets.x.label = "total eGenes",
            mainbar.y.label = "eGenes in common"
);
TF <- upset(fromList(o), 
            queries = list(
              list(query = intersects, params = list("Tfh female"), color = "orange", active = T),
              list(query = intersects, params = list("Tfh male"), color = "orange", active = T)
            ),
            order.by = c("freq"),
            nsets = 3, nintersects = NA,
            text.scale = 1.95,
            point.size = 3,
            line.size = 1.5,
            main.bar.color = "gray23",
            matrix.color = "gray23",
            sets.bar.color = "gray23",
            sets.x.label = "total eGenes",
            mainbar.y.label = "eGenes in common"
);
NB_c <- cowplot::plot_grid(NULL, NB$Main_bar, NB$Sizes, NB$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
GC_c <- cowplot::plot_grid(NULL, GC$Main_bar, GC$Sizes, GC$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
NT_c <- cowplot::plot_grid(NULL, NT$Main_bar, NT$Sizes, NT$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
TF_c <- cowplot::plot_grid(NULL, TF$Main_bar, TF$Sizes, TF$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
png("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/Upset_sex.png", width=15, height=15, unit="in", res=300);
plot_grid(NB_c, GC_c, NT_c, TF_c, labels = c("A", "B", "C", "D"), ncol = 2, label_size = 25)

dev.off()

#print gene overlaps
cycle <- list(list(l, "NaiveB"), list(m,"GCB"), list(n,"NaiveT"), list(o,"TFH"))

for (i in cycle){
  upsetlist = i[[1]]
  celltype = i[[2]]
  df2 <- data.frame(gene=unique(unlist(upsetlist)))
  df1 <- lapply(upsetlist,function(x){
    data.frame(gene = x)
  }) %>% 
    bind_rows(.id = "celltype")
  
  df_int <- lapply(df2$gene,function(x){
    # pull the name of the intersections
    intersection <- df1 %>% 
      dplyr::filter(gene==x) %>% 
      arrange(celltype) %>% 
      pull("celltype") %>% 
      paste0(collapse = "|")
    
    # build the dataframe
    data.frame(gene = x,int = intersection)
  }) %>% 
    bind_rows()
  
  df_summary <- df_int %>% 
    group_by(int) %>% 
    summarise(n=n()) %>% 
    arrange(desc(n))
  
  #add gene names
  genes <- as.data.frame(genes(EnsDb.Hsapiens.v86))
  
  df_int <- separate(df_int, col = "gene", into = c("ENSG", "transcript?"), sep = "[.]")
  
  df_int_named <- left_join(df_int, select(genes, gene_id, gene_name, gene_biotype, seqnames, start, end), c("ENSG" = "gene_id") )
  
  write.table(df_summary, file= paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/",celltype,"_sex_intersections_summary.txt",sep=''), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(df_int_named, file= paste("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/",celltype,"_sex_intersections_named.txt",sep=''), sep="\t", row.names=FALSE, quote=FALSE)
}
