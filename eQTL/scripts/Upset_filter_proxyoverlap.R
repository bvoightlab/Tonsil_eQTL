library(UpSetR)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

NaiveB1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveB.sig_egenes.txt", header = TRUE)
GCB1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/GCB.sig_egenes.txt", header = TRUE)
NaiveT1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/NaiveT.sig_egenes.txt", header = TRUE)
TFH1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/TFH.sig_egenes.txt", header = TRUE)

NaiveB2 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/DICE_data/B_CELL_NAIVE_filtered_genes.txt", header = TRUE)
NaiveT2 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/DICE_data/CD4_NAIVE_filtered_genes.txt", header = TRUE)
TFH2 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/DICE_data/TFH_filtered_genes.txt", header = TRUE)

eQTLgen1 <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/eQTLgen_filtered_genes.txt", header = TRUE)


NaiveB1 <- separate(NaiveB1, phenotype_id, c("gene", "transcript"), sep = "\\.")
GCB1 <- separate(GCB1, phenotype_id, c("gene", "transcript"), sep = "\\.")
NaiveT1 <- separate(NaiveT1, phenotype_id, c("gene", "transcript"), sep = "\\.")
TFH1 <- separate(TFH1, phenotype_id, c("gene", "transcript"), sep = "\\.")

leads_with_proxies <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tonsil_leads_proxies/leads_with_DICE_proxies.txt", header = FALSE)
leads_with_proxies_eQTLgen <- read.delim("/project/voightlab_01/lorenzk/Romberg_eQTL/tonsil_leads_proxies/leads_with_eQTLgen_proxies.txt", header = FALSE)

#screen our leads to only keep those with a proxy checked in DICE

NaiveB1_overlap <- filter(NaiveB1, rsid_dbSNP154 %in% leads_with_proxies$V1)
GCB1_overlap <- filter(GCB1, rsid_dbSNP154 %in% leads_with_proxies$V1)
NaiveT1_overlap <- filter(NaiveT1, rsid_dbSNP154 %in% leads_with_proxies$V1)
TFH1_overlap <- filter(TFH1, rsid_dbSNP154 %in% leads_with_proxies$V1)

NaiveB1_overlap_eQTLgen <- filter(NaiveB1, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)
GCB1_overlap_eQTLgen <- filter(GCB1, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)
NaiveT1_overlap_eQTLgen <- filter(NaiveT1, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)
TFH1_overlap_eQTLgen <- filter(TFH1, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)


l <- list(
  "Naive B" = (NaiveB1[[1]]),
  "GC B" = (GCB1[[1]]),
  "Naive T" = (NaiveT1[[1]]),
  Tfh = (TFH1[[1]])
)

png("/project/voightlab_01/lorenzk/Romberg_eQTL/Paper_Figures/Upset_ours.png", width=10, height=5.5, unit="in", res=300);
upset(fromList(l), 
      order.by = "freq",
      nsets = 4, nintersects = NA,
      text.scale = 1.85,
      point.size = 3,
      line.size = 1.5,
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
);
dev.off();

#vs eQTLgen
g <- list(
  NaiveB = (NaiveB1_overlap_eQTLgen[[1]]),
  GCB = (GCB1_overlap_eQTLgen[[1]]),
  NaiveT = (NaiveT1_overlap_eQTLgen[[1]]),
  TFH = (TFH1_overlap_eQTLgen[[1]]),
  eQTLgen = (eQTLgen1[[8]])
)

png("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/Upset_eQTLGen_proxy.png", width=15, height=5.5, unit="in", res=300);
upset(fromList(g), 
      queries = list(
        list(query = intersects, params = list("NaiveB"), color = "orange", active = T),
        list(query = intersects, params = list("NaiveT"), color = "orange", active = T),
        list(query = intersects, params = list("GCB"), color = "orange", active = T),
        list(query = intersects, params = list("TFH"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "NaiveB", "TFH", "NaiveT"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "TFH", "NaiveT"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "NaiveB", "NaiveT"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "NaiveB", "TFH"), color = "orange", active = T),
        list(query = intersects, params = list("NaiveT", "NaiveB", "TFH"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "NaiveB"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "TFH"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "NaiveT"), color = "orange", active = T),
        list(query = intersects, params = list("NaiveB", "TFH"), color = "orange", active = T),
        list(query = intersects, params = list("NaiveB", "NaiveT"), color = "orange", active = T),
        list(query = intersects, params = list("TFH", "NaiveT"), color = "orange", active = T)
      ),
      #query.legend = "bottom",
      order.by = "freq",
      nsets = 5, nintersects = NA,
      text.scale = 1.85,
      point.size = 3,
      line.size = 1.5,
      #show.numbers = F, 
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total egenes",
      mainbar.y.label = "egenes in common"
); 
dev.off();


#Tcell comparisons
t <- list(  
  NaiveT_DICE = (NaiveT2[[1]]),
  TFH_DICE = (TFH2[[1]]),
  NaiveT = (NaiveT1_overlap[[1]]),
  TFH = (TFH1_overlap[[1]])
)
png("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/Upset_DICE_Tcell_proxy.png", width=10, height=7, unit="in", res=300);
upset(fromList(t), 
      queries = list(
        list(query = intersects, params = list("NaiveT"), color = "orange", active = T),
        list(query = intersects, params = list("TFH"), color = "orange", active = T),
        list(query = intersects, params = list("TFH", "NaiveT"), color = "orange", active = T)
      ),
      order.by = "freq",
      nsets = 4, nintersects = NA,
      text.scale = 1.85,
      point.size = 3,
      line.size = 1.5,
      #show.numbers = F, 
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total egenes",
      mainbar.y.label = "egenes in common"
);
dev.off()

#Bcell comparisons
b <- list(
  NaiveB = (NaiveB1_overlap[[1]]),
  GCB = (GCB1_overlap[[1]]),
  NaiveB_DICE = (NaiveB2[[1]])
)

png("/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/Upset_DICE_Bcell_proxy.png", width=10, height=7, unit="in", res=300);
upset(fromList(b), 
      queries = list(
        list(query = intersects, params = list("NaiveB"), color = "orange", active = T),
        list(query = intersects, params = list("GCB"), color = "orange", active = T),
        list(query = intersects, params = list("GCB", "NaiveB"), color = "orange", active = T)
      ),
      order.by = "freq",
      nsets = 3, nintersects = NA,
      text.scale = 1.85,
      point.size = 3,
      line.size = 1.5,
      #show.numbers = F, 
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total egenes",
      mainbar.y.label = "egenes in common"
);
dev.off()

#one comparison with all the things. 

NaiveB1_overlap_both <- filter(NaiveB1_overlap, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)
GCB1_overlap_both <- filter(GCB1_overlap, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)
NaiveT1_overlap_both <- filter(NaiveT1_overlap, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)
TFH1_overlap_both <- filter(TFH1_overlap, rsid_dbSNP154 %in% leads_with_proxies_eQTLgen$V1)

#make one list of DICE stuff
alldice <- unique(c((NaiveT2[[1]]), (TFH2[[1]]), (NaiveB2[[1]])))
Tdice <- unique(c((NaiveT2[[1]]), (TFH2[[1]])))

TonsilB <- unique(c((NaiveB1_overlap_both[[1]]), (GCB1_overlap_both[[1]])))
TonsilT <- unique(c((NaiveT1_overlap_both[[1]]), (TFH1_overlap_both[[1]])))

E <- list(
  "Tonsil B" = (TonsilB),
  "Tonsil T" = (TonsilT),
  "DICE B" = (NaiveB2[[1]]),
  "DICE T" = (Tdice),
  eQTLgen = (eQTLgen1[[8]])
)

png("/project/voightlab_01/lorenzk/Romberg_eQTL/Paper_Figures/Upset_compare_all_proxy.png", width=17, height=5.5, unit="in", res=300);
upset(fromList(E), 
      queries = list(
        list(query = intersects, params = list("Tonsil B"), color = "orange", active = T),
        list(query = intersects, params = list("Tonsil B","Tonsil T"), color = "orange", active = T),
        list(query = intersects, params = list("Tonsil T"), color = "orange", active = T)
      ),
      #query.legend = "bottom",
      order.by = "freq",
      sets = c("eQTLgen", "DICE T", "DICE B", "Tonsil T", "Tonsil B"),
      keep.order=TRUE,
      nsets = 5, nintersects = NA,
      text.scale = 1.85,
      point.size = 3,
      line.size = 1.5,
      #show.numbers = F, 
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
); 
dev.off()

png("/project/voightlab_01/lorenzk/Romberg_eQTL/Paper_Figures/Upset_compare_all_proxy_20.png", width=17, height=5.5, unit="in", res=300);
upset(fromList(E), 
      queries = list(
        list(query = intersects, params = list("Tonsil B"), color = "orange", active = T),
        list(query = intersects, params = list("Tonsil B","Tonsil T"), color = "orange", active = T),
        list(query = intersects, params = list("Tonsil T"), color = "orange", active = T)
      ),
      #query.legend = "bottom",
      order.by = "freq",
      sets = c("eQTLgen", "DICE T", "DICE B", "Tonsil T", "Tonsil B"),
      keep.order=TRUE,
      nsets = 5, nintersects = 20,
      text.scale = 2.2,
      point.size = 3,
      line.size = 1.5,
      #show.numbers = F, 
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
); 
dev.off()

png("/project/voightlab_01/lorenzk/Romberg_eQTL/Paper_Figures/Upset_compare_all_proxy_15.png", width=17, height=5.5, unit="in", res=300);
upset(fromList(E), 
      queries = list(
        list(query = intersects, params = list("Tonsil B"), color = "orange", active = T),
        list(query = intersects, params = list("Tonsil B","Tonsil T"), color = "orange", active = T),
        list(query = intersects, params = list("Tonsil T"), color = "orange", active = T)
      ),
      #query.legend = "bottom",
      order.by = "freq",
      sets = c("eQTLgen", "DICE T", "DICE B", "Tonsil T", "Tonsil B"),
      keep.order=TRUE,
      nsets = 5, nintersects = 15,
      text.scale = 2.2,
      point.size = 3,
      line.size = 1.5,
      #show.numbers = F, 
      main.bar.color = "gray23",
      matrix.color = "gray23",
      sets.bar.color = "gray23",
      sets.x.label = "total eGenes",
      mainbar.y.label = "eGenes in common"
); 
dev.off()
#print gene overlaps

df2 <- data.frame(gene=unique(unlist(E)))
df1 <- lapply(E,function(x){
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
 

 df_int_named <- left_join(df_int, select(genes, gene_id, gene_name, gene_biotype, seqnames, start, end), c("gene" = "gene_id") )
 
 write.table(df_summary, file="/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/all_intersections_summary.txt", sep="\t", row.names=FALSE, quote=FALSE)
 write.table(df_int_named, file="/project/voightlab_01/lorenzk/Romberg_eQTL/tensorqtl/all_intersections_named.txt", sep="\t", row.names=FALSE, quote=FALSE)



