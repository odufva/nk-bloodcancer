
# GSEA of NK CRISPR screen hits

# load libraries
library(data.table)
library(parallel)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(fgsea)
library(ggplot2)
library(cowplot)
library(scales)

# load screen results
data <- fread("crispr_mageck_combined.txt", data.table = F)

# run fgsea
pathways_reactome <- gmtPathways("c2.cp.reactome.v7.0.symbols.gmt")
pathways_hallmark <- gmtPathways("h.all.v7.0.symbols.gmt")
pathways_biocarta <- gmtPathways("c2.cp.biocarta.v7.0.symbols.gmt")

dir.create("results_GSEA")

run_gsea <- function(DATA, NAME, CELLLINE, PATHWAYS){
  
  ranks <- DATA %>% 
    filter(cell_line == CELLLINE) %>%
    dplyr::select(Gene, lfc, p) %>%
    mutate(p = ifelse(p == 0, 0.00000001, p)) %>%
    mutate(signed_p = -log10(p)*sign(lfc)) %>%
    dplyr::select(Gene, signed_p) %>%
    deframe()
  
  # GSEA
  fgseaRes <- fgsea(PATHWAYS, ranks, maxSize=500)
  fgseaRes <- fgseaRes %>%
    mutate(pos_neg = ifelse(ES > 0, "pos", "neg")) %>%
    arrange(pos_neg, pval) %>%
    select(-pos_neg)
  
  # table of top pathways
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  pdf(paste0("results_GSEA/", NAME, "_", CELLLINE, ".pdf"), height = 7, width = 20)
  plotGseaTable(PATHWAYS[topPathways], ranks, fgseaRes, gseaParam = 0.5)
  dev.off()
  
  fwrite(fgseaRes, paste0("results_GSEA/", NAME, "_", CELLLINE, ".txt"), quote = F, row.names = F, sep = "\t")
}


celllines <- unique(data$cell_line)

lapply(celllines, run_gsea, DATA = data, NAME = "NK_CRISPR_HALLMARK_GSEA", PATHWAYS = pathways_hallmark)
lapply(celllines, run_gsea, DATA = data, NAME = "NK_CRISPR_REACTOME_GSEA", PATHWAYS = pathways_reactome)
lapply(celllines, run_gsea, DATA = data, NAME = "NK_CRISPR_BIOCARTA_GSEA", PATHWAYS = pathways_biocarta)




# Dot plot 

# read GSEA
read_gsea <- function(PATHWAYS){
  
  PATTERN = paste0(PATHWAYS, ".*", ".txt")
  
  files <- list.files(path = "results_GSEA", pattern = PATTERN)
  
  read_file <- function(FILE){
    fread(paste0("results_GSEA/", FILE), data.table = F) %>% mutate(cancer_type = gsub(".*GSEA_|.txt", "", FILE))
  }
  
  lapply(files, read_file) %>% bind_rows()
}

# REACTOME
gsea <- read_gsea("REACTOME")

gsea <- gsea %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))

signif <- gsea %>% filter(pval < 0.001) %>% select(pathway) %>% unique()

signif_sentenceCase <- signif %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))

# subset data to recurrent pathways with p value < 0.05
gsea_plot <- gsea %>%
  filter(pval < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cancer_type = factor(cancer_type, levels = unique(cancer_type)[col.ord])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>1)

# list of pathways with adjusted p value < 0.15
pathways_signif <- gsea %>%
  filter(padj < 0.15) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cancer_type = factor(cancer_type, levels = unique(cancer_type)[col.ord])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>1) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

gsea_plot <- gsea_plot %>% filter(pathway_sentenceCase %in% pathways_signif)

# order by NES
nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

gsea_plot <- gsea_plot %>%
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order)) %>% 
  mutate(cancer_type = gsub("_", " ", cancer_type)) %>% 
  mutate(cancer_type = factor(cancer_type, levels = c("K562", "MOLM14", "SUDHL4", "MM1S", "LP1", "KMS11", "KMS11 KHYG1", "NALM6", "MM1S GOF", "LP1 KHYG1 GOF", "KMS11 KHYG1 GOF")))

 
# Write supplement table (Table S4O)
fwrite(gsea_plot_reactome, "gsea_reactome_recurrent_supplement.txt", sep = "\t")

