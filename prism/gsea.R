
# Compute GSEA of gene expression correlations for NK PRISM sensitivity in CCLE dataset

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

# load correlations
cor_all <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_all.tsv", data.table = F)
cor_mm <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_MM.tsv", data.table = F)
cor_aml <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_AML.tsv", data.table = F)
cor_ball <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_B-ALL.tsv", data.table = F)
cor_tall <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_T-ALL.tsv", data.table = F)
cor_bcl <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_BCL.tsv", data.table = F)

dir.create("results_GSEA_20Q4_complete")

# run fgsea

pathways_reactome <- gmtPathways("c2.cp.reactome.v7.0.symbols.gmt")
pathways_hallmark <- gmtPathways("h.all.v7.0.symbols.gmt")
pathways_biocarta <- gmtPathways("c2.cp.biocarta.v7.0.symbols.gmt")


run_gsea <- function(RESULT, NAME, PATHWAYS){
  
  ranks <- RESULT %>% 
    filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:GEXP") %>%
    dplyr::select(featureB, cor, p) %>%
    mutate(p = ifelse(p == 0, 0.00000001, p)) %>%
    mutate(signed_p = -log10(p)*sign(cor),
           featureB = gsub("N:GEXP:", "", featureB)) %>%
    dplyr::select(featureB, signed_p) %>%
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

  pdf(paste0("results_GSEA_20Q4_complete/", NAME, "_GSEA.pdf"), height = 7, width = 20)
  plotGseaTable(PATHWAYS[topPathways], ranks, fgseaRes, gseaParam = 0.5)
  dev.off()
  
  fwrite(fgseaRes, paste0("results_GSEA_20Q4_complete/", NAME, "_GSEA.txt"), quote = F, row.names = F, sep = "\t")
}

run_gsea(cor_all, "NK_PRISM_all_HALLMARK", pathways_hallmark)
run_gsea(cor_mm, "NK_PRISM_MM_HALLMARK", pathways_hallmark)
run_gsea(cor_aml, "NK_PRISM_AML_HALLMARK", pathways_hallmark)
run_gsea(cor_ball, "NK_PRISM_BALL_HALLMARK", pathways_hallmark)
run_gsea(cor_tall, "NK_PRISM_TALL_HALLMARK", pathways_hallmark)
run_gsea(cor_bcl, "NK_PRISM_BCL_HALLMARK", pathways_hallmark)

run_gsea(cor_all, "NK_PRISM_all_REACTOME", pathways_reactome)
run_gsea(cor_mm, "NK_PRISM_MM_REACTOME", pathways_reactome)
run_gsea(cor_aml, "NK_PRISM_AML_REACTOME", pathways_reactome)
run_gsea(cor_ball, "NK_PRISM_BALL_REACTOME", pathways_reactome)
run_gsea(cor_tall, "NK_PRISM_TALL_REACTOME", pathways_reactome)
run_gsea(cor_bcl, "NK_PRISM_BCL_REACTOME", pathways_reactome)

run_gsea(cor_all, "NK_PRISM_all_BIOCARTA", pathways_biocarta)
run_gsea(cor_mm, "NK_PRISM_MM_BIOCARTA", pathways_biocarta)
run_gsea(cor_aml, "NK_PRISM_AML_BIOCARTA", pathways_biocarta)
run_gsea(cor_ball, "NK_PRISM_BALL_BIOCARTA", pathways_biocarta)
run_gsea(cor_tall, "NK_PRISM_TALL_BIOCARTA", pathways_biocarta)
run_gsea(cor_bcl, "NK_PRISM_BCL_BIOCARTA", pathways_biocarta)
