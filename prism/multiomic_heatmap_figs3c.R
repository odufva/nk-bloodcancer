
# Plot multi-omic heatmap of NK PRISM sensitivity correlates in CCLE dataset (Figure S3C)
# 20Q4 complete data
# normalized % viability-based AUC

# load scripts
source("compute.pairwise.R")
source("functions_statistics.R")
source("regulon.feats.R")

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
library(GSVA)

# load feature matrix
fm <- readRDS("data/CCLE.fm_20Q4_complete.rds")

# load correlations
results <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_all.tsv", data.table = F)

# add NK cell sensitivity PRISM data to feature matrix
nk_prism <- fread("data/NK_PRISM_heme_pctviability.txt")
nk_prism_fm <- nk_prism %>%
  dplyr::select(Condition, pctviability, DepMap_ID) %>%
  pivot_wider(names_from = DepMap_ID, values_from = pctviability) %>%
  mutate(Condition = paste0("N:PRSM:", Condition))

nk_prism_auc <- fread("data//NK_PRISM_heme_pctviability_auc_norm.txt")
nk_prism_fm_auc <- nk_prism_auc %>%
  filter(Condition == "0") %>%
  dplyr::select(auc_norm, DepMap_ID) %>%
  mutate(Condition = "N:PRSM:AUC") %>%
  pivot_wider(names_from = DepMap_ID, values_from = auc_norm)

nk_prism_fm <- rbind(nk_prism_fm[,colnames(nk_prism_fm_auc)], nk_prism_fm_auc)

nk_prism_fm <- as.data.frame(nk_prism_fm)
rownames(nk_prism_fm) <- nk_prism_fm$Condition
nk_prism_fm$Condition <- NULL
nk_prism_fm$rownames <- rownames(nk_prism_fm)

fm$row_names <- rownames(fm)

fm_full <- full_join(fm, nk_prism_fm)
rownames(fm_full) <- c(rownames(fm), rownames(nk_prism_fm))

# subset matrix to cell lines with PRISM data
fm_prism <- fm_full[,!is.na(fm_full["N:PRSM:D1_NK_5",])]
fm_prism$rownames <- rownames(fm_full)

## -----------------------------------

# heatmaps

# order fm by PRISM NK sensitivity
fm_prism_clean <- fm_prism[,!grepl("sum|var|rownames", colnames(fm_prism))]
order <- order(as.numeric(fm_prism_clean["N:PRSM:AUC",]))
fm_prism_clean <- fm_prism_clean[,order]

# load CCLE sample info
ccle_annot <- fread("../NK_PRISM/sample_info.csv")

# match subtype from sample info to feature matrix
Subtype <- ccle_annot$Subtype[match(colnames(fm_prism_clean), ccle_annot$DepMap_ID)]
subtype_short = ifelse(grepl("AML", Subtype), "AML",
                       ifelse(grepl("ALL), B-cell", Subtype), "B-ALL",
                              ifelse(grepl("ALL), T-cell", Subtype), "T-ALL",
                                     ifelse(grepl("Multiple", Subtype), "MM",
                                            ifelse(grepl("CLL", Subtype), "CLL",
                                                   ifelse(grepl("Mantle", Subtype), "MCL",
                                                          ifelse(grepl("B-cell, Hodgkins", Subtype), "CHL",
                                                                 ifelse(grepl("Burkitt", Subtype), "BL",
                                                                        ifelse(grepl("ALCL", Subtype), "ALCL",
                                                                               ifelse(grepl("DLBCL", Subtype), "DLBCL", 
                                                                                      ifelse(grepl("B-cell, Non", Subtype), "BCL other",
                                                                                             ifelse(grepl("Cutaneous", Subtype), "CTCL",
                                                                                                    ifelse(grepl("^T-cell", Subtype), "TCL other", Subtype)))))))))))))


# data frame with subtype for heatmap annotation
annot <- data.frame(cancertype = subtype_short)
annot$cancertype[annot$cancertype %in% c("DLBCL", "BCL other", "CHL", "BL", "MCL", "CLL")] <- "BCL"
annot$cancertype[annot$cancertype %in% c("CTCL", "ALCL", "ATL", "TCL other")] <- "TCL"

colnames(annot) <- "Cancer type"

# replace DepMap IDs with cell line names
colnames(fm_prism_clean) <- ccle_annot$stripped_cell_line_name[match(colnames(fm_prism_clean), ccle_annot$DepMap_ID)]

# add HLA I score
dat_a=fm_prism_clean[rownames(fm_prism_clean)%in%c("N:GEXP:B2M",
                                                   "N:GEXP:HLA-A",
                                                   "N:GEXP:HLA-B",
                                                   "N:GEXP:HLA-C"),]

dat=2^dat_a+0.01
gm1=log2(t(apply(dat, 2, gm_mean)))
rownames(gm1)="N:SAMP:HLAIScore"

fm_prism_clean <- rbind(fm_prism_clean, gm1)

# function to create scaled matrices with selected alterations
heatmap_scaled_matrix <- function(DATAPAIR, ETRATIO, PVAL = 1, ADJP = 1, SCALE = T, TOPN = 50){
  
  # filter features to plot
  feats_top <- results %>%
    filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
    arrange(p) %>%
    filter(cor > 0) %>%
    top_n(n = -TOPN, wt = p) %>%
    arrange(desc(cor)) %>%
    dplyr::select(featureB) %>%
    deframe()
  
  feats_bottom <- results %>%
    filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
    arrange(p) %>%
    filter(cor < 0) %>%
    top_n(n = -TOPN, wt = p) %>%
    arrange(cor) %>%
    dplyr::select(featureB) %>%
    deframe()
  
  # select features from matrix
  mat <- fm_prism_clean[c(feats_bottom, feats_top),]
  
  # Z-score
  if (SCALE == T) {
    mat_scaled <- t(apply(mat, 1, scale))
    colnames(mat_scaled) <- colnames(mat)
    rownames(mat_scaled) <- gsub("N:....:", "", rownames(mat_scaled)) # clean row names
    return(mat_scaled)
  }
  
  # no Z-score
  else
  { 
    rownames(mat) <- gsub("^.:....:", "", rownames(mat))
    return(mat)}
  
}

# apply function to all data types to get matrices for heatmaps
mat_gexp <- heatmap_scaled_matrix(DATAPAIR = "GEXP", ETRATIO = "AUC", TOPN = 31)
mat_rppa <- heatmap_scaled_matrix(DATAPAIR = "RPPA", ETRATIO = "AUC", TOPN = 5)
mat_mirn <- heatmap_scaled_matrix(DATAPAIR = "MIRN", ETRATIO = "AUC", TOPN = 5)
mat_lcms <- heatmap_scaled_matrix(DATAPAIR = "LCMS", ETRATIO = "AUC", TOPN = 5)

# function to create annotations of correlation direction
heatmap_annotation_correlation <- function(DATAPAIR, ETRATIO, PVAL = 1, ADJP = 1, TOPN = 50){
  
  # filter features to plot
  feats_top <- results %>%
    filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
    arrange(p) %>%
    filter(cor > 0) %>%
    top_n(n = -TOPN, wt = p) %>%
    arrange(desc(cor)) %>%
    dplyr::select(featureB) %>%
    deframe()
  
  feats_bottom <- results %>%
    filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
    arrange(p) %>%
    filter(cor < 0) %>%
    top_n(n = -TOPN, wt = p) %>%
    arrange(desc(cor)) %>%
    dplyr::select(featureB) %>%
    deframe()
  
  # select features from matrix
  vector <- c(rep("1", length(feats_top)), rep("2", length(feats_bottom)))
  
  return(vector)
}

# apply function to all data types to get annotations for heatmaps
dir_rppa <- heatmap_annotation_correlation(DATAPAIR = "GEXP", ETRATIO = "AUC", TOPN = 5)
dir_mirn <- heatmap_annotation_correlation(DATAPAIR = "GNAB", ETRATIO = "AUC", TOPN = 5)
dir_lcms <- heatmap_annotation_correlation(DATAPAIR = "GEXP", ETRATIO = "AUC", TOPN = 5)


# GSEA

gsea <- fread("results_gsea/NK_PRISM_all_HALLMARK_GSEA.txt", data.table = F)

pathways_hallmark <- fgsea::gmtPathways("h.all.v7.0.symbols.gmt")

fm_gsva <- fm_prism_clean[grepl("N:GEXP", rownames(fm_prism_clean)),]
rownames(fm_gsva) <- gsub("N:GEXP:", "", rownames(fm_gsva))

# run GSVA
gsva_results <- gsva(as.matrix(fm_gsva), pathways_hallmark)

genesets_top <- gsea %>% filter(padj < 0.05) %>% select(pathway) %>% tibble::deframe()
mat_gsea <- t(apply(gsva_results[genesets_top,], 1, scale))
rownames(mat_gsea) <- gsub("Il6 jak stat3", "IL6-JAK-STAT3",
                           gsub("Il2 stat5", "IL2-STAT5",
                                gsub("Dna", "DNA",
                                     gsub("Myc", "MYC",
                                          gsub("E2f", "E2F",
                                               gsub("Kras", "KRAS",
                                                    gsub("Tnfa", "TNFA",
                                                         gsub("G2m", "G2M",
                                                              gsub("nfkb|Nfkb", "NF-kB",
                                                                   gsub("Ifng", "IFNy",
                                                                        gsub("Nk response", "Core NK cell response",
                                                                             stringr::str_to_sentence(gsub("_", " ",
                                                                                                           sub(".*?_", "", rownames(mat_gsea)))))))))))))))
colnames(mat_gsea) <- colnames(gsva_results)

dir_gsea <- gsea %>% filter(padj < 0.05) %>% mutate(direction = ifelse(sign(NES)>0, "2", "1")) %>% select(direction) %>% tibble::deframe()

# read cancer types and colors
cols <- fread("nk_crispr_colors.txt", data.table = F) %>% dplyr::rename(cancer = cancer_type)
cols$cancer_main <- gsub("ALCL", "TCL", gsub("DLBCL", "BCL", as.character(cols$cancer)))

cols_vector <- cols$color
names(cols_vector) <- cols$cancer_main

# create heatmap annotations
ha1 <- HeatmapAnnotation(`PRISM AUC` = anno_barplot(as.numeric(fm_prism_clean["N:PRSM:AUC",]), 
                                                    bar_width = 0.75, 
                                                    border = FALSE, 
                                                    axis = TRUE,
                                                    axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                    gp = gpar(col = NA, fill = "grey30")),
                         df = annot, col = list(`Cancer type` = cols_vector),
                         gap = unit(0.75, "mm"),
                         annotation_legend_param = list(`Cancer type` = list(title = "Cancer type", title_gp = gpar(fontsize = 5), 
                                                                             labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))),
                         height = unit(1, "cm"),
                         simple_anno_size_adjust = T,
                         show_annotation_name = T,
                         annotation_name_gp = gpar(fontsize = 8)
)

# gexp vs meth correlations

# genelist
genelist = rownames(mat_gexp)

# choose feature types for correlations
extrafeatures=NULL

l.regulon.gene=regulon.feats(fm_prism, genelist)

# run correlations
results_gexp=pairwise.correlation(l.regulon.gene, fm_prism, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

gexp_meth_cor <- results_gexp %>%
  filter(datapairs == "GEXP:METH") %>%
  mutate(featureA = gsub("N:GEXP:", "", featureA))

rownames(gexp_meth_cor) <- gexp_meth_cor$featureA
gexp_meth_cor <- gexp_meth_cor[rownames(mat_gexp),]
gexp_meth_cor_df <- data.frame(R = as.numeric(gexp_meth_cor$cor),
                               FDR = gexp_meth_cor$signifCode)
gexp_meth_cor_df$FDR <- "<5%"
gexp_meth_cor_df$FDR[as.numeric(gexp_meth_cor$adj.p)>0.05] <- ">5%"
gexp_meth_cor_df$FDR[is.na(gexp_meth_cor$signifCode)] <- ">5%"

meth_genes <- as.character(na.omit(rownames(mat_gexp)[as.numeric(gexp_meth_cor$adj.p) < 0.002]))


mat_meth <- fm_prism_clean[paste0("N:METH:", meth_genes),]
rownames(mat_meth) <- gsub("N:METH:", "", rownames(mat_meth))

gexp_meth_cor_df[rownames(mat_gexp) %in% meth_genes,]

meth_cor <- results %>%
  filter(datapairs == "PRSM:GEXP", featureA == "N:PRSM:AUC", featureB %in% paste0("N:GEXP:", meth_genes)) %>%
  dplyr::select(featureB, cor)

rownames(meth_cor) <- gsub("N:GEXP:", "", meth_cor$featureB)

dir_meth <- ifelse(sign(meth_cor[rownames(mat_meth),]$cor)>0, "1", "2")

# clean row names
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

rownames(mat_rppa) <- gsub("_", " ", gsub("_Caution", "", rownames(mat_rppa)))

rownames(mat_lcms) <- gsub("sm", "SM", gsub("tag", "TAG", firstup(rownames(mat_lcms))))



# plot heatmaps

ht_gsea <- Heatmap(mat_gsea,
                   name = "gsea",
                   split = dir_gsea,
                   row_title = c("Sens", "Res"),
                   row_title_gp = gpar(fontsize = 6),
                   col = colorRamp2(seq(quantile(mat_gsea, 0.975)*(-1), quantile(mat_gsea, 0.975), length.out = 9), pals::ocean.thermal(9)),
                   top_annotation = ha1,
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   show_row_names = T,
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Gene set (GSVA)",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               at = c(-2,2),
                                               labels = c("Low", "High"),
                                               title_position = "topcenter",
                                               legend_direction = "horizontal",
                                               border = NA
                   )
)

ht_rppa <- Heatmap(mat_rppa,
                   name = "rppa",
                   split = dir_rppa,
                   row_title = c("Sens", "Res"),
                   row_title_gp = gpar(fontsize = 6),
                   col = colorRamp2(seq(quantile(mat_rppa, 0.975, na.rm = T)*(-1), quantile(mat_rppa, 0.975, na.rm = T), length.out = 11),  rev(brewer.pal(11, "PuOr"))),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Protein (RPPA)",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               title_position = "topcenter",
                                               legend_direction = "horizontal",
                                               at = c(-2,2),
                                               labels = c("Low", "High"),
                                               border = NA)
)

ht_mirn <- Heatmap(mat_mirn,
                   name = "mirn",
                   split = dir_mirn,
                   row_title = c("Sens", "Res"),
                   row_title_gp = gpar(fontsize = 6),
                   col = colorRamp2(seq(quantile(mat_mirn, 0.975, na.rm = T)*(-1), quantile(mat_mirn, 0.975, na.rm = T), length.out = 11),  rev(brewer.pal(11, "RdYlGn"))),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "miRNA",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               title_position = "topcenter",
                                               legend_direction = "horizontal",
                                               at = c(-2,2),
                                               labels = c("Low", "High"),
                                               border = NA)
)

ht_lcms <- Heatmap(mat_lcms,
                   name = "lcms",
                   split = dir_lcms,
                   row_title = c("Sens", "Res"),
                   row_title_gp = gpar(fontsize = 6),
                   col = colorRamp2(seq(quantile(mat_lcms, 0.975, na.rm = T)*(-1), quantile(mat_lcms, 0.975, na.rm = T), length.out = 11),  rev(brewer.pal(11, "PRGn"))),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Metabolite (LCMS)",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               title_position = "topcenter",
                                               legend_direction = "horizontal",
                                               at = c(-2,2),
                                               labels = c("Low", "High"),
                                               border = NA)
)


ht_meth <- Heatmap(mat_meth,
                   name = "meth",
                   split = dir_meth,
                   row_title = c("Sens", "Res"),
                   row_title_gp = gpar(fontsize = 6),
                   col = colorRamp2(seq(0, 1, length.out = 9), brewer.pal(9, "PuRd")),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Methylation",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               title_position = "topcenter",
                                               legend_direction = "horizontal",
                                               border = NA)
)


pdf("CCLE_NK_PRISM_heatmap_multiomic_manuscript_small_pctviability_supplement.pdf", height = 4.5, width = 6)
draw(ht_gsea %v% ht_meth %v% ht_rppa %v% ht_mirn %v% ht_lcms, heatmap_legend_side = "left",
     annotation_legend_side = "left", merge_legend = T, ht_gap = unit(c(2,2,2,2), "mm"))
dev.off()

