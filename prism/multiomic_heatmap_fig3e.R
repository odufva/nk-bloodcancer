
# Plot multi-omic heatmap of NK PRISM sensitivity correlates in CCLE dataset (Figure 3E)

# CCLE 20Q4 complete data
# normalized % viability-based AUC

# load scripts
source("compute.pairwise.R")
source("functions_generate_fm.R")
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
fm <- readRDS("CCLE.fm_20Q4_complete.rds")

# load correlations
results <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_all.tsv", data.table = F)

# add NK cell sensitivity PRISM data to feature matrix
nk_prism <- fread("NK_PRISM_heme_pctviability.txt")
nk_prism_fm <- nk_prism %>%
  dplyr::select(Condition, pctviability, DepMap_ID) %>%
  pivot_wider(names_from = DepMap_ID, values_from = pctviability) %>%
  mutate(Condition = paste0("N:PRSM:", Condition))

nk_prism_auc <- fread("NK_PRISM_heme_pctviability_auc_norm.txt")
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
ccle_annot <- fread("sample_info.csv")

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
mat_gdsc <- heatmap_scaled_matrix(DATAPAIR = "GDSC", ETRATIO = "AUC", TOPN = 5)
mat_mirn <- heatmap_scaled_matrix(DATAPAIR = "MIRN", ETRATIO = "AUC", TOPN = 5)
mat_lcms <- heatmap_scaled_matrix(DATAPAIR = "LCMS", ETRATIO = "AUC", TOPN = 5)
mat_gnab <- heatmap_scaled_matrix(DATAPAIR = "GNAB", ETRATIO = "AUC", SCALE = F, TOPN = 12)

# remove 13th resistance row
mat_gnab <- mat_gnab[1:24,]

mat_gnab[mat_gnab == "1"] <- "Yes"
mat_gnab[mat_gnab == "0"] <- "No"

mat_score <- fm_prism_clean["N:SAMP:HLAIScore",]
mat_score_scaled <- t(apply(mat_score, 1, scale))
colnames(mat_score_scaled) <- colnames(mat_score)
rownames(mat_score_scaled) <- "HLA I score"

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
dir_gexp <- heatmap_annotation_correlation(DATAPAIR = "GEXP", ETRATIO = "AUC", TOPN = 31)
dir_gnab <- c(rep("1", 12), rep("2", 12))


# GSEA

core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")

genesets = list(SCRNASEQ_NK_RESPONSE = core)#,

gsea <- fread("results_gsea/NK_PRISM_all_HALLMARK_GSEA.txt", data.table = F)

pathways_hallmark <- fgsea::gmtPathways("h.all.v7.0.symbols.gmt")

fm_gsva <- fm_prism_clean[grepl("N:GEXP", rownames(fm_prism_clean)),]
rownames(fm_gsva) <- gsub("N:GEXP:", "", rownames(fm_gsva))

# run GSVA
gsva_results <- gsva(as.matrix(fm_gsva), c(pathways_hallmark, genesets))

genesets_top <- gsea %>% filter(padj < 0.05) %>% select(pathway) %>% tibble::deframe()
mat_gsea <- t(apply(gsva_results[c("SCRNASEQ_NK_RESPONSE", genesets_top),], 1, scale))
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

# only NK response gene set
mat_gsea <- mat_gsea[1,, drop = F]

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
                                                    gp = gpar(col = NA, fill = "grey30")),#,
                         df = annot, col = list(`Cancer type` = cols_vector),
                         gap = unit(0.75, "mm"),
                         annotation_legend_param = list(`Cancer type` = list(title = "Cancer type", title_gp = gpar(fontsize = 5), 
                                                                             labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))),
                         height = unit(1, "cm"),
                         simple_anno_size_adjust = T,
                         show_annotation_name = T,
                         annotation_name_gp = gpar(fontsize = 10)
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

row_ha = rowAnnotation(df = gexp_meth_cor_df,
                       col = list(R = colorRamp2(seq(-1, 0, length.out = 9), rev(brewer.pal(9, "PuRd"))),
                                  FDR = structure(c("<5%" = "grey30", ">5%" = "grey90"))),
                       simple_anno_size_adjust = TRUE,
                       width = unit(0.35, "cm"),
                       annotation_legend_param = list(R = list(title = "Correlation with methylation",
                                                               title_gp = gpar(fontsize = 5), 
                                                               labels_gp = gpar(fontsize = 5),
                                                               grid_height = unit(0.2, "cm"),
                                                               grid_width = unit(2, "mm"),
                                                               title_position = "leftcenter-rot",
                                                               legend_direction = "vertical"),
                                                      FDR = list(title = "FDR",
                                                                 title_gp = gpar(fontsize = 5), 
                                                                 labels_gp = gpar(fontsize = 5),
                                                                 grid_height = unit(0.2, "cm"),
                                                                 grid_width = unit(2, "mm"),
                                                                 title_position = "leftcenter-rot",
                                                                 legend_direction = "vertical")),
                       annotation_name_gp = gpar(fontsize = 5)
)

# clean row names
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

rownames(mat_rppa) <- gsub("_", " ", gsub("_Caution", "", rownames(mat_rppa)))

rownames(mat_lcms) <- gsub("sm", "SM", gsub("tag", "TAG", gsub("\\.", " ", firstup(rownames(mat_lcms)))))


# plot heatmaps
ht_gexp <- Heatmap(mat_gexp,
                   name = "gexp",
                   split = dir_gexp,
                   row_title = c("Sensitivity", "Resistance"),
                   row_title_gp = gpar(fontsize = 7),
                   col = colorRamp2(seq(quantile(mat_gexp, 0.95)*(-1), quantile(mat_gexp, 0.95), length.out = 9), pals::ocean.deep(9)),
                   rect_gp = gpar(col = "white", lwd = unit(0.4, "mm")),
                   top_annotation = ha1, 
                   right_annotation = row_ha,
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_row_dend = FALSE,
                   show_heatmap_legend = T,
                   border = T,
                   border_gp = gpar(col= "black", lwd = unit(0.5, "mm")),
                   heatmap_legend_param = list(title = "Scaled expression",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               at = c(-2,2),
                                               labels = c("Low", "High"),
                                               title_position = "leftcenter-rot",
                                               legend_direction = "vertical",
                                               border = NA),
                   height = unit(10, "cm")
)

ht_score <- Heatmap(mat_score_scaled,
                    name = "score",
                    col = colorRamp2(seq(quantile(mat_score_scaled, 0.95)*(-1), quantile(mat_score_scaled, 0.95), length.out = 9), rev(brewer.pal(9, "YlGnBu"))),
                    rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                    row_names_side = "right",
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 5),
                    show_column_dend = FALSE,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    show_row_dend = FALSE,
                    row_title_gp = gpar(fontsize = 5),
                    show_heatmap_legend = T,
                    border = T,
                    border_gp = gpar(col= "black", lwd = unit(0.5, "mm")),
                    heatmap_legend_param = list(title = "Scaled expression",
                                                title_gp = gpar(fontsize = 5),
                                                labels_gp = gpar(fontsize = 5),
                                                grid_height = unit(0.2, "cm"),
                                                grid_width = unit(2, "mm"),
                                                at = c(-2,2),
                                                labels = c("Low", "High"),
                                                title_position = "leftcenter-rot",
                                                legend_direction = "vertical",
                                                border = NA),
                    height = unit(0.15, "cm")
)


ht_gnab <- Heatmap(mat_gnab,
                   name = "gnab",
                   split = dir_gnab,
                   row_title = c("Sensitivity", "Resistance"),
                   row_title_gp = gpar(fontsize = 7),
                   col = pals::ocean.matter(10)[c(10,6)],
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   show_heatmap_legend = T,
                   border = T,
                   border_gp = gpar(col= "black", lwd = unit(0.5, "mm")),
                   heatmap_legend_param = list(title = "Mutation",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               title_position = "leftcenter-rot",
                                               legend_direction = "vertical"),
                   height = unit(4, "cm")
)

ht_gsea <- Heatmap(mat_gsea,
                   name = "gsea",
                   col = colorRamp2(seq(quantile(mat_gsea, 0.975)*(-1), quantile(mat_gsea, 0.975), length.out = 9), pals::ocean.thermal(9)),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   show_row_names = T,
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   border = T,
                   border_gp = gpar(col= "black", lwd = unit(0.5, "mm")),
                   heatmap_legend_param = list(title = "GSVA Z-score",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               at = c(-2,2),
                                               labels = c("Low", "High"),
                                               title_position = "leftcenter-rot",
                                               legend_direction = "vertical",
                                               border = NA
                   ),
                   height = unit(0.15, "cm")
)




pdf("CCLE_NK_PRISM_heatmap_multiomic_manuscript_small_pctviability_mainfigure.pdf", height = 7, width = 6)
draw(ht_gexp %v% ht_gsea %v% ht_score %v% ht_gnab, heatmap_legend_side = "left",
     annotation_legend_side = "left", merge_legend = T, ht_gap = unit(c(5,2,2,2), "mm"))
dev.off()



# save supplement table (Table S3C)
results_supplement <- results %>% 
  filter(featureA == "N:PRSM:AUC") %>% 
  mutate(datapairs = factor(datapairs, levels = c("PRSM:GEXP", "PRSM:GNAB", "PRSM:METH",
                                                  "PRSM:CNVR", "PRSM:RPPA", "PRSM:MIRN",
                                                  "PRSM:LCMS", "PRSM:GDSC", "PRSM:GDEP"))) %>% 
  mutate(signed_p = sign(cor)*-log10(p)) %>% 
  arrange(datapairs, signed_p) %>% 
  dplyr::select(-test.group, -signed_p, -test.method)

fwrite(results_supplement, "CCLE_NK_PRISM_allcancers_correlations_supplement.txt")

