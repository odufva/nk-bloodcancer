
# Compute pairwise correlations for NK PRISM sensitivity in CCLE dataset (scripts by Petri Pölönen)


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

# add NK cell sensitivity PRISM data to feature matrix
nk_prism <- fread("data/NK_PRISM_heme_pctviability.txt")
nk_prism_fm <- nk_prism %>%
  dplyr::select(Condition, pctviability, DepMap_ID) %>%
  pivot_wider(names_from = DepMap_ID, values_from = pctviability) %>%
  mutate(Condition = paste0("N:PRSM:", Condition))

nk_prism_auc <- fread("data/NK_PRISM_heme_pctviability_auc_norm.txt")
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
fm_prism$rownames <- NULL


# "genelist" of NK PRISM features
genelist = c("D1_NK_0.625", "D1_NK_1.25", "D1_NK_2.5", "D1_NK_5", "AUC")

# choose feature types for correlations
extrafeatures=c(grep("^B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^N:RPPA|^N:LCMS|^N:DRUG|^N:METH|^N:GDSC|^N:MIRN|^N:GDEP|^N:GEXP", rownames(fm_prism), value=T))

l.regulon.gene=regulon.feats(fm_prism, genelist)#, cnv_annot)

# run correlations
results=pairwise.correlation(l.regulon.gene, fm_prism, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

dir.create("results_20Q4_complete")
fwrite(results, "results_20Q4_complete/NK_PRISM_heme_correlations_all.tsv", sep ="\t")


# filtered correlations

# get sums and variances for features and specify variance qunatile cutoffs
fm_prism$sum <- rowSums(fm_prism[,1:63])
fm_prism$var <- rowVars(as.matrix(fm_prism[,1:63]), na.rm = T)
gexp_vars <- fm_prism[grepl("N:GEXP", rownames(fm_prism)),"var"]
gexp_var_cutoff <- quantile(gexp_vars, 0.75, na.rm = T)
mirn_vars <- fm_prism[grepl("N:MIRN", rownames(fm_prism)),"var"]
mirn_var_cutoff <- quantile(mirn_vars, 0.75, na.rm = T)
meth_vars <- fm_prism[grepl("N:METH", rownames(fm_prism)),"var"]
meth_var_cutoff <- quantile(meth_vars, 0.75, na.rm = T)
gdep_vars <- fm_prism[grepl("N:GDEP", rownames(fm_prism)),"var"]
gdep_var_cutoff <- quantile(gdep_vars, 0.75, na.rm = T)


# filter feature matrix to remove low variance features
fm_prism_filtered <- fm_prism %>%
  filter(!(grepl("N:GEXP", rownames) & var < gexp_var_cutoff)) %>%
  filter(!(grepl("N:MIRN", rownames) & var < mirn_var_cutoff)) %>%
  filter(!(grepl("N:METH", rownames) & var < meth_var_cutoff)) %>%
  filter(!(grepl("N:GDEP", rownames) & var < gdep_var_cutoff))


extrafeatures=c(grep("^B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^N:RPPA|^N:LCMS|^N:DRUG|^N:METH|^N:GDSC|^N:MIRN|^N:GDEP|^N:GEXP", rownames(fm_prism_filtered), value=T))
l.regulon.gene=regulon.feats(fm_prism_filtered, genelist)
results_filtered=pairwise.correlation(l.regulon.gene, fm_prism_filtered, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

fwrite(results_filtered, "results_20Q4_complete/NK_PRISM_heme_correlations_filtered.tsv", sep ="\t")

## -----------------------------------

# heatmaps

# order fm by PRISM NK sensitivity
fm_prism_clean <- fm_prism[,!grepl("sum|var|rownames", colnames(fm_prism))]
order <- order(fm_prism_clean["N:PRSM:AUC",])
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
annot <- data.frame(subtype = subtype_short)

# replace DepMap IDs with cell line names
colnames(fm_prism_clean) <- ccle_annot$stripped_cell_line_name[match(colnames(fm_prism_clean), ccle_annot$DepMap_ID)]


# function to create scaled matrices with selected alterations
heatmap_scaled_matrix <- function(DATAPAIR, ETRATIO, PVAL = 1, ADJP = 1, SCALE = T){
  
  # filter features to plot
  feats_top <- results %>%
    filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
    arrange(p) %>%
    filter(cor > 0) %>%
    top_n(n = -50, wt = p) %>%
    arrange(desc(cor)) %>%
    dplyr::select(featureB) %>%
    deframe()
  
  feats_bottom <- results %>%
    filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
    arrange(p) %>%
    filter(cor < 0) %>%
    top_n(n = -50, wt = p) %>%
    arrange(desc(cor)) %>%
    dplyr::select(featureB) %>%
    deframe()
  
  # select features from matrix
  mat <- fm_prism_clean[c(feats_top, feats_bottom),]
  
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
mat_gexp <- heatmap_scaled_matrix(DATAPAIR = "GEXP", ETRATIO = "AUC")
mat_rppa <- heatmap_scaled_matrix(DATAPAIR = "RPPA", ETRATIO = "AUC")
mat_gdsc <- heatmap_scaled_matrix(DATAPAIR = "GDSC", ETRATIO = "AUC")
mat_mirn <- heatmap_scaled_matrix(DATAPAIR = "MIRN", ETRATIO = "AUC")
mat_lcms <- heatmap_scaled_matrix(DATAPAIR = "LCMS", ETRATIO = "AUC")
mat_gdep <- heatmap_scaled_matrix(DATAPAIR = "GDEP", ETRATIO = "AUC", SCALE = F)
mat_gnab <- heatmap_scaled_matrix(DATAPAIR = "GNAB", ETRATIO = "AUC", SCALE = F)
mat_cnvr <- heatmap_scaled_matrix(DATAPAIR = "CNVR", ETRATIO = "AUC", SCALE = F)
mat_meth <- heatmap_scaled_matrix(DATAPAIR = "METH", ETRATIO = "AUC", SCALE = F)
mat_samp <- heatmap_scaled_matrix(DATAPAIR = "SAMP", ETRATIO = "AUC", SCALE = F)


# create heatmap annotations
ha1 <- HeatmapAnnotation(df = annot, col = list(subtype = structure(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(annot$subtype))), 
                                                                    names = as.character(unique(annot$subtype)))),
                         
                         PRISM = anno_barplot(as.numeric(fm_prism_clean["N:PRSM:AUC",]), 
                                              bar_width = 0.75, 
                                              border = FALSE, 
                                              axis = TRUE,
                                              axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                              gp = gpar(col = NA, fill = "grey50")
                         ),
                         gap = unit(0.75, "mm"),
                         annotation_legend_param = list(subtype = list(title = "Subtype", title_gp = gpar(fontsize = 5), 
                                                                       labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))),
                         height = unit(1.5, "cm"),
                         simple_anno_size_adjust = T,
                         show_annotation_name = T
)


# plot heatmaps
ht_gexp <- Heatmap(mat_gexp,
                   name = "heatmap",
                   col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Expression (log2)\nZ-score",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_rppa <- Heatmap(mat_rppa,
                   name = "rppa",
                   col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Z-score",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_gdsc <- Heatmap(mat_gdsc,
                   name = "gdsc",
                   col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Z-score",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_mirn <- Heatmap(mat_mirn,
                   name = "mirn",
                   col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Z-score",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_lcms <- Heatmap(mat_lcms,
                   name = "lcms",
                   col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Z-score",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_gdep <- Heatmap(mat_gdep,
                   name = "gdep",
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Dependency",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_gnab <- Heatmap(mat_gnab,
                   name = "gnab",
                   col = c("1" = "black", "0" = "grey80"),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Mutation",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_cnvr <- Heatmap(mat_cnvr,
                   name = "cnvr",
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Copy number",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_meth <- Heatmap(mat_meth,
                   name = "meth",
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Methylation",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

ht_samp <- Heatmap(mat_samp,
                   name = "heatmap",
                   col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "GSVA score",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)



# print heatmaps
pdf("results_20Q4_complete/NK_PRISM_GEXP_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_gexp)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_RPPA_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_rppa)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_GDSC_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_gdsc)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_MIRN_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_mirn)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_LCMS_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_lcms)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_GDEP_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_gdep)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_GNAB_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_gnab)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_CNVR_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_cnvr)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_METH_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_meth)
dev.off()

pdf("results_20Q4_complete/NK_PRISM_SAMP_heatmap_AUC.pdf", height = 7, width = 7)
draw(ht_samp)
dev.off()

