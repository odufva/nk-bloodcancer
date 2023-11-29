
# Compute pairwise correlations for NK PRISM sensitivity in CCLE dataset by subtype

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
fm_prism$rownames <- rownames(fm_full)

# "genelist" of NK PRISM features
genelist = c("D1_NK_0.625", "D1_NK_1.25", "D1_NK_2.5", "D1_NK_5", "AUC")

# correlations by subtype

Subtype <- nk_prism$Subtype[match(colnames(fm_prism), nk_prism$DepMap_ID)]
subtype_short <- ifelse(grepl("AML", Subtype), "AML",
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
subtype_grouped <- ifelse(subtype_short %in% c("DLBCL", "BCL other", "BL", "CHL", "MCL"), "BCL", subtype_short)

correlations_subtype <- function(SUBTYPE){
  
  fm_prism_subtype <- fm_prism[,subtype_grouped %in% c(SUBTYPE)]
  
  extrafeatures=c(grep("^B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^N:RPPA|^N:LCMS|^N:DRUG|^N:METH|^N:GDSC|^N:MIRN|^N:GDEP|^N:GEXP", rownames(fm_prism_subtype), value=T))
  l.regulon.gene=regulon.feats(fm_prism_subtype, genelist)#, cnv_annot)
  results_subtype=pairwise.correlation(l.regulon.gene, fm_prism_subtype, extrafeatures, filter.val = 3, cores = 2, adjust.method = "BH")
  
  fwrite(results_subtype, paste0("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_", SUBTYPE, ".tsv"), sep ="\t")
  
  heatmap_subtype <- function(ET_RATIO, DATA, RESULTS){
    
    # order fm by PRISM NK sensitivity
    data <- DATA
    fm_prism_clean <- data[,!grepl("sum|var|rownames", colnames(data))]
    order <- order(fm_prism_clean[paste0("N:PRSM:", ET_RATIO),])
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
    
    DATAPAIR = "GEXP"
    ETRATIO = ET_RATIO
    PVAL = 1
    ADJP = 1
    # function to create scaled matrices with selected alterations
    heatmap_scaled_matrix <- function(DATAPAIR, ETRATIO, PVAL = 1, ADJP = 1, SCALE = T){
      
      # filter features to plot
      feats_top <- results_subtype %>%
        filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
        arrange(p) %>%
        filter(cor > 0) %>%
        top_n(n = -50, wt = p) %>%
        arrange(desc(cor)) %>%
        dplyr::select(featureB) %>%
        deframe()
      
      feats_bottom <- results_subtype %>%
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
    mat_gexp <- heatmap_scaled_matrix(DATAPAIR = "GEXP", ETRATIO = ET_RATIO)
    mat_rppa <- heatmap_scaled_matrix(DATAPAIR = "RPPA", ETRATIO = ET_RATIO)
    mat_gdsc <- heatmap_scaled_matrix(DATAPAIR = "GDSC", ETRATIO = ET_RATIO)
    mat_mirn <- heatmap_scaled_matrix(DATAPAIR = "MIRN", ETRATIO = ET_RATIO)
    mat_lcms <- heatmap_scaled_matrix(DATAPAIR = "LCMS", ETRATIO = ET_RATIO)
    mat_gdep <- heatmap_scaled_matrix(DATAPAIR = "GDEP", ETRATIO = ET_RATIO, SCALE = F)
    mat_gnab <- heatmap_scaled_matrix(DATAPAIR = "GNAB", ETRATIO = ET_RATIO, SCALE = F)
    mat_cnvr <- heatmap_scaled_matrix(DATAPAIR = "CNVR", ETRATIO = ET_RATIO, SCALE = F)
    mat_meth <- heatmap_scaled_matrix(DATAPAIR = "METH", ETRATIO = ET_RATIO, SCALE = F)
    
    
    
    # create heatmap annotations
    ha1 <- HeatmapAnnotation(df = annot, col = list(subtype = structure(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(annot$subtype))), 
                                                                        names = as.character(unique(annot$subtype)))),
                             
                             PRISM = anno_barplot(as.numeric(fm_prism_clean[paste0("N:PRSM:", ET_RATIO),]), 
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
    
    
    
    # print heatmaps
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_GEXP_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_gexp)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_RPPA_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_rppa)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_GDSC_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_gdsc)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_MIRN_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_mirn)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_LCMS_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_lcms)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_GDEP_heatmap_", SUBTYPE, "_",ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_gdep)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_GNAB_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_gnab)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_CNVR_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 14, width = 7)
    draw(ht_cnvr)
    dev.off()
    
    pdf(paste0("results_20Q4_complete_pctviability/NK_PRISM_METH_heatmap_", SUBTYPE, "_", ET_RATIO, ".pdf"), height = 10, width = 7)
    draw(ht_meth)
    dev.off()
    
  }
  
  lapply(c("D1_NK_0.625", "D1_NK_1.25", "D1_NK_2.5", "D1_NK_5", "AUC"), heatmap_subtype, DATA = fm_prism_subtype)
  
}

lapply(c("MM", "B-ALL", "AML", "T-ALL", "BCL"), correlations_subtype)
