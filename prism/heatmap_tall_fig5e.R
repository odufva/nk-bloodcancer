
# Plot selected correlates in T-ALL NK PRISM sensitivity in CCLE dataset (Figure 5E)

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
library(viridis)

# load feature matrix
fm <- readRDS("data/CCLE.fm_20Q4_complete.rds")

# load methylation data
meth <- fread("data/CCLE_RRBS_TSS_1kb_20180614.txt", data.table = F)
meth_cpg <- fread("data/CCLE_RRBS_TSS_CpG_clusters_20180614.txt", data.table = F)

# load CCLE sample info
ccle_annot <- fread("data/sample_info.csv")

# load correlations
results <- fread("NK_PRISM_heme_correlations_T-ALL.tsv", data.table = F)

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

## -----------------------------------

# heatmaps

# order fm by PRISM NK sensitivity
fm_prism_clean <- fm_prism[,!grepl("sum|var|rownames", colnames(fm_prism))]
order <- order(as.numeric(fm_prism_clean["N:PRSM:AUC",]))
fm_prism_clean <- fm_prism_clean[,order]

## add methylation data (CpG and TSS averaged to include more genes)
# combine TSS and CpG methylation data
celllines <- ccle_annot$CCLE_Name[match(colnames(fm_prism_clean), ccle_annot$DepMap_ID)]

# add NA columns for missing cell lines
missing <- celllines[!celllines%in%colnames(meth)]
meth[,missing] <- NA
missing <- celllines[!celllines%in%colnames(meth_cpg)]
meth_cpg[,missing] <- NA

colnames(meth) <- gsub("_name|TSS_|cluster_", "", colnames(meth_cpg))
colnames(meth_cpg) <- gsub("_name|TSS_|cluster_", "", colnames(meth_cpg))

colnames(meth[,c("id", "gene", celllines)])==colnames(meth_cpg[,c("id", "gene", celllines)]) # check that column names match before joining
meth_prism <- rbind(meth[,c("id", "gene", celllines)], meth_cpg[,c("id", "gene", celllines)])
meth_prism[,!colnames(meth_prism)%in%c("id", "gene")] <- sapply(meth_prism[,!colnames(meth_prism)%in%c("id", "gene")], as.numeric) # make numeric

mat_meth <- as.matrix(meth_prism[,celllines])
rownames(mat_meth) <- meth_prism$gene
mat_meth <- limma::avereps(mat_meth) # average over methylation values for each gene
mat_meth <- as.data.frame(mat_meth)
rownames(mat_meth) <- paste0("N:METH:", rownames(mat_meth))
colnames(mat_meth) <- colnames(fm_prism_clean)

# replace methylation data in feature matrix with new data
fm_prism_clean <- fm_prism_clean[!grepl("N:METH", rownames(fm_prism_clean)),]
fm_prism_clean <- rbind(fm_prism_clean, mat_meth)

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


fm_prism_clean <- fm_prism_clean[,subtype_short %in% c("T-ALL")]

# data frame with subtype for heatmap annotation
annot <- data.frame(cancertype = subtype_short[subtype_short %in% c("T-ALL")])
colnames(annot) <- "Cancer type"

# replace DepMap IDs with cell line names
colnames(fm_prism_clean) <- ccle_annot$stripped_cell_line_name[match(colnames(fm_prism_clean), ccle_annot$DepMap_ID)]


# function to create scaled matrices with selected alterations
heatmap_scaled_matrix <- function(DATAPAIR, ETRATIO, PVAL = 1, ADJP = 1, SCALE = T, TOPN = 50, GENES = NULL){
  
  if (is.null(GENES)) {
    
    # filter features to plot
    feats_top <- results %>%
      filter(datapairs == paste0("PRSM:", DATAPAIR) & featureA == paste0("N:PRSM:", ETRATIO) & p <= PVAL & adj.p <= ADJP) %>%
      arrange(p) %>%
      filter(cor > 0) %>%
      top_n(n = -TOPN, wt = p) %>%
      arrange(desc(cor)) %>%
      dplyr::select(featureB) %>%
      deframe()
    
    
    # select features from matrix
    mat <- fm_prism_clean[feats_top,]
    
  }
  
  else
  {
    if (DATAPAIR == "GEXP") {
      # selected genes
      mat <- fm_prism_clean[paste0("N:GEXP:", GENES),]
    }
    
    else {
      
      mat <- fm_prism_clean[paste0("N:METH:", GENES),]
    }
    
  }
  
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
mat_gexp <- heatmap_scaled_matrix(DATAPAIR = "GEXP", ETRATIO = "AUC", GENES = c("FAS", "PVR", "ULBP1", "CD44"))
mat_meth <- heatmap_scaled_matrix(DATAPAIR = "METH", ETRATIO = "AUC", SCALE = F, GENES = c("FAS", "PVR", "ULBP1", "CD44"))


# create heatmap annotations

prism_auc <- as.numeric(scale(as.numeric(fm_prism_clean["N:PRSM:AUC",])))

subtypes <- read_excel("../cell_line_subtypes.xlsx")
subtypes$subtype <- gsub("TAL1/LMO.", "TAL1",
                         gsub("TRA/D-", "", subtypes$subtype))
subtypes$subtype[is.na(subtypes$subtype)] <- "Other"
subtypes$subtype[subtypes$subtype %in% c("TCL1", "NA")] <- "Other"
annot <- subtypes[match(colnames(mat_gexp), subtypes$cell_line),]
annot <- data.frame(`PRISM AUC` = prism_auc,
                    Subtype = annot[,3])
colnames(annot) <- c("PRISM AUC", "Subtype")
annot$Subtype <- factor(annot$Subtype, levels = c("TAL1", "LMO2", "TLX1", "TLX3", "Other"))

ha1 <- HeatmapAnnotation(df = annot, col = list(`PRISM AUC` = colorRamp2(seq(-2, 2, length.out = 20), viridis_pal(direction = -1)(20)),
                                                Subtype = structure(c(colorRampPalette(brewer.pal(11, "Spectral")[c(1:4, 8:11)])(8)[c(1,3,5,6)], "grey80"), 
                                                                    names = c("TAL1", "LMO2", "TLX1", "TLX3", "Other"))
),
gap = unit(0.1, "mm"),
annotation_legend_param = list(Subtype = list(title = "Subtype", title_gp = gpar(fontsize = 4), 
                                              title_position = "topcenter",
                                              labels_gp = gpar(fontsize = 4), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"),
                                              legend_direction = "horizontal"),
                               `PRISM AUC` = list(title = "PRISM AUC\nZ-score", title_gp = gpar(fontsize = 4), 
                                                  title_position = "topcenter",
                                                  labels_gp = gpar(fontsize = 4), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"),
                                                  legend_direction = "horizontal")#,
),
annotation_name_gp = gpar(fontsize = 4),
height = unit(0.4, "cm"),
simple_anno_size_adjust = T,
show_annotation_name = T
)


# plot heatmaps
ht_gexp <- Heatmap(mat_gexp,
                   name = "heatmap",
                   col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   top_annotation = ha1, 
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                   column_names_gp = gpar(fontsize = 3),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 4),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Expression (log2)\nZ-score",
                                               title_gp = gpar(fontsize = 4),
                                               title_position = "topcenter",
                                               labels_gp = gpar(fontsize = 4),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               legend_direction = "horizontal")
)

ht_meth <- Heatmap(mat_meth,
                   name = "meth",
                   col = colorRamp2(seq(0, 1, length.out = 9), brewer.pal(9, "RdPu")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                   column_names_gp = gpar(fontsize = 3),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 4),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Methylation",
                                               title_gp = gpar(fontsize = 4),
                                               title_position = "topcenter",
                                               labels_gp = gpar(fontsize = 4),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               legend_direction = "horizontal")
)

ht_gnab <- Heatmap(mat_gnab,
                   name = "gnab",
                   col = c("1" = "black", "0" = "grey90"),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                   column_names_gp = gpar(fontsize = 3),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 4),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "Mutation",
                                               title_gp = gpar(fontsize = 4),
                                               title_position = "topcenter",
                                               labels_gp = gpar(fontsize = 4),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               legend_direction = "horizontal")
)

ht_list <- ht_gexp %v% ht_meth

pdf("CCLE_NK_PRISM_heatmap_multiomic_TALL.pdf", height = 3.25, width = 1.1)
draw(ht_list, ht_gap = unit(0.5, "mm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

