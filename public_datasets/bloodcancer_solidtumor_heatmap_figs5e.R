
# Expression and methylation of CRISPR screen hits in heme and solid tumors CCLE + TCGA heatmap (Figure S5E)

# load libraries
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(data.table)
library(ggsci)
library(ggrepel)
library(viridis)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(ggnewscale)
library(wesanderson)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# load CCLE feature matrix
fm <- readRDS("data/CCLE.fm_20Q4_complete.rds")

# load CCLE methylation data
meth <- fread("data/CCLE_RRBS_TSS_1kb_20180614.txt", data.table = F)
meth_cpg <- fread("data/CCLE_RRBS_TSS_CpG_clusters_20180614.txt", data.table = F)

# load CCLE annotations
ccle_annot <- fread("data/sample_info.csv")
ccle_annot <- ccle_annot[match(colnames(fm), ccle_annot$DepMap_ID),]

fm <- fm[!is.na(ccle_annot$DepMap_ID),]


ccle_annot <- ccle_annot[ccle_annot$DepMap_ID %in% colnames(fm),]

fm <- fm[,colnames(fm) %in% ccle_annot$DepMap_ID]

# add subtypes
ccle_annot$subtype_short = ifelse(grepl("AML", ccle_annot$Subtype), "AML",
                                  ifelse(grepl("CML", ccle_annot$Subtype), "CML",
                                         ifelse(grepl("ALL), B-cell", ccle_annot$Subtype), "B-ALL",
                              ifelse(grepl("ALL), T-cell", ccle_annot$Subtype), "T-ALL",
                                     ifelse(grepl("Multiple", ccle_annot$Subtype), "MM",
                                            ifelse(grepl("CLL", ccle_annot$Subtype), "CLL",
                                                   ifelse(grepl("Mantle", ccle_annot$Subtype), "MCL",
                                                          ifelse(grepl("B-cell, Hodgkins|^Hodgkins", ccle_annot$Subtype), "CHL",
                                                                 ifelse(grepl("Burkitt", ccle_annot$Subtype), "BL",
                                                                        ifelse(grepl("ALCL", ccle_annot$Subtype), "ALCL",
                                                                               ifelse(grepl("DLBCL", ccle_annot$Subtype), "DLBCL", 
                                                                                      ifelse(grepl("B-cell, Non|^B-cell$", ccle_annot$Subtype), "BCL other",
                                                                                             ifelse(grepl("Cutaneous", ccle_annot$Subtype), "CTCL",
                                                                                                    ifelse(grepl("^T-cell", ccle_annot$Subtype), "TCL other", ccle_annot$Subtype))))))))))))))

ccle_annot$subtype_main <- ccle_annot$subtype_short
ccle_annot$subtype_main[ccle_annot$subtype_short %in% c("DLBCL", "BCL other", "CHL", "BL", "MCL", "CLL")] <- "BCL"
ccle_annot$subtype_main[ccle_annot$subtype_short %in% c("CTCL", "ALCL", "ATL", "TCL other")] <- "TCL"

## add methylation data (CpG and TSS averaged to include more genes)
# combine TSS and CpG methylation data
celllines <- ccle_annot$CCLE_Name[match(colnames(fm), ccle_annot$DepMap_ID)]

# add NA columns for missing cell lines
missing <- celllines[!celllines%in%colnames(meth)]
meth[,missing] <- NA
missing <- celllines[!celllines%in%colnames(meth_cpg)]
meth_cpg[,missing] <- NA

colnames(meth) <- gsub("_name|TSS_|cluster_", "", colnames(meth))
colnames(meth_cpg) <- gsub("_name|TSS_|cluster_", "", colnames(meth_cpg))

colnames(meth[,c("id", "gene", celllines)])==colnames(meth_cpg[,c("id", "gene", celllines)]) # check that column names match before joining
meth_fm <- rbind(meth[,c("id", "gene", celllines)], meth_cpg[,c("id", "gene", celllines)])
meth_fm[,!colnames(meth_fm)%in%c("id", "gene")] <- sapply(meth_fm[,!colnames(meth_fm)%in%c("id", "gene")], as.numeric) # make numeric

mat_meth <- as.matrix(meth_fm[,celllines])
rownames(mat_meth) <- meth_fm$gene
mat_meth <- limma::avereps(mat_meth) # average over methylation values for each gene
mat_meth <- as.data.frame(mat_meth)
rownames(mat_meth) <- paste0("N:METH:", rownames(mat_meth))
colnames(mat_meth) <- colnames(fm)

# replace methylation data in feature matrix with new data
fm <- fm[!grepl("N:METH", rownames(fm)),]
fm <- rbind(fm, mat_meth)

# only cell lines with methylation data
fm <- fm[,!is.na(fm["N:METH:KRAS",])]
ccle_annot <- ccle_annot[match(colnames(fm), ccle_annot$DepMap_ID),]

# load TCGA data
annot <- fread("data/merged_sample_quality_annotations.tsv", data.table = FALSE)

color <-  fread("data/PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv", data.table = FALSE)

gexp_meth_df <- fread("data/tcga_pancancer_expression_methylation_1kb_tss_mean_df.txt", data.table = FALSE)
rownames(gexp_meth_df) <- gexp_meth_df$V1
gexp_meth_df$V1 <- NULL

gexp_meth_df$cancer_type_heme_solid <- "Solid tumor"
gexp_meth_df$cancer_type_heme_solid[gexp_meth_df$cancer_type == "LAML"] <- "AML"
gexp_meth_df$cancer_type_heme_solid[gexp_meth_df$cancer_type == "DLBC"] <- "DLBCL"

# select genes
genelist <- c("CD48", "SPN", "RHOH", "MYB", "SELPLG", "TNFRSF1B", "PVR", "ULBP3")

# matrix of CCLE gene expression medians 
ccle_gexp <- as.data.frame(t(fm[paste0("N:GEXP:", genelist),]))
ccle_gexp_mat <- ccle_gexp %>%
  mutate(cancer_type = ifelse(ccle_annot$primary_disease %in% c("Leukemia", "Lymphoma", "Myeloma"),
                                       ccle_annot$subtype_main, "Solid tumor")) %>% 
  group_by(cancer_type) %>% 
  summarize(across(starts_with("N:GEXP"), median, na.rm = T))

cancer_type <- ccle_gexp_mat$cancer_type
ccle_gexp_mat$cancer_type <- NULL
ccle_gexp_mat_scaled <- t(apply(ccle_gexp_mat, 2, scale))
colnames(ccle_gexp_mat_scaled) <- cancer_type
rownames(ccle_gexp_mat_scaled) <- gsub("N:GEXP:", "", colnames(ccle_gexp_mat))

ccle_gexp_mat_scaled <- ccle_gexp_mat_scaled[,c("AML", "CML", "TCL", "T-ALL", "B-ALL", "BCL", "MM", "Solid tumor")]

ht_ccle_gexp <- Heatmap(ccle_gexp_mat_scaled,
                   name = "Expression CCLE",
                   col = colorRamp2(seq(-2, 2, length.out = 11), rev(brewer.pal(11, "RdBu"))),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   column_names_side = "top",
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                   column_names_gp = gpar(fontsize = 5),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 5),
                   show_heatmap_legend = F,
                   heatmap_legend_param = list(title = "Scaled median\nexpression",
                                               title_gp = gpar(fontsize = 5),
                                               labels_gp = gpar(fontsize = 5),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"))
)

# matrix of CCLE methylation medians 
ccle_meth <- as.data.frame(t(fm[paste0("N:METH:", genelist),]))
ccle_meth_mat <- ccle_meth %>%
  mutate(cancer_type = ifelse(ccle_annot$primary_disease %in% c("Leukemia", "Lymphoma", "Myeloma"),
                              ccle_annot$subtype_main, "Solid tumor")) %>% 
  group_by(cancer_type) %>% 
  summarize(across(starts_with("N:METH"), median, na.rm = T))

cancer_type <- ccle_meth_mat$cancer_type
ccle_meth_mat$cancer_type <- NULL
ccle_meth_mat <- t(ccle_meth_mat)
colnames(ccle_meth_mat) <- cancer_type
rownames(ccle_meth_mat) <- gsub("N:METH:", "", rownames(ccle_meth_mat))

ccle_meth_mat <- ccle_meth_mat[,c("AML", "CML", "TCL", "T-ALL", "B-ALL", "BCL", "MM", "Solid tumor")]

ht_ccle_meth <- Heatmap(ccle_meth_mat,
                        name = "Methylation CCLE",
                        col = colorRamp2(seq(0, 1, length.out = 9), brewer.pal(9, "RdPu")),
                        rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                        column_names_side = "top",
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                        column_names_gp = gpar(fontsize = 5),
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        row_title_gp = gpar(fontsize = 5),
                        show_heatmap_legend = F,
                        heatmap_legend_param = list(title = "Median\nmethylation",
                                                    title_gp = gpar(fontsize = 5),
                                                    labels_gp = gpar(fontsize = 5),
                                                    grid_height = unit(0.2, "cm"),
                                                    grid_width = unit(2, "mm"))
)

# matrix of TCGA gene expression medians 
tcga_gexp <- as.data.frame(gexp_meth_df[,c(paste0(genelist, "_gexp"), "cancer_type_heme_solid")])
tcga_gexp_mat <- tcga_gexp %>%
  group_by(cancer_type_heme_solid) %>% 
  summarize(across(ends_with("_gexp"), median, na.rm = T))

cancer_type <- tcga_gexp_mat$cancer_type_heme_solid
tcga_gexp_mat$cancer_type_heme_solid <- NULL
tcga_gexp_mat_scaled <- t(apply(tcga_gexp_mat, 2, scale))
colnames(tcga_gexp_mat_scaled) <- cancer_type
rownames(tcga_gexp_mat_scaled) <- gsub("_gexp", "", colnames(tcga_gexp_mat))

ht_tcga_gexp <- Heatmap(tcga_gexp_mat_scaled,
                        name = "Expression TCGA",
                        col = colorRamp2(seq(-1, 1, length.out = 11), pals::ocean.deep(11)),
                        rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                        show_row_names = F,
                        column_names_side = "top",
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                        column_names_gp = gpar(fontsize = 5),
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        row_title_gp = gpar(fontsize = 5),
                        show_heatmap_legend = T,
                        heatmap_legend_param = list(title = "Scaled median\nexpression",
                                                    title_gp = gpar(fontsize = 5),
                                                    labels_gp = gpar(fontsize = 5),
                                                    grid_height = unit(0.2, "cm"),
                                                    grid_width = unit(2, "mm"))
)


# matrix of TCGA methylation medians 
tcga_meth <- as.data.frame(gexp_meth_df[,c(paste0(genelist, "_meth"), "cancer_type_heme_solid")])
tcga_meth_mat <- tcga_meth %>%
  group_by(cancer_type_heme_solid) %>% 
  summarize(across(ends_with("_meth"), median, na.rm = T))

cancer_type <- tcga_meth_mat$cancer_type_heme_solid
tcga_meth_mat$cancer_type_heme_solid <- NULL
tcga_meth_mat <- t(tcga_meth_mat)
colnames(tcga_meth_mat) <- cancer_type
rownames(tcga_meth_mat) <- gsub("_meth", "", rownames(tcga_meth_mat))

ht_tcga_meth <- Heatmap(tcga_meth_mat,
                        name = "Methylation TCGA",
                        col = colorRamp2(seq(0, 1, length.out = 9), brewer.pal(9, "RdPu")),
                        rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                        column_names_side = "top",
                        show_row_names = F,
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                        column_names_gp = gpar(fontsize = 5),
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        row_title_gp = gpar(fontsize = 5),
                        show_heatmap_legend = F,
                        heatmap_legend_param = list(title = "Median\nmethylation",
                                                    title_gp = gpar(fontsize = 5),
                                                    labels_gp = gpar(fontsize = 5),
                                                    grid_height = unit(0.2, "cm"),
                                                    grid_width = unit(2, "mm"))
)

## --------------------------------------------------------------------------------------

ht_list_ccle <- ht_ccle_gexp %v% ht_ccle_meth
ht_list_tcga <- ht_tcga_gexp %v% ht_tcga_meth

dir.create("heme_solid_comparison")

pdf("heme_solid_comparison/ccle_heatmap.pdf", height = 2, width = 1.3)
draw(ht_list_ccle)
dev.off()

pdf("heme_solid_comparison/tcga_heatmap.pdf", height = 2, width = 0.95)
draw(ht_list_tcga)
dev.off()

## --------------------------------------------------------------------------------------

# all solid tumor subtypes shown

# matrix of CCLE gene expression medians 
ccle_gexp <- as.data.frame(t(fm[paste0("N:GEXP:", genelist),]))
ccle_gexp_mat <- ccle_gexp %>%
  mutate(cancer_type = ifelse(ccle_annot$primary_disease %in% c("Leukemia", "Lymphoma", "Myeloma"),
                              ccle_annot$subtype_main, gsub(" Cancer", "", ccle_annot$primary_disease))) %>% 
  group_by(cancer_type) %>% 
  summarize(across(starts_with("N:GEXP"), median, na.rm = T))

cancer_type <- ccle_gexp_mat$cancer_type
ccle_gexp_mat$cancer_type <- NULL
ccle_gexp_mat_scaled <- t(apply(ccle_gexp_mat, 2, scale))
colnames(ccle_gexp_mat_scaled) <- cancer_type
rownames(ccle_gexp_mat_scaled) <- gsub("N:GEXP:", "", colnames(ccle_gexp_mat))

solid_tumors <- gsub(" Cancer", "", unique(ccle_annot$primary_disease)[grepl("Cancer", unique(ccle_annot$primary_disease))])
ccle_gexp_mat_scaled <- ccle_gexp_mat_scaled[,c("AML", "CML", "TCL", "T-ALL", "B-ALL", "BCL", "MM", solid_tumors)]

rownames(ccle_gexp_mat) <- cancer_type
ccle_gexp_mat <- t(ccle_gexp_mat)
ccle_gexp_mat <- ccle_gexp_mat[,c("AML", "CML", "TCL", "T-ALL", "B-ALL", "BCL", "MM", solid_tumors)]
rownames(ccle_gexp_mat) <- gsub("N:GEXP:", "", rownames(ccle_gexp_mat))

ht_ccle_gexp <- Heatmap(ccle_gexp_mat,
                        name = "Expression CCLE",
                        col = colorRamp2(seq(0, 8, length.out = 11), pals::ocean.deep(11)),
                        rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                        column_names_side = "top",
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                        column_names_gp = gpar(fontsize = 7),
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        row_title_gp = gpar(fontsize = 5),
                        heatmap_legend_param = list(title = "Median\nexpression\n(log2)",
                                                    title_gp = gpar(fontsize = 7),
                                                    labels_gp = gpar(fontsize = 7),
                                                    grid_height = unit(0.2, "cm"),
                                                    grid_width = unit(3, "mm"))
)

# matrix of CCLE methylation medians 
ccle_meth <- as.data.frame(t(fm[paste0("N:METH:", genelist),]))
ccle_meth_mat <- ccle_meth %>%
  mutate(cancer_type = ifelse(ccle_annot$primary_disease %in% c("Leukemia", "Lymphoma", "Myeloma"),
                              ccle_annot$subtype_main, gsub(" Cancer", "", ccle_annot$primary_disease))) %>% 
  group_by(cancer_type) %>% 
  summarize(across(starts_with("N:METH"), median, na.rm = T))

cancer_type <- ccle_meth_mat$cancer_type
ccle_meth_mat$cancer_type <- NULL
ccle_meth_mat <- t(ccle_meth_mat)
colnames(ccle_meth_mat) <- cancer_type
rownames(ccle_meth_mat) <- gsub("N:METH:", "", rownames(ccle_meth_mat))

ccle_meth_mat <- ccle_meth_mat[,c("AML", "CML", "TCL", "T-ALL", "B-ALL", "BCL", "MM", solid_tumors)]

ht_ccle_meth <- Heatmap(ccle_meth_mat,
                        name = "Methylation CCLE",
                        col = colorRamp2(seq(0, 1, length.out = 9), brewer.pal(9, "RdPu")),
                        rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                        column_names_side = "top",
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                        column_names_gp = gpar(fontsize = 7),
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        row_title_gp = gpar(fontsize = 5),
                        heatmap_legend_param = list(title = "Median\nmethylation",
                                                    title_gp = gpar(fontsize = 7),
                                                    labels_gp = gpar(fontsize = 7),
                                                    grid_height = unit(0.2, "cm"),
                                                    grid_width = unit(3, "mm"))
)


ht_list_ccle <- ht_ccle_gexp %v% ht_ccle_meth

pdf("heme_solid_comparison/ccle_heatmap_allsolidtumors.pdf", height = 2.75, width = 4.25)
draw(ht_list_ccle)
decorate_heatmap_body("Expression CCLE", {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.5))
  
})
decorate_heatmap_body("Methylation CCLE", {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.5))
  
})
dev.off()

