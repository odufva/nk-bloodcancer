
# Plot selected correlates in BCL NK PRISM sensitivity in CCLE dataset (Figure 5H)


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
library(GSVA)

# load feature matrix
fm <- readRDS("data/CCLE.fm_20Q4_complete.rds")

# load correlations
results <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_BCL.tsv", data.table = F)

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
ccle_annot <- fread("data/sample_info.csv")

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


fm_prism_clean <- fm_prism_clean[,subtype_short %in% c("MCL", "CHL", "BL", "DLBCL", "BCL other")]

# data frame with subtype for heatmap annotation
annot <- data.frame(cancertype = subtype_short[subtype_short %in% c("MCL", "CHL", "BL", "DLBCL", "BCL other")])
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
rownames(gm1)="N:GEXP:HLAIScore"

fm_prism_clean <- rbind(fm_prism_clean, gm1)

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
      mat <- fm_prism_clean[paste0("B:GNAB:", GENES),]
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
mat_gexp <- heatmap_scaled_matrix(DATAPAIR = "GEXP", ETRATIO = "AUC", GENES = c("NCR3LG1", "HLA-C", "NLRC5", "TAP1", "NFKBIA", "MYB", "HLAIScore"))
mat_gnab <- heatmap_scaled_matrix(DATAPAIR = "GNAB", ETRATIO = "AUC", SCALE = F, GENES = c("GNA13"))

rownames(mat_gexp) <- gsub("HLAIScore", "HLA I score", rownames(mat_gexp))

# GSEA

ifn <- fread("../../../cropseq_olli/results/combine/ifng_score.txt", data.table = F)$gene
nfkb <- fread("../../../cropseq_olli/results/combine/nfkb_score.txt", data.table = F)$gene

core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")


genesets = list(SCRNASEQ_NK_RESPONSE = core)

gsea <- fread("results_gsea/NK_PRISM_MM_HALLMARK_GSEA.txt", data.table = F)

pathways_hallmark <- fgsea::gmtPathways("/Users/odufva/Documents/Labra/Scripts_programs/msigdb/msigdb_v7.0_GMTs/h.all.v7.0.symbols.gmt")

fm_gsva <- fm_prism_clean[grepl("N:GEXP", rownames(fm_prism_clean)),]
rownames(fm_gsva) <- gsub("N:GEXP:", "", rownames(fm_gsva))

# run GSVA
gsva_results <- gsva(as.matrix(fm_gsva), c(pathways_hallmark, genesets))

genesets_top <- gsea %>% filter(grepl("E2F|G2M|MITOTIC|TNFA|GAMMA", pathway)) %>% select(pathway) %>% tibble::deframe()
mat_gsea <- t(apply(gsva_results[c("SCRNASEQ_NK_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),], 1, scale))
rownames(mat_gsea) <- gsub("Il6 jak stat3", "IL6-JAK-STAT3",
                           gsub("Il2 stat5", "IL2-STAT5",
                                gsub("Dna", "DNA",
                                     gsub("Myc", "MYC",
                                          gsub("E2f", "E2F",
                                               gsub("Kras", "KRAS",
                                                    gsub("Tnfa", "TNFA",
                                                         gsub("G2m", "G2M",
                                                              gsub("nfkb|Tnfa signaling via nfkb", "NF-kB",
                                                                   gsub("Interferon gamma response", "IFNy",
                                                                        gsub("Nk", "NK",
                                                                             stringr::str_to_sentence(gsub("_", " ",
                                                                                                           sub(".*?_", "", rownames(mat_gsea)))))))))))))))
colnames(mat_gsea) <- colnames(gsva_results)


# create heatmap annotations

prism_auc <- as.numeric(scale(as.numeric(fm_prism_clean["N:PRSM:AUC",])))

subtypes <- read_excel("../cell_line_subtypes.xlsx")
subtypes$cancer[subtypes$cancer %in% "DLBCL"] <- "DLBCL other"
subtypes$cancer[subtypes$subtype %in% "GCB"] <- "GCB DLBCL"
subtypes$cancer[subtypes$subtype %in% "ABC"] <- "ABC DLBCL"
subtypes$cancer[subtypes$cancer %in% "BCL other"] <- "Other"

annot_subtype <- subtypes[match(colnames(mat_gexp), subtypes$cell_line),]
annot <- data.frame(`PRISM AUC` = prism_auc,
                    Subtype = annot_subtype[,2])
colnames(annot) <- c("PRISM AUC", "Subtype")

ha1 <- HeatmapAnnotation(df = annot, col = list(`PRISM AUC` = colorRamp2(seq(-2, 2, length.out = 20), viridis_pal(direction = -1)(20)),
                                                Subtype = structure(c(brewer.pal(9, "Set1")[c(5,4,3)], "#CFA42D", "#B25C3F", "#999999", "grey80"),
                                                                    names = c("ABC DLBCL", "GCB DLBCL", "DLBCL other", "MCL", "BL", "CHL", "Other"))
),
gap = unit(0.1, "mm"),
annotation_legend_param = list(Subtype = list(title = "Subtype", title_gp = gpar(fontsize = 4), 
                                              title_position = "topcenter",
                                              labels_gp = gpar(fontsize = 4), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"),
                                              legend_direction = "horizontal"),
                               `PRISM AUC` = list(title = "PRISM AUC\nZ-score", title_gp = gpar(fontsize = 4), 
                                                  title_position = "topcenter",
                                                  labels_gp = gpar(fontsize = 4), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"),
                                                  legend_direction = "horizontal")
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
                   column_names_gp = gpar(fontsize = 4),
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

ht_gnab <- Heatmap(mat_gnab,
                   name = "gnab",
                   col = c("1" = "black", "0" = "grey90"),
                   rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 4, fontface = "italic"),
                   column_names_gp = gpar(fontsize = 4),
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

ht_gsea <- Heatmap(mat_gsea,
                   name = "gsea",
                   col = colorRamp2(seq(quantile(mat_gsea, 0.975)*(-1), quantile(mat_gsea, 0.975), length.out = 9), pals::ocean.thermal(9)),
                   show_row_names = T,
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize = 4),
                   column_names_gp = gpar(fontsize = 3),
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   show_row_dend = FALSE,
                   row_title_gp = gpar(fontsize = 4),
                   show_heatmap_legend = T,
                   heatmap_legend_param = list(title = "GSVA Z-score",
                                               title_gp = gpar(fontsize = 4),
                                               labels_gp = gpar(fontsize = 4),
                                               grid_height = unit(0.2, "cm"),
                                               grid_width = unit(2, "mm"),
                                               title_position = "topcenter",
                                               legend_direction = "horizontal"))

ht_list <- ht_gexp %v% ht_gnab %v% ht_gsea

pdf("CCLE_NK_PRISM_heatmap_multiomic_BCL_manuscript_selected.pdf", height = 4, width = 1.2)
draw(ht_list, ht_gap = unit(0.5, "mm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

