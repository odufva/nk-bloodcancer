
# Expression and methylation of genes correlated with PRISM AUC scatter plots (Figures 3F-I and S5I)

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
library(patchwork)
library(tidyr)
library(GSVA)

# load scripts
source("compute.pairwise.R")
source("functions_generate_fm.R")
source("functions_statistics.R")
source("regulon.feats.R")

# load feature matrix
fm <- readRDS("CCLE.fm_20Q4_complete.rds")

# load methylation data
meth <- fread("data/CCLE_RRBS_TSS_1kb_20180614.txt", data.table = F)
meth_cpg <- fread("data/CCLE_ligand_methylation/CCLE_RRBS_TSS_CpG_clusters_20180614.txt", data.table = F)

# load annotations
ccle_annot <- fread("../NK_PRISM/sample_info.csv")
ccle_annot <- ccle_annot[match(colnames(fm), ccle_annot$DepMap_ID),]

fm <- fm[,!is.na(ccle_annot$DepMap_ID),]
ccle_annot <- ccle_annot[ccle_annot$DepMap_ID %in% colnames(fm),]

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

# subset to PRISM cell lines

# add NK cell sensitivity PRISM data to feature matrix
nk_prism <- fread("../NK_PRISM/NK_PRISM_heme_pctviability.txt")
nk_prism_fm <- nk_prism %>%
  dplyr::select(Condition, pctviability, DepMap_ID) %>%
  pivot_wider(names_from = DepMap_ID, values_from = pctviability) %>%
  mutate(Condition = paste0("N:PRSM:", Condition))

nk_prism_auc <- fread("../NK_PRISM/NK_PRISM_heme_pctviability_auc_norm.txt")
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
prism_celllines <- colnames(fm_prism)

fm <- fm[, prism_celllines]

# GSVA
core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")

genesets = list(SCRNASEQ_NK_RESPONSE = core)

fm_gsva <- fm[grepl("N:GEXP", rownames(fm)),]
rownames(fm_gsva) <- gsub("N:GEXP:", "", rownames(fm_gsva))

# run GSVA
gsva_results <- gsva(as.matrix(fm_gsva), genesets)

mat_gsea <- t(apply(gsva_results["SCRNASEQ_NK_RESPONSE",, drop = F], 1, scale))
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



# match subtype from sample info to feature matrix
Subtype <- ccle_annot$Subtype[match(colnames(fm), ccle_annot$DepMap_ID)]
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

# read cancer types and colors
cols <- fread("nk_crispr_colors.txt", data.table = F) %>% dplyr::rename(cancer = cancer_type)
cols$cancer_main <- gsub("ALCL", "TCL", gsub("DLBCL", "BCL", as.character(cols$cancer)))

cols_vector <- cols$color
names(cols_vector) <- cols$cancer_main


# NCL3LG1 and ULBP1 scatter plots (Figures 3F-G)

plot_scatterplot <- function(gene){
  
  df <- data.frame(gexp = as.numeric(fm[paste0("N:GEXP:", gene),]), 
                   meth = as.numeric(fm[paste0("N:METH:", gene),]),
                   cancer_type = factor(annot$cancertype, levels = c("AML", "TCL", "T-ALL", "BCL", "MM", "B-ALL")),
                   auc = as.numeric(fm_prism["N:PRSM:AUC",]))
  
  df <- df %>% filter(!is.na(meth)) # remove samples without methylation data
  
  ggplot(df, aes(x = gexp, y = auc, fill = meth)) +
    geom_point(pch = 21, size = 3) +
    geom_smooth(method='lm', se = F, color = "black") +
    scale_fill_distiller("Methylation", palette = "RdPu", direction = 1, na.value = "grey70") +
    theme_cowplot() +
    theme(plot.title = element_text(face = "italic", hjust = 0.5),
      legend.justification = "right",
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(0.3, "cm"),
      plot.margin = margin(t=5,1,1,1, "lines"),
      legend.direction = "horizontal",
      legend.position = c(0.9, 1.5)) +
    xlab("Expression (log2)") +
    ylab("PRISM AUC") +
    ggtitle(gene) +
    guides(color = "none",
           fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    ggpubr::stat_cor(label.sep = "\n", label.x = 2.75)
  
}

p1 <- plot_scatterplot("NCR3LG1")
p2 <- plot_scatterplot("ULBP1") + guides(fill = "none") + theme(plot.margin = margin(t=1,1,1,1, "lines"))

p1 / p2

ggsave("genes_methylation_scatterplots.pdf", height = 7, width = 3)


# NK response gene set (Figure 3H)
cols_vector <- cols_vector[c("AML", "TCL", "T-ALL", "BCL", "MM", "B-ALL")]

df <- data.frame(nkresponse = as.numeric(mat_gsea),
                 cancer_type = factor(annot$cancertype, levels = c("AML", "TCL", "T-ALL", "BCL", "MM", "B-ALL")),
                 auc = as.numeric(fm_prism["N:PRSM:AUC",]))

ggpubr::ggscatter(df, x = "nkresponse", y = "auc", shape = 21, size = 3, fill = "cancer_type",
                  add = "reg.line", cor.coef = T, cor.method = "spearman", palette = cols_vector, cor.coeff.args = list(label.sep = "\n")) +
  ylab("PRISM AUC") +
  xlab("Core NK cell response\n(GSVA Z-score)") +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 15)) +
  guides(fill = "none")

ggsave("core_nk_response_auc_scatter.pdf", height = 3, width = 3)

# faceted by cancer type (Figure S5I)
ggpubr::ggscatter(df, x = "nkresponse", y = "auc", shape = 21, size = 3, fill = "cancer_type",
                  add = "reg.line", cor.coef = T, cor.method = "spearman", palette = cols_vector[levels(df$cancer_type)],
                  cor.coeff.args = list(label.sep = "\n")) +
  ylab("PRISM AUC") +
  xlab("Core NK cell response\n(GSVA Z-score)") +
  theme(legend.title = element_blank()) +
  facet_grid(. ~ cancer_type) +
  theme(legend.title = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

ggsave("core_nk_response_auc_scatter_bycancertype.pdf", height = 3, width = 10)


# HLA I score

# add HLA I score
dat_a=fm_prism[rownames(fm_prism)%in%c("N:GEXP:B2M",
                                       "N:GEXP:HLA-A",
                                       "N:GEXP:HLA-B",
                                       "N:GEXP:HLA-C"),]

dat=2^dat_a+0.01
gm1=log2(t(apply(dat, 2, gm_mean)))
rownames(gm1)="N:SAMP:HLAIScore"

fm_prism<- rbind(fm_prism, gm1)

# scatter plot (Figure 3I)
df <- data.frame(hlaiscore = as.numeric(fm_prism["N:SAMP:HLAIScore",]),
                 cancer_type = factor(annot$cancertype, levels = c("AML", "TCL", "T-ALL", "BCL", "MM", "B-ALL")),
                 auc = as.numeric(fm_prism["N:PRSM:AUC",]))

ggpubr::ggscatter(df, x = "hlaiscore", y = "auc", shape = 21, size = 3, fill = "cancer_type",
                  add = "reg.line", cor.coef = T, cor.method = "spearman", cor.coeff.args = list(label.sep = "\n"), palette = cols_vector[levels(df$cancer_type)]) +
  ylab("PRISM AUC") +
  xlab("HLA I score") +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 15)) +
  guides(fill = "none")

ggsave("hlaiscore_auc_scatter.pdf", height = 3, width = 3)

# faceted by cancer type (Figure S5I)
ggpubr::ggscatter(df, x = "hlaiscore", y = "auc", shape = 21, size = 3, fill = "cancer_type",
                  add = "reg.line", cor.coef = T, cor.method = "spearman", palette = cols_vector[levels(df$cancer_type)],
                  cor.coeff.args = list(label.sep = "\n")) +
  ylab("PRISM AUC") +
  xlab("HLA I score") +
  theme(legend.title = element_blank()) +
  facet_grid(. ~ cancer_type) +
  theme(legend.title = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_legend(nrow = 1))


ggsave("hlaiscore_auc_scatter_bycancertype.pdf", height = 3, width = 10)

