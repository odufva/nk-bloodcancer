
# Compare DEGs of NK vs no NK in CROP-seq to genes correlating with NK cell infiltration in patient data
# Reddy DLBCL and TCGA AML data

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(viridis)
library(GSVA)

# laod plotting scripts
source("public_datasets/plotting_functions_may2020.R")
source("public_datasets/compute.pairwise.R")
source("public_datasets/regulon.feats.R")
source("public_datasets/functions_statistics.R")

dir.create("results/celllinepanel_combined/patient_integration//")
dir.create("results/celllinepanel_combined/patient_integration/bulkrnaseq/")


## CROP-seq DEG NK vs no NK (no logfc threshold)

# load hashing DEGs and take median avg_log2FC and p_val_adj over cell lines
hashing <- fread("results/celllinepanel_combined/targets/deg/deg_mainclusters_all.txt", data.table = F)
hashing_median <- hashing %>%
  group_by(gene) %>% 
  summarize(avg_log2FC = median(avg_log2FC),
            p_val_adj = median(p_val_adj))

topgenes <- hashing %>% 
  filter(avg_log2FC > 0.5) %>% 
  group_by(gene) %>% 
  summarize(cell_line_count = n()) %>% 
  top_n(50, cell_line_count) %>%
  arrange(desc(cell_line_count)) %>% 
  select(gene) %>%
  tibble::deframe()

type1ifn <- c("MX1", "MX2", "OAS1", "OAS2", "IRF7")


## Pairwise correlations of gene expression with cytolytic score 

# Reddy et al. DLBCL

# load feature matrix
fm=get(load("data/REDDY_DLBCL_fm.Rdata"))

# Compute pairwise correlations for cytolytic score/NK infiltration
genelist = c("CytolyticScore")

# nested correlations:
extrafeatures=c(grep("N:GEXP:", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH", prettyNumbers = F)

# write result table
fwrite(results, "results/celllinepanel_combined/patient_integration/bulkrnaseq/reddy_cytscore_correlations.tsv", sep = "\t")


# load results
hashing_median <- hashing_median %>% 
  select(gene, avg_log2FC_scrnaseq = avg_log2FC, adj_p_scrnaseq = p_val_adj)

cytscore_correlation <- fread("results/celllinepanel_combined/patient_integration/bulkrnaseq/reddy_cytscore_correlations.tsv", data.table = F)
cytscore_correlation <- cytscore_correlation %>%
  select(gene = featureB, cor_patient = cor, p_patient = p, adj_p_patient = adj.p) %>%
  mutate(gene = gsub("N:GEXP:", "", gene)) %>% 
  mutate(gene = gsub("HLA_", "HLA-", gene)) %>% 
  mutate(gene = gsub("\\_", "\\.", gene))


# merge result tables
combined_de <- fread("results/celllinepanel_combined/targets/deg/targets_deg_combined.txt", data.table = F)
combined_de <- combined_de %>% 
  select(gene, avg_log2FC_scrnaseq = avg_log2FC, adj_p_scrnaseq = p_val_adj)
data <- merge(combined_de, cytscore_correlation, by = "gene") %>% arrange(desc(cor_patient))

fwrite(data, "results/celllinepanel_combined/patient_integration/bulkrnaseq/nk_hashing_reddy_cytscore.tsv", sep = "\t")

data <- fread("results/celllinepanel_combined/patient_integration/bulkrnaseq/nk_hashing_reddy_cytscore.tsv", data.table = F)

#Thresholds
threshold_cor_patient = 0.1
threshold_adj_p_patient = 0.05
threshold_avg_log2FC_scrnaseq = 0.3
threshold_adj_p_scrnaseq = 0.05

data_plot <- data %>% 
  mutate(label = ifelse(cor_patient > threshold_cor_patient &
                          adj_p_patient < threshold_adj_p_patient &
                          avg_log2FC_scrnaseq > threshold_avg_log2FC_scrnaseq &
                          adj_p_scrnaseq < threshold_adj_p_scrnaseq, "show_label", "hide_label")) %>% 
  mutate(p_patient = ifelse(p_patient== 0, 1e-150, p_patient)) %>% 
  mutate(label = ifelse(cor_patient > 0.8 | avg_log2FC_scrnaseq > 0.3, "show_label", label))

ggplot(data_plot, aes(y = cor_patient, x = avg_log2FC_scrnaseq)) +
  ggrastr::geom_point_rast(data = data_plot[data_plot$label!="show_label",], alpha = 0.1, color = "grey50") +
  geom_point(data = data_plot[data_plot$label=="show_label" & data_plot$avg_log2FC_scrnaseq > 0.3,], color = brewer.pal(9, "RdYlGn")[2], shape = 16, size = 2, alpha = 0.8) +
  geom_point(data = data_plot[data_plot$label=="show_label" & data_plot$avg_log2FC_scrnaseq <= 0.3,], color = brewer.pal(9, "RdYlGn")[8], shape = 16, size = 2, alpha = 0.8) +
  geom_text_repel(data = data_plot[data_plot$label=="show_label",], aes(label = ifelse(label == "show_label", gene, "")),
                  size = 2.5,
                  color = "black",
                  max.overlaps = 100,
                  point.padding = 0.5,
                  fontface = "italic") +
  theme_cowplot() +
  xlab("Fold change (log2)\nWith vs without NK co-culture") +
  ylab("Correlation with cytolytic score") +
  labs(size = "P value (-log10)\nscRNA-seq",
       color = "P value (-log10)\nPatients") +
  guides(alpha = "none", fill = "none") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50")

ggsave("results/celllinepanel_combined/patient_integration/bulkrnaseq/reddy_dlbcl_scatter.pdf", height = 4, width = 5)


# flipped axes (Figure 2G)
core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "STAT1", "IRF1", "IRF9", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "UBE2L6", "GNLY", "CCL5")

threshold_adj_p_patient = 0.05
threshold_avg_log2FC_scrnaseq = 0.2
threshold_adj_p_scrnaseq = 0.05

data_plot <- data %>% 
  mutate(label = ifelse(adj_p_patient < threshold_adj_p_patient &
                          avg_log2FC_scrnaseq > threshold_avg_log2FC_scrnaseq &
                          adj_p_scrnaseq < threshold_adj_p_scrnaseq, "show_label", "hide_label")) %>%
  mutate(p_patient = ifelse(p_patient== 0, 1e-150, p_patient))


ggplot(data_plot, aes(y = avg_log2FC_scrnaseq, x = cor_patient)) +
  ggrastr::geom_point_rast(data = data_plot[data_plot$label!="show_label",], alpha = 0.1, color = "grey50") +
  geom_point(data = data_plot[data_plot$label=="show_label",], color = brewer.pal(9, "RdYlGn")[2], shape = 16, size = 2, alpha = 0.8) +
  geom_point(data = data_plot[data_plot$label=="show_label" & data_plot$gene %in% core,], color = brewer.pal(9, "RdYlGn")[8], shape = 16, size = 2, alpha = 0.8) +
  geom_text_repel(data = data_plot[data_plot$label=="show_label",], aes(label = ifelse(label == "show_label", gene, "")),
                  size = 2.5,
                  color = "black",
                  max.overlaps = 100,
                  point.padding = 0.5,
                  fontface = "italic") +
  theme_cowplot() +
  ylab("Fold change (log2)\nWith vs without NK co-culture") +
  xlab("Correlation with cytolytic score") +
  labs(size = "P value (-log10)\nscRNA-seq",
       color = "P value (-log10)\nPatients") +
  guides(alpha = "none", fill = "none") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50")

ggsave("results/celllinepanel_combined/patient_integration/bulkrnaseq/reddy_dlbcl_scatter_flipped.pdf", height = 3.5, width = 4.5)


# Correlation of cytolytic score with NK-induced genes
core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")

fm_gsva <- fm[grepl("N:GEXP", rownames(fm)),]
rownames(fm_gsva) <- gsub("N:GEXP:", "", rownames(fm_gsva))
fm_gsva <- fm_gsva[,!is.na(fm_gsva["KRAS",])]

# run GSVA
gsva_results <- gsva(as.matrix(fm_gsva), list(NKinduced = core))

df <- data.frame(cytscore = as.numeric(fm["N:SAMP:CytolyticScore", !is.na(fm["N:GEXP:KRAS",])]),
                 nkinduced = as.numeric(gsva_results))

df <- data.frame(cytscore = as.numeric(fm["N:GEXP:HLA", !is.na(fm["N:GEXP:KRAS",])]),
                 nkinduced = as.numeric(gsva_results))

ggplot(df, aes(x = cytscore, y = nkinduced)) +
  geom_point() +
  geom_smooth(method = "lm", color = "darkred") +
  theme_cowplot() +
  ylab("NK-induced gene signature (in vitro)") +
  xlab("Cytolytic score (lymphoma patients)") +
  ggpubr::stat_cor(label.sep = "\n", label.y = 1)

ggsave("results/celllinepanel_combined/patient_integration/bulkrnaseq/reddy_dlbcl_nkinduced_gsva_cytscore.pdf", height = 3.75, width = 4)


# TCGA AML

# load feature matrix
fm=get(load("data/DUFVA_TCGA_AML_FM_meth.Rdata"))

# Compute pairwise correlations for cytolytic score/NK infiltration
genelist = c("DUFVA_CYTOLYTIC_SCORE")

# nested correlations:
extrafeatures=c(grep("N:GEXP:", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH", prettyNumbers = F)

# write result table
fwrite(results, "results/celllinepanel_combined/patient_integration/bulkrnaseq/tcga_aml_cytscore_correlations.tsv", sep = "\t")


# load results
cytscore_correlation <- fread("results/celllinepanel_combined/patient_integration/bulkrnaseq/tcga_aml_cytscore_correlations.tsv", data.table = F)
cytscore_correlation <- cytscore_correlation %>%
  select(gene = featureB, cor_patient = cor, p_patient = p, adj_p_patient = adj.p) %>%
  mutate(gene = gsub("N:GEXP:|:chr.*", "", gene))


# merge result tables
data <- merge(combined_de, cytscore_correlation, by = "gene") %>% arrange(desc(cor_patient))

fwrite(data, "results/celllinepanel_combined/patient_integration/bulkrnaseq/nk_hashing_tcga_aml_cytscore.tsv", sep = "\t")

data <- fread("results/celllinepanel_combined/patient_integration/bulkrnaseq/nk_hashing_tcga_aml_cytscore.tsv", data.table = F)

# Thresholds
threshold_cor_patient = 0.1
threshold_adj_p_patient = 0.05
threshold_avg_log2FC_scrnaseq = 0.3
threshold_adj_p_scrnaseq = 0.05

data_plot <- data %>% 
  mutate(label = ifelse(cor_patient > threshold_cor_patient &
                          adj_p_patient < threshold_adj_p_patient &
                          avg_log2FC_scrnaseq > threshold_avg_log2FC_scrnaseq &
                          adj_p_scrnaseq < threshold_adj_p_scrnaseq, "show_label", "hide_label")) %>% 
  mutate(p_patient = ifelse(p_patient== 0, 1e-150, p_patient)) %>% 
  mutate(label = ifelse(cor_patient > 0.85 | avg_log2FC_scrnaseq > 0.3, "show_label", label))

ggplot(data_plot, aes(y = cor_patient, x = avg_log2FC_scrnaseq)) +
  ggrastr::geom_point_rast(data = data_plot[data_plot$label!="show_label",], alpha = 0.1, color = "grey50") +
  geom_point(data = data_plot[data_plot$label=="show_label" & data_plot$avg_log2FC_scrnaseq > 0.3,], color = brewer.pal(9, "RdYlGn")[2], shape = 16, size = 2, alpha = 0.8) +
  geom_point(data = data_plot[data_plot$label=="show_label" & data_plot$avg_log2FC_scrnaseq <= 0.3,], color = brewer.pal(9, "RdYlGn")[8], shape = 16, size = 2, alpha = 0.8) +
  geom_text_repel(data = data_plot[data_plot$label=="show_label",], aes(label = ifelse(label == "show_label", gene, "")),
                  size = 2.5,
                  color = "black",
                  max.overlaps = 100,
                  point.padding = 0.5,
                  fontface = "italic") +
  theme_cowplot() +
  xlab("Fold change (log2)\nWith vs without NK cell co-culture") +
  ylab("Correlation with cytolytic score") +
  labs(size = "P value (-log10)\nscRNA-seq",
       color = "P value (-log10)\nPatients") +
  guides(alpha = "none", fill = "none") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50")

ggsave("results/celllinepanel_combined/patient_integration/bulkrnaseq/tcga_aml_scatter.pdf", height = 4, width = 5)


# flipped axes (Figure S2C)
core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")

threshold_adj_p_patient = 0.05
threshold_avg_log2FC_scrnaseq = 0.2
threshold_adj_p_scrnaseq = 0.05

data_plot <- data %>% 
  mutate(label = ifelse(adj_p_patient < threshold_adj_p_patient &
                          avg_log2FC_scrnaseq > threshold_avg_log2FC_scrnaseq &
                          adj_p_scrnaseq < threshold_adj_p_scrnaseq, "show_label", "hide_label")) %>%
  mutate(p_patient = ifelse(p_patient== 0, 1e-150, p_patient))

ggplot(data_plot, aes(y = avg_log2FC_scrnaseq, x = cor_patient)) +
  ggrastr::geom_point_rast(data = data_plot[data_plot$label!="show_label",], alpha = 0.1, color = "grey50") +
  geom_point(data = data_plot[data_plot$label=="show_label" & data_plot$avg_log2FC_scrnaseq > 0.2,], color = brewer.pal(9, "RdYlGn")[2], shape = 16, size = 2, alpha = 0.8) +
  geom_point(data = data_plot[data_plot$label=="show_label" & data_plot$gene %in% core,], color = brewer.pal(9, "RdYlGn")[8], shape = 16, size = 2, alpha = 0.8) +
  geom_text_repel(data = data_plot[data_plot$label=="show_label",], aes(label = ifelse(label == "show_label", gene, "")),
                  size = 2.5,
                  color = "black",
                  max.overlaps = 100,
                  point.padding = 0.5,
                  fontface = "italic") +
  theme_cowplot() +
  ylab("Fold change (log2)\nWith vs without NK co-culture") +
  xlab("Correlation with cytolytic score") +
  labs(size = "P value (-log10)\nscRNA-seq",
       color = "P value (-log10)\nPatients") +
  guides(alpha = "none", fill = "none") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50")

ggsave("results/celllinepanel_combined/patient_integration/bulkrnaseq/tcga_aml_scatter_flipped.pdf", height = 3.1, width = 4)



