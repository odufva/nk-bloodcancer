
# NK cell co-culture with panel of 26 cell lines, hashing scRNA-seq analysis
# Run GSEA

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(data.table)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(fgsea)

pathways_reactome <- gmtPathways("c2.cp.reactome.v7.0.symbols.gmt")
pathways_hallmark <- gmtPathways("h.all.v7.0.symbols.gmt")
pathways_biocarta <- gmtPathways("c2.cp.biocarta.v7.0.symbols.gmt")

dir.create("results/celllinepanel_combined/gsea/")

## D14 NK treated

# load DE results
data <- fread("results/celllinepanel_combined/targets/deg/deg_mainclusters_all.txt", data.table = F)

# GSEA function
run_gsea <- function(DATA, NAME, CELLLINE, PATHWAYS){
  
  ranks <- DATA %>% 
    filter(cell_line == CELLLINE) %>%
    dplyr::select(gene, avg_log2FC, p_val) %>%
    dplyr::select(gene, avg_log2FC) %>%
    tibble::deframe()
  
  # GSEA
  fgseaRes <- fgsea(PATHWAYS, ranks, maxSize=500, eps = 0)
  fgseaRes <- fgseaRes %>%
    mutate(pos_neg = ifelse(ES > 0, "pos", "neg")) %>%
    arrange(pos_neg, pval) %>%
    select(-pos_neg) %>% 
    mutate(cell_line = CELLLINE)
  
  # table of top pathways
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  pdf(paste0("results/celllinepanel_combined/gsea/", NAME, "_", CELLLINE, ".pdf"), height = 7, width = 20)
  plotGseaTable(PATHWAYS[topPathways], ranks, fgseaRes, gseaParam = 0.5)
  dev.off()
  
  return(fgseaRes)
}

# run GSEA with 3 different gene set collections
celllines <- unique(data$cell_line)

result_hallmark <- lapply(celllines, run_gsea, DATA = data, NAME = "HALLMARK", PATHWAYS = pathways_hallmark) %>% bind_rows()
fwrite(result_hallmark, "results/celllinepanel_combined/gsea/hallmark.txt", quote = F, row.names = F, sep = "\t")

result_reactome <- lapply(celllines, run_gsea, DATA = data, NAME = "REACTOME", PATHWAYS = pathways_reactome) %>% bind_rows()
fwrite(result_reactome, "results/celllinepanel_combined/gsea/reactome.txt", quote = F, row.names = F, sep = "\t")

result_biocarta <- lapply(celllines, run_gsea, DATA = data, NAME = "BIOCARTA", PATHWAYS = pathways_biocarta) %>% bind_rows()
fwrite(result_biocarta, "results/celllinepanel_combined/gsea/biocarta.txt", quote = F, row.names = F, sep = "\t")

## ---------------------------

# Dot plot 

# HALLMARK

# read GSEA results
gsea <- fread("results/celllinepanel_combined/gsea/hallmark.txt", data.table = F)

gsea <- gsea %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))
signif <- gsea %>% filter(padj < 0.05) %>% select(pathway) %>% unique()

signif_sentenceCase <- signif %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))

gsea_clust <- gsea %>%
  filter(pathway %in% signif$pathway) %>%
  select(cell_line, pathway_sentenceCase, NES) %>%
  tidyr::pivot_wider(names_from = cell_line, values_from = NES) %>%
  as.data.frame()

rownames(gsea_clust) <- gsea_clust$pathway_sentenceCase
gsea_clust$pathway_sentenceCase <- NULL
gsea_clust <- gsea_clust[signif_sentenceCase$pathway_sentenceCase,]

# clustering with hclust
dd.row <- as.dendrogram(hclust(dist(gsea_clust), method = "ward.D2"))
dd.col <- as.dendrogram(hclust(dist(t(gsea_clust)), method = "ward.D2"))

# ordering based on clustering
row.ord <- order.dendrogram(dd.row)
col.ord <- order.dendrogram(dd.col)


# dot plot
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = signif_sentenceCase$pathway_sentenceCase[row.ord])) %>%
  mutate(celL_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))]))


ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_x_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/hallmark_dotplot.pdf", height = 30, width = 15)

# recurrent pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(all(NES > 0) | all(NES < 0)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

nes_order_celllines <- gsea_plot %>% 
  filter(pathway_sentenceCase == "INTERFERON GAMMA RESPONSE") %>% 
  arrange(padj) %>% 
  select(cell_line) %>% 
  tibble::deframe()

gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order)) %>% 
  mutate(cell_line = factor(cell_line, levels = nes_order_celllines))

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  #scale_size(limits = c(0,6)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  #plot.margin = unit(c(1.75,1.75,0,0), "cm")) +
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/hallmark_dotplot_recurrent.pdf", height = 5, width = 12)

# top pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(any(NES > 2) | any(NES < -2)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

nes_order_celllines <- gsea_plot %>%
  filter(pathway_sentenceCase == "INTERFERON GAMMA RESPONSE") %>%
  arrange(padj) %>%
  select(cell_line) %>%
  tibble::deframe()

gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order)) %>% 
  mutate(cell_line = factor(cell_line, levels = nes_order_celllines))

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/hallmark_dotplot_top.pdf", height = 6, width = 12)


## ---------------------------------------

# REACTOME

# read GSEA results
gsea <- fread("results/celllinepanel_combined/gsea/reactome.txt", data.table = F)

gsea <- gsea %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))
signif <- gsea %>% filter(padj < 0.05) %>% select(pathway) %>% unique()

signif_sentenceCase <- signif %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))

gsea_clust <- gsea %>%
  filter(pathway %in% signif$pathway) %>%
  select(cell_line, pathway_sentenceCase, NES) %>%
  tidyr::pivot_wider(names_from = cell_line, values_from = NES) %>%
  as.data.frame()

rownames(gsea_clust) <- gsea_clust$pathway_sentenceCase
gsea_clust$pathway_sentenceCase <- NULL
gsea_clust <- gsea_clust[signif_sentenceCase$pathway_sentenceCase,]

# clustering with hclust
dd.row <- as.dendrogram(hclust(dist(gsea_clust), method = "ward.D2"))
dd.col <- as.dendrogram(hclust(dist(t(gsea_clust)), method = "ward.D2"))

# ordering based on clustering
row.ord <- order.dendrogram(dd.row)
col.ord <- order.dendrogram(dd.col)


# dot plot

gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = signif_sentenceCase$pathway_sentenceCase[row.ord])) %>%
  mutate(celL_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))]))


ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_x_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/reactome_dotplot.pdf", height = 8, width = 18)

# recurrent pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(all(NES > 0) | all(NES < 0)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order))

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/reactome_dotplot_recurrent.pdf", height = 20, width = 20)


# all pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(all(NES > 0) | all(NES < 0)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()


gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order))

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/reactome_dotplot_top.pdf", height = 30, width = 20)

# top pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(any(NES > 2.75) | any(NES < -2.75)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

nes_order_celllines <- gsea_plot %>%
  filter(pathway_sentenceCase == "CYTOKINE SIGNALING IN IMMUNE SYSTEM") %>%
  arrange(padj) %>%
  select(cell_line) %>%
  tibble::deframe()

nes_order_celllines <- c(as.character(nes_order_celllines), celllines[!celllines %in% as.character(nes_order_celllines)])

gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order)) %>% 
  mutate(cell_line = factor(cell_line, levels = nes_order_celllines))

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/reactome_dotplot_top.pdf", height = 8, width = 20)


## ---------------------------------------

# BIOCARTA

# read GSEA results
gsea <- fread("results/celllinepanel_combined/gsea/biocarta.txt", data.table = F)

gsea <- gsea %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))

signif <- gsea %>% filter(padj < 0.05) %>% select(pathway) %>% unique()

signif_sentenceCase <- signif %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway)))

gsea_clust <- gsea %>%
  filter(pathway %in% signif$pathway) %>%
  select(cell_line, pathway_sentenceCase, NES) %>%
  tidyr::pivot_wider(names_from = cell_line, values_from = NES) %>%
  as.data.frame()

rownames(gsea_clust) <- gsea_clust$pathway_sentenceCase
gsea_clust$pathway_sentenceCase <- NULL
gsea_clust <- gsea_clust[signif_sentenceCase$pathway_sentenceCase,]

# clustering with hclust
dd.row <- as.dendrogram(hclust(dist(gsea_clust), method = "ward.D2"))
dd.col <- as.dendrogram(hclust(dist(t(gsea_clust)), method = "ward.D2"))

# ordering based on clustering
row.ord <- order.dendrogram(dd.row)
col.ord <- order.dendrogram(dd.col)


# dot plot

gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = signif_sentenceCase$pathway_sentenceCase[row.ord])) %>%
  mutate(celL_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))]))


ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_x_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/biocarta_dotplot.pdf", height = 8, width = 18)

# recurrent pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(all(NES > 0) | all(NES < 0)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order)) 

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/biocarta_dotplot_recurrent.pdf", height = 20, width = 20)


# all pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(all(NES > 0) | all(NES < 0)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order))

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/biocarta_dotplot_top.pdf", height = 30, width = 20)

# top pathways
gsea_plot <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(cell_line = factor(cell_line, levels = unique(cell_line)[order(unique(cell_line))])) %>%
  group_by(pathway_sentenceCase) %>%
  filter(n()>2)

nes_order <- gsea_plot %>% 
  group_by(pathway_sentenceCase) %>%
  filter(any(NES > 0.5) | any(NES < -0.5)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

nes_order_celllines <- gsea_plot %>%
  filter(pathway_sentenceCase == "MHC PATHWAY") %>%
  arrange(padj) %>%
  select(cell_line) %>%
  tibble::deframe()

gsea_plot <- gsea_plot %>% 
  filter(pathway_sentenceCase %in% nes_order) %>% 
  mutate(pathway_sentenceCase = factor(pathway_sentenceCase, levels = nes_order)) %>% 
  mutate(cell_line = factor(cell_line, levels = nes_order_celllines))

ggplot(gsea_plot, aes(y = pathway_sentenceCase, x = cell_line, color = NES, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "BrBG", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(gsea_plot$NES)) * c(-1, 1)) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_cowplot() +
  theme(plot.title = element_text(face = "plain", hjust = 0.5)) +#,
  labs(color = "NES", size = "FDR (-log10)")

ggsave("results/celllinepanel_combined/gsea/biocarta_dotplot_top.pdf", height = 8, width = 20)

