
# NK cell co-culture with panel of 26 cell lines, hashing scRNA-seq analysis
# NK cell analysis
# Correct batch by regressing out in ScaleData

# load libraries
library(dplyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(data.table)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(SingleR)

source("scripts/fun_getGenes.R")

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))

dir.create("results/celllinepanel/hto_singlets/nk/batchcorrected")

hash_seurat_singlet <- readRDS("results/celllinepanel/celllinepanel_nk1_seurat.rds")

# subset to expanded NK cells
hash_seurat_nk <- subset(hash_seurat_singlet, seurat_clusters %in% c("1", "28"))
hash_seurat_nk <- subset(hash_seurat_nk, orig.ident %in% c("Batch1_NK", "Batch2_NK"))

# Select the top 1000 most variable features
hash_seurat_nk <- FindVariableFeatures(hash_seurat_nk, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_nk <- CellCycleScoring(hash_seurat_nk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_nk <- ScaleData(hash_seurat_nk, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_nk))

# Run PCA
hash_seurat_nk <- RunPCA(hash_seurat_nk, features = VariableFeatures(hash_seurat_nk))

# Plot PCA with batch

DimPlot(hash_seurat_nk, reduction = "pca", dims = c(1,2), group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(5))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/pca_batch_beforecorrection_dims12.png", width = 5, height = 4)

DimPlot(hash_seurat_nk, reduction = "pca", dims = c(2,3), group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(5))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/pca_batch_beforecorrection_dims23.png", width = 5, height = 4)

# Regress out batch
hash_seurat_nk <- ScaleData(hash_seurat_nk, vars.to.regress = c("S.Score", "G2M.Score", "orig.ident"), features = VariableFeatures(hash_seurat_nk))

# Run PCA
hash_seurat_nk <- RunPCA(hash_seurat_nk, features = VariableFeatures(hash_seurat_nk))

# Select the top 20 PCs for clustering and UMAP
hash_seurat_nk <- FindNeighbors(hash_seurat_nk, reduction = "pca", dims = c(1:20))
hash_seurat_nk <- FindClusters(hash_seurat_nk, resolution = 0.3, verbose = FALSE)
hash_seurat_nk <- RunUMAP(hash_seurat_nk, reduction = "pca", dims = c(1:20))

hash_seurat_nk$cluster <- factor(hash_seurat_nk$seurat_clusters,
                                 levels = c(0:4),
                                 labels = c("Resting (0)",
                                            "Adaptive (1)",
                                            "Activated (2)",
                                            "Type I IFN (3)",
                                            "Cytokine (4)"))


## Visualize
DimPlot(hash_seurat_nk, group.by = "cluster", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(6))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/umap.png", width = 7, height = 4)

DimPlot(hash_seurat_nk, group.by = "orig.ident", label = F, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(5))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/umap_batch.png", width = 7, height = 4)

DimPlot(hash_seurat_nk, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/umap_cycle.png", width = 6, height = 4)

# UMAPs split by cell lines
DimPlot(hash_seurat_nk, group.by = "cluster", split.by = "hash.ID", label = F, repel = T, ncol = 7) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(5))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/umap_celllines.png", width = 14, height = 8)


# DEGs of clusters
Idents(hash_seurat_nk) <- "cluster"
all_deg <- FindAllMarkers(hash_seurat_nk, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/celllinepanel/hto_singlets/nk/batchcorrected/clusters_deg.txt", sep = "\t", quote = F, row.names = F)
xlsx::write.xlsx(all_deg, "../NK_resistance Heme CollabPaper/Manuscript/Supplement/TableS1.xlsx", sheetName = "B. Cluster DEG (expanded NK)", row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_nk,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/cluster_deg_dotplot.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_nk,
               features = top25$gene,
               raster = F,
               angle = 90) + NoLegend()
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/cluster_deg_heatmap.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_nk, "results/celllinepanel/hto_singlets/nk/batchcorrected/celllinepanel_nk1_nk_batchcorrected_seurat.rds")

## ----------------------------

# HTO singlet cell numbers

Idents(hash_seurat_nk) <- "hash.ID"

cell_numbers <- data.frame(table(Idents(hash_seurat_nk))) %>%
  mutate(Var1 = factor(Var1, levels = unique(hash_seurat_nk$hash.ID)))

ggplot(cell_numbers, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_col() +
  scale_fill_manual(values = getPalette3(27)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0.2, 0.2, 0.2, 2), "cm")) +
  xlab("") +
  ylab("Number of singlets") +
  guides(fill = F)

ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/cell_numbers.pdf", height = 5, width = 7)


# cluster fractions

cluster_freq <- data.frame(table(Idents(hash_seurat_nk), hash_seurat_nk$cluster))
cluster_freq_pert <- cluster_freq %>%
  group_by(Var1, Var2) %>%
  summarize(count = sum(Freq)) %>%
  mutate(freq = count / sum(count)) %>% 
  mutate(Var1 = gsub("-NK-expanded", "", Var1))

sample_order <- cluster_freq_pert %>%
  filter(Var2 %in% c("Resting (0)", "Adaptive (1)")) %>%
  group_by(Var1) %>% 
  summarize(freq_sum = sum(freq)) %>%
  arrange(freq_sum) %>%
  dplyr::select(Var1) %>%
  tibble::deframe()

cluster_freq_pert <- cluster_freq_pert %>%
  mutate(Var1 = factor(Var1, levels = sample_order))

ggplot(cluster_freq_pert, aes(x = Var1, y = freq, fill = Var2)) +
  geom_col() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = rev(getPalette(5))) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") +
  ylab("Cluster fraction") +
  labs(fill = "Cluster")

ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/cluster_fraction.pdf", height = 4, width = 9)

## -----------------------------------------

## PBMC NK

dir.create("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc")

# subset to expanded NK cells
hash_seurat_nk <- subset(hash_seurat_singlet, seurat_clusters %in% c("5"))
hash_seurat_nk <- subset(hash_seurat_nk, orig.ident %in% c("Batch1_NK2587_D0", "Batch2_NK2587_D0"))

# Select the top 1000 most variable features
hash_seurat_nk <- FindVariableFeatures(hash_seurat_nk, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_nk <- CellCycleScoring(hash_seurat_nk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_nk <- ScaleData(hash_seurat_nk, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_nk))

# Plot PCA with batch

DimPlot(hash_seurat_nk, reduction = "pca", dims = c(1,2), group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(5))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/pca_batch_beforecorrection_dims12.png", width = 5, height = 4)

DimPlot(hash_seurat_nk, reduction = "pca", dims = c(2,3), group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(5))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/pca_batch_beforecorrection_dims23.png", width = 5, height = 4)

# Regress out batch
hash_seurat_nk <- ScaleData(hash_seurat_nk, vars.to.regress = c("S.Score", "G2M.Score", "orig.ident"), features = VariableFeatures(hash_seurat_nk))

# Run PCA
hash_seurat_nk <- RunPCA(hash_seurat_nk, features = VariableFeatures(hash_seurat_nk))

# Select the top 20 PCs for clustering and UMAP
hash_seurat_nk <- FindNeighbors(hash_seurat_nk, reduction = "pca", dims = c(1:20))
hash_seurat_nk <- FindClusters(hash_seurat_nk, resolution = 0.5, verbose = FALSE)
hash_seurat_nk <- RunUMAP(hash_seurat_nk, reduction = "pca", dims = c(1:20))

hash_seurat_nk$cluster <- factor(hash_seurat_nk$seurat_clusters,
                                 levels = c(0:5),
                                 labels = c("CD56dim (0)",
                                            "Adaptive (1)",
                                            "Type I IFN (2)",
                                            "Cytokine (3)",
                                            "Activated (4)",
                                            "CD56bright (5)"))

## Visualize
DimPlot(hash_seurat_nk, group.by = "cluster", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(6))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/umap.png", width = 7, height = 4)

DimPlot(hash_seurat_nk, group.by = "orig.ident", label = F, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(6))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/umap_batch.png", width = 7, height = 4)

DimPlot(hash_seurat_nk, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/umap_cycle.png", width = 6, height = 4)

DimPlot(hash_seurat_nk, group.by = "cluster", split.by = "hash.ID", label = F, repel = T, ncol = 7) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(6))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/umap_celllines.png", width = 14, height = 8)

# DEGs of clusters
Idents(hash_seurat_nk) <- "cluster"
all_deg <- FindAllMarkers(hash_seurat_nk, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/clusters_deg.txt", sep = "\t", quote = F, row.names = F)
xlsx::write.xlsx(all_deg, "../NK_resistance Heme CollabPaper/Manuscript/Supplement/TableS1.xlsx", sheetName = "C. Cluster DEG (PBMC NK)", row.names = F, append = T)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_nk,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/cluster_deg_dotplot.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_nk,
               features = top25$gene,
               raster = F,
               angle = 90) + NoLegend()
ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/cluster_deg_heatmap.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_nk, "results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/celllinepanel_nk1_batchcorrected_seurat.rds")

## ----------------------------

# HTO singlet cell numbers

Idents(hash_seurat_nk) <- "hash.ID"

cell_numbers <- data.frame(table(Idents(hash_seurat_nk))) %>%
  mutate(Var1 = factor(Var1, levels = unique(hash_seurat_nk$hash.ID)))

ggplot(cell_numbers, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_col() +
  scale_fill_manual(values = getPalette3(27)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,2), "cm")) +
  xlab("") +
  ylab("Number of singlets") +
  guides(fill = F)

ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/cell_numbers.pdf", height = 5, width = 12)


# cluster fractions

cluster_freq <- data.frame(table(Idents(hash_seurat_nk), hash_seurat_nk$cluster))
cluster_freq_pert <- cluster_freq %>%
  #filter(Var2 %in% c(0:3)) %>%
  group_by(Var1, Var2) %>%
  summarize(count = sum(Freq)) %>%
  mutate(freq = count / sum(count)) %>% 
  mutate(Var1 = gsub("-NK-PBMC", "", Var1))

sample_order <- cluster_freq_pert %>%
  filter(Var2 %in% c("CD56dim (0)", "Adaptive (1)", "CD56bright (5)")) %>%
  group_by(Var1) %>% 
  summarize(freq_sum = sum(freq)) %>%
  arrange(freq_sum) %>%
  dplyr::select(Var1) %>%
  tibble::deframe()

cluster_freq_pert <- cluster_freq_pert %>%
  mutate(Var1 = factor(Var1, levels = sample_order))

ggplot(cluster_freq_pert, aes(x = Var1, y = freq, fill = Var2)) +
  geom_col() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = rev(getPalette(6))) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") +
  ylab("Cluster fraction") +
  labs(fill = "Cluster")

ggsave("results/celllinepanel/hto_singlets/nk/batchcorrected/nk_pbmc/cluster_fraction.pdf", height = 4, width = 9)

