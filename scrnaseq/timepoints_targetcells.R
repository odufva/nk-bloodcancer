
# NK cell co-culture time course experiment, hashing scRNA-seq analysis
# Target cell analysis


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
library(scales)
library(patchwork)

dir.create("results/timepoints/hto_singlets/targets")

hash_seurat_singlet <- readRDS("results/timepoints/timepoints_seurat.rds")


## K562

hash_seurat_k562 <- subset(hash_seurat_singlet, seurat_clusters %in% c("2", "3", "16", "18"))

# Select the top 1000 most variable features
hash_seurat_k562 <- FindVariableFeatures(hash_seurat_k562, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_k562 <- CellCycleScoring(hash_seurat_k562, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_k562 <- ScaleData(hash_seurat_k562, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_k562))

# Run PCA
hash_seurat_k562 <- RunPCA(hash_seurat_k562, features = VariableFeatures(hash_seurat_k562))

ElbowPlot(hash_seurat_k562, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_k562 <- FindNeighbors(hash_seurat_k562, reduction = "pca", dims = 1:20)
hash_seurat_k562 <- FindClusters(hash_seurat_k562, resolution = 0.4, verbose = FALSE)
hash_seurat_k562 <- RunUMAP(hash_seurat_k562, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_k562, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_k562.png", width = 6, height = 4)

DimPlot(hash_seurat_k562, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_cycle_k562.png", width = 6, height = 4)

DimPlot(hash_seurat_k562, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_k562.png", width = 14, height = 12)

DimPlot(hash_seurat_k562, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_batch_k562.png", width = 6, height = 4)


# selected conditions

hash_seurat_k562_selected <- subset(hash_seurat_k562, hash.ID %in% c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

hash_seurat_k562_selected$hash.ID <- factor(hash_seurat_k562_selected$hash.ID, levels =  c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

DimPlot(hash_seurat_k562_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) &
  xlab("") &
  ylab("")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_order_k562.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_k562) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_k562, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/clusters_deg_k562.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_k562,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_dotplot_k562.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_k562,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_heatmap_k562.pdf", p, height = 15, width = 15)


# remove outlier clusters
hash_seurat_k562 <- subset(hash_seurat_k562, seurat_clusters %in% c("0", "1", "2", "4", "6"))

# Select the top 1000 most variable features
hash_seurat_k562 <- FindVariableFeatures(hash_seurat_k562, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_k562 <- CellCycleScoring(hash_seurat_k562, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_k562 <- ScaleData(hash_seurat_k562, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_k562))

# Run PCA
hash_seurat_k562 <- RunPCA(hash_seurat_k562, features = VariableFeatures(hash_seurat_k562))

ElbowPlot(hash_seurat_k562, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_k562 <- FindNeighbors(hash_seurat_k562, reduction = "pca", dims = 1:20)
hash_seurat_k562 <- FindClusters(hash_seurat_k562, resolution = 0.4, verbose = FALSE)
hash_seurat_k562 <- RunUMAP(hash_seurat_k562, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_k562, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_k562_2.png", width = 6, height = 4)

DimPlot(hash_seurat_k562, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_cycle_k562_2.png", width = 6, height = 4)

DimPlot(hash_seurat_k562, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_k562_2.png", width = 14, height = 12)

DimPlot(hash_seurat_k562, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_batch_k562_2.png", width = 6, height = 4)


# selected conditions

hash_seurat_k562_selected <- subset(hash_seurat_k562, hash.ID %in% c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

hash_seurat_k562_selected$hash.ID <- factor(hash_seurat_k562_selected$hash.ID, levels =  c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

DimPlot(hash_seurat_k562_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) &
  xlab("") &
  ylab("")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_order_k562_2.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_k562) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_k562, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/clusters_deg_k562_2.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_k562,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_dotplot_k562_2.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_k562,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_heatmap_k562_2.pdf", p, height = 15, width = 15)


p1 <- VlnPlot(hash_seurat_k562, features = c("IFNB1"), slot = "data")
p2 <- VlnPlot(hash_seurat_k562, features = c("IFNB1"), slot = "counts")

p1 | p2
ggsave("results/timepoints/hto_singlets/targets/ifnb1_vlnplot_k562.pdf", height = 4, width = 10)


# timepoint DEG
Idents(hash_seurat_k562_selected) <- "hash.ID"
all_deg <- FindAllMarkers(hash_seurat_k562_selected, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/timepoints_deg_k562.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top timepoint DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_k562_selected,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_dotplot_k562.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_k562_selected,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_heatmap_k562.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_k562, "results/timepoints/hto_singlets/targets/timepoints_k562_seurat.rds")

## ----------------------------

## GDM1

hash_seurat_gdm1 <- subset(hash_seurat_singlet, seurat_clusters %in% c("1", "5", "10", "11", "13"))

# Select the top 1000 most variable features
hash_seurat_gdm1 <- FindVariableFeatures(hash_seurat_gdm1, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_gdm1 <- CellCycleScoring(hash_seurat_gdm1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_gdm1 <- ScaleData(hash_seurat_gdm1, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_gdm1))

# Run PCA
hash_seurat_gdm1 <- RunPCA(hash_seurat_gdm1, features = VariableFeatures(hash_seurat_gdm1))

ElbowPlot(hash_seurat_gdm1, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_gdm1 <- FindNeighbors(hash_seurat_gdm1, reduction = "pca", dims = 1:20)
hash_seurat_gdm1 <- FindClusters(hash_seurat_gdm1, resolution = 0.4, verbose = FALSE)
hash_seurat_gdm1 <- RunUMAP(hash_seurat_gdm1, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_gdm1, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_gdm1.png", width = 6, height = 4)


# Remove tiny outlier cluster and re-normalize 

hash_seurat_gdm1 <- subset(hash_seurat_gdm1, seurat_clusters %in% c("0", "1", "2", "3"))

# Select the top 1000 most variable features
hash_seurat_gdm1 <- FindVariableFeatures(hash_seurat_gdm1, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_gdm1 <- CellCycleScoring(hash_seurat_gdm1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_gdm1 <- ScaleData(hash_seurat_gdm1, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_gdm1))

# Run PCA
hash_seurat_gdm1 <- RunPCA(hash_seurat_gdm1, features = VariableFeatures(hash_seurat_gdm1))

ElbowPlot(hash_seurat_gdm1, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_gdm1 <- FindNeighbors(hash_seurat_gdm1, reduction = "pca", dims = 1:20)
hash_seurat_gdm1 <- FindClusters(hash_seurat_gdm1, resolution = 0.4, verbose = FALSE)
hash_seurat_gdm1 <- RunUMAP(hash_seurat_gdm1, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_gdm1, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_gdm1_2.png", width = 6, height = 4)

DimPlot(hash_seurat_gdm1, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_cycle_gdm1.png", width = 6, height = 4)

DimPlot(hash_seurat_gdm1, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_gdm1.png", width = 14, height = 12)

DimPlot(hash_seurat_gdm1, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_batch_gdm1.png", width = 6, height = 4)


# selected conditions

hash_seurat_gdm1_selected <- subset(hash_seurat_gdm1, hash.ID %in%  c("GDM1-0h", "GDM1-1h", "GDM1-3h", "GDM1-6h", "GDM1-12h", "GDM1-24h"))

hash_seurat_gdm1_selected$hash.ID <- factor(hash_seurat_gdm1_selected$hash.ID, levels = c("GDM1-0h", "GDM1-1h", "GDM1-3h", "GDM1-6h", "GDM1-12h", "GDM1-24h"))

DimPlot(hash_seurat_gdm1_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) &
  xlab("") &
  ylab("")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_order_gdm1.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_gdm1) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_gdm1, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/clusters_deg_gdm1.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_gdm1,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_dotplot_gdm1.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_gdm1,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_heatmap_gdm1.pdf", p, height = 15, width = 15)


# timepoint DEG
Idents(hash_seurat_gdm1_selected) <- "hash.ID"
all_deg <- FindAllMarkers(hash_seurat_gdm1_selected, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/timepoints_deg_gdm1.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top timepoint DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_gdm1_selected,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_dotplot_gdm1.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_gdm1_selected,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_heatmap_gdm1.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_gdm1, "results/timepoints/hto_singlets/targets/timepoints_gdm1_seurat.rds")

p1 <- VlnPlot(hash_seurat_gdm1, features = c("IFNB1"), slot = "data")
p2 <- VlnPlot(hash_seurat_gdm1, features = c("IFNB1"), slot = "counts")

p1 | p2
ggsave("results/timepoints/hto_singlets/targets/ifnb1_vlnplot_gdm1.pdf", height = 4, width = 10)

p1 <- VlnPlot(hash_seurat_gdm1_selected, features = c("IFNB1"), slot = "data") + NoLegend() + xlab("") + theme(plot.title = element_text(face = "italic"))
p2 <- VlnPlot(hash_seurat_gdm1_selected, features = c("IFNB1"), slot = "counts") + NoLegend() + xlab("") + theme(plot.title = element_text(face = "italic"))

p1 | p2
ggsave("results/timepoints/hto_singlets/targets/ifnb1_vlnplot_timepoints_gdm1.pdf", height = 4, width = 10)

## ----------------------------



## RI1

hash_seurat_ri1 <- subset(hash_seurat_singlet, seurat_clusters %in% c("4", "8", "9", "17"))

# Select the top 1000 most variable features
hash_seurat_ri1 <- FindVariableFeatures(hash_seurat_ri1, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_ri1 <- CellCycleScoring(hash_seurat_ri1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_ri1 <- ScaleData(hash_seurat_ri1, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_ri1))

# Run PCA
hash_seurat_ri1 <- RunPCA(hash_seurat_ri1, features = VariableFeatures(hash_seurat_ri1))

ElbowPlot(hash_seurat_ri1, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_ri1 <- FindNeighbors(hash_seurat_ri1, reduction = "pca", dims = 1:20)
hash_seurat_ri1 <- FindClusters(hash_seurat_ri1, resolution = 0.4, verbose = FALSE)
hash_seurat_ri1 <- RunUMAP(hash_seurat_ri1, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_ri1, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_ri1.png", width = 6, height = 4)

DimPlot(hash_seurat_ri1, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_cycle_ri1.png", width = 6, height = 4)

DimPlot(hash_seurat_ri1, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_ri1.png", width = 14, height = 12)

DimPlot(hash_seurat_ri1, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_batch_ri1.png", width = 6, height = 4)


# selected conditions

hash_seurat_ri1_selected <- subset(hash_seurat_ri1, hash.ID %in% c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

hash_seurat_ri1_selected$hash.ID <- factor(hash_seurat_ri1_selected$hash.ID, levels = c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

DimPlot(hash_seurat_ri1_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) &
  xlab("") &
  ylab("")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_order_ri1.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_ri1) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_ri1, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/clusters_deg_ri1.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_ri1,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_dotplot_ri1.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_ri1,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_heatmap_ri1.pdf", p, height = 15, width = 15)


# remove outlie clusters

hash_seurat_ri1 <- subset(hash_seurat_ri1, seurat_clusters %in% c("0", "1", "2", "3", "5"))

# Select the top 1000 most variable features
hash_seurat_ri1 <- FindVariableFeatures(hash_seurat_ri1, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_ri1 <- CellCycleScoring(hash_seurat_ri1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_ri1 <- ScaleData(hash_seurat_ri1, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_ri1))

# Run PCA
hash_seurat_ri1 <- RunPCA(hash_seurat_ri1, features = VariableFeatures(hash_seurat_ri1))

ElbowPlot(hash_seurat_ri1, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_ri1 <- FindNeighbors(hash_seurat_ri1, reduction = "pca", dims = 1:20)
hash_seurat_ri1 <- FindClusters(hash_seurat_ri1, resolution = 0.4, verbose = FALSE)
hash_seurat_ri1 <- RunUMAP(hash_seurat_ri1, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_ri1, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_ri1.png", width = 6, height = 4)

DimPlot(hash_seurat_ri1, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_cycle_ri1.png", width = 6, height = 4)

DimPlot(hash_seurat_ri1, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_ri1.png", width = 14, height = 12)

DimPlot(hash_seurat_ri1, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_batch_ri1.png", width = 6, height = 4)


# selected conditions

hash_seurat_ri1_selected <- subset(hash_seurat_ri1, hash.ID %in% c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

hash_seurat_ri1_selected$hash.ID <- factor(hash_seurat_ri1_selected$hash.ID, levels = c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

DimPlot(hash_seurat_ri1_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) &
  xlab("") &
  ylab("")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_order_ri1_2.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_ri1) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_ri1, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/clusters_deg_ri1_2.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_ri1,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_dotplot_ri1_2.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_ri1,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_heatmap_ri1_2.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_ri1, "results/timepoints/hto_singlets/targets/timepoints_ri1_seurat.rds")

p1 <- VlnPlot(hash_seurat_ri1, features = c("IFNB1"), slot = "data")
p2 <- VlnPlot(hash_seurat_ri1, features = c("IFNB1"), slot = "counts")

p1 | p2
ggsave("results/timepoints/hto_singlets/targets/ifnb1_vlnplot_ri1.pdf", height = 4, width = 10)


# timepoint DEG
Idents(hash_seurat_ri1_selected) <- "hash.ID"
all_deg <- FindAllMarkers(hash_seurat_ri1_selected, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/timepoints_deg_ri1.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top timepoint DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_ri1_selected,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_dotplot_ri1.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_ri1_selected,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_heatmap_ri1.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_ri1, "results/timepoints/hto_singlets/targets/hash_seurat_ri1.rds")


## ----------------------------



## JURKAT

hash_seurat_jurkat <- subset(hash_seurat_singlet, seurat_clusters %in% c("7", "14"))

# Select the top 1000 most variable features
hash_seurat_jurkat <- FindVariableFeatures(hash_seurat_jurkat, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_jurkat <- CellCycleScoring(hash_seurat_jurkat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_jurkat <- ScaleData(hash_seurat_jurkat, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_jurkat))

# Run PCA
hash_seurat_jurkat <- RunPCA(hash_seurat_jurkat, features = VariableFeatures(hash_seurat_jurkat))

ElbowPlot(hash_seurat_jurkat, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_jurkat <- FindNeighbors(hash_seurat_jurkat, reduction = "pca", dims = 1:20)
hash_seurat_jurkat <- FindClusters(hash_seurat_jurkat, resolution = 0.4, verbose = FALSE)
hash_seurat_jurkat <- RunUMAP(hash_seurat_jurkat, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_jurkat, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_jurkat.png", width = 6, height = 4)

DimPlot(hash_seurat_jurkat, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_cycle_jurkat.png", width = 6, height = 4)

DimPlot(hash_seurat_jurkat, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_jurkat.png", width = 14, height = 12)

DimPlot(hash_seurat_jurkat, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_batch_jurkat.png", width = 6, height = 4)


# selected conditions

hash_seurat_jurkat_selected <- subset(hash_seurat_jurkat, hash.ID %in% c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

hash_seurat_jurkat_selected$hash.ID <- factor(hash_seurat_jurkat_selected$hash.ID, levels = c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

DimPlot(hash_seurat_jurkat_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) &
  xlab("") &
  ylab("")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_order_jurkat.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_jurkat) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_jurkat, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/clusters_deg_jurkat.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_jurkat,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_dotplot_jurkat.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_jurkat,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_heatmap_jurkat.pdf", p, height = 15, width = 15)



# remove outlier clusters


hash_seurat_jurkat <- subset(hash_seurat_jurkat, seurat_clusters %in% c("0", "1", "2"))

# Select the top 1000 most variable features
hash_seurat_jurkat <- FindVariableFeatures(hash_seurat_jurkat, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_jurkat <- CellCycleScoring(hash_seurat_jurkat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_jurkat <- ScaleData(hash_seurat_jurkat, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_jurkat))

# Run PCA
hash_seurat_jurkat <- RunPCA(hash_seurat_jurkat, features = VariableFeatures(hash_seurat_jurkat))

ElbowPlot(hash_seurat_jurkat, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_jurkat <- FindNeighbors(hash_seurat_jurkat, reduction = "pca", dims = 1:20)
hash_seurat_jurkat <- FindClusters(hash_seurat_jurkat, resolution = 0.4, verbose = FALSE)
hash_seurat_jurkat <- RunUMAP(hash_seurat_jurkat, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_jurkat, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_jurkat_2.png", width = 6, height = 4)

DimPlot(hash_seurat_jurkat, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_cycle_jurkat_2.png", width = 6, height = 4)

DimPlot(hash_seurat_jurkat, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_jurkat_2.png", width = 14, height = 12)

DimPlot(hash_seurat_jurkat, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/targets/umap_batch_jurkat_2.png", width = 6, height = 4)


# selected conditions

hash_seurat_jurkat_selected <- subset(hash_seurat_jurkat, hash.ID %in% c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

hash_seurat_jurkat_selected$hash.ID <- factor(hash_seurat_jurkat_selected$hash.ID, levels = c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

DimPlot(hash_seurat_jurkat_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(8))) &
  xlab("") &
  ylab("")
ggsave("results/timepoints/hto_singlets/targets/umap_conditions_order_jurkat_2.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_jurkat) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_jurkat, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/clusters_deg_jurkat_2.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_jurkat,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_dotplot_jurkat_2.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_jurkat,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/cluster_deg_heatmap_jurkat_2.pdf", p, height = 15, width = 15)


# timepoint DEG
Idents(hash_seurat_jurkat_selected) <- "hash.ID"
all_deg <- FindAllMarkers(hash_seurat_jurkat_selected, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/targets/timepoints_deg_jurkat.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top timepoint DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_jurkat_selected,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_dotplot_jurkat.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_jurkat_selected,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_heatmap_jurkat.pdf", p, height = 15, width = 15)


# Save object
saveRDS(hash_seurat_jurkat, "results/timepoints/hto_singlets/targets/timepoints_jurkat_seurat.rds")


## ----------------------------------------------

## timepoint DEG vs 0h

# function
de_vs0h <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0.05, LOGFC=0.1){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST)
  
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  
  result$timepoint <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  result$p_adj <- p.adjust(result$p_val, method = "fdr")
  return(result)
}

timepoints <- c("1h", "3h", "6h", "12h", "24h")

# K562
hash_seurat_k562 <- readRDS("results/timepoints/hto_singlets/targets/timepoints_k562_seurat.rds")

hash_seurat_k562_selected <- subset(hash_seurat_k562, hash.ID %in% c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

hash_seurat_k562_selected$hash.ID <- factor(hash_seurat_k562_selected$hash.ID, levels =  c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

Idents(hash_seurat_k562_selected) <- "hash.ID"
result_target <- lapply(paste0("K562-", timepoints), de_vs0h, DATA = hash_seurat_k562_selected, IDENT2 = "K562-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_k562.txt", row.names = F, quote = F, sep = "\t")


# JURKAT
hash_seurat_jurkat <- readRDS("results/timepoints/hto_singlets/targets/timepoints_jurkatseurat.rds")

hash_seurat_jurkat_selected <- subset(hash_seurat_jurkat, hash.ID %in% c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

hash_seurat_jurkat_selected$hash.ID <- factor(hash_seurat_jurkat_selected$hash.ID, levels =  c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

Idents(hash_seurat_jurkat_selected) <- "hash.ID"
result_target <- lapply(paste0("JURKAT-", timepoints), de_vs0h, DATA = hash_seurat_jurkat_selected, IDENT2 = "JURKAT-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_jurkat.txt", row.names = F, quote = F, sep = "\t")


# RI1
hash_seurat_ri1 <- readRDS("results/timepoints/hto_singlets/targets/timepoints_ri1_seurat.rds")

hash_seurat_ri1_selected <- subset(hash_seurat_ri1, hash.ID %in% c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

hash_seurat_ri1_selected$hash.ID <- factor(hash_seurat_ri1_selected$hash.ID, levels =  c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

Idents(hash_seurat_ri1_selected) <- "hash.ID"
result_target <- lapply(paste0("RI1-", timepoints), de_vs0h, DATA = hash_seurat_ri1_selected, IDENT2 = "RI1-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_ri1.txt", row.names = F, quote = F, sep = "\t")


# GDM1
hash_seurat_gdm1 <- readRDS("results/timepoints/hto_singlets/targets/timepoints_gdm1_seurat.rds")

hash_seurat_gdm1_selected <- subset(hash_seurat_gdm1, hash.ID %in% c("GDM1-0h", "GDM1-1h", "GDM1-3h", "GDM1-6h", "GDM1-12h", "GDM1-24h"))

hash_seurat_gdm1_selected$hash.ID <- factor(hash_seurat_gdm1_selected$hash.ID, levels =  c("GDM1-0h", "GDM1-1h", "GDM1-3h", "GDM1-6h", "GDM1-12h", "GDM1-24h"))

Idents(hash_seurat_gdm1_selected) <- "hash.ID"
result_target <- lapply(paste0("GDM1-", timepoints), de_vs0h, DATA = hash_seurat_gdm1_selected, IDENT2 = "GDM1-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_gdm1.txt", row.names = F, quote = F, sep = "\t")


## ----------------------------------------------

## timepoint DEG vs 0h (all genes without thresholds)

# function
de_vs0h <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0.00, LOGFC=0){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST,
                        min.pct = 0)
  
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  
  result$timepoint <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  result$p_adj <- p.adjust(result$p_val, method = "fdr")
  return(result)
}

timepoints <- c("1h", "3h", "6h", "12h", "24h")

# K562
hash_seurat_k562 <- readRDS("results/timepoints/hto_singlets/targets/timepoints_k562_seurat.rds")

hash_seurat_k562_selected <- subset(hash_seurat_k562, hash.ID %in% c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

hash_seurat_k562_selected$hash.ID <- factor(hash_seurat_k562_selected$hash.ID, levels =  c("K562-0h", "K562-1h", "K562-3h", "K562-6h", "K562-12h", "K562-24h"))

Idents(hash_seurat_k562_selected) <- "hash.ID"
result_target <- lapply(paste0("K562-", timepoints), de_vs0h, DATA = hash_seurat_k562_selected, IDENT2 = "K562-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_k562.txt", row.names = F, quote = F, sep = "\t")


# JURKAT
hash_seurat_jurkat <- readRDS("results/timepoints/hto_singlets/targets/timepoints_jurkat_seurat.rds")

hash_seurat_jurkat_selected <- subset(hash_seurat_jurkat, hash.ID %in% c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

hash_seurat_jurkat_selected$hash.ID <- factor(hash_seurat_jurkat_selected$hash.ID, levels =  c("JURKAT-0h", "JURKAT-1h", "JURKAT-3h", "JURKAT-6h", "JURKAT-12h", "JURKAT-24h"))

Idents(hash_seurat_jurkat_selected) <- "hash.ID"
result_target <- lapply(paste0("JURKAT-", timepoints), de_vs0h, DATA = hash_seurat_jurkat_selected, IDENT2 = "JURKAT-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_jurkat.txt", row.names = F, quote = F, sep = "\t")


# RI1
hash_seurat_ri1 <- readRDS("results/timepoints/hto_singlets/targets/timepoints_ri1_seurat.rds")

hash_seurat_ri1_selected <- subset(hash_seurat_ri1, hash.ID %in% c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

hash_seurat_ri1_selected$hash.ID <- factor(hash_seurat_ri1_selected$hash.ID, levels =  c("RI1-0h", "RI1-1h", "RI1-3h", "RI1-6h", "RI1-12h", "RI1-24h"))

Idents(hash_seurat_ri1_selected) <- "hash.ID"
result_target <- lapply(paste0("RI1-", timepoints), de_vs0h, DATA = hash_seurat_ri1_selected, IDENT2 = "RI1-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_ri1.txt", row.names = F, quote = F, sep = "\t")


# GDM1
hash_seurat_gdm1 <- readRDS("results/timepoints/hto_singlets/targets/timepoints_gdm1_seurat.rds")

hash_seurat_gdm1_selected <- subset(hash_seurat_gdm1, hash.ID %in% c("GDM1-0h", "GDM1-1h", "GDM1-3h", "GDM1-6h", "GDM1-12h", "GDM1-24h"))

hash_seurat_gdm1_selected$hash.ID <- factor(hash_seurat_gdm1_selected$hash.ID, levels =  c("GDM1-0h", "GDM1-1h", "GDM1-3h", "GDM1-6h", "GDM1-12h", "GDM1-24h"))

Idents(hash_seurat_gdm1_selected) <- "hash.ID"
result_target <- lapply(paste0("GDM1-", timepoints), de_vs0h, DATA = hash_seurat_gdm1_selected, IDENT2 = "GDM1-0h") %>% bind_rows()
fwrite(result_target, "results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_gdm1.txt", row.names = F, quote = F, sep = "\t")


## ----------------------------------------

# Top DEG in all cell lines + selected type I IFN genes

k562 <- fread("results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_k562.txt") %>% mutate(cell_line = "K562")
jurkat <- fread("results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_jurkat.txt") %>% mutate(cell_line = "JURKAT")
gdm1 <- fread("results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_gdm1.txt") %>% mutate(cell_line = "GDM1")
ri1 <- fread("results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_ri1.txt") %>% mutate(cell_line = "RI1")

k562 <- fread("../../nk_hashing/results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_k562.txt") %>% mutate(cell_line = "K562")
jurkat <- fread("../../nk_hashing/results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_jurkat.txt") %>% mutate(cell_line = "JURKAT")
gdm1 <- fread("../../nk_hashing/results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_gdm1.txt") %>% mutate(cell_line = "GDM1")
ri1 <- fread("../../nk_hashing/results/timepoints/hto_singlets/targets/timepoints_deg_vs0h_all_ri1.txt") %>% mutate(cell_line = "RI1")

data_vs0h <- rbind(k562, jurkat, gdm1, ri1) %>% 
  mutate(timepoint = gsub(".*-", "", timepoint)) %>% 
  mutate(timepoint = factor(timepoint, levels = c("1h", "3h", "6h", "12h", "24h"))) 

topgenes_all <- data_vs0h %>% 
  group_by(timepoint, gene) %>% 
  filter(p_val_adj < 0.05) %>% 
  summarize(count = n(), mean_avg_log2FC = mean(avg_log2FC)) %>% 
  group_by(timepoint, gene) %>% 
  filter(count > 3) %>% 
  group_by(timepoint) %>% 
  slice_max(mean_avg_log2FC, n = 10)

toptimepoint <- topgenes_all %>% 
  group_by(gene) %>% 
  slice_max(n = 1, mean_avg_log2FC) %>% 
  arrange(timepoint, desc(mean_avg_log2FC)) %>% 
  select(gene, timepoint) %>% 
  mutate(column = "Peak timepoint")

ifn1genes <- c("IFNB1", "MX2", "OAS1", "OAS2", "IRF7", "IFI44L")

toptimepoint_ifn1 <-  data_vs0h %>% 
  filter(gene %in% ifn1genes) %>% 
  group_by(timepoint, gene) %>% 
  summarize(count = n(), mean_avg_log2FC = mean(avg_log2FC)) %>% 
  group_by(gene) %>% 
  slice_max(n = 1, mean_avg_log2FC) %>% 
  arrange(timepoint, desc(mean_avg_log2FC)) %>% 
  select(gene, timepoint) %>% 
  mutate(column = "Peak timepoint")


genelist <- toptimepoint %>% 
  select(gene) %>%
  unique() %>%
  tibble::deframe()

genelist_final <- c(genelist, toptimepoint_ifn1$gene)

df <- data_vs0h %>% 
  filter(gene %in% genelist_final) %>% 
  mutate(gene = factor(gene, levels = genelist_final)) %>% 
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-310, p_val_adj),
         p_val = ifelse(p_val == 0, 1e-310, p_val)) %>% 
  mutate(cell_line = factor(cell_line, levels = c("JURKAT", "K562", "GDM1", "RI1")))


p_dotplot <- ggplot(df, aes(x = cell_line, y = gene, size = -log10(p_val_adj), fill = avg_log2FC)) +
  geom_point(data = df[!grepl("Cancer|HALLMARK", df$gene),], pch = 21, color = "white") +
  scale_fill_distiller("Fold\nchange\n(log2)", palette = "RdBu", 
                       type = "div", limits = quantile(abs(as.numeric(df$avg_log2FC[!grepl("HALLMARK", df$gene)])), 0.95) * c(-1, 1), oob = squish,
                       guide = guide_colorbar(title.position = "top")) +
  geom_point(data = df[!grepl("Cancer|HALLMARK", df$gene),], pch = 1, color = ifelse(df[!grepl("Cancer|HALLMARK", df$gene),]$p_val_adj < 0.05, "grey50", "white")) +
  ggnewscale::new_scale_fill() +
  geom_point(data = df[grepl("IFNB1", df$gene),], aes(fill = avg_log2FC), pch = 21, color = "white") +
  scale_fill_distiller("IFNB1\nfold\nchange\n(log2)", palette = "YlGn", values = seq(0, 1, length.out = 11), direction = 1,
                       type = "div", limits = c(0, max(abs(as.numeric(df$avg_log2FC[grepl("IFNB1", df$gene)])))),
                       guide = guide_colorbar(title.position = "top")) +
  scale_size("Adjusted\np value\n(-log10)", range = c(2, 4.5), guide = guide_legend(title.position = "top",
                                                                                    override.aes = list(pch = 16, color = "black"))) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "bottom", drop = F) +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(2.5,0,0,0), "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  ylab("") +
  xlab("") +
  facet_grid(cols = vars(timepoint))

df_toptimepoint <- rbind(toptimepoint, toptimepoint_ifn1) %>% 
  filter(gene %in% genelist_final) %>% 
  mutate(gene = factor(gene, levels = rev(genelist_final)))

p_toptimepoint <- ggplot(df_toptimepoint, aes(x = column, y = gene, color = timepoint)) +
  geom_point() +
  scale_color_manual("Peak timepoint", values = rev(pals::ocean.matter(6)[2:6])) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "vertical")


p_toptimepoint | p_dotplot

ggsave("results/timepoints/hto_singlets/targets/timepoints_deg_dotplot_ifn1_small_fig2e.pdf", height = 7.75, width = 8.5)



