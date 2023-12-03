
# NK cell co-culture media transfer experiment
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
library(slingshot)
library(patchwork)


getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))

dir.create("results/media_transfer/hto_singlets/targets")

hash_seurat_singlet <- readRDS("results/media_transfer/mediatransfer_seurat.rds")


## K562

hash_seurat_k562 <- subset(hash_seurat_singlet, seurat_clusters %in% c("0", "7", "11", "12"))

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
ggsave("results/media_transfer/hto_singlets/targets/umap_k562.png", width = 6, height = 4)

DimPlot(hash_seurat_k562, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/media_transfer/hto_singlets/targets/umap_cycle_k562.png", width = 6, height = 4)

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
ggsave("results/media_transfer/hto_singlets/targets/umap_conditions_k562.png", width = 14, height = 12)

DimPlot(hash_seurat_k562, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/media_transfer/hto_singlets/targets/umap_batch_k562.png", width = 6, height = 4)


# selected conditions

hash_seurat_k562_selected <- subset(hash_seurat_k562, hash.ID %in% c("K562", "K562-NK", "K562-medium-from-K562-NK",
                                                                     "K562-medium-from-NK", "K562-medium-from-697-NK",
                                                                     "K562-medium-from-GDM1-NK"))

hash_seurat_k562_selected$hash.ID <- factor(hash_seurat_k562_selected$hash.ID, levels = c("K562", "K562-NK", "K562-medium-from-K562-NK",
                                                                                          "K562-medium-from-NK", "K562-medium-from-697-NK",
                                                                                          "K562-medium-from-GDM1-NK"))

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
ggsave("results/media_transfer/hto_singlets/targets/umap_conditions_order_k562.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_k562) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_k562, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/media_transfer/hto_singlets/targets/clusters_deg_k562.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_k562,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/media_transfer/hto_singlets/targets/cluster_deg_dotplot_k562.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_k562,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/media_transfer/hto_singlets/targets/cluster_deg_heatmap_k562.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_k562, "results/media_transfer/hto_singlets/targets/mediatransfer_k562_seurat.rds")

## ----------------------------

## GDM1

hash_seurat_gdm1 <- subset(hash_seurat_singlet, seurat_clusters %in% c("2", "5", "8", "13"))

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
ggsave("results/media_transfer/hto_singlets/targets/umap_gdm1.png", width = 6, height = 4)


# Remove tiny outlier cluster and re-normalize 

hash_seurat_gdm1 <- subset(hash_seurat_gdm1, seurat_clusters %in% c("0", "1", "2", "3", "4"))

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
ggsave("results/media_transfer/hto_singlets/targets/umap_gdm1_2.png", width = 6, height = 4)

DimPlot(hash_seurat_gdm1, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/media_transfer/hto_singlets/targets/umap_cycle_gdm1.png", width = 6, height = 4)

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
ggsave("results/media_transfer/hto_singlets/targets/umap_conditions_gdm1.png", width = 14, height = 12)

DimPlot(hash_seurat_gdm1, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/media_transfer/hto_singlets/targets/umap_batch_gdm1.png", width = 6, height = 4)


# selected conditions

hash_seurat_gdm1_selected <- subset(hash_seurat_gdm1, hash.ID %in% c("GDM1", "GDM1-NK", "GDM1-medium-from-GDM1-NK",
                                                                     "GDM1-medium-from-NK", "GDM1-medium-from-697-NK",
                                                                     "GDM1-medium-from-K562-NK"))

hash_seurat_gdm1_selected$hash.ID <- factor(hash_seurat_gdm1_selected$hash.ID, levels = c("GDM1", "GDM1-NK", "GDM1-medium-from-GDM1-NK",
                                                                                          "GDM1-medium-from-NK", "GDM1-medium-from-697-NK",
                                                                                          "GDM1-medium-from-K562-NK"))

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
ggsave("results/media_transfer/hto_singlets/targets/umap_conditions_order_gdm1.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_gdm1) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_gdm1, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/media_transfer/hto_singlets/targets/clusters_deg_gdm1.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_gdm1,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/media_transfer/hto_singlets/targets/cluster_deg_dotplot_gdm1.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_gdm1,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/media_transfer/hto_singlets/targets/cluster_deg_heatmap_gdm1.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_gdm1, "results/media_transfer/hto_singlets/targets/hash_seurat_gdm1.rds")

hash_seurat_gdm1_selected <- subset(hash_seurat_gdm1, hash.ID %in% c("GDM1", "GDM1-NK", "GDM1-medium-from-GDM1-NK",
                                                                     "GDM1-medium-from-NK", "GDM1-medium-from-697-NK",
                                                                     "GDM1-medium-from-K562-NK"))

hash_seurat_gdm1_selected$hash.ID <- factor(hash_seurat_gdm1_selected$hash.ID, levels = c("GDM1", "GDM1-NK", "GDM1-medium-from-GDM1-NK",
                                                                                          "GDM1-medium-from-NK", "GDM1-medium-from-697-NK",
                                                                                          "GDM1-medium-from-K562-NK"))

# IFNB1 violin plot
hash_seurat_gdm1 <- readRDS("results/media_transfer/hto_singlets/targets/mediatransfer_gdm1_seurat.rds")

p1 <- VlnPlot(hash_seurat_gdm1, features = c("IFNB1"), slot = "data")
p2 <- VlnPlot(hash_seurat_gdm1, features = c("IFNB1"), slot = "counts")

p1 | p2
ggsave("results/media_transfer/hto_singlets/targets/ifnb1_vlnplot_gdm1.pdf", height = 4, width = 10)


# IFNB1 violin plot split by condition
Idents(hash_seurat_gdm1_selected) <- "hash.ID"
p1 <- VlnPlot(hash_seurat_gdm1_selected, features = c("IFNB1"), slot = "data") + NoLegend() + xlab("") + theme(plot.title = element_text(face = "italic"))
p2 <- VlnPlot(hash_seurat_gdm1_selected, features = c("IFNB1"), slot = "counts") + NoLegend() + xlab("") + theme(plot.title = element_text(face = "italic"))

p1 | p2
ggsave("results/media_transfer/hto_singlets/targets/ifnb1_vlnplot_conditions_gdm1.pdf", height = 4, width = 10)


## ----------------------------



## 697

hash_seurat_697 <- subset(hash_seurat_singlet, seurat_clusters %in% c("1", "4", "9"))

# Select the top 1000 most variable features
hash_seurat_697 <- FindVariableFeatures(hash_seurat_697, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_697 <- CellCycleScoring(hash_seurat_697, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_697 <- ScaleData(hash_seurat_697, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_697))

# Run PCA
hash_seurat_697 <- RunPCA(hash_seurat_697, features = VariableFeatures(hash_seurat_697))

ElbowPlot(hash_seurat_697, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_697 <- FindNeighbors(hash_seurat_697, reduction = "pca", dims = 1:20)
hash_seurat_697 <- FindClusters(hash_seurat_697, resolution = 0.4, verbose = FALSE)
hash_seurat_697 <- RunUMAP(hash_seurat_697, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_697, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/media_transfer/hto_singlets/targets/umap_697.png", width = 6, height = 4)

DimPlot(hash_seurat_697, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/media_transfer/hto_singlets/targets/umap_cycle_697.png", width = 6, height = 4)

DimPlot(hash_seurat_697, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/media_transfer/hto_singlets/targets/umap_conditions_697.png", width = 14, height = 12)

DimPlot(hash_seurat_697, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/media_transfer/hto_singlets/targets/umap_batch_697.png", width = 6, height = 4)


# selected conditions

hash_seurat_697_selected <- subset(hash_seurat_697, hash.ID %in% c("697", "697-NK", "697-medium-from-697-NK",
                                                                     "697-medium-from-NK", "697-medium-from-K562-NK",
                                                                     "697-medium-from-GDM1-NK"))

hash_seurat_697_selected$hash.ID <- factor(hash_seurat_697_selected$hash.ID, levels = c("697", "697-NK", "697-medium-from-697-NK",
                                                                                          "697-medium-from-NK", "697-medium-from-K562-NK",
                                                                                          "697-medium-from-GDM1-NK"))

DimPlot(hash_seurat_697_selected, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 8) +
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
ggsave("results/media_transfer/hto_singlets/targets/umap_conditions_order_697.png", width = 12, height = 2.5)


# DEGs of clusters
Idents(hash_seurat_697) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_697, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/media_transfer/hto_singlets/targets/clusters_deg_697.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_697,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/media_transfer/hto_singlets/targets/cluster_deg_dotplot_697.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_697,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/media_transfer/hto_singlets/targets/cluster_deg_heatmap_697.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_697, "results/media_transfer/hto_singlets/targets/mediatransfer_697_seurat.rds")

## ------------------------------------


# Plot ith all UMAPs with conditions colored

hash_seurat_k562 <- readRDS("results/media_transfer/hto_singlets/targets/mediatransfer_k562_seurat.rds")
hash_seurat_gdm1 <- readRDS("results/media_transfer/hto_singlets/targets/mediatransfer_gdm1_seurat.rds")
hash_seurat_697 <- readRDS("results/media_transfer/hto_singlets/targets/mediatransfer_697_seurat.rds")


cols = c("grey50", LaCroixColoR::lacroix_palette("PeachPear", 6)[c(1,2,3,4,5)])
names(cols) <- c("Untreated", "NK-treated", "Medium from K562-NK co-culture",
                 "Medium from NK cells",
                 "Medium from GDM1-NK co-culture",
                 "Medium from 697-NK co-culture")

core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "STAT1", "IRF1", "IRF9", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "UBE2L6", "GNLY", "CCL5")




# K562

hash_seurat_k562_selected <- subset(hash_seurat_k562, hash.ID %in% c("K562", "K562-NK", "K562-medium-from-K562-NK",
                                                                     "K562-medium-from-NK", 
                                                                     "K562-medium-from-GDM1-NK",
                                                                     "K562-medium-from-697-NK"))

hash_seurat_k562_selected$hash.ID <- factor(hash_seurat_k562_selected$hash.ID, levels = c("K562", "K562-NK", "K562-medium-from-K562-NK",
                                                                                          "K562-medium-from-NK", 
                                                                                          "K562-medium-from-GDM1-NK",
                                                                                          "K562-medium-from-697-NK"),
                                            labels = c("Untreated", "NK-treated", "Medium from K562-NK co-culture",
                                                       "Medium from NK cells",
                                                       "Medium from GDM1-NK co-culture",
                                                       "Medium from 697-NK co-culture"))


hash_seurat_k562_selected <- AddModuleScore(hash_seurat_k562_selected, features = list(core), name = c("nk_response_score"))


p1 <- DimPlot(hash_seurat_k562_selected, group.by = "hash.ID", label = F, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = cols) +
  xlab("") +
  ylab("") +
  ggtitle("K562 (CML)") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5))

p4 <- FeaturePlot(hash_seurat_k562_selected, features = "nk_response_score1", label = F, max.cutoff = 0.8) +
  theme_bw(base_size = 12) +
  scale_color_distiller(palette = "RdBu") +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom", 
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5))
  



# GDM1

hash_seurat_gdm1_selected <- subset(hash_seurat_gdm1, hash.ID %in% c("GDM1", "GDM1-NK", "GDM1-medium-from-GDM1-NK",
                                                                     "GDM1-medium-from-NK", "GDM1-medium-from-697-NK",
                                                                     "GDM1-medium-from-K562-NK"))

hash_seurat_gdm1_selected$hash.ID <- factor(hash_seurat_gdm1_selected$hash.ID, levels = c("GDM1", "GDM1-NK", "GDM1-medium-from-K562-NK",
                                                                                          "GDM1-medium-from-NK", "GDM1-medium-from-GDM1-NK",
                                                                                          "GDM1-medium-from-697-NK"),
                                            labels = c("Untreated", "NK-treated", "Medium from K562-NK co-culture",
                                                                                                   "Medium from NK cells",
                                                                                                   "Medium from GDM1-NK co-culture",
                                                                                                   "Medium from 697-NK co-culture"))

hash_seurat_gdm1_selected <- AddModuleScore(hash_seurat_gdm1_selected, features = list(core), name = c("nk_response_score"))

p2 <- DimPlot(hash_seurat_gdm1_selected, group.by = "hash.ID", label = F, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = cols) +
  xlab("") +
  ylab("") +
  ggtitle("GDM1 (AML)") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5))


p5 <- FeaturePlot(hash_seurat_gdm1_selected, features = "nk_response_score1", label = F, max.cutoff = 0.8) +
  theme_bw(base_size = 12) +
  scale_color_distiller(palette = "RdBu") +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom", 
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5))


# 697

hash_seurat_697_selected <- subset(hash_seurat_697, hash.ID %in% c("697", "697-NK", "697-medium-from-697-NK",
                                                                   "697-medium-from-NK", "697-medium-from-K562-NK",
                                                                   "697-medium-from-GDM1-NK"))

hash_seurat_697_selected$hash.ID <- factor(hash_seurat_697_selected$hash.ID, levels = c("697", "697-NK", "697-medium-from-K562-NK",
                                                                                        "697-medium-from-NK", "697-medium-from-GDM1-NK",
                                                                                        "697-medium-from-697-NK"),
                                              labels = c("Untreated", "NK-treated", "Medium from K562-NK co-culture",
                                                       "Medium from NK cells",
                                                       "Medium from GDM1-NK co-culture",
                                                       "Medium from 697-NK co-culture"))

hash_seurat_697_selected <- AddModuleScore(hash_seurat_697_selected, features = list(core), name = c("nk_response_score"))


p3 <- DimPlot(hash_seurat_697_selected, group.by = "hash.ID", label = F, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = cols) +
  xlab("") +
  ylab("") +
  ggtitle("697 (B-ALL)") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))

p6 <- FeaturePlot(hash_seurat_697_selected, features = "nk_response_score1", label = F, max.cutoff = 0.7, min.cutoff = -0.1) +
  theme_bw(base_size = 12) +
  scale_color_distiller(palette = "RdBu") +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom", 
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5))


(p1 + p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")) /
  (p4 + p5 + p6 + plot_layout(ncol = 3))

ggsave("results/media_transfer/hto_singlets/targets/umaps_all_celllines.pdf", height = 5.5, width = 7)



