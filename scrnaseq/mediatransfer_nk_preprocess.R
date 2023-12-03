
# NK cell co-culture media transfer experiment, hashing scRNA-seq analysis
# NK cell UMAPs


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

dir.create("results/media_transfer/hto_singlets/nk")

hash_seurat_singlet <- readRDS("results/media_transfer/mediatransfer_seurat.rds")

hash_seurat_nk <- subset(hash_seurat_singlet, seurat_clusters %in% c("3", "10", "14"))

# Select the top 1000 most variable features
hash_seurat_nk <- FindVariableFeatures(hash_seurat_nk, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_nk <- CellCycleScoring(hash_seurat_nk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_nk <- ScaleData(hash_seurat_nk, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_nk))

# Run PCA
hash_seurat_nk <- RunPCA(hash_seurat_nk, features = VariableFeatures(hash_seurat_nk))

ElbowPlot(hash_seurat_nk, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_nk <- FindNeighbors(hash_seurat_nk, reduction = "pca", dims = 1:20)
hash_seurat_nk <- FindClusters(hash_seurat_nk, resolution = 0.4, verbose = FALSE)
hash_seurat_nk <- RunUMAP(hash_seurat_nk, reduction = "pca", dims = 1:20)


## Visualize
DimPlot(hash_seurat_nk, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/media_transfer/hto_singlets/nk/umap.png", width = 6, height = 4)

DimPlot(hash_seurat_nk, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/media_transfer/hto_singlets/nk/umap_cycle.png", width = 6, height = 4)

DimPlot(hash_seurat_nk, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/media_transfer/hto_singlets/nk/umap_celllines.png", width = 14, height = 12)

DimPlot(hash_seurat_nk, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/media_transfer/hto_singlets/nk/umap_batch.png", width = 6, height = 4)

# DEGs of clusters
Idents(hash_seurat_nk) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_nk, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/media_transfer/hto_singlets/nk/clusters_deg.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_nk,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/media_transfer/hto_singlets/nk/cluster_deg_dotplot.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_nk,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/media_transfer/hto_singlets/nk/cluster_deg_heatmap.pdf", p, height = 15, width = 15)

# Save object
saveRDS(hash_seurat_nk, "results/media_transfer/hto_singlets/nk/mediatransfer_nk_seurat.rds")


