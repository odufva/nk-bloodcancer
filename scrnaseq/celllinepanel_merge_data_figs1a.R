
# NK cell co-culture with panel of 26 cell lines, hashing scRNA-seq analysis
# merge all data from original and revision experiments

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
library(limma)

## Load functions
me = "user"
source("scripts/fun_helper.R")
source("scripts/fun_getGenes.R")

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))


dir.create("results/celllinepanel_combined/")

"%ni%" <- Negate("%in%")

# original (donor NK1)
hash_seurat_singlet1 <- readRDS("results/celllinepanel_nk1/celllinepanel_nk1_seurat.rds")

# revision (donors NK2-NK6)
hash_seurat_singlet2 <- readRDS("results/celllinepanel_nk2nk3/celllinepanel_nk2nk3_seurat.rds")


hash_seurat_singlet <- merge(hash_seurat_singlet1, hash_seurat_singlet2)


# Select the top 1000 most variable features
hash_seurat_singlet <- FindVariableFeatures(hash_seurat_singlet, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_singlet <- CellCycleScoring(hash_seurat_singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_singlet <- ScaleData(hash_seurat_singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_singlet))

# Run PCA
hash_seurat_singlet <- RunPCA(hash_seurat_singlet, features = VariableFeatures(hash_seurat_singlet))

ElbowPlot(hash_seurat_singlet, ndims = 50, reduction = "pca")

# Select the top 20 PCs for clustering and UMAP
hash_seurat_singlet <- FindNeighbors(hash_seurat_singlet, reduction = "pca", dims = 1:20)
hash_seurat_singlet <- FindClusters(hash_seurat_singlet, resolution = 0.6, verbose = FALSE)
hash_seurat_singlet <- RunUMAP(hash_seurat_singlet, reduction = "pca", dims = 1:20)

## Visualize
DimPlot(hash_seurat_singlet, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(91))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/celllinepanel_combined/umap.png", width = 7, height = 5)

DimPlot(hash_seurat_singlet, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel_combined/umap_cycle.png", width = 6, height = 4)

hash_seurat_singlet$hash.ID <- gsub("NK-expanded", "NK1", gsub("NK-PBMC", "NK1-PBMC", hash_seurat_singlet$hash.ID))
hash_seurat_singlet$orig.ident <- gsub("NK_", "", gsub("_D14", "", gsub("D0", "PBMC", gsub("NK2587", "NK1", gsub("NK2770|FM2770", "NK2", gsub("NK2763|FM2763", "NK3",
                                  gsub("NK2768|FM2768", "NK4", gsub("NK2771|FM2771", "NK5", gsub("NK2931|FM2931", "NK6", hash_seurat_singlet$orig.ident)))))))))


hash_seurat_singlet$barcode[is.na(hash_seurat_singlet$hash.ID)]
cells <- hash_seurat_singlet$barcode[!is.na(hash_seurat_singlet$hash.ID)]

hash_seurat_singlet <- subset(hash_seurat_singlet, cells = cells)


# Save object
saveRDS(hash_seurat_singlet, "results/celllinepanel_combined/hash_seurat_singlet.rds")

## ---------------------------------------

# Plots 

DimPlot(hash_seurat_singlet, group.by = "hash.ID", label = T, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(155))) +
  labs("UMAP 1", "UMAP 2") +
  NoLegend()
ggsave("results/celllinepanel_combined/hto_class_umap_colored.pdf", height = 10, width = 16)

# Experiment
DimPlot(hash_seurat_singlet, split.by = "orig.ident", ncol = 3, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(41))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel_combined/umap_experiment.png", height = 16, width = 24)

# Experiment with sample labels
DimPlot(hash_seurat_singlet, split.by = "orig.ident", ncol = 3, group.by = "hash.ID", label = T, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(155))) +
  labs("UMAP 1", "UMAP 2") +
  NoLegend()
ggsave("results/celllinepanel_combined/umap_experiment_samples.png", height = 12, width = 16)

# SingleR
DimPlot(hash_seurat_singlet, group.by = "singler_blueprint_pred", label = T, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(23))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel_combined/singler_blueprint_pred.pdf", height = 7, width = 12)

DimPlot(hash_seurat_singlet, group.by = "singler_hpca_pred", label = T, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(41))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/celllinepanel_combined/singler_hpca_pred.pdf", height = 7, width = 12)

## ---------------------------------------


## UMAP of all cells in the experiment

hash_seurat_singlet <- readRDS("results/celllinepanel_combined/hash_seurat_singlet.rds")

hash_seurat_singlet$plotlabels <- gsub("-NK.*", "", gsub("^NK.*", " ", hash_seurat_singlet$hash.ID))
hash_seurat_singlet$plotlabels[hash_seurat_singlet$seurat_clusters %in% c("0", "23", "28")] <- "NK cells"

hash_seurat_singlet$treatment <- "Targets only"
hash_seurat_singlet$treatment[grepl("NK", hash_seurat_singlet$orig.ident)] <- "NK-treated"
hash_seurat_singlet$treatment <- factor(hash_seurat_singlet$treatment, levels = c("Targets only","NK-treated"))

p1 <- DimPlot(hash_seurat_singlet, group.by = "plotlabels", label = T, repel = T, raster = T, pt.size = 0.01) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(length(unique(hash_seurat_singlet$plotlabels))))) &
  ggtitle("Cell identity") &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  NoLegend()

# Plot Figure S1A)
p1
ggsave("results/celllinepanel_combined/umap_labels_manuscript.pdf", height = 6, width = 9)

p2 <- DimPlot(hash_seurat_singlet, group.by = "treatment", label = F, cols = c("grey70", viridis::inferno(3, begin = 0.25, end = 0.75, direction = -1)[2]), raster = T, pt.size = 0.01) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "bottom") +
  ggtitle("Treatment condition") &
  xlab("UMAP 1") &
  ylab("UMAP 2")

p1 + p2
ggsave("results/celllinepanel_combined/umap_celllines_conditions_manuscript.pdf", height = 6, width = 15)



