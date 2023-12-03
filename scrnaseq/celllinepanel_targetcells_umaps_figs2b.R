
# NK cell co-culture with panel of 26 cell lines, hashing scRNA-seq analysis
# Combine target cells from original and revision experiments and plot UMAPs (Figure S2B)

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
library(viridis)


## Load functions
me = "user"
source("scripts/fun_helper.R")
source("scripts/fun_getGenes.R")

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))


## Get cell cycle scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


dir.create("results/celllinepanel_combined/targets")
dir.create("results/celllinepanel_combined/targets/all")

"%ni%" <- Negate("%in%")

hash_seurat_singlet <- readRDS("results/celllinepanel_combined/hash_seurat_singlet.rds")

# remove NK cell clusters
hash_seurat_targets <- subset(hash_seurat_singlet, seurat_clusters %ni% c("0", "23", "28"))

# function to extract one cell line, normalize and compute and plot UMAP, save as separate object
plotTargetUMAPs <- function(CELLLINE){
  
  # select cell line
  hash_seurat_target <- subset(hash_seurat_targets, hash.ID %in% c(CELLLINE, paste0(CELLLINE, "-NK1"), paste0(CELLLINE, "-NK1-PBMC"), paste0(CELLLINE, "-NK2"),
                                                                            paste0(CELLLINE, "-NK3"), paste0(CELLLINE, "-NK4"), paste0(CELLLINE, "-NK5"), paste0(CELLLINE, "-NK6")))
  
  # Select the top 1000 most variable features
  hash_seurat_target <- FindVariableFeatures(hash_seurat_target, selection.method = "mean.var.plot")
  
  # Scale RNA data and regress out cell cycle
  hash_seurat_target <- CellCycleScoring(hash_seurat_target, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  hash_seurat_target <- ScaleData(hash_seurat_target, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_target))
  
  # Run PCA
  hash_seurat_target <- RunPCA(hash_seurat_target, features = VariableFeatures(hash_seurat_target))
  
  # Select the top 20 PCs for clustering and UMAP
  hash_seurat_target <- FindNeighbors(hash_seurat_target, reduction = "pca", dims = 1:20)
  hash_seurat_target <- FindClusters(hash_seurat_target, resolution = 0.4, verbose = FALSE)
  hash_seurat_target <- RunUMAP(hash_seurat_target, reduction = "pca", dims = 1:20)
  
  ## Visualize
  DimPlot(hash_seurat_target, group.by = "seurat_clusters", label = T, repel = T) +
    theme_bw(base_size = 12) +
    scale_color_manual(values = rev(getPalette(length(unique(hash_seurat_target$seurat_clusters))))) +
    labs("UMAP 1", "UMAP 2") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())
  ggsave(paste0("results/celllinepanel_combined/targets/all/umap_", CELLLINE, ".png"), width = 6, height = 4)
  
  p_condition <- DimPlot(hash_seurat_target, group.by = "hash.ID", label = T, repel = T) +
    theme_bw(base_size = 12) +
    scale_color_manual(values = rev(getPalette(8))) +
    labs("UMAP 1", "UMAP 2") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())
  ggsave(paste0("results/celllinepanel_combined/targets/all/umap_conditioncolor_", CELLLINE, ".png"), p_condition, width = 6, height = 4)
  
  DimPlot(hash_seurat_target, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
  ggsave(paste0("results/celllinepanel_combined/targets/all/umap_cycle_", CELLLINE, ".png"), width = 6, height = 4)
  
  DimPlot(hash_seurat_target, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          strip.background = element_blank()) +
    scale_color_manual(values = rev(getPalette(length(unique(hash_seurat_target$seurat_clusters))))) +
    labs("UMAP 1", "UMAP 2")
  ggsave(paste0("results/celllinepanel_combined/targets/all/umap_conditions_", CELLLINE, ".png"), width = 10, height = 3)
  
  # DEGs of clusters
  Idents(hash_seurat_target) <- "seurat_clusters"
  all_deg <- FindAllMarkers(hash_seurat_target, test.use = "t", only.pos = T, return.thresh = 0.05)
  fwrite(all_deg, paste0("results/celllinepanel_combined/targets/all/clusters_deg_", CELLLINE, ".txt"), sep = "\t", quote = F, row.names = F)
  
  # Dot plot of top cluster DEGs
  top10 <- all_deg %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  DotPlot(hash_seurat_target,
          features = unique(top10$gene),
          cols = "RdBu") + RotatedAxis()
  ggsave(paste0("results/celllinepanel_combined/targets/all/cluster_deg_dotplot_", CELLLINE, ".pdf"), height = 4, width = 18)
  
  # Heatmap of top cluster DEGs
  top25 <- all_deg %>%
    group_by(cluster) %>%
    top_n(n = 25, wt = avg_log2FC)
  
  p <- DoHeatmap(hash_seurat_target,
                 features = top25$gene,
                 raster = F,
                 angle = 0) + NoLegend()
  ggsave(paste0("results/celllinepanel_combined/targets/all/cluster_deg_heatmap_", CELLLINE, ".pdf"), p, height = 15, width = 15)
  
  # Save object
  saveRDS(hash_seurat_target, paste0("results/celllinepanel_combined/targets/all/hash_seurat_target_", CELLLINE, ".rds"))
  
  p_condition
  
}

# run function on all cell lines
celllines <- unique(gsub("\\-.*", "", hash_seurat_singlet$hash.ID))
celllines <- celllines[!grepl("^NK", celllines)]
celllines <- celllines[!is.na(celllines)]


p1 <- lapply(celllines, plotTargetUMAPs)
m1 <- marrangeGrob(p1, ncol = 7, nrow = 4, top = NULL)
ggsave("results/celllinepanel_combined/targets/all/umap_conditioncolor_all.pdf", m1, height = 10, width = 35)


## --------------------------------------------------


## Remove small outlier clusters and repeat 

dir.create("results/celllinepanel_combined/targets/mainclusters")

plotTargetUMAPsMainClusters <- function(CELLLINE){
  
  # select cell line
  hash_seurat_target <- subset(hash_seurat_targets, hash.ID %in% c(CELLLINE, paste0(CELLLINE, "-NK1"), paste0(CELLLINE, "-NK1-PBMC"), paste0(CELLLINE, "-NK2"),
                                                                            paste0(CELLLINE, "-NK3"), paste0(CELLLINE, "-NK4"), paste0(CELLLINE, "-NK5"), paste0(CELLLINE, "-NK6")))
  
  # subset data to include only clusters comprising > 5% of all target cells of one cell line
  cluster_cellcounts <- table(hash_seurat_target$seurat_clusters)
  maxclusters <- names(cluster_cellcounts[cluster_cellcounts>0.05*sum(cluster_cellcounts)])
  hash_seurat_target <- subset(hash_seurat_target, seurat_clusters %in% maxclusters)
  
  # Select the top 1000 most variable features
  hash_seurat_target <- FindVariableFeatures(hash_seurat_target, selection.method = "mean.var.plot")
  
  # Scale RNA data and regress out cell cycle
  hash_seurat_target <- CellCycleScoring(hash_seurat_target, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  hash_seurat_target <- ScaleData(hash_seurat_target, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_target))
  
  # Run PCA
  hash_seurat_target <- RunPCA(hash_seurat_target, features = VariableFeatures(hash_seurat_target))
  
  # Select the top 20 PCs for clustering and UMAP
  hash_seurat_target <- FindNeighbors(hash_seurat_target, reduction = "pca", dims = 1:20)
  hash_seurat_target <- FindClusters(hash_seurat_target, resolution = 0.4, verbose = FALSE)
  hash_seurat_target <- RunUMAP(hash_seurat_target, reduction = "pca", dims = 1:20)
  
  ## Visualize
  DimPlot(hash_seurat_target, group.by = "seurat_clusters", label = T, repel = T) +
    theme_bw(base_size = 12) +
    scale_color_manual(values = rev(getPalette(length(unique(hash_seurat_target$seurat_clusters))))) +
    labs("UMAP 1", "UMAP 2") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())
  ggsave(paste0("results/celllinepanel_combined/targets/mainclusters/umap_", CELLLINE, ".png"), width = 6, height = 4)
  
  p_condition <- DimPlot(hash_seurat_target, group.by = "hash.ID", label = T, repel = T) +
    theme_bw(base_size = 12) +
    scale_color_manual(values = rev(getPalette(8))) +
    labs("UMAP 1", "UMAP 2") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())
  ggsave(paste0("results/celllinepanel_combined/targets/mainclusters/umap_conditioncolor_", CELLLINE, ".png"), p_condition, width = 6, height = 4)
  
  DimPlot(hash_seurat_target, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
  ggsave(paste0("results/celllinepanel_combined/targets/mainclusters/umap_cycle_", CELLLINE, ".png"), width = 6, height = 4)
  
  DimPlot(hash_seurat_target, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          strip.background = element_blank()) +
    scale_color_manual(values = rev(getPalette(length(unique(hash_seurat_target$seurat_clusters))))) +
    labs("UMAP 1", "UMAP 2")
  ggsave(paste0("results/celllinepanel_combined/targets/mainclusters/umap_conditions_", CELLLINE, ".png"), width = 10, height = 3)
  
  # DEGs of clusters
  Idents(hash_seurat_target) <- "seurat_clusters"
  all_deg <- FindAllMarkers(hash_seurat_target, test.use = "t", only.pos = T, return.thresh = 0.05)
  fwrite(all_deg, paste0("results/celllinepanel_combined/targets/mainclusters/clusters_deg_", CELLLINE, ".txt"), sep = "\t", quote = F, row.names = F)
  
  # Dot plot of top cluster DEGs
  top10 <- all_deg %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  DotPlot(hash_seurat_target,
          features = unique(top10$gene),
          cols = "RdBu") + RotatedAxis()
  ggsave(paste0("results/celllinepanel_combined/targets/mainclusters/cluster_deg_dotplot_", CELLLINE, ".pdf"), height = 4, width = 18)
  
  # Heatmap of top cluster DEGs
  top25 <- all_deg %>%
    group_by(cluster) %>%
    top_n(n = 25, wt = avg_log2FC)
  
  p <- DoHeatmap(hash_seurat_target,
                 features = top25$gene,
                 raster = F,
                 angle = 0) + NoLegend()
  ggsave(paste0("results/celllinepanel_combined/targets/mainclusters/cluster_deg_heatmap_", CELLLINE, ".pdf"), p, height = 15, width = 15)
  
  # Save object
  saveRDS(hash_seurat_target, paste0("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_", CELLLINE, ".rds"))
  
  ## ----------------------------
  
  p_condition
  
}

# run function on all cell lines
p1 <- lapply(celllines, plotTargetUMAPsMainClusters)
m1 <- marrangeGrob(p1, ncol = 7, nrow = 4, top = NULL)
ggsave("results/celllinepanel_combined/targets/mainclusters/umap_conditioncolor_all.pdf", m1, height = 10, width = 35)


## ---------------------------------------------


## Plot target cell UMAPs

x697 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_697.rds")
kasumi2 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_KASUMI2.rds")
rchacv <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_RCHACV.rds")
lp1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_LP1.rds")
jjn3 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_JJN3.rds")
pl21 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_PL21.rds")
skm1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_SKM1.rds")
mm1s <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_MM1S.rds")
monomac1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_MONOMAC1.rds")
molm13 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_MOLM13.rds")
allsil <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_ALLSIL.rds")
molm14 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_MOLM14.rds")
granta519 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_GRANTA519.rds")
nudhl1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_NUDHL1.rds")
supt11 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_SUPT11.rds")
nalm6 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_NALM6.rds")
amo1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_AMO1.rds")
l363 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_L363.rds")
dnd41 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_DND41.rds")
thp1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_THP1.rds")
sudhl4 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_SUDHL4.rds")
ocim1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_OCIM1.rds")
jurkat <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_JURKAT.rds")
ri1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_RI1.rds")
k562 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_K562.rds")
gdm1 <- readRDS("results/celllinepanel_combined/targets/mainclusters/hash_seurat_target_GDM1.rds")

# merge and save object for pseudobulk DE analysis
objects <- list(x697, kasumi2, lp1, jjn3, rchacv, monomac1, pl21, mm1s, nalm6, skm1, allsil, molm13, molm14,
                supt11, thp1, nudhl1, amo1, l363, dnd41, granta519, sudhl4, ri1, jurkat, ocim1, k562, gdm1)

cell.ids = c("697", "KASUMI2", "LP1", "JJN3", "RCHACV", "MONOMAC1", "PL21", "MM1S", "NALM6", "SKM1", "ALLSIL", "MOLM13",
             "MOLM14", "SUPT11", "THP1", "NUDHL1", "AMO1", "L363", "DND41", "GRANTA519", "SUDHL4", "RI1", "JURKAT", "OCIM1", "K562", "GDM1")

merge_seurat <- merge(objects[[1]], objects[-1], add.cell.ids = cell.ids)

# Save merged object with all target cells (small outlier clusters removed)
saveRDS(merge_seurat, "results/celllinepanel_combined/targets/hash_seurat_targets_combined_mainclusters.rds")

plot_DimPlot <- function(OBJECT){
  
  seurat <- OBJECT
  seurat <- subset(seurat, cells = seurat$barcode[!grepl("NK4|NK5|NK6", seurat$hash.ID)])
  
  DimPlot(seurat, group.by = "hash.ID", label = F, raster = T, pt.size = 5) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "black"),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    scale_color_manual(values = c("grey70", viridis(4, option = "magma", begin = 0.25, end = 0.75, direction = -1))) +
    ggtitle(unique(gsub("-.*", "", seurat$hash.ID))) +
    xlab("") +
    ylab("") +
    NoLegend()
}

p_filler <-  DimPlot(x697, group.by = "hash.ID", label = F, cols = rep("white", 8)) +
  theme_void() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_blank()) +
  xlab("") +
  ylab("") +
  NoLegend()

p <- lapply(c(x697, kasumi2, lp1, jjn3, rchacv, monomac1, pl21, mm1s, nalm6, skm1, allsil, molm13, molm14,
              supt11, thp1, nudhl1, amo1, l363, dnd41, granta519, sudhl4, ri1, jurkat, ocim1, k562, gdm1), plot_DimPlot)

p <- list(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]], p[[17]], p[[18]],
          p[[19]], p[[20]], p[[21]], p[[22]], p[[23]], p[[24]], p[[25]], p[[26]])

m <- marrangeGrob(p, top = NULL, ncol = 7, nrow = 4, layout_matrix = matrix(1:28, 4, 7, TRUE))

ggsave("results/celllinepanel_combined/targets/mainclusters/umap_conditions_manuscript.pdf", m, width = 14, height = 8)
ggsave("results/celllinepanel_combined/targets/mainclusters/umap_conditions_manuscript.png", m, width = 13, height = 7)

## target cell UMAP examples for Figure 2A)
(p[[26]] + ggtitle("GDM1 (AML)")) + (p[[25]] + ggtitle("K562 (CML)")) + (p[[1]] + ggtitle("697 (B-ALL)"))
ggsave("results/celllinepanel_combined/targets/mainclusters/umap_conditions_examples_manuscript.pdf", width = 6, height = 2)
