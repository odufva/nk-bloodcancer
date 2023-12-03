
# NK cell co-culture with edited K562/SUDHL4 hashing scRNA-seq analysis
# HTO data processed using new CITE-seq-count version

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
library(cluster)

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))


## Load functions
me = "user"
source("scripts/fun_helper.R")
source("scripts/fun_getGenes.R")

## Get data
files <- list.files(recursive = T, include.dirs = T) %>% 
  grep(pattern = "K562_SUDHL4_CRISPR", value = T) %>% 
  grep(pattern = "HTO", value = T, invert = T)

# get files for each sample/lane
files1 <- files %>% grep(pattern = "K562_SUDHL4_CRISPR", value = T)

# Create Seurat objects of each sample/lane and get sample/lane names from filenames
extractBatchName <- function(filenames){gsub("data/|_matrix.mtx.gz", "", filenames[grepl("mtx", filenames)])}

scrnaseq_files1 <- ReadMtx(files1[grepl("mtx", files1)], files1[grepl("barcodes", files1)], files1[grepl("features", files1)]) %>% CreateSeuratObject(project = extractBatchName(files1), min.cells = 3, min.features = 200)

hash_seurat    <- scrnaseq_files1


## Basic QC
dir.create("./results/k562_sudhl4_new/qc/", recursive = T, showWarnings = F)
dir.create("./results/k562_sudhl4_new/qc/before_qc1/", recursive = T, showWarnings = F)
dir.create("./results/k562_sudhl4_new/qc/after_qc1/", recursive = T, showWarnings = F)
dir.create("./results/k562_sudhl4_new/qc/hto/", recursive = T, showWarnings = F)

clonality_genes <- getClonalityGenes(hash_seurat)
unwanted_genes  <- getUnwantedGenes(hash_seurat)
hash_seurat     <- PercentageFeatureSet(hash_seurat, pattern = "^MT-", col.name = "percent.mt")
hash_seurat     <- PercentageFeatureSet(hash_seurat, pattern = "^RP", col.name = "percent.ribo")
hash_seurat     <- PercentageFeatureSet(hash_seurat, features = cycle.genes, col.name = "percent.cycle")
hash_seurat@meta.data$barcode <- colnames(hash_seurat)


## Plots and cluster DEG before QC
hash_seurat %>% plotQC(folder = "results/k562_sudhl4_new/qc/before_qc1/")

hash_seurat_beforeQC <- hash_seurat %>% preprocessSeurat(cells.to.use = colnames(hash_seurat))

## Get optimal clustering
hash_seurat_beforeQC <- hash_seurat_beforeQC %>% getClustering()
hash_seurat_beforeQC %>% plotClustering()
ggsave("results/k562_sudhl4_new/qc/before_qc1/scatter_clustering.png", width = 5, height = 4)


Idents(hash_seurat_beforeQC) <- hash_seurat_beforeQC$RNA_snn_res.0.3
hash_seurat_beforeQC$cluster <- Idents(hash_seurat_beforeQC)

## Visualize
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))

DimPlot(hash_seurat_beforeQC, group.by = "RNA_snn_res.0.3", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(15)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/k562_sudhl4_new/qc/before_qc1/umap.png", width = 5, height = 4)

all_deg <- FindAllMarkers(hash_seurat_beforeQC, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/k562_sudhl4_new/qc/before_qc1/all_deg.txt", sep = "\t", quote = F, row.names = F)

## --------------------------

## QC (no doublet removal, will be done using HTOs)

hash_seurat <- hash_seurat %>% getQC()

## Get SingleR predictions
hash_seurat <- hash_seurat %>% getSingler()
relevant_hpca_clusters <- hash_seurat@meta.data %>% group_by(singler_hpca_pred) %>% summarise(n=n()) %>% filter(n>20) %>% pull(singler_hpca_pred)
relevant_blue_clusters <- hash_seurat@meta.data %>% group_by(singler_blueprint_pred) %>% summarise(n=n()) %>% filter(n>20) %>% pull(singler_blueprint_pred)

hash_seurat$singler_hpca_pred      <- ifelse(hash_seurat$singler_hpca_pred %in% relevant_hpca_clusters, hash_seurat$singler_hpca_pred, "rare")
hash_seurat$singler_blueprint_pred <- ifelse(hash_seurat$singler_blueprint_pred %in% relevant_blue_clusters, hash_seurat$singler_blueprint_pred, "rare")


## Get doublets
hash_seurat <- hash_seurat %>% getDoublets()
hash_seurat$hybrid_doublet_score

hash_seurat@meta.data %>% 
  ggplot(aes(orig.ident, hybrid_doublet_score/2)) + geom_violin(draw_quantiles = 0.5)

hash_seurat@meta.data %>% 
  ggplot(aes(hybrid_doublet_score/2)) + geom_histogram()

table(hash_seurat$hybrid_doublet_score/2 > 0.90)

## Get Seurat object
hash_seurat <- hash_seurat %>% preprocessSeurat(cells.to.use = colnames(hash_seurat))

## Get optimal clustering
hash_seurat <- hash_seurat %>% getClustering()
hash_seurat %>% plotClustering()
ggsave("results/k562_sudhl4_new/qc/after_qc1/scatter_clustering.png", width = 5, height = 4)

Idents(hash_seurat) <- hash_seurat$RNA_snn_res.0.4
hash_seurat$cluster <- Idents(hash_seurat)

## Visualize
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))

DimPlot(hash_seurat, group.by = "RNA_snn_res.0.4", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(11)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/k562_sudhl4_new/qc/after_qc1/umap.png", width = 5, height = 4)


## Get cell cycle scores and visualize
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
hash_seurat_cycle <- CellCycleScoring(hash_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_cycle <- ScaleData(hash_seurat_cycle, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_cycle))

DimPlot(hash_seurat_cycle, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/k562_sudhl4_new/qc/after_qc1/umap_cycle.png", width = 5, height = 4)

DimPlot(hash_seurat, group.by = "singler_hpca_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(5)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/k562_sudhl4_new/qc/after_qc1/umap_singler_hpca.png", width = 7, height = 4)

DimPlot(hash_seurat, group.by = "singler_blueprint_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(7)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/k562_sudhl4_new/qc/after_qc1/umap_blue_blueprint.png", width = 7, height = 4)

FeaturePlot(hash_seurat, features = "hybrid_doublet_score", label = T, repel = T, cols = c("lightgrey", "tomato"), min.cutoff = 0) + theme_bw(base_size = 12)  + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/k562_sudhl4_new/qc/after_qc1/umap_dbl.png", width = 5, height = 4)

hash_seurat %>% plotQC(folder = "results/k562_sudhl4_new/qc/after_qc1/")

all_deg <- FindAllMarkers(hash_seurat, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/k562_sudhl4_new/qc/after_qc1/all_deg.txt", sep = "\t", quote = F, row.names = F)

## --------------------------

## Analysis using HTOs

## Get data
files_hto <- list.files(recursive = T, include.dirs = T) %>% 
  grep(pattern = "K562_SUDHL4_CRISPR", value = T) %>% 
  grep(pattern = "HTO", value = T)

# get files for each sample/lane
files1_hto <- files_hto %>% grep(pattern = "K562_SUDHL4_CRISPR", value = T) 

# Create Seurat objects of each sample/lane and get sample/lane names from filenames
extractBatchName <- function(filenames){gsub("data/|_matrix.mtx.gz", "", filenames[grepl("mtx", filenames)])}


## Load HTO count matrix
htos <- ReadMtx(files1_hto[grepl("mtx", files1_hto)], files1_hto[grepl("barcodes", files1_hto)], files1_hto[grepl("features", files1_hto)], feature.column = 1)
rowSums(htos)
htos <- htos[1:14,-1]

## Load HTO ids
hto_id <- fread("data/NK1_K562_SUDHL4_CRISPR_HTO_sample_list.txt", data.table = F)

# Replace HTO matrix rownames with ids
rownames(htos) <- hto_id$sample[match(gsub("^.*\\-", "", rownames(htos)), hto_id$sequence)]

# Modify column names to match expression matrix
colnames(htos) <- paste0(colnames(htos), "-1")

# Subset RNA and HTO counts by joint cell barcodes
joint.bcs <- intersect(colnames(hash_seurat), colnames(htos))
hash_seurat <- subset(hash_seurat, barcode %in% joint.bcs)
htos <- as.matrix(htos[, hash_seurat$barcode])

# Add HTO data as a new assay independent from RNA
hash_seurat[["HTO"]] <- CreateAssayObject(counts = htos)

# Normalize HTO data using centered log-ratio (CLR) transformation
hash_seurat <- NormalizeData(hash_seurat, assay = "HTO", normalization.method = "CLR")

# Save object
saveRDS(hash_seurat, "results/k562_sudhl4_new/hash_seurat.rds")


# Demultiplex cells based on HTO enrichment
hash_seurat_demux <- HTODemux(hash_seurat, assay = "HTO", positive.quantile = 0.999)


# Visualize demultiplexing results

# Global classification results
table(hash_seurat_demux$HTO_classification.global)

# Group cells based on the max HTO signal
hash_seurat_demux$HTO_maxID <- factor(hash_seurat_demux$HTO_maxID, levels = rev(rownames(hash_seurat_demux[["HTO"]])))
Idents(hash_seurat_demux) <- "HTO_maxID"

RidgePlot(hash_seurat_demux, assay = "HTO", features = rownames(hash_seurat_demux[["HTO"]]), ncol = 7) &
  scale_fill_manual(values = getPalette3(14)) &
  ylab("") &
  theme(plot.title = element_text(face = "plain", hjust = 0.5))
ggsave("results/k562_sudhl4_new/qc/hto/ridgeplot.pdf", height = 10, width = 30)

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
p1 <- FeatureScatter(subset(hash_seurat_demux, HTO_classification.global == "Singlet"),
                     feature1 = "hto_K562-sgCtrl.1", feature2 = "hto_K562-sgNCR3LG1.1",
                     cols = getPalette3(14))

p2 <- FeatureScatter(subset(hash_seurat_demux, HTO_classification.global == "Singlet"),
                     feature1 = "hto_K562-sgCtrl.1", feature2 = "hto_NK-only",
                     cols = getPalette3(14))

p1 + p2
ggsave("results/k562_sudhl4_new/qc/hto/scatterplots.pdf", height = 5, width = 12)

# Compare number of UMIs for singlets, doublets, and negative cells
Idents(hash_seurat_demux) <- "HTO_classification.global"
VlnPlot(hash_seurat_demux, features = "nCount_RNA", pt.size = 0.0001, log = TRUE) +
  scale_fill_manual(values = getPalette3(3)) &
  xlab("") &
  theme(plot.title = element_text(face = "plain", hjust = 0.5))
ggsave("results/k562_sudhl4_new/qc/hto/ncount_rna_violinplot.pdf", height = 4, width = 4)


# Generate a two dimensional tSNE embedding for HTOs

# Remove negative cells from the object
hash_seurat_subset <- subset(hash_seurat_demux, idents = "Negative", invert = TRUE)

dir.create("results/k562_sudhl4_new/hto_singlets/", recursive = T, showWarnings = F)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(hash_seurat_subset) <- "HTO"
hash_seurat_subset <- ScaleData(hash_seurat_subset, features = rownames(hash_seurat_subset), 
                                verbose = FALSE)
hash_seurat_subset <- RunPCA(hash_seurat_subset, features = rownames(hash_seurat_subset), approx = FALSE)
hash_seurat_subset <- RunTSNE(hash_seurat_subset, dims = 1:14, perplexity = 100)
DimPlot(hash_seurat_subset) +
  scale_color_manual(values = brewer.pal(9, "RdGy")[c(2,9)])
ggsave("results/k562_sudhl4_new/hto_singlets/tsne.pdf", height = 5, width = 6)


# Create an HTO heatmap
hash_seurat_demux$HTO_classification <- factor(hash_seurat_demux$HTO_classification, levels = c(rownames(hash_seurat_demux[["HTO"]]), "NA"))
HTOHeatmap(hash_seurat_demux, assay = "HTO", raster = F, singlet.names = hash_seurat_demux$HTO_classification)
ggsave("results/k562_sudhl4_new/hto_singlets/heatmap.pdf", height = 4, width = 8)


# Cluster and visualize cells
dir.create("results/k562_sudhl4_new/hto_singlets/", recursive = T, showWarnings = F)

# Extract singlets and remove cancer cell cluster
hash_seurat_singlet <- subset(hash_seurat_demux, idents = "Singlet")
hash_seurat_singlet <- subset(hash_seurat_singlet, cluster != "6") # remove target cell cluster
hash_seurat_singlet$HTO_classification <- factor(hash_seurat_singlet$HTO_classification, levels = rownames(hash_seurat_demux[["HTO"]]))

# Select the top 1000 most variable features
hash_seurat_singlet <- FindVariableFeatures(hash_seurat_singlet, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_singlet <- CellCycleScoring(hash_seurat_singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_singlet <- ScaleData(hash_seurat_singlet, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_singlet))

# Run PCA
hash_seurat_singlet <- RunPCA(hash_seurat_singlet, features = VariableFeatures(hash_seurat_singlet))

# Select the top 10 PCs for clustering and UMAP
hash_seurat_singlet <- FindNeighbors(hash_seurat_singlet, reduction = "pca", dims = 1:10)
hash_seurat_singlet <- FindClusters(hash_seurat_singlet, resolution = 0.4, verbose = FALSE)
hash_seurat_singlet <- RunUMAP(hash_seurat_singlet, reduction = "pca", dims = 1:10)

hash_seurat_singlet$cluster <- factor(hash_seurat_singlet$seurat_clusters,
                                      levels = c(0:4),
                                      labels = c("CD56bright (0)",
                                                 "Activated (1)",
                                                 "Proliferating (2)",
                                                 "Adaptive-like (3)",
                                                 "CX3CR1+ (4)"))

# group biological replicates (cell lines)
hash_seurat_singlet$HTO_classification_grouped <- gsub("sg|\\..*|\\-.*", "", hash_seurat_singlet$HTO_classification)
hash_seurat_singlet$HTO_classification_grouped <- factor(hash_seurat_singlet$HTO_classification_grouped,
                                                         levels = c("NK", "K562", "SUDHL4"),
                                                         labels = c("NK cells only", "NK cells + K562", "NK cells + SUDHL4"))

## Visualize
DimPlot(hash_seurat_singlet, group.by = "cluster", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(8))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/k562_sudhl4_new/hto_singlets/umap.png", width = 7, height = 4)

DimPlot(hash_seurat_singlet, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/k562_sudhl4_new/hto_singlets/umap_cycle.png", width = 6, height = 4)

# Save object
saveRDS(hash_seurat_singlet, "results/k562_sudhl4_new/hash_seurat_singlet.rds")

