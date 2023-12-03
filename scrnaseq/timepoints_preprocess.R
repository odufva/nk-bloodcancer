
# NK cell co-culture time course experiment, hashing scRNA-seq analysis

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

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))


## Load functions
me = "user"
source("scrnaseq/fun_helper.R")
source("scrnaseq/fun_getGenes.R")

## Get data
files <- list.files(recursive = T, include.dirs = T) %>% 
  grep(pattern = "Timepoints", value = T) %>% 
  grep(pattern = "HTO", value = T, invert = T)

# get files for each sample/lane
files1 <- files %>% grep(pattern = "plate1", value = T)
files2 <- files %>% grep(pattern = "plate2", value = T)

# Create Seurat objects of each sample/lane and get sample/lane names from filenames
extractBatchName <- function(filenames){gsub("data/|_matrix.mtx.gz", "", filenames[grepl("mtx", filenames)])}

scrnaseq_files1 <- ReadMtx(files1[grepl("mtx", files1)], files1[grepl("barcodes", files1)], files1[grepl("features", files1)]) %>% CreateSeuratObject(project = extractBatchName(files1), min.cells = 3, min.features = 200)
scrnaseq_files2 <- ReadMtx(files2[grepl("mtx", files2)], files2[grepl("barcodes", files2)], files2[grepl("features", files2)]) %>% CreateSeuratObject(project = extractBatchName(files2), min.cells = 3, min.features = 200)

hash_seurat    <- merge(scrnaseq_files1, scrnaseq_files2, add.cell.ids = extractBatchName(c(files1, files2)))


## Basic QC
dir.create("./results/timepoints/qc/", recursive = T, showWarnings = F)
dir.create("./results/timepoints/qc/before_qc1/", recursive = T, showWarnings = F)
dir.create("./results/timepoints/qc/after_qc1/", recursive = T, showWarnings = F)
dir.create("./results/timepoints/qc/hto/", recursive = T, showWarnings = F)

clonality_genes <- getClonalityGenes(hash_seurat)
unwanted_genes  <- getUnwantedGenes(hash_seurat)
hash_seurat     <- PercentageFeatureSet(hash_seurat, pattern = "^MT-", col.name = "percent.mt")
hash_seurat     <- PercentageFeatureSet(hash_seurat, pattern = "^RP", col.name = "percent.ribo")
hash_seurat     <- PercentageFeatureSet(hash_seurat, features = cycle.genes, col.name = "percent.cycle")
hash_seurat@meta.data$barcode <- colnames(hash_seurat)


## Plots and cluster DEG before QC
hash_seurat %>% plotQC(folder = "results/timepoints/qc/before_qc1/", min_mito = 0, max_mito = 15, min_ribo = 5, max_ribo = 50, min_features = 300, max_features = 10e3,
                       min_counts = 700, max_counts = 50e4)

hash_seurat_beforeQC <- hash_seurat %>% preprocessSeurat(cells.to.use = colnames(hash_seurat))

## Get optimal clustering
hash_seurat_beforeQC <- hash_seurat_beforeQC %>% getClustering()
hash_seurat_beforeQC %>% plotClustering()
ggsave("results/timepoints/qc/before_qc1/scatter_clustering.png", width = 5, height = 4)


Idents(hash_seurat_beforeQC) <- hash_seurat_beforeQC$RNA_snn_res.0.4
hash_seurat_beforeQC$cluster <- Idents(hash_seurat_beforeQC)

## Visualize
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))

DimPlot(hash_seurat_beforeQC, group.by = "RNA_snn_res.0.4", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(20)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/timepoints/qc/before_qc1/umap.png", width = 8, height = 6)

# Plot cell type marker genes
FeaturePlot(hash_seurat_beforeQC, features = c("GZMA", "CD79A", "CD3E", "JCHAIN", "MPO", "HLA-DRA", "CD14"), label = T, repel = T, cols = c("lightgrey", "tomato"), min.cutoff = 0) & theme_bw(base_size = 12)  & labs("UMAP 1", "UMAP 2") & theme(legend.position = "right")
ggsave("results/timepoints/qc/before_qc1/umap_genes.png", width = 20, height = 15)

# Plot experiments
DimPlot(hash_seurat_beforeQC, group.by = "RNA_snn_res.0.4", split.by = "orig.ident", label = T, repel = T, ncol = 3) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(77)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/timepoints/qc/before_qc1/umap_experiments.png", width = 24, height = 12)

#all_deg <- FindAllMarkers(hash_seurat_beforeQC, test.use = "t", only.pos = T, return.thresh = 0.05)
#fwrite(all_deg, "results/timepoints/qc/before_qc1/all_deg.txt", sep = "\t", quote = F, row.names = F)

## --------------------------

## QC (no doublet removal, will be done using HTOs)

hash_seurat <- hash_seurat %>% getQC(min_mito = 0, max_mito = 15, min_ribo = 5, max_ribo = 50, min_features = 300, max_features = 10e3,
                                     min_counts = 700, max_counts = 50e4)

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
ggsave("results/timepoints/qc/after_qc1/scatter_clustering.png", width = 5, height = 4)

Idents(hash_seurat) <- hash_seurat$RNA_snn_res.0.5
hash_seurat$cluster <- Idents(hash_seurat)

## Visualize
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))

DimPlot(hash_seurat, group.by = "RNA_snn_res.0.5", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(44)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/timepoints/qc/after_qc1/umap.png", width = 5, height = 4)

# Plot experiments
DimPlot(hash_seurat_beforeQC, group.by = "RNA_snn_res.0.5", split.by = "orig.ident", label = T, repel = T, ncol = 3) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(60)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/timepoints/qc/after_qc1/umap_experiments.png", width = 24, height = 12)


## Get cell cycle scores and visualize
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
hash_seurat_cycle <- CellCycleScoring(hash_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_cycle <- ScaleData(hash_seurat_cycle, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_cycle))

DimPlot(hash_seurat_cycle, group.by = "Phase", repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/qc/after_qc1/umap_cycle.png", width = 5, height = 4)

DimPlot(hash_seurat, group.by = "singler_hpca_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(28)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/timepoints/qc/after_qc1/umap_singler_hpca.png", width = 14, height = 7)

DimPlot(hash_seurat, group.by = "singler_blueprint_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(23)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/timepoints/qc/after_qc1/umap_singler_blueprint.png", width = 10, height = 6)

FeaturePlot(hash_seurat, features = "hybrid_doublet_score", label = T, repel = T, cols = c("lightgrey", "tomato"), min.cutoff = 0) + theme_bw(base_size = 12)  + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/timepoints/qc/after_qc1/umap_dbl.png", width = 5, height = 4)

hash_seurat %>% plotQC(folder = "results/timepoints/qc/after_qc1/")

#all_deg <- FindAllMarkers(hash_seurat, test.use = "t", only.pos = T, return.thresh = 0.05)
#fwrite(all_deg, "results/timepoints/qc/after_qc1/all_deg.txt", sep = "\t", quote = F, row.names = F)

## --------------------------

## Analysis using HTOs

## Get data
files_hto <- list.files(recursive = T, include.dirs = T) %>% 
  grep(pattern = "Timepoints", value = T) %>% 
  grep(pattern = "HTO", value = T)

# get files for each sample/lane
files1_hto <- files_hto %>% grep(pattern = "plate1", value = T)
files2_hto <- files_hto %>% grep(pattern = "plate2", value = T)

# Create Seurat objects of each sample/lane and get sample/lane names from filenames
extractBatchName <- function(filenames){gsub("data/|_matrix.mtx.gz", "", filenames[grepl("mtx", filenames)])}


# check that rowSums look like what's expected for the used HTOs
htos_plate1 <- ReadMtx(files1_hto[grepl("mtx", files1_hto)], files1_hto[grepl("barcodes", files1_hto)], files1_hto[grepl("features", files1_hto)], feature.column = 1)
rowSums(htos_plate1)
htos_plate1 <- htos_plate1[1:13,-1]

htos_plate2 <- ReadMtx(files2_hto[grepl("mtx", files2_hto)], files2_hto[grepl("barcodes", files2_hto)], files2_hto[grepl("features", files2_hto)], feature.column = 1)
rowSums(htos_plate2)
htos_plate2 <- htos_plate2[1:13,-1]


## Load HTO ids
hto_id <- as.data.frame(
  rbind(fread("data/Timepoints_plate1_HTO_sample_list.txt"),
        fread("data/Timepoints_plate2_HTO_sample_list.txt")))

## ---------------------------------

## Run demultiplxing for each experiment separately


# function for demultiplexing an experiment
runDemultiplex <- function(hto, prefix, experiment){
  
  htos <- hto
  
  # Replace HTO matrix rownames with ids
  hto_id_subset <- hto_id[hto_id$batch_10x_lane == gsub("_$", "", prefix),]
  rownames(htos) <- hto_id_subset$sample[match(gsub("^.*\\-", "", rownames(htos)), hto_id_subset$sequence)]
  
  # Modify column names to match expression matrix
  colnames(htos) <- paste0(prefix, colnames(htos), "-1")
  
  # Subset RNA and HTO counts by joint cell barcodes
  joint.bcs <- intersect(colnames(hash_seurat), colnames(htos))
  hash_seurat_hto <- subset(hash_seurat, barcode %in% joint.bcs)
  htos <- as.matrix(htos[, hash_seurat_hto$barcode])
  
  # Add HTO data as a new assay independent from RNA
  hash_seurat_hto[["HTO"]] <- CreateAssayObject(counts = htos)
  
  # Normalize HTO data using centered log-ratio (CLR) transformation
  hash_seurat_hto <- NormalizeData(hash_seurat_hto, assay = "HTO", normalization.method = "CLR")
  
  # Demultiplex cells based on HTO enrichment
  hash_seurat_demux <- HTODemux(hash_seurat_hto, assay = "HTO", positive.quantile = 0.99)
  
  
  # Visualize demultiplexing results
  
  # Global classification results
  table(hash_seurat_demux$HTO_classification.global)
  
  # Group cells based on the max HTO signal
  hash_seurat_demux$HTO_maxID <- factor(hash_seurat_demux$HTO_maxID, levels = rev(rownames(hash_seurat_demux[["HTO"]])))
  Idents(hash_seurat_demux) <- "HTO_maxID"
  
  dir.create(paste0("results/timepoints/qc/hto/", experiment))
  
  RidgePlot(hash_seurat_demux, assay = "HTO", features = rownames(hash_seurat_demux[["HTO"]]), ncol = 7) &
    scale_fill_manual(values = getPalette3(14)) &
    ylab("") &
    theme(plot.title = element_text(face = "plain", hjust = 0.5))
  ggsave(paste0("results/timepoints/qc/hto/", experiment, "/ridgeplot.pdf"), height = 10, width = 30)
  
  # Compare number of UMIs for singlets, doublets, and negative cells
  Idents(hash_seurat_demux) <- "HTO_classification.global"
  VlnPlot(hash_seurat_demux, features = "nCount_RNA", pt.size = 0.1, log = TRUE) +
    scale_fill_manual(values = getPalette3(3)) &
    xlab("") &
    theme(plot.title = element_text(face = "plain", hjust = 0.5))
  ggsave(paste0("results/timepoints/qc/hto/", experiment, "/ncount_rna_violinplot.pdf"), height = 4, width = 4)
  
  
  # Generate a two dimensional tSNE embedding for HTOs
  
  # Remove negative cells from the object
  hash_seurat_subset <- subset(hash_seurat_demux, idents = "Negative", invert = TRUE)
  
  dir.create(paste0("results/timepoints/hto_singlets/", experiment), recursive = T, showWarnings = F)
  
  # Calculate a tSNE embedding of the HTO data
  DefaultAssay(hash_seurat_subset) <- "HTO"
  hash_seurat_subset <- ScaleData(hash_seurat_subset, features = rownames(hash_seurat_subset), 
                                  verbose = FALSE)
  hash_seurat_subset <- RunPCA(hash_seurat_subset, features = rownames(hash_seurat_subset), approx = FALSE)
  hash_seurat_subset <- RunTSNE(hash_seurat_subset, dims = 1:14, perplexity = 100)
  DimPlot(hash_seurat_subset) +
    scale_color_manual(values = brewer.pal(9, "RdGy")[c(9,2)])
  ggsave(paste0("results/timepoints/hto_singlets/", experiment, "/tsne.pdf"), height = 5, width = 6)
  
  
  # Create an HTO heatmap
  hash_seurat_demux$HTO_classification <- factor(hash_seurat_demux$HTO_classification, levels = c(rownames(hash_seurat_demux[["HTO"]]), "NA"))
  HTOHeatmap(hash_seurat_demux, assay = "HTO", raster = F, singlet.names = hash_seurat_demux$HTO_classification)
  ggsave(paste0("results/timepoints/hto_singlets/", experiment, "/heatmap.pdf"), height = 4, width = 8)
  
  return(hash_seurat_demux)
  
}

hash_seurat_demux_plate1 <- runDemultiplex(htos_plate1, "Timepoints_plate1_", "timepoints_plate1")
hash_seurat_demux_plate2 <- runDemultiplex(htos_plate2, "Timepoints_plate2_", "timepoints_plate2")

# extract sample identities and merge with original object
labels <- data.frame(barcode = c(names(hash_seurat_demux_plate1$hash.ID), names(hash_seurat_demux_plate2$hash.ID)),
                     hash.ID =  c(as.character(hash_seurat_demux_plate1$hash.ID), as.character(hash_seurat_demux_plate2$hash.ID)))

df <- hash_seurat@meta.data %>% left_join(labels, by = "barcode")
rownames(df) <- df$barcode
hash_seurat@meta.data <- df


# Cluster and visualize cells

# Extract singlets and remove cancer cell cluster
"%ni%" <- Negate("%in%")
Idents(hash_seurat) <- "hash.ID"
hash_seurat_singlet <- subset(hash_seurat, hash.ID %ni% c("Doublet", "Negative"))
#hash_seurat_singlet <- subset(hash_seurat_singlet, cluster != "10")
#hash_seurat_singlet$HTO_classification <- factor(hash_seurat_singlet$HTO_classification, levels = rownames(hash_seurat_demux[["HTO"]]))

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
ggsave("results/timepoints/hto_singlets/umap.png", width = 7, height = 5)

DimPlot(hash_seurat_singlet, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/umap_cycle.png", width = 6, height = 4)

# Save object
saveRDS(hash_seurat_singlet, "results/timepoints/hash_seurat_singlet.rds")

## ----------------------------

# Plots + DE analysis

# Projecting singlet identities on UMAP visualization
DimPlot(hash_seurat_singlet, split.by = "hash.ID", ncol = 9, group.by = "seurat_clusters") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(41))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/hto_class_umap.pdf", height = 16, width = 16)

# Larger plot
DimPlot(hash_seurat_singlet, split.by = "hash.ID", ncol = 9, group.by = "seurat_clusters") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(41))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/hto_class_umap_large.pdf", height = 32, width = 32)

# Larger plot with singleR colors
DimPlot(hash_seurat_singlet, split.by = "hash.ID", ncol = 9, group.by = "singler_blueprint_pred") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(41))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/hto_class_umap_large_singler.pdf", height = 32, width = 32)

DimPlot(hash_seurat_singlet, group.by = "hash.ID", label = T, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(101))) +
  labs("UMAP 1", "UMAP 2") +
  NoLegend()
ggsave("results/timepoints/hto_singlets/hto_class_umap_colored.pdf", height = 10, width = 16)

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
ggsave("results/timepoints/hto_singlets/umap_experiment.png", height = 16, width = 24)

# Experiment with sample labels
DimPlot(hash_seurat_singlet, split.by = "orig.ident", ncol = 3, group.by = "hash.ID", label = T, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(101))) +
  labs("UMAP 1", "UMAP 2") +
  NoLegend()
ggsave("results/timepoints/hto_singlets/umap_experiment_samples.png", height = 12, width = 16)

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
ggsave("results/timepoints/hto_singlets/singler_blueprint_pred.pdf", height = 7, width = 12)

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
ggsave("results/timepoints/hto_singlets/singler_hpca_pred.pdf", height = 7, width = 12)

## ---------------------------------------

# Save processed object
hash_seurat_singlet <- readRDS("results/timepoints/timepoints_nk_seurat.rds")


# HTO singlet cell numbers
Idents(hash_seurat_singlet) <- "hash.ID"

cell_numbers <- data.frame(table(Idents(hash_seurat_singlet))) #%>%
#mutate(Var1 = factor(Var1, levels = rownames(hash_seurat_singlet[["HTO"]])))

ggplot(cell_numbers, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_col() +
  scale_fill_manual(values = getPalette3(101)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,2), "cm")) +
  xlab("") +
  ylab("Number of singlets") +
  guides(fill = F)

ggsave("results/timepoints/hto_singlets/hto_cell_numbers.pdf", height = 5, width = 25)

# HTO singlet cell numbers stacked barplot with NK and target cells colored

cell_numbers <- data.frame(sample = hash_seurat_singlet$hash.ID, cluster = hash_seurat_singlet$seurat_clusters) %>% 
  mutate(celltype = ifelse(cluster %in% c(0, 28), "NK cell", "Target cell")) %>% 
  select(-cluster)

cell_numbers <- data.frame(table(cell_numbers)) %>%
  mutate(treatment = ifelse(grepl("NK2", sample), "NK2",
                            ifelse(grepl("NK3", sample), "NK3",
                                   ifelse(grepl("NK4", sample), "NK4",
                                          ifelse(grepl("NK5", sample), "NK5",
                                                 ifelse(grepl("NK6", sample), "NK6", "Targets only")))))) %>% 
  mutate(cell_line = gsub("\\-.*", "", sample)) %>% 
  arrange(treatment, cell_line)

cell_numbers$cell_line <- factor(cell_numbers$cell_line, levels = unique(cell_numbers$cell_line))


ggplot(cell_numbers, aes(x = cell_line, y = Freq, fill = celltype)) +
  geom_col() +
  scale_fill_manual(values = c("darkred", "grey70")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,2), "cm")) +
  xlab("") +
  ylab("Number of singlets") +
  facet_wrap(. ~ treatment, ncol = 1)

ggsave("results/timepoints/hto_singlets/hto_cell_numbers_celltype.pdf", height = 7, width = 10)



