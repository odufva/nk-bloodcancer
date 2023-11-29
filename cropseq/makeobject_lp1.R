
# Make LP1 CROP-seq object

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(reshape)
library(scds)


theme_set(theme_classic(base_size = 12))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))


## Load functions
me = "user"
source("scrnaseq/fun_helper.R")
source("scrnaseq/fun_getGenes.R")

## Get data
files <- list.files(recursive = T, include.dirs = T) %>% 
  grep(pattern = "CROPseq_LP1", value = T) %>% 
  grep(pattern = "sgRNA", value = T, invert = T)

# get files for each sample/lane
files1 <- files %>% grep(pattern = "noNK", value = T)
files2 <- files %>% grep(pattern = "NK1_1_16", value = T)
files3 <- files %>% grep(pattern = "NK1_1_4", value = T)

# Create Seurat objects of each sample/lane and get sample/lane names from filenames
extractBatchName <- function(filenames){gsub("data/|_matrix.mtx.gz", "", filenames[grepl("mtx", filenames)])}

scrnaseq_files1 <- ReadMtx(files1[grepl("mtx", files1)], files1[grepl("barcodes", files1)], files1[grepl("features", files1)]) %>% CreateSeuratObject(project = extractBatchName(files1), min.cells = 3, min.features = 200)
scrnaseq_files2 <- ReadMtx(files2[grepl("mtx", files2)], files2[grepl("barcodes", files2)], files2[grepl("features", files2)]) %>% CreateSeuratObject(project = extractBatchName(files2), min.cells = 3, min.features = 200)
scrnaseq_files3 <- ReadMtx(files3[grepl("mtx", files3)], files3[grepl("barcodes", files3)], files3[grepl("features", files3)]) %>% CreateSeuratObject(project = extractBatchName(files3), min.cells = 3, min.features = 200)

crop_seurat    <- merge(scrnaseq_files1, c(scrnaseq_files2, scrnaseq_files3), add.cell.ids = extractBatchName(c(files1, files2, files3)))


## Basic QC
dir.create("results", showWarnings = T)
dir.create("results/lp1/", showWarnings = T)
dir.create("results/lp1/qc/", showWarnings = F)
dir.create("results/lp1/qc/before_qc1/", showWarnings = F)
dir.create("results/lp1/qc/after_qc1/", showWarnings = F)

clonality_genes <- getClonalityGenes(crop_seurat)
unwanted_genes  <- getUnwantedGenes(crop_seurat)
crop_seurat     <- PercentageFeatureSet(crop_seurat, pattern = "^MT-", col.name = "percent.mt")
crop_seurat     <- PercentageFeatureSet(crop_seurat, pattern = "^RP", col.name = "percent.ribo")
crop_seurat     <- PercentageFeatureSet(crop_seurat, features = cycle.genes, col.name = "percent.cycle")
crop_seurat@meta.data$barcode <- colnames(crop_seurat)

crop_seurat %>% plotQC(folder = "results/lp1/qc/before_qc1/")
crop_seurat <- crop_seurat %>% getQC()



## Get SingleR predictions
crop_seurat <- crop_seurat %>% getSingler()
relevant_hpca_clusters <- crop_seurat@meta.data %>% group_by(singler_hpca_pred) %>% summarise(n=n()) %>% filter(n>20) %>% pull(singler_hpca_pred)
relevant_blue_clusters <- crop_seurat@meta.data %>% group_by(singler_blueprint_pred) %>% summarise(n=n()) %>% filter(n>20) %>% pull(singler_blueprint_pred)

crop_seurat$singler_hpca_pred      <- ifelse(crop_seurat$singler_hpca_pred %in% relevant_hpca_clusters, crop_seurat$singler_hpca_pred, "rare")
crop_seurat$singler_blueprint_pred <- ifelse(crop_seurat$singler_blueprint_pred %in% relevant_blue_clusters, crop_seurat$singler_blueprint_pred, "rare")


## Get doublets
crop_seurat <- crop_seurat %>% getDoublets()
crop_seurat$hybrid_doublet_score

crop_seurat@meta.data %>% 
  ggplot(aes(orig.ident, hybrid_doublet_score/2)) + geom_violin(draw_quantiles = 0.5)

crop_seurat@meta.data %>% 
  ggplot(aes(hybrid_doublet_score/2)) + geom_histogram()

table(crop_seurat$hybrid_doublet_score/2 > 0.90)

## Get Seurat object
crop_seurat <- crop_seurat %>% preprocessSeurat(cells.to.use = colnames(crop_seurat))

## Get optimal clustering
crop_seurat <- crop_seurat %>% getClustering()

## Decide on clustering
crop_seurat %>% plotClustering()
ggsave("results/lp1/qc/after_qc1/scatter_clustering.png", width = 5, height = 4)

Idents(crop_seurat) <- crop_seurat$RNA_snn_res.0.3
crop_seurat$cluster <- Idents(crop_seurat)

## Visualize
DimPlot(crop_seurat, group.by = "RNA_snn_res.0.3", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(10)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/lp1/qc/after_qc1/umap.png", width = 5, height = 4)

DimPlot(crop_seurat, group.by = "orig.ident", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(6)[c(1,4,5)]) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "top")
ggsave("results/lp1/qc/after_qc1/umap_ident.png", width = 5, height = 4)

DimPlot(crop_seurat, group.by = "singler_hpca_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette2(10)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc1/umap_singler_hpca.png", width = 7, height = 4)

DimPlot(crop_seurat, group.by = "singler_blueprint_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette2(11)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc1/umap_blue_blueprint.png", width = 7, height = 4)

FeaturePlot(crop_seurat, features = "hybrid_doublet_score", label = T, repel = T, cols = c("lightgrey", "tomato"), min.cutoff = 0) + theme_bw(base_size = 12)  + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc1/umap_dbl.png", width = 5, height = 4)

crop_seurat %>% plotQC(folder = "results/qc/after_qc1/")
dev.off()

all_deg <- FindAllMarkers(crop_seurat, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/lp1/qc/after_qc1/all_deg.txt", sep = "\t", quote = F, row.names = F)

saveRDS(crop_seurat, "results/lp1/qc/after_qc1/lp1_cropseq_seurat.rds")
getwd()
setwd("/Users/sgandalf85/Dropbox (HRUH)/NK_resistance Heme CollabPaper/Analysis/CROP_seq_LP1/R/")

######## Add the sgRNA information
dir.create("results/lp1/qc/sgrna_qc/", showWarnings = F)

cropseq_no_nk_data <- fread("data/CROPseq_LP1_noNK_sgRNA.txt") %>% dplyr::rename(sequence = barcode, barcode = cell) %>% mutate(barcode = paste0("CROPseq_LP1_noNK_", barcode))
cropseq_nk_1_16_data    <- fread("data/CROPseq_LP1_NK1_1_16_sgRNA.txt") %>% dplyr::rename(sequence = barcode, barcode = cell) %>% mutate(barcode = paste0("CROPseq_LP1_NK1_1_16_", barcode))
cropseq_nk_1_4_data    <- fread("data/CROPseq_LP1_NK1_1_4_sgRNA.txt") %>% dplyr::rename(sequence = barcode, barcode = cell) %>% mutate(barcode = paste0("CROPseq_LP1_NK1_1_4_", barcode))
cropseq_library    <- fread("data/CROPseq_LP1_sgRNA_library.txt")

cropseq_data       <- rbind(cropseq_no_nk_data, cropseq_nk_1_16_data, cropseq_nk_1_4_data) %>% left_join(cropseq_library, by = "sequence")
cropseq_data       <- cropseq_data %>% filter(!is.na(sgrna))

## Look for barcodes that are found in the seurat object
table(crop_seurat$barcode %in% cropseq_data$barcode)
table(cropseq_data$barcode %in% crop_seurat$barcode)
cropseq_data <- cropseq_data %>% filter(barcode %in% crop_seurat$barcode)
fwrite(cropseq_data, "results/lp1/qc/sgrna_qc/cropseq_data.txt", sep = "\t", quote = F, row.names = F)


## Look what is up with duplicated barcodes (ie same cells with at least 2 different sgRNAs)
duplicated_barcodes <- cropseq_data[duplicated(cropseq_data$barcode), ]$barcode %>% unique()
cropseq_data_dup    <- cropseq_data %>% filter(barcode %in% duplicated_barcodes)

plotdata <- cropseq_data_dup %>% group_by(barcode, gene) %>% summarise(total = sum(umi_count)) %>% 
  mutate(prop=total/sum(total)) %>% group_by(barcode) %>% top_n(n = 1, wt = prop)

ggplot(plotdata, aes(prop,log10(total))) + geom_point(alpha=0.5, size = 0.4) +
  labs(x = "top sgRNA \nUMI fraction", y = "top sgRNA \nUMI counts (log10)")
ggsave("results/lp1/qc/sgrna_qc/scatter_top_fraction_counts.png", width = 5, height = 4)

ggplot(plotdata, aes(prop)) + geom_histogram(bins = 100) +
  labs(x = "top sgRNA \nUMI fraction", y = "count")
ggsave("results/lp1/qc/sgrna_qc/histogram_top_fraction.png", width = 5, height = 4)

rank_df <- cropseq_data_dup %>% group_by(barcode, gene) %>% summarise(total = sum(umi_count)) %>% mutate(prop=total/sum(total)) %>% group_by(barcode) %>% top_n(n = 2, wt = prop) %>% group_by(barcode) %>% mutate(rank = dense_rank(prop)) # %>% filter(total >= 3)
rank_1  <- rank_df %>% filter(rank == 1)
rank_2  <- rank_df %>% filter(rank == 2)
rank_df <- dplyr::select(rank_1, barcode, total, prop) %>% left_join(dplyr::select(rank_2, barcode, total, prop), by = "barcode") 
rank_df <- rank_df[!duplicated(rank_df), ]

label_df <- rank_df %>% 
  mutate(status = ifelse(as.numeric(prop.y) > 0.5 & log10(total.y) > 1, "single", "artefact")) %>% 
  mutate(status = ifelse(as.numeric(prop.x) > 0.1 & log10(total.x) > 1, "doublet", status))

table(label_df$status)

rank_df %>%
  left_join(label_df) %>%
  ggplot(aes(prop.y, log10(total.x), color = status)) +
  geom_point(size = 0.4, alpha = 0.5) +
  labs(x = "top 1 sgRNA \nUMI fraction", y = "top 2 sgRNA \nUMI counts") +
  add_guide +
  scale_color_manual(values = c("lightgrey", "salmon", "dodgerblue"))
ggsave("results/lp1/qc/sgrna_qc/scatter_top1_fraction_top2_counts.png", width = 5, height = 4)

rank_df %>%
  left_join(label_df) %>%
  ggplot(aes(prop.y, log10(total.y), color = status)) +
  geom_point(size = 0.4, alpha = 0.5) +
  labs(x = "top 1 sgRNA \nUMI fraction", y = "top 1 sgRNA \nUMI counts") +
  add_guide +
  scale_color_manual(values = c("lightgrey", "salmon", "dodgerblue"))
ggsave("results/lp1/qc/sgrna_qc/scatter_top1_fraction_top1_counts.png", width = 5, height = 4)

rank_df %>%
  left_join(label_df) %>%
  ggplot(aes(prop.x, log10(total.x), color = status)) +
  geom_point(size = 0.4, alpha = 0.5) +
  labs(x = "top 2 sgRNA \nUMI fraction", y = "top 2 sgRNA \nUMI counts") +
  add_guide +
  scale_color_manual(values = c("lightgrey", "salmon", "dodgerblue"))
ggsave("results/lp1/qc/sgrna_qc/scatter_top2_fraction_top2_counts.png", width = 5, height = 4)

rank_df %>%
  left_join(label_df) %>%
  ggplot(aes(prop.x, log10(total.y), color = status)) +
  geom_point(size = 0.4, alpha = 0.5) +
  labs(x = "top 2 sgRNA \nUMI fraction", y = "top 1 sgRNA \nUMI counts") +
  add_guide +
  scale_color_manual(values = c("lightgrey", "salmon", "dodgerblue"))
ggsave("results/lp1/qc/sgrna_qc/scatter_top2_fraction_top1_counts.png", width = 5, height = 4)

fwrite(cropseq_data_dup, "results/lp1/qc/sgrna_qc/cropseq_data_duplicated_barcodes.txt", sep = "\t", quote = F, row.names = F)
fwrite(rank_df, "results/lp1/qc/sgrna_qc/rank_df.txt", sep = "\t", quote = F, row.names = F)
fwrite(label_df, "results/lp1/qc/sgrna_qc/label_df.txt", sep = "\t", quote = F, row.names = F)


## Add labels to cropseq duplicates
cropseq_data_dup <- cropseq_data_dup %>% left_join(dplyr::select(label_df, barcode, status))

singlet_lookup_df <- cropseq_data_dup %>%
  filter(status == "single") %>%
  group_by(barcode) %>% top_n(n = 1, wt = umi_count) %>%
  dplyr::select(barcode,gene, sgrna_name, status) %>%
  group_by(barcode) %>%
  mutate(n = n()) %>%
  filter(n == 1) %>%
  dplyr::select(-n) %>%
  mutate(gene2 = NA, sgrna_name2 = NA) %>%
  dplyr::select(barcode, status, gene, sgrna_name, gene2, sgrna_name2)

doublet_lookup_df <- cropseq_data_dup %>%
  filter(status == "doublet") %>%
  group_by(barcode) %>% top_n(n = 2, wt = umi_count) %>%
  dplyr::select(barcode, gene, sgrna_name, status, umi_count)

doublet_lookup_df1 <- doublet_lookup_df %>%
  group_by(barcode) %>% top_n(n = 1, wt = umi_count) %>%
  dplyr::select(barcode,gene, sgrna_name, status)

doublet_lookup_df2 <- doublet_lookup_df %>%
  group_by(barcode) %>% top_n(n = 1, wt = dplyr::desc(umi_count)) %>%
  dplyr::select(barcode,gene, sgrna_name, status) %>%
  dplyr::rename(gene2 = gene, sgrna_name2 = sgrna_name)


doublet_lookup_df <- merge(doublet_lookup_df1, doublet_lookup_df2, by = c("barcode", "status"), suffixes = c("1", "2")) %>%
  group_by(barcode) %>%
  mutate(n = n()) %>%
  filter(n == 1) %>%
  dplyr::select(-n)

cropseq_labels <- rbind(singlet_lookup_df, doublet_lookup_df)

cropseq_labels[duplicated(cropseq_labels$barcode), ]$barcode

cropseq_labels <- cropseq_labels[!duplicated(cropseq_labels), ]
fwrite(cropseq_labels, "results/lp1/qc/sgrna_qc/cropseq_labels_in_seurat.txt", sep = "\t", quote = F, row.names = F)


## Add to seurat
df <- crop_seurat@meta.data %>% left_join(cropseq_labels, by = "barcode")
rownames(df) <- df$barcode
crop_seurat@meta.data <- df
crop_seurat$gene <- ifelse(grepl(pattern = "Control", crop_seurat$gene), "Control", crop_seurat$gene)
crop_seurat$gene2 <- ifelse(grepl(pattern = "Control", crop_seurat$gene2), "Control", crop_seurat$gene2)
saveRDS(crop_seurat, "results/lp1/qc/after_qc1/lp1_crop_seurat_qc1.rds")

DimPlot(crop_seurat, split.by = "gene", ncol = 5, label = F, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = getPalette3(25)) +
  labs("UMAP 1", "UMAP 2") +
  theme(legend.position = "none")
ggsave("results/lp1/qc/after_qc1/umap_pert.png", width = 10, height = 8)

## ------------------------------------------------------------

## QC2: Remove clusters 5-n (NK cells/low quality?)

dir.create("results/lp1/qc/after_qc2/", showWarnings = F)
crop_seurat <- subset(crop_seurat, idents = c(0:3,6))
crop_seurat <- crop_seurat %>% preprocessSeurat(cells.to.use = colnames(crop_seurat)) %>% getClustering()

## Decide on clustering
crop_seurat %>% plotClustering()
ggsave("results/lp1/qc/after_qc2/scatter_clustering.png", width = 5, height = 4)

Idents(crop_seurat) <- crop_seurat$RNA_snn_res.0.5
crop_seurat$cluster <- Idents(crop_seurat)

## Visualize
DimPlot(crop_seurat, group.by = "RNA_snn_res.0.5", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(10)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/lp1/qc/after_qc2/umap.png", width = 5, height = 4)

DimPlot(crop_seurat, group.by = "orig.ident", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(6)[c(1,4,5)]) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "top")
ggsave("results/lp1/qc/after_qc2/umap_ident.png", width = 5, height = 4)

DimPlot(crop_seurat, group.by = "singler_hpca_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette2(10)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc2/umap_singler_hpca.png", width = 7, height = 4)

DimPlot(crop_seurat, group.by = "singler_blueprint_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette2(11)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc2/umap_blue_blueprint.png", width = 7, height = 4)

FeaturePlot(crop_seurat, features = "hybrid_doublet_score", label = T, repel = T, cols = c("lightgrey", "tomato"), min.cutoff = 0) + theme_bw(base_size = 12)  + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc2/umap_dbl.png", width = 5, height = 4)

crop_seurat %>% plotQC(folder = "results/qc/after_qc2/")
dev.off()

all_deg <- FindAllMarkers(crop_seurat, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/lp1/qc/after_qc2/all_deg.txt", sep = "\t", quote = F, row.names = F)

saveRDS(crop_seurat, "results/lp1/qc/after_qc2/lp1_cropseq_seurat.rds")

## QC3: Remove small NK cell cluster 
dir.create("results/lp1/qc/after_qc3/", showWarnings = F)
crop_seurat <- subset(crop_seurat, idents = c(0:6,8))
crop_seurat <- crop_seurat %>% preprocessSeurat(cells.to.use = colnames(crop_seurat)) %>% getClustering()

## Decide on clustering
crop_seurat %>% plotClustering()
ggsave("results/lp1/qc/after_qc3/scatter_clustering.png", width = 5, height = 4)

Idents(crop_seurat) <- crop_seurat$RNA_snn_res.0.6
crop_seurat$cluster <- Idents(crop_seurat)

# Remove cell doublets
VlnPlot(crop_seurat, features = "nCount_RNA", split.by = "status")
ggsave("results/lp1/qc/after_qc3/sgrna_doublet_nCount_RNA.png", width = 5, height = 4)

VlnPlot(crop_seurat, features = "bcds_doublet_score_norm", split.by = "status")
ggsave("results/lp1/qc/after_qc3/sgrna_doublet_bcds_score.png", width = 5, height = 4)

crop_seurat$cell_doublet <- "singlet"
crop_seurat$cell_doublet[crop_seurat$status%in%"doublet" & crop_seurat$bcds_doublet_score>0.1] <- "doublet"

VlnPlot(crop_seurat, features = "nFeature_RNA", split.by = "cell_doublet")
ggsave("results/lp1/qc/after_qc3/cell_doublet.png", width = 5, height = 4)

VlnPlot(subset(crop_seurat, cell_doublet == "singlet"), features = "nCount_RNA", split.by = "status")
ggsave("results/lp1/qc/after_qc2/sgrna_doublet_nCount_RNA_filtered.png", width = 5, height = 4)

crop_seurat <- subset(crop_seurat, cell_doublet == "singlet" & status %in% c("single", "doublet"))

# save object including cells with doublet sgRNAs
saveRDS(crop_seurat, "results/lp1/lp1_crop_seurat_withdoublets_2.rds")

# subset only to singlets (which are used in the analysis)
crop_seurat <- subset(crop_seurat, cell_doublet == "singlet" & status == "single")

## Add cell-cycle info
exp.mat     <- read.table(file = "data/genes/cell_cycle/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)
s.genes     <- cc.genes$s.genes
g2m.genes   <- cc.genes$g2m.genes
crop_seurat <- CellCycleScoring(crop_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(crop_seurat) <- crop_seurat$RNA_snn_res.0.6
crop_seurat$cluster <- Idents(crop_seurat)

DimPlot(crop_seurat, group.by = "RNA_snn_res.0.6", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(10)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "none")
ggsave("results/lp1/qc/after_qc3/umap.png", width = 5, height = 4)

DimPlot(crop_seurat, group.by = "orig.ident", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(6)[c(1,4,5)]) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc3/umap_ident.png", width = 7, height = 4)

DimPlot(crop_seurat, group.by = "singler_hpca_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(10)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc3/umap_singler_hpca.png", width = 7, height = 4)

DimPlot(crop_seurat, group.by = "singler_blueprint_pred", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(11)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc3/umap_blue_blueprint.png", width = 7, height = 4)

DimPlot(crop_seurat, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(4)) + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc3/umap_phase.png", width = 5.5, height = 4)

FeaturePlot(crop_seurat, features = "hybrid_doublet_score", label = T, repel = T, cols = c("lightgrey", "tomato"), min.cutoff = 0) + theme_bw(base_size = 12)  + labs("UMAP 1", "UMAP 2") + theme(legend.position = "right")
ggsave("results/lp1/qc/after_qc3/umap_dbl.png", width = 5, height = 4)

all_deg <- FindAllMarkers(crop_seurat, test.use = "t", only.pos = T, return.thresh = 0.05, logfc.threshold = 0.05)
fwrite(all_deg, "results/lp1/qc/after_qc3/all_deg.txt", sep = "\t", quote = F, row.names = F)
all_deg %>% group_by(cluster) %>% top_n(20, dplyr::desc(avg_log2FC)) %>% View

saveRDS(crop_seurat, "results/lp1/lp1_crop_seurat_singlet_2.rds")

DimPlot(crop_seurat, split.by = "gene", ncol = 5, label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = getPalette3(25)) +
  labs("UMAP 1", "UMAP 2") +
  theme(legend.position = "none")
ggsave("results/lp1/qc/after_qc3/umap_pert.png", width = 10, height = 8)



