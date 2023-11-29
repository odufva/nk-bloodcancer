
# CROP-seq mixscape analysis, SUDHL4 cell line as example

# load libraries
library(Seurat)
library(patchwork)
library(data.table)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)

dir.create("results/sudhl4/mixscape/singlet/")

data <- readRDS("results/sudhl4/sudhl4_crop_seurat_singlet.rds")

run_analyze_mixscape <- function(ORIGIDENT, CONDITION){
  
  # select only one condition
  data <- subset(data, orig.ident==ORIGIDENT)
  
  # subset to singlets
  data <- subset(data, status == "single")
  
  # new vector for NT vs perturbed
  data$NT <- ifelse(data$gene == "Control", "NT", "Perturbed")
  data$gene <- gsub("Control", "NT", data$gene)
  data$gene <- data$gene
  
  # theme for plots
  custom_theme <- theme(
    plot.title = element_text(size=16, hjust = 0.5), 
    legend.key.size = unit(0.7, "cm"), 
    legend.text = element_text(size = 14))
  
  # Normalize data, find variable features and scale data
  DefaultAssay(object = data) <- 'RNA'
  data <- NormalizeData(object = data) %>% FindVariableFeatures() %>% ScaleData()
  
  # PCA
  data <- RunPCA(object = data)
  
  # UMAP
  data <- RunUMAP(object = data, dims = 1:40)
  
  # Check if clustering is driven by perturbation status
  DimPlot(
    object = data, 
    group.by = 'NT', 
    pt.size = 0.2, 
    reduction = "umap", 
    split.by = "NT", 
    ncol = 1, 
    cols = c("grey39","goldenrod3")
  ) + 
    ggtitle("Perturbation Status") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_umap_pert_status.pdf"), height = 5, width = 6)
  
  DimPlot(
    object = data, 
    group.by = 'Phase', 
    pt.size = 0.2, 
    reduction = "umap", 
    ncol = 1, 
  ) + 
    ggtitle("Cell cycle phase") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_umap_phase.pdf"), height = 5, width = 6)
  
  DimPlot(
    object = data, 
    group.by = 'gene', 
    pt.size = 0.2, 
    reduction = "umap", 
    ncol = 1, 
  ) + 
    ggtitle("Target gene") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_umap_pert.pdf"), height = 5, width = 6)
  
  # Calculate perturbation signature (PRTB)
  data <- CalcPerturbSig(
    object = data, 
    assay = "RNA", 
    slot = "data", 
    gd.class ="gene", 
    nt.cell.class = "NT", 
    reduction = "pca", 
    ndims = 15, 
    num.neighbors = 20, 
    new.assay.name = "PRTB")
  
  # Prepare PRTB assay for dimensionality reduction: 
  # Normalize data, find variable features and center data
  DefaultAssay(object = data) <- 'PRTB'
  
  # Use variable features from RNA assay.
  VariableFeatures(object = data) <- VariableFeatures(object = data[["RNA"]])
  data <- ScaleData(object = data, do.scale = F, do.center = T)
  
  # PCA
  data <- RunPCA(object = data, reduction.key = 'prtbpca', reduction.name = 'prtbpca')
  
  # UMAP
  data <- RunUMAP(
    object = data, 
    dims = 1:40, 
    reduction = 'prtbpca', 
    reduction.key = 'prtbumap', 
    reduction.name = 'prtbumap')
  
  # Check if clustering is driven by perturbation status
  DimPlot(
    object = data, 
    group.by = 'NT', 
    pt.size = 0.2, 
    reduction = 'prtbumap', 
    split.by = "NT", 
    ncol = 1, 
    cols = c("grey39","goldenrod3")
  ) + 
    ggtitle("Perturbation Status") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
  
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_umap_pert_status_prtbsig.pdf"), height = 5, width = 6)
  
  DimPlot(
    object = data, 
    group.by = 'Phase', 
    pt.size = 0.2, 
    reduction = 'prtbumap', 
    ncol = 1, 
  ) + 
    ggtitle("Cell cycle phase") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
  
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_umap_phase_prtbsig.pdf"), height = 5, width = 6)
  
  DimPlot(
    object = data, 
    group.by = 'gene', 
    pt.size = 0.2, 
    reduction = 'prtbumap', 
    ncol = 1, 
  ) + 
    ggtitle("Target gene") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
  
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_umap_pert_prtbsig.pdf"), height = 5, width = 6)
  
  # Run mixscape
  # Set lower logfc threshold than default (e.g. 0.025) to detect smaller perturbations
  Idents(data) <- "NT"
  data <- RunMixscape(
    object = data,
    assay = "PRTB", 
    slot = "scale.data", 
    labels = "gene", 
    nt.class.name = "NT", 
    min.de.genes = 5, 
    iter.num = 10, 
    de.assay = "RNA", 
    verbose = F,
    logfc.threshold = 0.025)
  
  # Calculate percentage of KO cells for all target gene classes
  df <- prop.table(table(data$mixscape_class.global, data$sgrna_name),2)
  
  df2 <- reshape2::melt(df)
  df2$Var2 <- as.character(df2$Var2)
  test <- df2[which(df2$Var1 == "KO"),]
  test <- test[order(test$value, decreasing = T),]
  new.levels <- test$Var2
  df2$Var2 <- factor(df2$Var2, levels = new.levels)
  df2$Var1 <- factor(df2$Var1, levels = c("NT", "NP", "KO"))
  df2$gene <- gsub("sg|.1$|.2$|.3$|.4$|.5$|.6$", "", df2$Var2)
  df2$guide_number <- gsub("^.*\\.", "", df2$Var2)
  df3 <- df2[-c(which(df2$gene == "Ctrl")),]
  
  ggplot(df3, aes(x = guide_number, y = value*100, fill= Var1)) +
    geom_bar(stat= "identity") +
    theme_classic()+
    scale_fill_manual(values = c("grey49", "grey79","coral1")) + 
    ylab("% of cells") +
    xlab("sgRNA") +
    theme(axis.text.x = element_text(size = 18, hjust = 1), 
          axis.text.y = element_text(size = 18), 
          axis.title = element_text(size = 16), 
          strip.text = element_text(size=16, face = "bold")) + 
    facet_wrap(vars(gene),ncol = 5, scales = "free") +
    labs(fill = "mixscape class") +
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_mixscape_class_probability.pdf"), height = 7, width = 10)
  

  plot_perturbscore_density <- function(TARGETGENE){
    
    PlotPerturbScore(object = data, 
                     target.gene.ident = TARGETGENE, 
                     group.by = "mixscape_class", 
                     col = "coral2") +
      labs(fill = "mixscape class")
  }
  
  # select pertubed targets
  targets <- df3 %>% filter(Var1 == "KO", value > 0) %>% select(gene) %>% tibble::deframe() %>% unique()
  targets <- targets[!grepl("doublet", targets)]
  
  p1 <- lapply(targets, plot_perturbscore_density)
  m1 <- marrangeGrob(p1, nrow = 2, ncol = 3, top = NULL)
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_mixscape_prtbscore_density.pdf"), m1, height = 7, width = 14)
  
  
  # Run DE analysis and visualize results on a heatmap ordering cells by their posterior probability values
  Idents(object = data) <- "gene"
  
  plot_mixscape_heatmap <- function(TARGETGENE){
    
    MixscapeHeatmap(object = data, 
                    ident.1 = "NT", 
                    ident.2 = TARGETGENE, 
                    balanced = F, 
                    assay = "RNA", 
                    max.genes = 20, angle = 0, 
                    group.by = "mixscape_class", 
                    max.cells.group = 300, 
                    size=3.5,
                    logfc.threshold = 0.025) + NoLegend()
  }
  
  p1 <- lapply(targets, plot_mixscape_heatmap)
  
  m1 <- marrangeGrob(p1, nrow = 2, ncol = 3, top = NULL)
  
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_mixscape_heatmap.pdf"), m1, height = 7, width = 14)
  
  saveRDS(data, paste0("results/sudhl4/mixscape/singlet/", CONDITION, ".rds"))
  
  # Remove non-perturbed cells and run LDA
  Idents(data) <- "mixscape_class.global"
  sub <- subset(data, idents = c("KO", "NT"))
  Idents(object = sub) <- "gene"
  
  # Run LDA.
  sub <- MixscapeLDA(
    object = sub, 
    assay = "RNA", 
    pc.assay = "PRTB", 
    labels = "gene", 
    nt.label = "NT", 
    npcs = 10, 
    logfc.threshold = 0.025, 
    verbose = F)
  
  # Use LDA results to run UMAP
  # Number of dimensions = number of labels minus one (NT)
  sub <- RunUMAP(
    object = sub,
    dims = 1:length(targets),
    reduction = 'lda',
    reduction.key = 'ldaumap',
    reduction.name = 'ldaumap')
  
  saveRDS(sub, paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_ko_nt.rds"))
  
  # Visualize UMAP clustering results
  Idents(sub) <- "mixscape_class"
  sub$mixscape_class <- as.factor(sub$mixscape_class)
  
  # Set colors for each perturbation
  col = setNames(object = scales::hue_pal()(length(targets)+1), nm = levels(sub$mixscape_class))
  col["NT"] <- "grey39"
  
  DimPlot(object = sub, 
          reduction = "ldaumap", 
          repel = T, 
          label.size = 5, 
          label = F) +
    scale_color_manual(values=col, drop=FALSE) + 
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_mixscape_ldaumap_prtb.pdf"), height = 4.5, width = 6.5)
  
  DefaultAssay(sub) <- "RNA"
  FeaturePlot(object = sub, 
              features = c("HLA-A", "HLA-C", "NFKBIA", "PSMB9"),
              reduction = "ldaumap")
  
  ggsave(paste0("results/sudhl4/mixscape/singlet/", CONDITION, "_mixscape_ldaumap_feats.pdf"), height = 5, width = 7)
  
}

run_analyze_mixscape("CROPseq_SUDHL4_NK1_1_16", "nk_1_16")
run_analyze_mixscape("CROPseq_SUDHL4_NK1_1_4", "nk_1_4")
run_analyze_mixscape("CROPseq_SUDHL4_noNK", "nonk")

