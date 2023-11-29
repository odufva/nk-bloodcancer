
# Analyse core NK responses score in SUDHL4 CROP-seq data

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(viridis)

theme_set(theme_classic(base_size = 12))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

# load object
crop_seurat <- readRDS("results/sudhl4/sudhl4_crop_seurat_singlet.rds")

data <- fread("results/sudhl4/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "SUDHL4")

# Gene set scores
pathways <- fgsea::gmtPathways("h.all.v7.0.symbols.gmt")


# nkresponse <- data %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>% select(gene) %>% tibble::deframe() %>% unique()
nkresponse <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "STAT1", "IRF1", "IRF9", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "UBE2L6", "GNLY", "CCL5")
ifn <- fread("results/combine/ifng_score.txt", data.table = F)$gene
nfkb <- fread("results/combine/nfkb_score.txt", data.table = F)$gene

crop_seurat <- AddModuleScore(crop_seurat, features = list(nkresponse), name = c("nk_response_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(ifn), name = c("ifng_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(nfkb), name = c("nfkb_score"))

scores <- t(data.frame(nkresponse = crop_seurat$nk_response_score1, ifng = crop_seurat$ifng_score1, nfkb = crop_seurat$nfkb_score1))
crop_seurat[["SCR"]] <- CreateAssayObject(data = scores)

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_SUDHL4_NK1_1_16" & status == "single")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_SUDHL4_NK1_1_4" & status == "single")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_SUDHL4_noNK" & status == "single")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/sudhl4/deg/modulescore/")

# function
de_perturbation <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0.05, LOGFC=0){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST,
                        assay = "SCR")
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  result$perturbation <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  return(result)
}

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows()
result_nk_1_16$condition  <- "NK 1:16"

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows()
result_nk_1_4$condition  <- "NK 1:4"

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows()
result_nonk$condition  <- "no NK"

result <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result, "results/sudhl4/deg/modulescore/pert_scores.txt", row.names = F, quote = F, sep = "\t")


# dot plots
dir.create("results/sudhl4/deg/modulescore/dotplots")

plot_dotplot <- function(DATA, PLOTNAME, NGENES=10){
  genelist <- DATA %>%
    filter(!grepl("doublet", perturbation)) %>%
    group_by(perturbation) %>%
    top_n(NGENES, dplyr::desc(p_val)) %>%
    ungroup() %>%
    select(gene) %>%
    unique() %>%
    tibble::deframe()
  
  plotdata <- DATA %>%
    filter(!grepl("doublet", perturbation) & gene %in% genelist) %>%
    filter(p_val_adj < 0.05)
  
  
  plotdata_wide <- dcast(plotdata, gene ~ perturbation, value.var = "avg_log2FC")
  plotdata_wide[is.na(plotdata_wide)] <- 0
  rownames(plotdata_wide) <- plotdata_wide$gene
  plotdata_wide$gene <- NULL
  
  # clustering with hclust on row and on column
  dd.col <- as.dendrogram(hclust(dist(1-cor(t(plotdata_wide), method="spearman"))))
  dd.row <- as.dendrogram(hclust(dist(1-cor(plotdata_wide, method="spearman"))))
  
  # ordering based on clustering
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  
  p <- plotdata %>%
    ungroup() %>%
    mutate(perturbation = factor(perturbation, levels = sort(unique(perturbation))[row.ord]),
           gene = factor(gene, levels = sort(unique(gene))[col.ord])) %>%
    
    ggplot(aes(y = perturbation, x = gene, size = -log10(p_val_adj), color = avg_log2FC)) +
    geom_point() +
    scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                          type = "div", limits = max(abs(plotdata$avg_log2FC)) * c(-1, 1)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",
         color = "Fold change (log2)") +
    ylab("") +
    xlab("") 
  ggsave(paste0("results/sudhl4/deg/modulescore/dotplots/", PLOTNAME, "_dotplot_top", NGENES, ".pdf"), p, height = 4, width = 4)
}

plot_dotplot(result_nk_1_16, "nk_1_16")
plot_dotplot(result_nk_1_4, "nk_1_4")
plot_dotplot(result_nonk, "nonk")


# test DE genes between NK and no nK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                          ident.1 = "CROPseq_SUDHL4_NK1_1_16",
                                          ident.2 = "CROPseq_SUDHL4_noNK",
                                          group.by = "orig.ident",
                                          test.use = "t",
                                          logfc.threshold = 0.025,
                                          assay = "SCR")
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$condition <- "NK 1:16"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                         ident.1 = "CROPseq_SUDHL4_NK1_1_4",
                                         ident.2 = "CROPseq_SUDHL4_noNK",
                                         group.by = "orig.ident",
                                         test.use = "t",
                                         logfc.threshold = 0.025,
                                         assay = "SCR")
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$condition <- "NK 1:4"


nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control)
fwrite(nk_de_markers, "results/sudhl4/deg/modulescore/condition_scores.txt", row.names = F, quote = F, sep = "\t")


## ----------------------------------------


# Analyse core NK cell responses score in K562 CROP-seq data

# load object
crop_seurat <- readRDS("results/k562/k562_crop_seurat_singlet.rds")

data <- fread("results/k562/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "K562")

# Gene set scores
pathways <- fgsea::gmtPathways("h.all.v7.0.symbols.gmt")

nkresponse <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "STAT1", "IRF1", "IRF9", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "UBE2L6", "GNLY", "CCL5")
ifn <- fread("results/combine/ifng_score.txt", data.table = F)$gene
nfkb <- fread("results/combine/nfkb_score.txt", data.table = F)$gene

crop_seurat <- AddModuleScore(crop_seurat, features = list(nkresponse), name = c("nk_response_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(ifn), name = c("ifng_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(nfkb), name = c("nfkb_score"))

scores <- t(data.frame(nkresponse = crop_seurat$nk_response_score1, ifng = crop_seurat$ifng_score1, nfkb = crop_seurat$nfkb_score1))
crop_seurat[["SCR"]] <- CreateAssayObject(data = scores)

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_K562_NK1_1_16" & status == "single")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_K562_noNK" & status == "single")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/k562/deg/modulescore/")

# function
de_perturbation <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0.05, LOGFC=0){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST,
                        assay = "SCR")
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  result$perturbation <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  return(result)
}

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows()
result_nk_1_16$condition  <- "NK 1:16"

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows()
result_nonk$condition  <- "no NK"

result <- rbind(result_nk_1_16, result_nonk)
fwrite(result, "results/k562/deg/modulescore/pert_scores.txt", row.names = F, quote = F, sep = "\t")


# dot plots
dir.create("results/k562/deg/modulescore/dotplots")

plot_dotplot <- function(DATA, PLOTNAME, NGENES=10){
  genelist <- DATA %>%
    filter(!grepl("doublet", perturbation)) %>%
    group_by(perturbation) %>%
    top_n(NGENES, dplyr::desc(p_val)) %>%
    ungroup() %>%
    select(gene) %>%
    unique() %>%
    tibble::deframe()
  
  plotdata <- DATA %>%
    filter(!grepl("doublet", perturbation) & gene %in% genelist) %>%
    filter(p_val_adj < 0.05)
  
  
  plotdata_wide <- dcast(plotdata, gene ~ perturbation, value.var = "avg_log2FC")
  plotdata_wide[is.na(plotdata_wide)] <- 0
  rownames(plotdata_wide) <- plotdata_wide$gene
  plotdata_wide$gene <- NULL
  
  # clustering with hclust on row and on column
  dd.col <- as.dendrogram(hclust(dist(1-cor(t(plotdata_wide), method="spearman"))))
  dd.row <- as.dendrogram(hclust(dist(1-cor(plotdata_wide, method="spearman"))))
  
  # ordering based on clustering
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  
  p <- plotdata %>%
    ungroup() %>%
    mutate(perturbation = factor(perturbation, levels = sort(unique(perturbation))[row.ord]),
           gene = factor(gene, levels = sort(unique(gene))[col.ord])) %>%
    
    ggplot(aes(y = perturbation, x = gene, size = -log10(p_val_adj), color = avg_log2FC)) +
    geom_point() +
    scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                          type = "div", limits = max(abs(plotdata$avg_log2FC)) * c(-1, 1)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",
         color = "Fold change (log2)") +
    ylab("") +
    xlab("") 
  ggsave(paste0("results/k562/deg/modulescore/dotplots/", PLOTNAME, "_dotplot_top", NGENES, ".pdf"), p, height = 4, width = 4)
}

plot_dotplot(result_nk_1_16, "nk_1_16")
plot_dotplot(result_nonk, "nonk")


# test DE genes between NK and no nK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                          ident.1 = "CROPseq_K562_NK1_1_16",
                                          ident.2 = "CROPseq_K562_noNK",
                                          group.by = "orig.ident",
                                          test.use = "t",
                                          logfc.threshold = 0.025,
                                          assay = "SCR")
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$condition <- "NK 1:16"


nk_de_markers <- rbind(nk_de_markers_1_16_control)
fwrite(nk_de_markers, "results/k562/deg/modulescore/condition_scores.txt", row.names = F, quote = F, sep = "\t")


## -------------------------------------


# Analyse core NK cell response score in MM1S CROP-seq data


# load object
crop_seurat <- readRDS("results/mm1s/mm1s_crop_seurat_singlet.rds")

data <- fread("results/mm1s/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "MM1S")

# Gene set scores
pathways <- fgsea::gmtPathways("h.all.v7.0.symbols.gmt")

nkresponse <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "STAT1", "IRF1", "IRF9", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "UBE2L6", "GNLY", "CCL5")
ifn <- fread("results/combine/ifng_score.txt", data.table = F)$gene
nfkb <- fread("results/combine/nfkb_score.txt", data.table = F)$gene

crop_seurat <- AddModuleScore(crop_seurat, features = list(nkresponse), name = c("nk_response_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(ifn), name = c("ifng_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(nfkb), name = c("nfkb_score"))

scores <- t(data.frame(nkresponse = crop_seurat$nk_response_score1, ifng = crop_seurat$ifng_score1, nfkb = crop_seurat$nfkb_score1))
crop_seurat[["SCR"]] <- CreateAssayObject(data = scores)

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_MM1S_NK1_1_16" & status == "single")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_MM1S_NK1_1_4" & status == "single")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_MM1S_noNK" & status == "single")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/mm1s/")
dir.create("results/mm1s/deg/")
dir.create("results/mm1s/deg/modulescore/")

# function
de_perturbation <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0.05, LOGFC=0){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST,
                        assay = "SCR")
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  result$perturbation <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  return(result)
}

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows()
result_nk_1_16$condition  <- "NK 1:16"

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows()
result_nk_1_4$condition  <- "NK 1:4"

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows()
result_nonk$condition  <- "no NK"

result <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result, "results/mm1s/deg/modulescore/pert_scores.txt", row.names = F, quote = F, sep = "\t")


# dot plots
dir.create("results/mm1s/deg/modulescore/dotplots")

plot_dotplot <- function(DATA, PLOTNAME, NGENES=10){
  genelist <- DATA %>%
    filter(!grepl("doublet", perturbation)) %>%
    group_by(perturbation) %>%
    top_n(NGENES, dplyr::desc(p_val)) %>%
    ungroup() %>%
    select(gene) %>%
    unique() %>%
    tibble::deframe()
  
  plotdata <- DATA %>%
    filter(!grepl("doublet", perturbation) & gene %in% genelist) %>%
    filter(p_val_adj < 0.05)
  
  
  plotdata_wide <- dcast(plotdata, gene ~ perturbation, value.var = "avg_log2FC")
  plotdata_wide[is.na(plotdata_wide)] <- 0
  rownames(plotdata_wide) <- plotdata_wide$gene
  plotdata_wide$gene <- NULL
  
  # clustering with hclust on row and on column
  dd.col <- as.dendrogram(hclust(dist(1-cor(t(plotdata_wide), method="spearman"))))
  dd.row <- as.dendrogram(hclust(dist(1-cor(plotdata_wide, method="spearman"))))
  
  # ordering based on clustering
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  
  p <- plotdata %>%
    ungroup() %>%
    mutate(perturbation = factor(perturbation, levels = sort(unique(perturbation))[row.ord]),
           gene = factor(gene, levels = sort(unique(gene))[col.ord])) %>%
    
    ggplot(aes(y = perturbation, x = gene, size = -log10(p_val_adj), color = avg_log2FC)) +
    geom_point() +
    scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                          type = "div", limits = max(abs(plotdata$avg_log2FC)) * c(-1, 1)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",
         color = "Fold change (log2)") +
    ylab("") +
    xlab("") 
  ggsave(paste0("results/mm1s/deg/modulescore/dotplots/", PLOTNAME, "_dotplot_top", NGENES, ".pdf"), p, height = 4, width = 4)
}

plot_dotplot(result_nk_1_16, "nk_1_16")
plot_dotplot(result_nk_1_4, "nk_1_4")
plot_dotplot(result_nonk, "nonk")


# test DE genes between NK and no nK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                          ident.1 = "CROPseq_MM1S_NK1_1_16",
                                          ident.2 = "CROPseq_MM1S_noNK",
                                          group.by = "orig.ident",
                                          test.use = "t",
                                          logfc.threshold = 0.025,
                                          assay = "SCR")
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$condition <- "NK 1:16"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                         ident.1 = "CROPseq_MM1S_NK1_1_4",
                                         ident.2 = "CROPseq_MM1S_noNK",
                                         group.by = "orig.ident",
                                         test.use = "t",
                                         logfc.threshold = 0.025,
                                         assay = "SCR")
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$condition <- "NK 1:4"


nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control)
fwrite(nk_de_markers, "results/mm1s/deg/modulescore/condition_scores.txt", row.names = F, quote = F, sep = "\t")

## ------------------------------------


# Analyse core NK cell response score in LP1 CROP-seq data

# load object
crop_seurat <- readRDS("results/lp1/lp1_crop_seurat_singlet_2.rds")

data <- fread("results/lp1/deg/singlet/nk_nonk_deg_qc3.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "LP1")

# Gene set scores
pathways <- fgsea::gmtPathways("h.all.v7.0.symbols.gmt")

nkresponse <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "STAT1", "IRF1", "IRF9", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "UBE2L6", "GNLY", "CCL5")
ifn <- fread("results/combine/ifng_score.txt", data.table = F)$gene
nfkb <- fread("results/combine/nfkb_score.txt", data.table = F)$gene

crop_seurat <- AddModuleScore(crop_seurat, features = list(nkresponse), name = c("nk_response_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(ifn), name = c("ifng_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(nfkb), name = c("nfkb_score"))

scores <- t(data.frame(nkresponse = crop_seurat$nk_response_score1, ifng = crop_seurat$ifng_score1, nfkb = crop_seurat$nfkb_score1))
crop_seurat[["SCR"]] <- CreateAssayObject(data = scores)

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_LP1_NK1_1_16" & status == "single")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_LP1_NK1_1_4" & status == "single")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_LP1_noNK" & status == "single")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/lp1/")
dir.create("results/lp1/deg/")
dir.create("results/lp1/deg/modulescore/")

# function
de_perturbation <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0.05, LOGFC=0){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST,
                        assay = "SCR")
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  result$perturbation <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  return(result)
}

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows()
result_nk_1_16$condition  <- "NK 1:16"

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows()
result_nk_1_4$condition  <- "NK 1:4"

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows()
result_nonk$condition  <- "no NK"

result <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result, "results/lp1/deg/modulescore/pert_scores.txt", row.names = F, quote = F, sep = "\t")


# dot plots
dir.create("results/lp1/deg/modulescore/dotplots")

plot_dotplot <- function(DATA, PLOTNAME, NGENES=10){
  genelist <- DATA %>%
    filter(!grepl("doublet", perturbation)) %>%
    group_by(perturbation) %>%
    top_n(NGENES, dplyr::desc(p_val)) %>%
    ungroup() %>%
    select(gene) %>%
    unique() %>%
    tibble::deframe()
  
  plotdata <- DATA %>%
    filter(!grepl("doublet", perturbation) & gene %in% genelist) %>%
    filter(p_val_adj < 0.05)
  
  
  plotdata_wide <- dcast(plotdata, gene ~ perturbation, value.var = "avg_log2FC")
  plotdata_wide[is.na(plotdata_wide)] <- 0
  rownames(plotdata_wide) <- plotdata_wide$gene
  plotdata_wide$gene <- NULL
  
  # clustering with hclust on row and on column
  dd.col <- as.dendrogram(hclust(dist(1-cor(t(plotdata_wide), method="spearman"))))
  dd.row <- as.dendrogram(hclust(dist(1-cor(plotdata_wide, method="spearman"))))
  
  # ordering based on clustering
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  
  p <- plotdata %>%
    ungroup() %>%
    mutate(perturbation = factor(perturbation, levels = sort(unique(perturbation))[row.ord]),
           gene = factor(gene, levels = sort(unique(gene))[col.ord])) %>%
    
    ggplot(aes(y = perturbation, x = gene, size = -log10(p_val_adj), color = avg_log2FC)) +
    geom_point() +
    scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                          type = "div", limits = max(abs(plotdata$avg_log2FC)) * c(-1, 1)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",
         color = "Fold change (log2)") +
    ylab("") +
    xlab("") 
  ggsave(paste0("results/lp1/deg/modulescore/dotplots/", PLOTNAME, "_dotplot_top", NGENES, ".pdf"), p, height = 4, width = 4)
}

plot_dotplot(result_nk_1_16, "nk_1_16")
plot_dotplot(result_nk_1_4, "nk_1_4")
plot_dotplot(result_nonk, "nonk")


# test DE genes between NK and no nK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                          ident.1 = "CROPseq_LP1_NK1_1_16",
                                          ident.2 = "CROPseq_LP1_noNK",
                                          group.by = "orig.ident",
                                          test.use = "t",
                                          logfc.threshold = 0.025,
                                          assay = "SCR")
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$condition <- "NK 1:16"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                         ident.1 = "CROPseq_LP1_NK1_1_4",
                                         ident.2 = "CROPseq_LP1_noNK",
                                         group.by = "orig.ident",
                                         test.use = "t",
                                         logfc.threshold = 0.025,
                                         assay = "SCR")
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$condition <- "NK 1:4"


nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control)
fwrite(nk_de_markers, "results/lp1/deg/modulescore/condition_scores.txt", row.names = F, quote = F, sep = "\t")


## --------------------------------------



# Analyse core NK cell response score in NALM6 CROP-seq data


# load object
crop_seurat <- readRDS("results/nalm6/nalm6_crop_seurat_singlet.rds")

data <- fread("results/nalm6/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "NALM6")

# Gene set scores
pathways <- fgsea::gmtPathways("h.all.v7.0.symbols.gmt")

nkresponse <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "STAT1", "IRF1", "IRF9", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "UBE2L6", "GNLY", "CCL5")
ifn <- fread("results/combine/ifng_score.txt", data.table = F)$gene
nfkb <- fread("results/combine/nfkb_score.txt", data.table = F)$gene

crop_seurat <- AddModuleScore(crop_seurat, features = list(nkresponse), name = c("nk_response_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(ifn), name = c("ifng_score"))
crop_seurat <- AddModuleScore(crop_seurat, features = list(nfkb), name = c("nfkb_score"))

scores <- t(data.frame(nkresponse = crop_seurat$nk_response_score1, ifng = crop_seurat$ifng_score1, nfkb = crop_seurat$nfkb_score1))
crop_seurat[["SCR"]] <- CreateAssayObject(data = scores)

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_NALM6_NK1_1_16" & status == "single")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_NALM6_NK1_1_4" & status == "single")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_NALM6_noNK" & status == "single")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/nalm6/deg/modulescore/")

# function
de_perturbation <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0.05, LOGFC=0){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST,
                        assay = "SCR")
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  result$perturbation <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  return(result)
}

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows()
result_nk_1_16$condition  <- "NK 1:16"

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows()
result_nk_1_4$condition  <- "NK 1:4"

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows()
result_nonk$condition  <- "no NK"

result <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result, "results/nalm6/deg/modulescore/pert_scores.txt", row.names = F, quote = F, sep = "\t")


# dot plots
dir.create("results/nalm6/deg/modulescore/dotplots")

plot_dotplot <- function(DATA, PLOTNAME, NGENES=10){
  genelist <- DATA %>%
    filter(!grepl("doublet", perturbation)) %>%
    group_by(perturbation) %>%
    top_n(NGENES, dplyr::desc(p_val)) %>%
    ungroup() %>%
    select(gene) %>%
    unique() %>%
    tibble::deframe()
  
  plotdata <- DATA %>%
    filter(!grepl("doublet", perturbation) & gene %in% genelist) %>%
    filter(p_val_adj < 0.05)
  
  
  plotdata_wide <- dcast(plotdata, gene ~ perturbation, value.var = "avg_log2FC")
  plotdata_wide[is.na(plotdata_wide)] <- 0
  rownames(plotdata_wide) <- plotdata_wide$gene
  plotdata_wide$gene <- NULL
  
  # clustering with hclust on row and on column
  dd.col <- as.dendrogram(hclust(dist(1-cor(t(plotdata_wide), method="spearman"))))
  dd.row <- as.dendrogram(hclust(dist(1-cor(plotdata_wide, method="spearman"))))
  
  # ordering based on clustering
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  
  p <- plotdata %>%
    ungroup() %>%
    mutate(perturbation = factor(perturbation, levels = sort(unique(perturbation))[row.ord]),
           gene = factor(gene, levels = sort(unique(gene))[col.ord])) %>%
    
    ggplot(aes(y = perturbation, x = gene, size = -log10(p_val_adj), color = avg_log2FC)) +
    geom_point() +
    scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                          type = "div", limits = max(abs(plotdata$avg_log2FC)) * c(-1, 1)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",
         color = "Fold change (log2)") +
    ylab("") +
    xlab("") 
  ggsave(paste0("results/nalm6/deg/modulescore/dotplots/", PLOTNAME, "_dotplot_top", NGENES, ".pdf"), p, height = 4, width = 4)
}

plot_dotplot(result_nk_1_16, "nk_1_16")
plot_dotplot(result_nk_1_4, "nk_1_4")
plot_dotplot(result_nonk, "nonk")


# test DE genes between NK and no nK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                          ident.1 = "CROPseq_NALM6_NK1_1_16",
                                          ident.2 = "CROPseq_NALM6_noNK",
                                          group.by = "orig.ident",
                                          test.use = "t",
                                          logfc.threshold = 0,
                                          assay = "SCR")
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$condition <- "NK 1:16"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"),
                                         ident.1 = "CROPseq_NALM6_NK1_1_4",
                                         ident.2 = "CROPseq_NALM6_noNK",
                                         group.by = "orig.ident",
                                         test.use = "t",
                                         logfc.threshold = 0.025,
                                         assay = "SCR")
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$condition <- "NK 1:4"


nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control)
fwrite(nk_de_markers, "results/nalm6/deg/modulescore/condition_scores.txt", row.names = F, quote = F, sep = "\t")


