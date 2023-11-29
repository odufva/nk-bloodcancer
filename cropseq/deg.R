
# Analyse CROP-seq data
# All perturbation DEGs without logFC cutoff

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(viridis)

# function
de_perturbation <- function(DATA, IDENT1, IDENT2, GENELIST=NULL, RETURNTHRESH=0, LOGFC=0){
  
  result <- FindMarkers(DATA,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        ident.1 = IDENT1,
                        ident.2 = IDENT2,
                        features = GENELIST)
  if(length(result) == 3) {
    result <- data.frame(pval = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  result$perturbation <- IDENT1
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  #result$p_adj <- p.adjust(result$p_val, method = "fdr")
  return(result)
}


# SUDHL4

# load object
crop_seurat <- readRDS("results/sudhl4/sudhl4_crop_seurat_singlet.rds")

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_SUDHL4_NK1_1_16")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_SUDHL4_NK1_1_4")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_SUDHL4_noNK")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/sudhl4/deg/singlet/")

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:16")

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:4")

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="no NK")

result_all <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result_all, "results/sudhl4/deg/singlet/pert_deg_all.txt", row.names = F, quote = F, sep = "\t")


# test DE genes between NK and no NK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_SUDHL4_NK1_1_16", ident.2 = "CROPseq_SUDHL4_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$comparison <- "NK 1:16 vs no NK"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_SUDHL4_NK1_1_4", ident.2 = "CROPseq_SUDHL4_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$comparison <- "NK 1:4 vs no NK"

# 1:4 vs 1:16
nk_de_markers_1_4_vs_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_SUDHL4_NK1_1_4", ident.2 = "CROPseq_SUDHL4_NK1_1_16", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_vs_1_16_control$gene <- rownames(nk_de_markers_1_4_vs_1_16_control)
nk_de_markers_1_4_vs_1_16_control <- nk_de_markers_1_4_vs_1_16_control[,c(6,1:5)]
nk_de_markers_1_4_vs_1_16_control$comparison <- "NK 1:4 vs NK 1:16"

nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control, nk_de_markers_1_4_vs_1_16_control)
fwrite(nk_de_markers, "results/sudhl4/deg/singlet/nk_nonk_deg_all.txt", quote = F, row.names = F, sep = "\t")




# NALM6

# load object
crop_seurat <- readRDS("results/nalm6/nalm6_crop_seurat_singlet.rds")

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_NALM6_NK1_1_16")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_NALM6_NK1_1_4")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_NALM6_noNK")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/nalm6/deg/singlet/")

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:16")

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:4")

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="no NK")

result_all <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result_all, "results/nalm6/deg/singlet/pert_deg_all.txt", row.names = F, quote = F, sep = "\t")


# test DE genes between NK and no NK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_NALM6_NK1_1_16", ident.2 = "CROPseq_NALM6_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$comparison <- "NK 1:16 vs no NK"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_NALM6_NK1_1_4", ident.2 = "CROPseq_NALM6_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$comparison <- "NK 1:4 vs no NK"

# 1:4 vs 1:16
nk_de_markers_1_4_vs_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_NALM6_NK1_1_4", ident.2 = "CROPseq_NALM6_NK1_1_16", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_vs_1_16_control$gene <- rownames(nk_de_markers_1_4_vs_1_16_control)
nk_de_markers_1_4_vs_1_16_control <- nk_de_markers_1_4_vs_1_16_control[,c(6,1:5)]
nk_de_markers_1_4_vs_1_16_control$comparison <- "NK 1:4 vs NK 1:16"

nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control, nk_de_markers_1_4_vs_1_16_control)
fwrite(nk_de_markers, "results/nalm6/deg/singlet/nk_nonk_deg_all.txt", quote = F, row.names = F, sep = "\t")




# MM1S

# load object
crop_seurat <- readRDS("results/mm1s/mm1s_crop_seurat_singlet.rds")

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_MM1S_NK1_1_16")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_MM1S_NK1_1_4")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_MM1S_noNK")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/mm1s/deg/singlet/")

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:16")

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:4")

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="no NK")

result_all <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result_all, "results/mm1s/deg/singlet/pert_deg_all.txt", row.names = F, quote = F, sep = "\t")


# test DE genes between NK and no NK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_MM1S_NK1_1_16", ident.2 = "CROPseq_MM1S_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$comparison <- "NK 1:16 vs no NK"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_MM1S_NK1_1_4", ident.2 = "CROPseq_MM1S_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$comparison <- "NK 1:4 vs no NK"

# 1:4 vs 1:16
nk_de_markers_1_4_vs_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_MM1S_NK1_1_4", ident.2 = "CROPseq_MM1S_NK1_1_16", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_vs_1_16_control$gene <- rownames(nk_de_markers_1_4_vs_1_16_control)
nk_de_markers_1_4_vs_1_16_control <- nk_de_markers_1_4_vs_1_16_control[,c(6,1:5)]
nk_de_markers_1_4_vs_1_16_control$comparison <- "NK 1:4 vs NK 1:16"

nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control, nk_de_markers_1_4_vs_1_16_control)
fwrite(nk_de_markers, "results/mm1s/deg/singlet/nk_nonk_deg_all.txt", quote = F, row.names = F, sep = "\t")


# MM1S

# load object
crop_seurat <- readRDS("results/mm1sa/mm1sa_crop_seurat_singlet_2.rds")

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "NK2587_1_16_MM1S_CROPseq_activ")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "NK2587_1_4_MM1S_CROPseq_activ")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "MM1S_CROPseq_activ")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/mm1sa/")
dir.create("results/mm1sa/deg/")
dir.create("results/mm1sa/deg/singlet/")

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:16")

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:4")

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="no NK")

result_all <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result_all, "results/mm1sa/deg/singlet/pert_deg_all.txt", row.names = F, quote = F, sep = "\t")


# LP1

# load object
crop_seurat <- readRDS("results/lp1/lp1_crop_seurat_singlet_2.rds")

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_LP1_NK1_1_16")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nk_1_4 <- subset(crop_seurat, orig.ident == "CROPseq_LP1_NK1_1_4")
table(crop_seurat_nk_1_4$gene)
Idents(crop_seurat_nk_1_4) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_LP1_noNK")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/lp1/deg/singlet/")

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:16")

# NK 1:4
result_nk_1_4 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_4, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:4")

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="no NK")

result_all <- rbind(result_nk_1_16, result_nk_1_4, result_nonk)
fwrite(result_all, "results/lp1/deg/singlet/pert_deg_all.txt", row.names = F, quote = F, sep = "\t")


# test DE genes between NK and no NK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_LP1_NK1_1_16", ident.2 = "CROPseq_LP1_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$comparison <- "NK 1:16 vs no NK"

# 1:4
nk_de_markers_1_4_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_LP1_NK1_1_4", ident.2 = "CROPseq_LP1_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_control$gene <- rownames(nk_de_markers_1_4_control)
nk_de_markers_1_4_control <- nk_de_markers_1_4_control[,c(6,1:5)]
nk_de_markers_1_4_control$comparison <- "NK 1:4 vs no NK"

# 1:4 vs 1:16
nk_de_markers_1_4_vs_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_LP1_NK1_1_4", ident.2 = "CROPseq_LP1_NK1_1_16", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_4_vs_1_16_control$gene <- rownames(nk_de_markers_1_4_vs_1_16_control)
nk_de_markers_1_4_vs_1_16_control <- nk_de_markers_1_4_vs_1_16_control[,c(6,1:5)]
nk_de_markers_1_4_vs_1_16_control$comparison <- "NK 1:4 vs NK 1:16"

nk_de_markers <- rbind(nk_de_markers_1_16_control, nk_de_markers_1_4_control, nk_de_markers_1_4_vs_1_16_control)
fwrite(nk_de_markers, "results/lp1/deg/singlet/nk_nonk_deg_all.txt", quote = F, row.names = F, sep = "\t")



# K562

# load object
crop_seurat <- readRDS("results/k562/k562_crop_seurat_singlet.rds")

# subset to conditions
crop_seurat_nk_1_16 <- subset(crop_seurat, orig.ident == "CROPseq_K562_NK1_1_16")
table(crop_seurat_nk_1_16$gene)
Idents(crop_seurat_nk_1_16) <- "gene"

crop_seurat_nonk <- subset(crop_seurat, orig.ident == "CROPseq_K562_noNK")
table(crop_seurat_nonk$gene)
Idents(crop_seurat_nonk) <- "gene"

# DEGs perturbation vs control
dir.create("results/k562/deg/singlet/")

# vector of all perturbations
all_pert <- as.character(unique(crop_seurat_nk_1_16$gene))
all_pert <- all_pert[!grepl("Control", all_pert)]
all_pert <- all_pert[!is.na(all_pert)]


## Get DEG perturbations vs ctrl

# NK 1:16
result_nk_1_16 <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nk_1_16, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="NK 1:16")

# no NK
result_nonk <- lapply(all_pert, de_perturbation, DATA = crop_seurat_nonk, IDENT2 = "Control") %>% bind_rows() %>%
  mutate(condition ="no NK")

result_all <- rbind(result_nk_1_16, result_nonk)
fwrite(result_all, "results/k562/deg/singlet/pert_deg_all.txt", row.names = F, quote = F, sep = "\t")


# test DE genes between NK and no nK

# 1:16
nk_de_markers_1_16_control <- FindMarkers(subset(crop_seurat, gene %in% "Control"), ident.1 = "CROPseq_K562_NK1_1_16", ident.2 = "CROPseq_K562_noNK", group.by = "orig.ident", test.use = "t", logfc.threshold = 0)
nk_de_markers_1_16_control$gene <- rownames(nk_de_markers_1_16_control)
nk_de_markers_1_16_control <- nk_de_markers_1_16_control[,c(6,1:5)]
nk_de_markers_1_16_control$comparison <- "NK 1:16 vs no NK"

nk_de_markers <- nk_de_markers_1_16_control
fwrite(nk_de_markers, "results/k562/deg/singlet/nk_nonk_deg_all.txt", quote = F, row.names = F, sep = "\t")

