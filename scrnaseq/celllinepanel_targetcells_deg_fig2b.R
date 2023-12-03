
# NK cell co-culture with panel of 26 cell lines, hashing scRNA-seq analysis
# DEGs with small clusters removed
# 3 NK cell donors (NK1, NK2, NK3) from original and revision experiments

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
library(tidyr)



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


dir.create("results/celllinepanel_combined/targets/deg")


"%ni%" <- Negate("%in%")

# load merged data
hash_seurat_singlet <- readRDS("results/celllinepanel_combined/hash_seurat_singlet.rds")

# remove NK cell clusters
hash_seurat_targets <- subset(hash_seurat_singlet, seurat_clusters %ni% c("0", "23", "28"))


cluster_cellcounts <- table(hash_seurat_targets$hash.ID, hash_seurat_targets$seurat_clusters)


# DEGs in NK-treated vs control target cells

# function
de_target <- function(DATA, IDENT2, GENELIST=NULL, RETURNTHRESH=0.05, LOGFC=0.1, MINPCT = 0.1, ASSAY = "RNA"){
  
  # subset data to include only clusters comprising > 10% of all target cells of one cell line
  data_subset <- subset(DATA, hash.ID %in% c(IDENT2, paste0(IDENT2, "-NK1"), paste0(IDENT2, "-NK2"),
                                             paste0(IDENT2, "-NK3")))
  data_subset$hash.ID <- gsub("NK1|NK2|NK3", "NK", data_subset$hash.ID)
  
  cluster_cellcounts <- table(data_subset$seurat_clusters)
  maxclusters <- names(cluster_cellcounts[cluster_cellcounts>0.05*sum(cluster_cellcounts)])
  
  data_subset <- subset(data_subset, seurat_clusters %in% maxclusters)
  
  Idents(data_subset) <- "hash.ID"
  
  result <- FindMarkers(data_subset,
                        test.use = "t",
                        return.thresh = RETURNTHRESH, 
                        logfc.threshold = LOGFC,
                        min.pct = MINPCT,
                        ident.1 = paste0(IDENT2, "-NK"),
                        ident.2 = IDENT2,
                        features = GENELIST,
                        assay = ASSAY)
  
  if(length(result) == 3) {
    result <- data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)
  }
  
  result$cell_line <- IDENT2
  result$gene <- rownames(result)
  result <- result[,c(7,1:6)]
  result$p_adj <- p.adjust(result$p_val, method = "fdr")
  return(result)
}

# vector of all perturbations (no NK-only)
nk_treated <- as.character(unique(hash_seurat_targets$hash.ID))
nk_treated <- nk_treated[grepl("NK2", nk_treated) & !grepl("^NK", nk_treated)]
untreated <- gsub("-NK2", "", nk_treated)

## Get DEG

Idents(hash_seurat_targets) <- "hash.ID"

# targets NK-treated vs resting
result_target <- lapply(untreated, de_target, DATA = hash_seurat_targets) %>% bind_rows()
fwrite(result_target, "results/celllinepanel_combined/targets/deg/deg_mainclusters.txt", row.names = F, quote = F, sep = "\t")


# all genes without threshold

# targets NK-treated vs resting
result_target <- lapply(untreated, de_target, DATA = hash_seurat_targets_combined, RETURNTHRESH = 0, LOGFC = 0, MINPCT = 0) %>% bind_rows()
fwrite(result_target, "results/celllinepanel_combined/targets/deg/deg_mainclusters_all.txt", row.names = F, quote = F, sep = "\t")


## Volcano plots

# volcano plot function
volcanoplot <- function(CELLLINE, DATA, CUTOFF_FDR=1, CUTOFF_PVAL=1, GENES = NULL, TOPGENES=F, NTOPGENES=30){
  
  plotdata <- DATA %>%
    filter(cell_line == CELLLINE)
  
  if (TOPGENES == T) {
    labelgenes <- plotdata %>%
      top_n(-NTOPGENES, wt = p_val) %>%
      dplyr::select(gene) %>%
      tibble::deframe()
  }
  
  else
    
  { labelgenes <- plotdata$gene }
  
  ggplot(plotdata, aes(x=avg_log2FC, y=-log10(p_val), label = gene)) +
    geom_point(aes(color = ifelse(avg_log2FC > 0, "up", "down"))) +
    geom_point(data = plotdata[plotdata$p_val_adj>0.05,], color = "grey80") +
    geom_text_repel(aes(label=ifelse(p_val_adj <= CUTOFF_FDR & p_val <= CUTOFF_PVAL & gene %in% labelgenes, as.character(gene),'')),
                    point.padding = 0.2,
                    max.overlaps = 100) +
    ylab("p value (-log10)") +
    xlab("Average fold change (log2)") +
    scale_color_manual(values = c(brewer.pal(11, "RdBu")[10], brewer.pal(11, "RdBu")[2])) +
    ggtitle(CELLLINE) +
    guides(color = F) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0, face = "plain"))
  
}

# targets
p1 <- lapply(untreated, volcanoplot, DATA = result_target, TOPGENES = T, NTOPGENES=50)
m1 <- marrangeGrob(p1, nrow = 2, ncol = 3, top = NULL)
ggsave("results/celllinepanel2/hto_singlets/volcanoplots_targets_mainclusters.pdf", m1, height = 10, width = 20)

# targets (d0)
p1 <- lapply(untreated, volcanoplot, DATA = result_target_d0, TOPGENES = T, NTOPGENES=50)
m1 <- marrangeGrob(p1, nrow = 2, ncol = 3, top = NULL)
ggsave("results/celllinepanel2/hto_singlets/volcanoplots_targets_d0_mainclusters.pdf", m1, height = 10, width = 20)

## -----------------------------


# Module score DE analysis

core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")

hash_seurat_targets_combined <- AddModuleScore(hash_seurat_targets_combined, features = list(core), name = c("core_nk_score"))

scores <- t(data.frame(core_nk = hash_seurat_targets_combined$core_nk_score1))
hash_seurat_targets_combined[["SCR"]] <- CreateAssayObject(data = scores)


# Differential module score in NK-treated vs control target cells

# targets NK-treated vs resting
result_target <- lapply(untreated, de_target, DATA = hash_seurat_targets_combined, ASSAY = "SCR", LOGFC = 0) %>% bind_rows()
fwrite(result_target, "results/celllinepanel_combined/targets/targets_differential_modulescore_mainclusters.txt", row.names = F, quote = F, sep = "\t")



## ---------------------------------------------------

## DE analysis of all target cell lines combined

hash_seurat_target <- readRDS("results/celllinepanel_combined/targets/hash_seurat_targets_combined_mainclusters.rds")

cells <- colnames(hash_seurat_target)[!grepl("PBMC", hash_seurat_target$hash.ID)]

hash_seurat_target <- subset(hash_seurat_target, cells = cells)

hash_seurat_target$treatment <- "Targets only"
hash_seurat_target$treatment[grepl("NK", hash_seurat_target$hash.ID)] <- "NK-treated"

sort(table(hash_seurat_target$hash.ID))
Idents(hash_seurat_target) <- "hash.ID"

# downsample to 25 cell per cell line
hash_seurat_target_downsampled <- subset(hash_seurat_target, downsample = 25)

Idents(hash_seurat_target_downsampled) <- "treatment"

# Expanded NK-treated vs untreated
result <- FindMarkers(hash_seurat_target_downsampled,
                      test.use = "t",
                      logfc.threshold = 0,
                      ident.1 = "NK-treated",
                      ident.2 = "Targets only")

result$gene <- rownames(result)
result <- result[,c(6,1:5)]

fwrite(result, "results/celllinepanel_combined/targets/deg/targets_deg_combined.txt", row.names = F, quote = F, sep = "\t")


# volcano plot of expanded NK-treatment DEG (Figure 2B)
plotdata <- fread("results/celllinepanel_combined/targets/deg/targets_deg_combined.txt", data.table = F)

CUTOFF_FDR = 0.05
CUTOFF_PVAL = 0.05

labelgenes <- plotdata %>%
  top_n(-20, wt = p_val) %>%
  dplyr::select(gene) %>%
  tibble::deframe()

core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")


ggplot(plotdata, aes(x=avg_log2FC, y=-log10(p_val_adj), label = gene)) +
  geom_point(data = plotdata[plotdata$p_val_adj>=0.05,], color = "grey80") +
  geom_point(data = plotdata[plotdata$p_val_adj<0.05,], aes(color = ifelse(avg_log2FC > 0, "up", "down"))) +
  geom_point(data = plotdata[plotdata$gene %in% core,], color = brewer.pal(11, "RdBu")[4]) +
  geom_text_repel(aes(label=ifelse(p_val_adj <= CUTOFF_FDR & p_val <= CUTOFF_PVAL & gene %in% labelgenes, as.character(gene),'')),
                  point.padding = 0.2,
                  max.overlaps = 100,
                  fontface = "italic") +
  ylab("Adjusted p value (-log10)") +
  xlab("Average fold change (log2)") +
  scale_color_manual(values = c(brewer.pal(11, "RdBu")[10], brewer.pal(11, "RdBu")[2])) +
  guides(color = "none") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0, face = "plain")) +
  scale_x_continuous(limits = c(-0.7, 0.7), breaks = c(-0.5, 0, 0.5))

ggsave("results/celllinepanel_combined/targets/deg/volcanoplot_all.pdf", height = 4, width = 5.5)








