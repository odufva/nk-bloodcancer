
# Plot CROP-seq UMAPs with core NK response module scores (Figure 6C)

# load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(reshape)

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

cols <- fread("perturbation_colors.txt", data.table = F)
cols_vector <- as.character(cols$`Hex Colors`)
names(cols_vector) <- cols$Perturbation

# load Seurat objects and add module scores

nkresponse <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")

k562_prtb_nonk <- readRDS("results/k562/mixscape/singlet/nonk_ko_nt.rds")
k562_prtb_nonk <- AddModuleScore(k562_prtb_nonk, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

k562_prtb_1_16 <- readRDS("results/k562/mixscape/singlet/nk_1_16_ko_nt.rds")
k562_prtb_1_16 <- AddModuleScore(k562_prtb_1_16, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

sudhl4_prtb_nonk <- readRDS("results/sudhl4/mixscape/singlet/nonk_ko_nt.rds")
sudhl4_prtb_nonk <- AddModuleScore(sudhl4_prtb_nonk, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

sudhl4_prtb_1_16 <- readRDS("results/sudhl4/mixscape/singlet/nk_1_16_ko_nt.rds")
sudhl4_prtb_1_16 <- AddModuleScore(sudhl4_prtb_1_16, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

sudhl4_prtb_1_4 <- readRDS("results/sudhl4/mixscape/singlet/nk_1_4_ko_nt.rds")
sudhl4_prtb_1_4 <- AddModuleScore(sudhl4_prtb_1_4, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

nalm6_prtb_nonk <- readRDS("results/nalm6/mixscape/singlet/nonk_ko_nt.rds")
nalm6_prtb_nonk <- AddModuleScore(nalm6_prtb_nonk, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

nalm6_prtb_1_16 <- readRDS("results/nalm6/mixscape/singlet/nk_1_16_ko_nt.rds")
nalm6_prtb_1_16 <- AddModuleScore(nalm6_prtb_1_16, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

nalm6_prtb_1_4 <- readRDS("results/nalm6/mixscape/singlet/nk_1_4_ko_nt.rds")
nalm6_prtb_1_4 <- AddModuleScore(nalm6_prtb_1_4, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

mm1s_prtb_nonk <- readRDS("results/mm1s/mixscape/singlet/nonk_ko_nt.rds")
mm1s_prtb_nonk <- AddModuleScore(mm1s_prtb_nonk, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

mm1s_prtb_1_16 <- readRDS("results/mm1s/mixscape/singlet/nk_1_16_ko_nt.rds")
mm1s_prtb_1_16 <- AddModuleScore(mm1s_prtb_1_16, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

mm1s_prtb_1_4 <- readRDS("results/mm1s/mixscape/singlet/nk_1_4_ko_nt.rds")
mm1s_prtb_1_4 <- AddModuleScore(mm1s_prtb_1_4, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

lp1_prtb_nonk <- readRDS("../NK_resistance Heme CollabPaper/Analysis/CROP_seq_LP1_b/R/results/lp1/mixscape/singlet/nonk_ko_nt.rds")
lp1_prtb_nonk <- AddModuleScore(lp1_prtb_nonk, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

lp1_prtb_1_16 <- readRDS("../NK_resistance Heme CollabPaper/Analysis/CROP_seq_LP1_b/R/results/lp1/mixscape/singlet/nk_1_16_ko_nt.rds")
lp1_prtb_1_16 <- AddModuleScore(lp1_prtb_1_16, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

lp1_prtb_1_4 <- readRDS("../NK_resistance Heme CollabPaper/Analysis/CROP_seq_LP1_b/R/results/lp1/mixscape/singlet/nk_1_4_ko_nt.rds")
lp1_prtb_1_4 <- AddModuleScore(lp1_prtb_1_4, features = list(nkresponse), name = c("nk_response_score"), assay = "RNA")

# order perturbations
k562_perturbations <- c("IFNGR2", "JAK1", "JAK2", "STAT1", "PTPN2", "GFI1B", "Control")
k562_prtb_nonk$gene <- gsub("NT", "Control", k562_prtb_nonk$gene)
k562_prtb_nonk$gene <- factor(k562_prtb_nonk$gene, levels = k562_perturbations)
k562_prtb_1_16$gene <- gsub("NT", "Control", k562_prtb_1_16$gene)
k562_prtb_1_16$gene <- factor(k562_prtb_1_16$gene, levels = k562_perturbations)

sudhl4_perturbations <- c("IFNGR2", "JAK1", "JAK2", "STAT1", "RFXAP", "STAG2", "PTEN", "METTL17", "PCGF1", "CHD7", "RUNX1", "YTHDF2",  "BID", "FADD", "CASP8", "Control")
sudhl4_prtb_nonk$gene <- gsub("NT", "Control", sudhl4_prtb_nonk$gene)
sudhl4_prtb_nonk$gene <- factor(sudhl4_prtb_nonk$gene, levels = sudhl4_perturbations)
sudhl4_prtb_1_16$gene <- gsub("NT", "Control", sudhl4_prtb_1_16$gene)
sudhl4_prtb_1_16$gene <- factor(sudhl4_prtb_1_16$gene, levels = sudhl4_perturbations)
sudhl4_prtb_1_4$gene <- gsub("NT", "Control", sudhl4_prtb_1_4$gene)
sudhl4_prtb_1_4$gene <- factor(sudhl4_prtb_1_4$gene, levels = sudhl4_perturbations)

nalm6_perturbations <- c("RFXAP", "CHD7", "KLF16", "KIAA0922", "CMIP", "Control")
nalm6_prtb_nonk$gene <- gsub("NT", "Control", nalm6_prtb_nonk$gene)
nalm6_prtb_nonk$gene <- factor(nalm6_prtb_nonk$gene, levels = nalm6_perturbations)
nalm6_prtb_1_16$gene <- gsub("NT", "Control", nalm6_prtb_1_16$gene)
nalm6_prtb_1_16$gene <- factor(nalm6_prtb_1_16$gene, levels = nalm6_perturbations)
nalm6_prtb_1_4$gene <- gsub("NT", "Control", nalm6_prtb_1_4$gene)
nalm6_prtb_1_4$gene <- factor(nalm6_prtb_1_4$gene, levels = nalm6_perturbations)

mm1s_perturbations <- c("IFNGR2", "JAK1", "JAK2", "STAT1", "NLRC5", "RFXAP", "NFKBIA", "NFKBIB", "TRAF2", "PTEN", "GNA13", "ARHGAP1", "PCGF5", "Control")
mm1s_prtb_nonk$gene <- gsub("NT", "Control", mm1s_prtb_nonk$gene)
mm1s_prtb_nonk$gene <- factor(mm1s_prtb_nonk$gene, levels = mm1s_perturbations)
mm1s_prtb_1_16$gene <- gsub("NT", "Control", mm1s_prtb_1_16$gene)
mm1s_prtb_1_16$gene <- factor(mm1s_prtb_1_16$gene, levels = mm1s_perturbations)
mm1s_prtb_1_4$gene <- gsub("NT", "Control", mm1s_prtb_1_4$gene)
mm1s_prtb_1_4$gene <- factor(mm1s_prtb_1_4$gene, levels = mm1s_perturbations)

lp1_perturbations <- c("IFNGR2", "JAK1", "JAK2","STAT1", "NLRC5", "RFXAP", "TRAF2", "PTEN", "GSK3B", "MYB", "MSI2", "CHD7", "Control")
lp1_prtb_nonk$gene <- gsub("NT", "Control", lp1_prtb_nonk$gene)
lp1_prtb_nonk$gene <- factor(lp1_prtb_nonk$gene, levels = lp1_perturbations)
lp1_prtb_1_16$gene <- gsub("NT", "Control", lp1_prtb_1_16$gene)
lp1_prtb_1_16$gene <- factor(lp1_prtb_1_16$gene, levels = lp1_perturbations)
lp1_prtb_1_4$gene <- gsub("NT", "Control", lp1_prtb_1_4$gene)
lp1_prtb_1_4$gene <- factor(lp1_prtb_1_4$gene, levels = lp1_perturbations)


# mixscape UMAPs

p_k562_prtb_nonk <- DimPlot(k562_prtb_nonk, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("K562") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_k562_prtb_1_16 <- DimPlot(k562_prtb_1_16, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[k562_perturbations[k562_perturbations %in% unique(k562_prtb_1_16$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("K562") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")


p_sudhl4_prtb_nonk <- DimPlot(sudhl4_prtb_nonk, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[sudhl4_perturbations[sudhl4_perturbations %in% unique(sudhl4_prtb_nonk$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("SUDHL4") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_sudhl4_prtb_1_16 <- DimPlot(sudhl4_prtb_1_16, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[sudhl4_perturbations[sudhl4_perturbations %in% unique(sudhl4_prtb_1_16$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("SUDHL4") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_sudhl4_prtb_1_4 <- DimPlot(sudhl4_prtb_1_4, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[sudhl4_perturbations[sudhl4_perturbations %in% unique(sudhl4_prtb_1_4$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("SUDHL4") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")


p_nalm6_prtb_nonk <- DimPlot(nalm6_prtb_nonk, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[nalm6_perturbations[nalm6_perturbations %in% unique(nalm6_prtb_nonk$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("NALM6") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_nalm6_prtb_1_16 <- DimPlot(nalm6_prtb_1_16, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[nalm6_perturbations[nalm6_perturbations %in% unique(nalm6_prtb_1_16$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("NALM6") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_nalm6_prtb_1_4 <- DimPlot(nalm6_prtb_1_4, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[nalm6_perturbations[nalm6_perturbations %in% unique(nalm6_prtb_1_4$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("NALM6") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")


# MM1S

p_mm1s_prtb_nonk <- DimPlot(mm1s_prtb_nonk, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[mm1s_perturbations[mm1s_perturbations %in% unique(mm1s_prtb_nonk$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("MM1S") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_mm1s_prtb_1_16 <- DimPlot(mm1s_prtb_1_16, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[mm1s_perturbations[mm1s_perturbations %in% unique(mm1s_prtb_1_16$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("MM1S") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_mm1s_prtb_1_4 <- DimPlot(mm1s_prtb_1_4, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[mm1s_perturbations[mm1s_perturbations %in% unique(mm1s_prtb_1_4$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("MM1S") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")


# LP1

p_lp1_prtb_nonk <- DimPlot(lp1_prtb_nonk, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[lp1_perturbations[lp1_perturbations %in% unique(lp1_prtb_nonk$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("LP1") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_lp1_prtb_1_16 <- DimPlot(lp1_prtb_1_16, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[lp1_perturbations[lp1_perturbations %in% unique(lp1_prtb_1_16$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("LP1") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")

p_lp1_prtb_1_4 <- DimPlot(lp1_prtb_1_4, group.by = "gene", reduction = "ldaumap", raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_manual(values = cols_vector[lp1_perturbations[lp1_perturbations %in% unique(lp1_prtb_1_4$gene)]]) &
  guides(color = guide_legend(nrow = 7, override.aes = list(size = 3))) &
  ggtitle("LP1") &
  xlab("LDAUMAP 1") &
  ylab("LDAUMAP 2")


# features 
DefaultAssay(k562_prtb_1_16) <- "RNA"

p_k562_1_16_feat1 <- FeaturePlot(k562_prtb_1_16, features = "nk_response_score1", reduction = "ldaumap", min.cutoff = -0.2, max.cutoff = 0.4, raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_k562_1_16_feat2 <- FeaturePlot(k562_prtb_1_16, features = "B2M", reduction = "ldaumap", min.cutoff = 1.5, max.cutoff = 3.5) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

DefaultAssay(sudhl4_prtb_1_16) <- "RNA"

p_sudhl4_1_16_feat1 <- FeaturePlot(sudhl4_prtb_1_16, features = "nk_response_score1", reduction = "ldaumap", max.cutoff = 0.4, min.cutoff = -0.2, raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_sudhl4_1_16_feat2 <- FeaturePlot(sudhl4_prtb_1_16, features = "NFKBIA", reduction = "ldaumap", max.cutoff = 1, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

DefaultAssay(sudhl4_prtb_1_4) <- "RNA"

p_sudhl4_1_4_feat1 <- FeaturePlot(sudhl4_prtb_1_4, features = "nk_response_score1", reduction = "ldaumap", max.cutoff = 0.4, min.cutoff = -0.2, raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_sudhl4_1_4_feat2 <- FeaturePlot(sudhl4_prtb_1_4, features = "NFKBIA", reduction = "ldaumap", max.cutoff = 1, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")


DefaultAssay(nalm6_prtb_1_16) <- "RNA"

p_nalm6_1_16_feat1 <- FeaturePlot(nalm6_prtb_1_16, features = "nk_response_score1", reduction = "ldaumap", max.cutoff = 0.3, min.cutoff = -0.1, raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_nalm6_1_16_feat2 <- FeaturePlot(nalm6_prtb_1_16, features = "DNTT", reduction = "ldaumap", max.cutoff = 2.5, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")


DefaultAssay(nalm6_prtb_nonk) <- "RNA"

p_nalm6_nonk_feat1 <- FeaturePlot(nalm6_prtb_nonk, features = "nk_response_score1", reduction = "ldaumap", max.cutoff = 0.3, min.cutoff = -0.1, raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_nalm6_nonk_feat2 <- FeaturePlot(nalm6_prtb_nonk, features = "DNTT", reduction = "ldaumap", max.cutoff = 2.5, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")



DefaultAssay(mm1s_prtb_1_16) <- "RNA"

p_mm1s_1_16_feat1 <- FeaturePlot(mm1s_prtb_1_16, features = "nk_response_score1", reduction = "ldaumap", max.cutoff = 0.6, min.cutoff = -0.2, raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu", breaks = c(-0.2, 0.6), labels = c("Low", "High"))# +
guides(color = "none")

p_mm1s_1_16_feat2 <- FeaturePlot(mm1s_prtb_1_16, features = "FAS", reduction = "ldaumap", max.cutoff = 0.75, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu")


DefaultAssay(mm1s_prtb_1_4) <- "RNA"

p_mm1s_1_4_feat1 <- FeaturePlot(mm1s_prtb_1_4, features = "HLA-B", reduction = "ldaumap", max.cutoff = 3.5, min.cutoff = 1) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_mm1s_1_4_feat2 <- FeaturePlot(mm1s_prtb_1_4, features = "FAS", reduction = "ldaumap", max.cutoff = 0.75, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu")



# LP1

DefaultAssay(lp1_prtb_1_16) <- "RNA"

p_lp1_1_16_feat1 <- FeaturePlot(lp1_prtb_1_16, features = "nk_response_score1", reduction = "ldaumap", max.cutoff = 0.6, min.cutoff = -0.2, raster = T, pt.size = 4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_lp1_1_16_feat2 <- FeaturePlot(lp1_prtb_1_16, features = "CFLAR", reduction = "ldaumap", max.cutoff = 2, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")


DefaultAssay(lp1_prtb_1_4) <- "RNA"

p_lp1_1_4_feat1 <- FeaturePlot(lp1_prtb_1_4, features = "HLA-B", reduction = "ldaumap", max.cutoff = 3, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")

p_lp1_1_4_feat2 <- FeaturePlot(lp1_prtb_1_4, features = "FAS", reduction = "ldaumap", max.cutoff = 0.75, min.cutoff = 0) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        plot.title = element_text(face = "italic", hjust = 0.04, margin = margin(t = 10, b = -17.5)),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,-5,-5)) +
  scale_color_distiller(palette = "RdBu") +
  guides(color = "none")


((p_k562_prtb_1_16 / p_k562_1_16_feat1) |
    (p_sudhl4_prtb_1_16 / p_sudhl4_1_16_feat1) |
    (p_mm1s_prtb_1_16 / p_mm1s_1_16_feat1) |
    (p_lp1_prtb_1_16 / p_lp1_1_16_feat1) |
    (p_nalm6_prtb_1_16 / p_nalm6_1_16_feat1)) & theme(axis.title = element_blank())

ggsave("results/visualization/cropseq_umaps_modulescore_onerow.pdf", height = 6, width = 12)

