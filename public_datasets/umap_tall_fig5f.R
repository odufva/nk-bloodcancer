
# UMAP of Liu et al T-ALL data (Figure 5F)

# load libraries
library(data.table)
library(parallel)
library(readxl)
library(dplyr)
library(GSVA)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggfortify)
library(Rtsne)
library(matrixStats)
library(cowplot)
library(umap)
library(viridis)
library(patchwork)
library(scales)

# load data
load("data/TALL_Liu_fm.Rdata")

gexp <- fm[grepl("N:GEXP", rownames(fm)), !is.na(fm["N:GEXP:KRAS",])]

gexp_15pct <- gexp[rowVars(as.matrix(gexp)) > quantile(rowVars(as.matrix(gexp)), 0.85),]


maturation <- rep("NA", length(fm))
maturation[fm["B:CLIN:Maturation_stage_Cortical",]==1] <- "Cortical"
maturation[fm["B:CLIN:Maturation_stage_Pre-cortical",]==1] <- "Pre-cortical"
maturation[fm["B:CLIN:Maturation_stage_Post-cortical",]==1] <- "Post-cortical"

subtypes <- rep("NA", length(fm))
subtypes[fm["B:CLIN:group_TLX1",]==1] <- "TLX1"
subtypes[fm["B:CLIN:group_TAL2",]==1] <- "TAL2"
subtypes[fm["B:CLIN:group_HOXA",]==1] <- "HOXA"
subtypes[fm["B:CLIN:group_TAL1",]==1] <- "TAL1"
subtypes[fm["B:CLIN:group_NKX2_1",]==1] <- "NKX2.1"
subtypes[fm["B:CLIN:group_TLX3",]==1] <- "TLX3"
subtypes[fm["B:CLIN:group_LMO2_LYL1",]==1] <- "LMO2/LYL1"
subtypes[fm["B:CLIN:group_LMO1/2",]==1] <- "LMO1/2"

subtypes <- factor(subtypes, levels = c("TAL1", "TAL2", "LMO1/2", "NKX2.1", "TLX1", "TLX3", "HOXA", "LMO2/LYL1", "NA"))
annot <- data.frame(subtype = subtypes)

colors_subtype <- c(colorRampPalette(brewer.pal(11, "Spectral")[c(1:4, 8:11)])(8), "grey70")
names(colors_subtype) <- c("TAL1", "TAL2", "LMO1/2", "NKX2.1", "TLX1", "TLX3", "HOXA", "LMO2/LYL1", "NA")

gexp_df <- apply(gexp, 1, scale)


# UMAP

set.seed(42)
umap <- umap(t(gexp))
umap <- umap(t(gexp_15pct))

umap_df <- as.data.frame(umap$layout)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# data frame for plotting
umap_df <- cbind(annot, umap_df, gexp_df)

point_size = 1

p1 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = subtype)) +
  geom_point(size = point_size) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 4, override.aes = list(size = 3))) +
  scale_color_manual(values = colors_subtype) +
  labs(color = "") +
  ggtitle("Subtype") +
  xlab("UMAP 1") +
  ylab("UMAP 2")
  
p2 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = `N:GEXP:PVR`)) +
  geom_point(size = point_size) +
  geom_point(data = umap_df[umap_df$`N:GEXP:PVR`>1,], size = point_size) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(3, "mm"),
        plot.title = element_text(hjust = 0.5, face = "italic")) +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  scale_color_distiller(palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                        type = "div", limits = max(abs(umap_df$`N:GEXP:PVR`)) * c(-1, 1)) +
  labs(color = "Scaled expression") +
  ggtitle("PVR") +
  xlab("UMAP 1") +
  ylab("UMAP 2")

p3 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = `N:GEXP:ULBP1`)) +
  geom_point(size = point_size) +
  geom_point(data = umap_df[umap_df$`N:GEXP:ULBP1`>1,], size = point_size) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(3, "mm"),
        plot.title = element_text(hjust = 0.5, face = "italic")) +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  scale_color_distiller(palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                        type = "div", limits = max(abs(umap_df$`N:GEXP:ULBP1`)) * c(-1, 1)) +
  labs(color = "Scaled expression") +
  ggtitle("ULBP1") +
  xlab("UMAP 1") +
  ylab("UMAP 2")

p1 + p2 + p3


ggsave("TALL_Liu_UMAP.pdf", height = 4.5, width = 8)
