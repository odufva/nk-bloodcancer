
# UMAP of Chapuy et al. DLBCL data (Figure 5I)

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
library(gridExtra)

# load data
load("data/GSE98588_fm.Rdata")

# subset to samples with gene expression data
gexp <- fm[grepl("N:GEXP", rownames(fm)), !is.na(fm["N:GEXP:KRAS",])]
fm_gexp <- fm[, !is.na(fm["N:GEXP:KRAS",])]

# take 15% most variable genes
gexp_15pct <- gexp[rowVars(as.matrix(gexp)) > quantile(rowVars(as.matrix(gexp)), 0.85),]

# simplify subtype annotations for plotting
subtype <- rep("NA", length(fm_gexp))
subtype[fm_gexp["B:SAMP:COO_byGEP_ABC",]==1] <- "ABC"
subtype[fm_gexp["B:SAMP:COO_byGEP_GCB",]==1] <- "GCB"
subtype[fm_gexp["B:SAMP:COO_byGEP_Unclassified",]==1] <- "Unclassified"

subtype <- factor(subtype, levels = c("ABC", "GCB", "Unclassified", "NA"))

cluster <- factor(as.character(fm_gexp["N:SAMP:CLUSTER",]))

annot <- data.frame(subtype = subtype, cluster = cluster)

colors_subtype <- c(colorRampPalette(brewer.pal(4, "Set1"))(3), "grey70")
names(colors_subtype) <- c("ABC", "GCB", "Unclassified", "NA")

# scale gene expression
gexp_df <- apply(gexp, 1, scale)

# mutations
gnab_df <- t(fm_gexp[grepl("B:GNAB", rownames(fm_gexp)),])
class(gnab_df) <- "character"

# clinical
clin_df <- t(fm_gexp[grepl("B:CLIN", rownames(fm_gexp)),])

# scaled scores
scores_df <- apply(fm_gexp[grepl("Score$", rownames(fm_gexp)),], 1, scale)


# UMAP

set.seed(42)

umap <- umap(t(gexp_15pct), spread = 5)

umap_df <- as.data.frame(umap$layout)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# data frame for plotting
umap_df <- cbind(annot, umap_df, gexp_df, gnab_df, clin_df, scores_df)

point_size = 1

# plot subtype
p1 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = subtype)) +
  geom_point(size = point_size) +
  geom_point(data = umap_df[umap_df$subtype!="NA",], size = point_size) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 3, override.aes = list(size = 3))) +
  scale_color_manual(values = brewer.pal(9, "Set1")[c(4,5,3)]) +
  labs(color = "") +
  ggtitle("Subtype") +
  xlab("UMAP 1") +
  ylab("UMAP 2")


# plot GNA13 mutations
p2 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = `B:GNAB:GNA13`)) +
  geom_point(size = point_size) +
  geom_point(data = umap_df[umap_df$`B:GNAB:GNA13`==1,], size = point_size) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(3, "mm"),
        plot.title = element_text(hjust = 0.5, face = "italic")) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 3, override.aes = list(size = 3))) +
  scale_color_manual(values = c("grey70", "black"), labels = c("", "GNA13 mut")) +
  labs(color = "") +
  ggtitle("GNA13") +
  xlab("UMAP 1") +
  ylab("UMAP 2")

# plot HLA I score
p3 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = `N:SAMP:HLAIScore`)) +
  geom_point(size = point_size) +
  geom_point(data = umap_df[umap_df$`N:SAMP:HLAIScore`<(-2),], size = point_size, color = brewer.pal(11, "RdBu")[11]) +
  geom_point(data = umap_df[umap_df$`N:SAMP:HLAIScore`>1,], size = point_size) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(3, "mm"),
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  scale_color_distiller(palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                        type = "div", limits = c(-2, 2)) +
  labs(color = "Scaled expression") +
  ggtitle("HLA I score") +
  xlab("UMAP 1") +
  ylab("UMAP 2")



p1 + p2 + p3

ggsave("Chapuy_DLBCL_UMAP.pdf", height = 4.25, width = 8)


