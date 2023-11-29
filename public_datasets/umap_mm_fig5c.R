
# UMAP of CoMMpass data

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

# load data
load("data/MM_COMPASS_FM.Rdata")

# subset to samples with gene expression data
gexp <- fm[grepl("N:GEXP", rownames(fm)), !is.na(fm["N:GEXP:KRAS",])]
fm_gexp <- fm[, !is.na(fm["N:GEXP:KRAS",])]
  
# take 15% most variable genes
gexp_15pct <- gexp[rowVars(as.matrix(gexp)) > quantile(rowVars(as.matrix(gexp)), 0.85),]

# simplify subtype and genetics annotations for plotting
subtype <- rep("NA", length(fm_gexp))
subtype[fm_gexp["B:SAMP:cancermap_subtypes_CCND1_Ig",]==1] <- "CCND1"
subtype[fm_gexp["B:SAMP:cancermap_subtypes_WHSC1_FGFR3_Ig",]==1] <- "WHSC1"
subtype[fm_gexp["B:SAMP:cancermap_subtypes_Hyperdiploid",]==1] <- "Hyperdiploid"
subtype[fm_gexp["B:SAMP:cancermap_subtypes_Hyperdiploid_amp1q",]==1] <- "Hyperdiploid (amp 1q)"
subtype[fm_gexp["B:SAMP:cancermap_subtypes_MAF_Ig",]==1] <- "MAF"
subtype[fm_gexp["B:SAMP:cancermap_subtypes_TRAF3_Aberrated",]==1] <- "TRAF3"

genetics <- rep("NA", length(fm_gexp))
genetics[fm_gexp["B:CNVR:Hyperdiploid_Call",]==1] <- "Hyperdiploid"
genetics[fm_gexp["B:CNVR:SeqWGS_WHSC1_Ig_translocation",]==1] <- "WHSC1"
genetics[fm_gexp["B:CNVR:SeqWGS_CCND1_Ig_translocation",]==1] <- "CCND1"
genetics[fm_gexp["B:CNVR:SeqWGS_CCND2_Ig_translocation",]==1] <- "CCND2"
genetics[fm_gexp["B:CNVR:SeqWGS_CCND3_Ig_translocation",]==1] <- "CCND3"
genetics[fm_gexp["B:CNVR:SeqWGS_MAF_Ig_translocation",]==1 |
           fm_gexp["B:CNVR:SeqWGS_MAFA_Ig_translocation",]==1 |
           fm_gexp["B:CNVR:SeqWGS_MAFB_Ig_translocation",]==1] <- "MAF"

traf3 <- rep("NA", length(fm_gexp))
traf3[fm_gexp["B:GNAB:TRAF3",]==1] <- "TRAF3 mut"

# discretize CNVs
traf3_cnv <- rep("NA", length(fm_gexp))
traf3_cnv[fm_gexp["N:CNVR:TRAF3",] < -1.5] <- "del"   
traf3_cnv[fm_gexp["N:CNVR:TRAF3",] < -0.5 & fm_gexp["N:CNVR:TRAF3",] >= -1.5] <- "loss"   
traf3_cnv[fm_gexp["N:CNVR:TRAF3",] > 0.25  & fm_gexp["N:CNVR:TRAF3",] <= 0.75] <- "gain"   
traf3_cnv[fm_gexp["N:CNVR:TRAF3",] > 0.75] <- "amp"   

traf3_mut <- fm_gexp["B:GNAB:TRAF3",]==1

# combine TRAF3 mut and del/loss
traf3[traf3_mut] <- "TRAF3 mut"
traf3[traf3_cnv%in%c("del", "loss")] <- "TRAF3 del/loss"
traf3[traf3_mut & traf3_cnv%in%c("del", "loss")] <- "TRAF3 mut + del/loss"

genetics <- factor(genetics, levels = c("WHSC1", "CCND1", "CCND2", "CCND3", "MAF", "Hyperdiploid", "NA"))
traf3 <- factor(traf3, levels = c("TRAF3 mut", "TRAF3 del/loss", "TRAF3 mut + del/loss", "NA"))
subtype <- factor(subtype, levels = c("WHSC1", "CCND1", "MAF", "Hyperdiploid", "Hyperdiploid (amp 1q)", "TRAF3", "NA"))

annot <- data.frame(genetics = genetics, traf3 = traf3, subtype = subtype)

# define colors
colors_genetics <- c(colorRampPalette(brewer.pal(11, "BrBG")[2:10])(6), "grey70")
names(colors_genetics) <- c("WHSC1", "CCND1", "CCND2", "CCND3", "MAF", "Hyperdiploid", "NA")

colors_subtype <- c(colorRampPalette(brewer.pal(11, "PuOr")[2:10])(6), "grey70")
names(colors_subtype) <- c("WHSC1", "CCND1", "MAF", "Hyperdiploid", "Hyperdiploid (amp 1q)", "TRAF3", "NA")

colors_traf3 <- c(brewer.pal(9, "YlOrRd")[c(5,7,9)], "grey70")
names(colors_traf3) <- c("TRAF3 del/loss", "TRAF3 mut", "TRAF3 mut + del/loss", "NA")

# scale gene expression
gexp_df <- apply(gexp, 1, scale)

# mutations
gnab_df <- t(fm_gexp[grepl("B:GNAB", rownames(fm_gexp)),])

# scaled scores
scores_df <- apply(fm_gexp[grepl("Score$", rownames(fm_gexp)),], 1, scale)

# UMAP

set.seed(42)

umap <- umap(t(gexp_15pct), spread = 10)

umap_df <- as.data.frame(umap$layout)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# data frame for plotting
umap_df <- cbind(annot, umap_df, gexp_df, gnab_df, scores_df)

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
  scale_color_manual(values = colors_subtype) +
  labs(color = "") +
  ggtitle("Subtype") +
  xlab("UMAP 1") +
  ylab("UMAP 2")

# plot TRAF3 alterations
p2 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = traf3)) +
  geom_point(size = point_size) +
  geom_point(data = umap_df[umap_df$traf3=="TRAF3 del/loss",], size = point_size) +
  geom_point(data = umap_df[umap_df$traf3=="TRAF3 mut",], size = point_size) +
  geom_point(data = umap_df[umap_df$traf3=="TRAF3 mut + del/loss",], size = point_size) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "italic")) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 3, override.aes = list(size = 3))) +
  scale_color_manual(values = colors_traf3) +
  labs(color = "") +
  ggtitle("TRAF3") +
  xlab("UMAP 1") +
  ylab("UMAP 2")

# plot CFLAR expression
p3 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = `N:GEXP:CFLAR`)) +
  geom_point(size = point_size) +
  geom_point(data = umap_df[umap_df$`N:GEXP:CFLAR`>1,], size = point_size) +
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
                        type = "div", limits = max(abs(umap_df$`N:GEXP:CFLAR`)) * c(-1, 1)) +
  labs(color = "Scaled expression") +
  ggtitle("CFLAR") +
  xlab("UMAP 1") +
  ylab("UMAP 2")


p1 + p2 + p3

ggsave("CoMMpass_UMAP.pdf", height = 4.25, width = 8)








