
# Map media transfer NK cell data to original 26 cell line panel UMAP

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
library(cluster)
library(ggridges)
library(awtools)
library(patchwork)

dir.create("results/media_transfer/integrated")


mediatransfer_seurat <- readRDS("results/media_transfer/hto_singlets/nk/mediatransfer_nk_seurat.rds")
nk1_seurat <- readRDS("results/celllinepanel/hto_singlets/nk/batchcorrected/celllinepanel_nk1_seurat.rds")

nk1_seurat$cluster <- factor(nk1_seurat$seurat_clusters,
                             levels = c(0:4),
                             labels = c("Resting (0)",
                                        "Adaptive (1)",
                                        "Activated (2)",
                                        "Type I IFN (3)",
                                        "Cytokine (4)"))




nk1_seurat <- ScaleData(nk1_seurat, assay = 'RNA')

# Run UMAP on cell line panel data with returbn.model = T to be able to map other data to it
nk1_seurat <- RunUMAP(nk1_seurat, reduction = "pca", dims = 1:20, return.model = TRUE)


# Map K562, SUDHL4, NALM6 data to cell line panel reference
anchors <- FindTransferAnchors(
  reference = nk1_seurat,
  query = mediatransfer_seurat,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:20,
  k.anchor = 20,
  k.score = 60,
  k.filter = 400
)

mediatransfer_seurat <- MapQuery(
  anchorset = anchors,
  query = mediatransfer_seurat,
  reference = nk1_seurat,
  refdata = list(
    cluster = "cluster"),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

mediatransfer_seurat$predicted.cluster <- factor(mediatransfer_seurat$predicted.cluster, levels = c("Resting (0)",
                                                                                          "Adaptive (1)",
                                                                                          "Activated (2)",
                                                                                          "Type I IFN (3)",
                                                                                          "Cytokine (4)"))

saveRDS(mediatransfer_seurat, "results/media_transfer/integrated/mediatransfer_nk_integrated_seurat.rds")



# UMAPs
DimPlot(mediatransfer_seurat, reduction = "ref.umap",  group.by = "predicted.cluster", label = T, repel = T) &
  theme_bw(base_size = 12) &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  scale_color_manual(values = a_palette[c(1,2,6,8,5)]) &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

sampleorder <- c("697", "KASUMI2", "RCHACV", "LP1", "JJN3", "PL21",
                 "SKM1", "MM1S", "MONOMAC1", "MOLM13", "ALLSIL", "MOLM14", "GRANTA519",
                 "NUDHL1", "SUPT11", "NALM6", "AMO1", "L363", "DND41", "THP1",
                 "SUDHL4", "OCIM1", "JURKAT", "RI1", "K562", "GDM1")

nkdonors <- c("NK1", "NK2", "NK3")

nkdonors_all <- c("NK1", "NK2", "NK3", "NK4", "NK5", "NK6")

grid <- expand.grid(sampleorder, nkdonors) %>% arrange(Var1)
plotorder < -apply(grid, 1, paste, collapse="-")
plotorder <- c("NK1", "NK2", "NK3", plotorder)

combined_seurat_plot <- combined_seurat
combined_seurat_plot$hash.ID <- gsub("NK-expanded", "NK1", combined_seurat_plot$hash.ID)
combined_seurat_plot <- subset(combined_seurat_plot, hash.ID %in% plotorder)
combined_seurat_plot$hash.ID <- factor(combined_seurat_plot$hash.ID, levels = plotorder)


mediatransfer_seurat_plot <- mediatransfer_seurat
mediatransfer_seurat_plot$hash.ID <- gsub("\\-1|\\-2", "", mediatransfer_seurat_plot$hash.ID)
mediatransfer_seurat_plot$hash.ID <- factor(mediatransfer_seurat_plot$hash.ID, levels = c("NK-only", "K562-NK", "GDM1-NK", "697-NK",
                                                                                          "NK-medium-from-K562-NK", "NK-medium-from-GDM1-NK", "NK-medium-from-697-NK",
                                                                                          "NK-medium-from-K562", "NK-medium-from-GDM1", "NK-medium-from-697"))


# Point UMAPs of mapped data

DimPlot(mediatransfer_seurat_plot, reduction = "ref.umap", group.by = "predicted.cluster", split.by = "hash.ID", ncol = 12) &
  theme_bw(base_size = 12) &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  scale_color_manual(values = a_palette[c(1,2,6,8,5)]) &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_blank(),
        strip.placement = "outside")

ggsave("results/media_transfer/integrated/umap_celllines_order.pdf", height = 2, width = 22)


# Density UMAP of mapped data
plotdata <- FetchData(mediatransfer_seurat, vars = c("refUMAP_1", "refUMAP_2", "hash.ID", "predicted.cluster"))
plotdata_selected <- plotdata %>% 
  mutate(treatment = ifelse(hash.ID %in% c("697-NK", "GDM1-NK", "K562-NK"), "Co-culture",
                            ifelse(hash.ID %in% c("NK-medium-from-697-NK", "NK-medium-from-GDM1-NK", "NK-medium-from-K562-NK"), "Medium from co-culture",
                                   ifelse(hash.ID %in% c("NK-medium-from-697", "NK-medium-from-GDM1", "NK-medium-from-K562"), "Medium from cell line", "Untreated")))) %>% 
  filter(treatment != "Untreated") %>% 
  mutate(cell_line = ifelse(grepl("K562", hash.ID), "K562",
                            ifelse(grepl("GDM1", hash.ID), "GDM1",
                                   ifelse(grepl("697", hash.ID), "697", "")))) %>% 
  mutate(treatment = factor(treatment, levels = c("Co-culture",  "Medium from cell line", "Medium from co-culture")),
         cell_line = factor(cell_line, levels = c("697", "K562", "GDM1")))
  
plotdata_background <- FetchData(nk1_seurat, vars = c("UMAP_1", "UMAP_2"))

p <- ggplot(plotdata_selected, aes(x = refUMAP_1, y = refUMAP_2)) +
  theme_bw() +
  stat_density_2d(data = plotdata_background, aes(x = UMAP_1, y = UMAP_2, fill = ..density..), contour = F, geom = "raster") +
  stat_density_2d(aes(x = refUMAP_1, y = refUMAP_2, color = ..level..), contour_var = "ndensity", h = 3) +
  scale_fill_gradientn(colors = brewer.pal(9, "Greys")[1:6]) +
  scale_color_viridis_c(name = "Density", option = "magma", begin = 0.1, end = 1) +
  scale_x_continuous(limits = c(-7.5, 7.5)) +
  scale_y_continuous(limits = c(-7, 5.5)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  guides(color = "none",
         fill = "none") +
  xlab("") +
  ylab("") +
  facet_grid(cols = vars(treatment), rows = vars(cell_line), switch = "y")

plotdata_nkonly <-  plotdata %>% 
  mutate(treatment = ifelse(hash.ID %in% c("697-NK", "GDM1-NK", "K562-NK"), "Co-culture",
                            ifelse(hash.ID %in% c("NK-medium-from-697-NK", "NK-medium-from-GDM1-NK", "NK-medium-from-K562-NK"), "Medium from co-culture",
                                   ifelse(hash.ID %in% c("NK-medium-from-697", "NK-medium-from-GDM1", "NK-medium-from-K562"), "Medium from cell line", "Untreated")))) %>% 
  filter(treatment == "Untreated")

p_nkonly <- ggplot(plotdata_nkonly, aes(x = refUMAP_1, y = refUMAP_2)) +
  theme_bw() +
  stat_density_2d(data = plotdata_background, aes(x = UMAP_1, y = UMAP_2, fill = ..density..), contour = F, geom = "raster") +
  stat_density_2d(aes(x = refUMAP_1, y = refUMAP_2, color = ..level..), contour_var = "ndensity", h = 3) +
  scale_fill_gradientn(colors = brewer.pal(9, "Greys")[1:6]) +
  scale_color_viridis_c(name = "Density", option = "magma", begin = 0.1, end = 1) +
  scale_x_continuous(limits = c(-7.5, 7.5)) +
  scale_y_continuous(limits = c(-7, 5.5)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5),
        strip.background = element_blank(),
        ) +
  guides(color = "none",
         fill = "none") +
  xlab("") +
  ylab("") +
  ggtitle("Untreated NK cells")

p_nk1 <- DimPlot(nk1_seurat, group.by = "cluster", label = F, repel = T) &
  theme_bw(base_size = 12) &
  #scale_color_manual(values = rev(getPalette(41))) +
  ylab("") &
  scale_color_manual(values = a_palette[c(1,2,6,8,5)]) &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

layout <- "
AACCCCC
BBCCCCC
"

p_nk1 + p_nkonly + p + plot_layout(design = layout, guides = "collect") & theme(legend.position = "top") 


ggsave("results/media_transfer/integrated/umap_celllines_density_grid.pdf", height = 6, width = 8.5)

## -----------------------------------------

# Plot UMAPs without density

mediatransfer_seurat <- readRDS("results/media_transfer/integrated/mediatransfer_nk_integrated_seurat.rds")


plotdata <- FetchData(mediatransfer_seurat, vars = c("UMAP_1", "UMAP_2", "hash.ID", "predicted.cluster"))
plotdata_selected <- plotdata %>%
  mutate(condition = gsub(".*-", "", hash.ID),
         cell_line = gsub("-.*", "", hash.ID)) %>% 
  mutate(condition = factor(timepoint, levels = timepoints),
         cell_line = factor(cell_line, levels = sampleorder))


p3 <- ggplot(plotdata, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual(values = a_palette[c(1,5,3,6,8)]) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  guides(color = "none",
         fill = "none") +
  xlab("") +
  ylab("") +
  facet_grid(rows = vars(cell_line), cols = vars(timepoint), switch = "y")

layout <- "
ACCC
BCCC
"

p1 + p2 + p3 + plot_layout(design = layout, guides = "collect") & theme(legend.position = "top") 
ggsave("results/timepoints/hto_singlets/nk/umap_celllines_order_labels_grid_clusters.pdf", width = 7.5, height = 5)


## -----------------------------

# cluster fractions barplot

Idents(mediatransfer_seurat) <- "hash.ID"

cluster_freq <- data.frame(table(Idents(mediatransfer_seurat), mediatransfer_seurat$predicted.cluster))
cluster_freq <- cluster_freq %>%
  group_by(Var1, Var2) %>%
  summarize(count = sum(Freq)) %>%
  mutate(freq = count / sum(count))

sample_order <-  c("NK-only", "697-NK", "K562-NK", "GDM1-NK",
                   "NK-medium-from-697-NK", "NK-medium-from-K562-NK", "NK-medium-from-GDM1-NK",
                   "NK-medium-from-697", "NK-medium-from-K562", "NK-medium-from-GDM1")

cluster_freq_pert <- cluster_freq %>%
  mutate(condition = factor(Var1, levels = sample_order)) %>% 
  mutate(freq = 100*freq)


# barplot of mean cluster fractions
p_cluster_barplot <- ggplot(cluster_freq_pert, aes(x = condition, y = freq, fill = Var2)) +
  geom_col() +
  theme_cowplot() +
  scale_fill_manual(values = a_palette[c(1,2,6,8,5)]) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") +
  ylab("% of NK cells in cluster") +
  labs(fill = "Cluster") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 10)) +
  guides(fill = "none")

