
# Map NK cells from all donors to original (donor NK1) 26 cell line panel UMAP

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
library(ComplexHeatmap)
library(circlize)


dir.create("results/celllinepanel_integrated")

# load object with donors NK2-NK6
nk_new_seurat <- readRDS("results/celllinepanel_nk2nk3/celllinepanel_nk2nk3_seurat.rds")

# load object with donor NK1
nk1_seurat <- readRDS("results/celllinepanel_nk1/hto_singlets/nk/batchcorrected/celllinepanel_nk1_batchcorrected_seurat.rds")

nk1_seurat$cluster <- factor(nk1_seurat$seurat_clusters,
                                       levels = c(0:4),
                                       labels = c("Resting (0)",
                                                  "Adaptive (1)",
                                                  "Activated (2)",
                                                  "Type I IFN (3)",
                                                  "Cytokine (4)"))


combined_seurat <- merge(nk1_seurat, nk_new_seurat)


nk1_seurat <- ScaleData(nk1_seurat, assay = 'RNA')

# Run UMAP on cell line panel data with returbn.model = T to be able to map other data to it
nk1_seurat <- RunUMAP(nk1_seurat, reduction = "pca", dims = 1:20, return.model = TRUE)


# Map K562, SUDHL4, NALM6 data to cell line panel reference
anchors <- FindTransferAnchors(
  reference = nk1_seurat,
  query = combined_seurat,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:20
)

combined_seurat <- MapQuery(
  anchorset = anchors,
  query = combined_seurat,
  reference = nk1_seurat,
  refdata = list(
    cluster = "cluster"),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

combined_seurat$predicted.cluster <- factor(combined_seurat$predicted.cluster, levels = c("Resting (0)",
                                                                                          "Adaptive (1)",
                                                                                          "Activated (2)",
                                                                                          "Type I IFN (3)",
                                                                                          "Cytokine (4)"))

saveRDS(combined_seurat, "results/celllinepanel_integrated/celllinepanel_nk.rds")

combined_seurat <- readRDS("results/celllinepanel_integrated/celllinepanel_nk.rds")

table(combined_seurat$predicted.cluster)

# UMAPs
DimPlot(combined_seurat, reduction = "ref.umap",  group.by = "predicted.cluster", label = T, repel = T) &
  theme_bw(base_size = 12) &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  scale_color_manual(values = a_palette[c(1,2,6,8,5)]) &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

# DEGs of clusters
Idents(combined_seurat) <- "predicted.cluster"
all_deg <- FindAllMarkers(combined_seurat, test.use = "t", only.pos = T, return.thresh = 0.05, logfc.threshold = 0.1)
fwrite(all_deg, "results/celllinepanel_integrated/clusters_deg.txt", sep = "\t", quote = F, row.names = F)

clu0_markers <- c("KLRC1", "GZMA", "KLRB1", "NCAM1", "GZMK")
clu1_markers <- c("KLRC2", "GZMH", "LAG3", "HLA-DRA")
clu2_markers <- c("TNFRSF18", "TNFRSF9", "CRTAM", "ENTPD1", "HAVCR2", "TIGIT", "TNFSF10", "BCL2L11")
clu3_markers <- c("ISG15", "MX1" ,"MX2", "OAS1", "OAS2", "OAS3", "IRF7", "IRF9", "EIF2AK2")
clu4_markers <- c("IFNG", "CCL4", "CCL3", "TNF", "CSF2", "CD69")

markers <- c(clu0_markers, clu1_markers, clu2_markers, clu3_markers, clu4_markers)


# Cluster marker gene dot plot (Figure S1B)
DotPlot(combined_seurat,
        group.by = "predicted.cluster",
        features = markers,
        cols = "RdBu") +
  coord_flip() +
  scale_x_discrete(limits = rev, drop = F) +
  scale_y_discrete(position = "right", drop = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black")) +
  xlab("") +
  ylab("") +
  guides(color = guide_colorbar(title = "Average\nexpression"),
         size = guide_legend(title = "Percent\nexpressed"))

ggsave("results/celllinepanel_combined/cluster_deg_dotplot_manuscript.pdf", height = 7.5, width = 4)


sampleorder <- c("697", "KASUMI2", "RCHACV", "LP1", "JJN3", "PL21",
                 "SKM1", "MM1S", "MONOMAC1", "MOLM13", "ALLSIL", "MOLM14", "GRANTA519",
                 "NUDHL1", "SUPT11", "NALM6", "AMO1", "L363", "DND41", "THP1",
                 "SUDHL4", "OCIM1", "JURKAT", "RI1", "K562", "GDM1")

nkdonors <- c("NK1", "NK2", "NK3")

nkdonors_all <- c("NK1", "NK2", "NK3", "NK4", "NK5", "NK6")

grid <- expand.grid(sampleorder, nkdonors) %>% arrange(Var1)
plotorder <- apply(grid, 1, paste, collapse="-")
plotorder <- c("NK1", "NK2", "NK3", plotorder)

combined_seurat_3donors <- combined_seurat
combined_seurat_3donors$hash.ID <- gsub("NK-expanded", "NK1", combined_seurat_3donors$hash.ID)
combined_seurat_3donors <- subset(combined_seurat_3donors, hash.ID %in% plotorder)
combined_seurat_3donors$hash.ID <- factor(combined_seurat_3donors$hash.ID, levels = plotorder)



p1 <- DimPlot(combined_seurat_plot, reduction = "ref.umap", group.by = "predicted.cluster", split.by = "hash.ID", ncol = 12) &
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


ggsave("results/celllinepanel_integrated/umap_celllines_3donors.pdf", p1, height = 9, width = 18)


# 6 donors

sampleorder <- c("697", "ALLSIL", "JURKAT", "K562", "RI1", "GDM1")
nkdonors <- c("NK1", "NK2", "NK3", "NK4", "NK5", "NK6")

grid <- expand.grid(sampleorder, nkdonors) %>% arrange(Var1)
plotorder <- apply(grid, 1, paste, collapse="-")
plotorder <- c("NK1", "NK2", "NK3", "NK4", "NK5", "NK6", plotorder)

combined_seurat_plot <- combined_seurat
combined_seurat_plot$hash.ID <- gsub("NK-expanded", "NK1", combined_seurat_plot$hash.ID)
combined_seurat_plot <- subset(combined_seurat_plot, hash.ID %in% plotorder)
combined_seurat_plot$hash.ID <- factor(combined_seurat_plot$hash.ID, levels = plotorder)


p1 <- DimPlot(combined_seurat_plot, reduction = "ref.umap", group.by = "predicted.cluster", split.by = "hash.ID", ncol = 6) &
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

ggsave("results/celllinepanel_integrated/umap_celllines_6donors.pdf", p1, height = 8, width = 9)


# 3 donors merged

min(table(combined_seurat_3donors$hash.ID))
Idents(combined_seurat_3donors) <- "hash.ID"

combined_seurat_3donors$hash.ID.grouped <- gsub("[0-9]$", "", combined_seurat_3donors$hash.ID)



p1 <- DimPlot(combined_seurat_3donors, reduction = "ref.umap", group.by = "predicted.cluster", split.by = "hash.ID.grouped", ncol = 12) &
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



sampleorder <- c("697", "KASUMI2", "RCHACV", "LP1", "JJN3", "PL21",
                 "SKM1", "MM1S", "MONOMAC1", "MOLM13", "ALLSIL", "MOLM14", "GRANTA519",
                 "NUDHL1", "SUPT11", "NALM6", "AMO1", "L363", "DND41", "THP1",
                 "SUDHL4", "OCIM1", "JURKAT", "RI1", "K562", "GDM1")


plotorder <- c("NK", paste0(sampleorder, "-NK"))

plotdata <- FetchData(combined_seurat_3donors, vars = c("refUMAP_1", "refUMAP_2", "hash.ID.grouped", "predicted.cluster"))
plotdata_selected <- plotdata %>%
  mutate(hash.ID.grouped = factor(hash.ID.grouped, levels = plotorder))
plotdata_background <- FetchData(combined_seurat_3donors, vars = c("refUMAP_1", "refUMAP_2"))

ggplot(plotdata_selected, aes(x = refUMAP_1, y = refUMAP_2)) +
  theme_bw() +
  stat_density_2d(data = plotdata_background, aes(x = refUMAP_1, y = refUMAP_2, fill = ..density..), contour = F, geom = "raster") +
  stat_density_2d(aes(x = refUMAP_1, y = refUMAP_2, color = ..level..), contour_var = "ndensity", h = 3) +
  scale_fill_gradientn(colors = brewer.pal(9, "Greys")[1:6]) +
  scale_color_viridis_c(name = "Density", option = "magma", begin = 0.1, end = 1) +
  scale_x_continuous(limits = c(-7, 7)) +
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
  facet_wrap(. ~ hash.ID.grouped, ncol = 12)



# Dotplot
clu0_markers <- c("KLRC1", "GZMA", "KLRB1", "NCAM1", "GZMK")
clu1_markers <- c("KLRC2", "GZMH", "LAG3", "HLA-DRA")
clu2_markers <- c("TNFRSF18", "TNFRSF9", "CRTAM", "ENTPD1", "HAVCR2", "TIGIT", "TNFSF10", "BCL2L11")
clu3_markers <- c("ISG15", "MX1" ,"MX2", "OAS1", "OAS2", "OAS3", "IRF7", "IRF9", "EIF2AK2")
clu4_markers <- c("IFNG", "CCL4", "CCL3", "TNF", "CD69")

markers <- c(clu0_markers, clu1_markers, clu2_markers, clu3_markers, clu4_markers)

DotPlot(combined_seurat_3donors,
        group.by = "predicted.cluster",
        features = markers,
        cols = "RdBu") +
  coord_flip() +
  scale_x_discrete(limits = rev, drop = F) +
  scale_y_discrete(position = "right", drop = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black")) +
  xlab("") +
  ylab("") +
  guides(color = guide_colorbar(title = "Average\nexpression"),
         size = guide_legend(title = "Percent\nexpressed"))


# cluster fractions

Idents(combined_seurat_3donors) <- "hash.ID"

saveRDS(combined_seurat_3donors, "results/celllinepanel_integrated/combined_seurat_integrated_3donors.rds")
     
combined_seurat_3donors <- readRDS("results/celllinepanel_integrated/combined_seurat_integrated_3donors.rds")

cluster_freq <- data.frame(table(Idents(combined_seurat_3donors), combined_seurat_3donors$predicted.cluster))
cluster_freq <- cluster_freq %>%
  filter(Var1 != "697") %>% 
  group_by(Var1, Var2) %>%
  summarize(count = sum(Freq)) %>%
  mutate(freq = count / sum(count))

cluster_freq_mean <- cluster_freq %>% 
  mutate(condition = gsub("^NK.*", "NK cells only", gsub("-NK.", "", Var1))) %>% 
  group_by(condition, Var2) %>% 
  summarize(mean_freq = mean(freq))

  
sample_order <- cluster_freq_mean %>%
  filter(Var2 %in% c("Resting (0)", "Adaptive (1)")) %>%
  group_by(condition) %>% 
  summarize(freq_sum = sum(mean_freq)) %>%
  arrange(freq_sum) %>%
  dplyr::select(condition) %>%
  tibble::deframe()

cluster_freq_pert <- cluster_freq_mean %>%
  mutate(condition = factor(condition, levels = rev(sample_order))) %>% 
  mutate(mean_freq = 100*mean_freq)


# Barplot of mean cluster fractions (Figure 1E)
p_cluster_barplot <- ggplot(cluster_freq_pert, aes(x = condition, y = mean_freq, fill = Var2)) +
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
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 10)) +
  guides(fill = "none")


# cancer type color dots

# read cancer types and colors
cols <- fread("data/nk_crispr_colors.txt", data.table = F) %>% dplyr::rename(cancer = cancer_type)
cols <- rbind(cols, data.frame(color = c("lightblue", "grey70"), cancer = c("CML", "NK cells only")))

cancertypes <- readxl::read_excel("data/cell_line_subtypes.xlsx")

cancertype_color <- merge(cancertypes, cols)

sample_order <- gsub("-NK*", "", sample_order)
  

cancertype_color_df <- cancertype_color %>% 
  filter(cell_line %in% sample_order) %>% 
  mutate(gene = "Cancer type") %>% 
  dplyr::select(gene, cell_line, cancer) %>% 
  mutate(cell_line = factor(cell_line, levels = rev(sample_order))) %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  mutate(cancer_main = factor(cancer_main, levels = c("T-ALL", "AML/CML", "BCL", "B-ALL", "MM"))) %>% 
  mutate(cancer_color = gsub("BCL", "DLBCL", cancer_main))


cancertype_color_vector <- cancertype_color %>% filter(cell_line %in% sample_order) %>% dplyr::select(cancer, color) %>% unique() %>% dplyr::select(color) %>% tibble::deframe()
names(cancertype_color_vector) <- cancertype_color %>% filter(cell_line %in% sample_order) %>% dplyr::select(cancer, color) %>% unique() %>% dplyr::select(cancer) %>% tibble::deframe()

cancertype_color_vector <- cancertype_color_vector[c("AML", "T-ALL", "DLBCL", "MCL", "MM", "B-ALL", "NK cells only")]
names(cancertype_color_vector) <- c("AML/CML", "T-ALL", "DLBCL", "MCL", "MM", "B-ALL", "NK cells only")

p_cancertypedots <- ggplot(cancertype_color_df, aes(x = cell_line, y = gene, color = cancer_color)) +
  geom_point(size = 3) +
  scale_color_manual("Cancer type", values = cancertype_color_vector,
                     guide = guide_legend(title.position = "top")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))

p_cluster_barplot / p_cancertypedots + plot_layout(heights = c(8,1))

ggsave("results/celllinepanel_integrated/cluster_fraction_manuscript.pdf", height = 3.5, width = 5)


# Barplot of cluster fractions with each donor (Figure S1D)
combined_seurat$hash.ID <- gsub("NK-expanded", "NK1", combined_seurat$hash.ID)
Idents(combined_seurat) <- "hash.ID"

cluster_freq <- data.frame(table(Idents(combined_seurat), combined_seurat$predicted.cluster))

cluster_freq <- cluster_freq %>%
  filter(Var1 != "697") %>% 
  group_by(Var1, Var2) %>%
  summarize(count = sum(Freq)) %>%
  mutate(freq = count / sum(count))

sample_order_6donors <- cluster_freq_mean %>%
  filter(Var2 %in% c("Resting (0)", "Adaptive (1)")) %>%
  group_by(condition) %>% 
  summarize(freq_sum = sum(mean_freq)) %>%
  arrange(freq_sum) %>%
  dplyr::select(condition) %>%
  tibble::deframe()

cluster_freq_donors <- cluster_freq %>%
  mutate(condition = gsub("^NK.*", "NK cells only", gsub("-NK.", "", Var1))) %>% 
  mutate(condition = factor(condition, levels = rev(sample_order_6donors))) %>% 
  mutate(Var1 = factor(Var1, levels = sort(as.character(unique(cluster_freq$Var1))))) %>% 
  mutate(freq = 100*freq)

p_cluster_barplot_donors <- ggplot(cluster_freq_donors, aes(x = Var1, y = freq, fill = Var2)) +
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
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.line.x = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 10),
        strip.text = element_text(angle = 90, hjust = 0, size = 10),
        strip.background = element_blank()) +
  guides(fill = "none") +
  facet_grid(. ~ condition, scales = "free_x", space = "free")

sample_order_nk <- gsub(" cells only", "", sample_order)

plotdata_6donors <- FetchData(combined_seurat, vars = c("refUMAP_1", "refUMAP_2", "hash.ID", "predicted.cluster"))
plotdata_selected_6donors <- plotdata_6donors %>%
  mutate(donor = gsub(".*-", "", hash.ID),
         cell_line = gsub("-.*", "", gsub("NK[0-9]", "NK", hash.ID))) %>% 
  mutate(donor = factor(donor, levels = c("NK1", "NK2", "NK3", "NK4", "NK5", "NK6")),
         cell_line = factor(cell_line, levels = rev(sample_order_nk)))


cancertype_color_df_donors <- plotdata_selected_6donors %>% select(condition = hash.ID, donor, cell_line) %>% unique() %>% left_join(cancertype_color_df) %>% 
  filter(!is.na(donor)) %>% 
  mutate(cell_line = gsub("NK", "NK cells only", cell_line)) %>% 
  mutate(gene = "NK cell donor") %>% 
  mutate(cell_line = factor(cell_line, levels = rev(sample_order))) %>% 
  mutate(donor = factor(donor, levels = c("NK1", "NK2", "NK3", "NK4", "NK5", "NK6")))

p_cancertypedots <- ggplot(cancertype_color_df_donors, aes(x = condition, y = gene, color = donor)) +
  geom_point(size = 3) +
  scale_color_manual("NK cell donor", values = brewer.pal(8, "Greys")[3:8],
                     guide = guide_legend(title.position = "top")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "right",
        plot.margin = unit(c(0,0,0,0), "cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  facet_grid(. ~ cell_line, scales = "free_x", space = "free", drop = F)


p_cluster_barplot_donors / p_cancertypedots + plot_layout(heights = c(8,1), guides = "collect")
ggsave("results/celllinepanel_integrated/cluster_fraction_donors_manuscript.pdf", height = 3.25, width = 18)


# Boxplot of clusters stratified by cancer type (Figure 1E)
cluster_freq_pert_cancertypes <- merge(cancertype_color_df, cluster_freq_pert, by.x = "cell_line", by.y = "condition")

cluster_freq_pert_maincancertypes <- cluster_freq_pert_cancertypes %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  mutate(cancer_main = factor(cancer_main, levels = c("BCL", "T-ALL", "AML/CML", "MM", "B-ALL", "NK cells only"))) %>% 
  mutate(cancer_color = gsub("AML/CML", "AML", gsub("BCL", "DLBCL", cancer_main))) %>% 
  filter(Var2 %in% c("Activated (2)", "Type I IFN (3)")) %>% 
  group_by(cancer_main, cancer_color, cell_line) %>% 
  summarize(cluster23sum = sum(mean_freq))


# test significance
comparisons = list(c("BCL", "T-ALL"), c("BCL", "AML/CML"), c("BCL", "MM"), c("BCL", "B-ALL"),
                   c("T-ALL", "AML/CML"), c("T-ALL", "MM"), c("T-ALL", "B-ALL"),
                   c("AML/CML", "MM"),  c("AML/CML", "B-ALL"),
                   c("MM", "B-ALL"))


test_wilcox_pairwise <- ggpubr::compare_means(cluster23sum ~ cancer_main, comparisons = comparisons, p.adjust.method = "BH", method = 'wilcox.test', data = cluster_freq_pert_maincancertypes) %>% arrange(p)
test_wilcox_pairwise <- test_wilcox_pairwise %>% mutate(y.position = c(120, 100, 110, rep(NA, 12)),
                                                        p.signif = ifelse(p < 0.05, "*", "ns"))

names(cancertype_color_vector) <- gsub("DLBCL", "BCL", names(cancertype_color_vector))

ggplot(cluster_freq_pert_maincancertypes, aes(x = cancer_main, y = cluster23sum)) +
  geom_boxplot(aes(fill = cancer_main), outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  theme_cowplot() +
  scale_fill_manual("Cancer type", values = cancertype_color_vector,
                    guide = guide_legend(title.position = "top")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,125), breaks = c(0,25,50,75)) +
  xlab("") +
  ylab("% of NK cells in activated/\ntype I IFN clusters") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.y = 50, label.x = "MM") +
  ggpubr::stat_pvalue_manual(test_wilcox_pairwise, label = "p.signif")


ggsave("results/celllinepanel_integrated/cluster_fraction_boxplot_manuscript.pdf", width = 2.75, height = 3.5)


# Heatmap of cluster fractions by cancer type (Figure 1E)
cluster_freq_pert_maincancertypes_wide <- cluster_freq_pert_cancertypes %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  mutate(cancer_main = factor(cancer_main, levels = c("BCL", "T-ALL", "AML/CML", "MM", "B-ALL", "NK cells only"))) %>% 
  mutate(cancer_color = gsub("AML/CML", "AML", gsub("BCL", "DLBCL", cancer_main))) %>% 
  group_by(cancer_main, Var2) %>% 
  summarize(mean_cluster_freq = mean(mean_freq)) %>% 
  tidyr::pivot_wider(names_from = cancer_main, values_from = mean_cluster_freq)


cluster_freq_pert_maincancertypes_mat <- as.matrix(cluster_freq_pert_maincancertypes_wide[,-1])
rownames(cluster_freq_pert_maincancertypes_mat) <- gsub("\\(.*", "", cluster_freq_pert_maincancertypes_wide$Var2)
cluster_freq_pert_maincancertypes_mat <- cluster_freq_pert_maincancertypes_mat[c(3,4),]


ht <- Heatmap(cluster_freq_pert_maincancertypes_mat,
        cluster_rows = F,
        cluster_columns = F,
        col = brewer.pal(name = "RdPu", n = 9),
        use_raster = T,
        show_column_names = F,
        show_row_names = T,
        row_names_side = "left",
        show_heatmap_legend = T,
        raster_quality = 3,
        heatmap_legend_param = list(title = "% NK cells in cluster",
                                    title_gp = gpar(fontsize = 10),
                                    labels_gp = gpar(fontsize = 10),
                                    grid_height = unit(0.2, "cm"),
                                    grid_width = unit(2, "mm"),
                                    title_position = "topcenter",
                                    legend_direction = "horizontal"))

pdf("results/celllinepanel_integrated/cluster_fraction_heatmap_manuscript.pdf", width = 2.75, height = 1.2)
draw(ht, heatmap_legend_side = "bottom")
dev.off()


# cluster fraction correlation with PRISM AUC
nk_prism_auc <- fread("data/NK_PRISM_heme_pctviability_auc_norm.txt")
nk_prism_fm_auc <- nk_prism_auc %>%
  filter(Condition == "0") %>%
  dplyr::select(auc_norm, cell_line) %>% 
  as.data.frame()

cluster_freq_pert_prism <- cluster_freq_pert %>%
  dplyr::rename(cell_line = condition, cluster = Var2) %>%
  filter(cluster %in% c("Activated (2)", "Type I IFN (3)")) %>%
  group_by(cell_line) %>% 
  summarize(freq_sum = sum(mean_freq)) %>% 
  left_join(nk_prism_fm_auc) %>%
  left_join(cancertype_color_df) #%>% 

# save df 
cluster_freq_pert_prism_d14 <- cluster_freq_pert_prism

cancertype_color_vector <- cancertype_color_vector[c("AML/CML", "T-ALL", "BCL", "MM", "B-ALL")]


# Scatter plot (Figure 3D)
ggpubr::ggscatter(cluster_freq_pert_prism, y = "auc_norm", x = "freq_sum", shape = 21, size = 3, fill = "cancer_main",
                  add = "reg.line", cor.coef = T, cor.coeff.args = list(label.sep = "\n"), cor.method = "spearman", palette = cancertype_color_vector,#, label = "cell_line", repel = T,
                  font.label = c(10, "plain"), cor.coef.coord = c(50, 1.1)) +
  xlab("% of NK cells in activated/\ntype I IFN clusters") +
  ylab("PRISM AUC") +
  theme(legend.title = element_blank(),
        plot.margin = unit(c(0,0.5,0,0), "cm")) +
  scale_y_continuous(breaks = c(0.5, 1)) +
  guides(fill="none")
ggsave("results/celllinepanel_integrated/cluster_fraction_prism_manuscript.pdf", height = 3, width = 3.3)


# Scatter plots (Figure S3B)
cluster_freq_pert_prism <- cluster_freq_pert %>%
  dplyr::rename(cell_line = condition, cluster = Var2) %>%
  group_by(cell_line) %>% 
  left_join(nk_prism_fm_auc) %>%
  left_join(cancertype_color_df)


ggpubr::ggscatter(cluster_freq_pert_prism, y = "auc_norm", x = "mean_freq", shape = 21, size = 3, fill = "cancer_main",
                  add = "reg.line", cor.coef = T, cor.method = "spearman", cor.coeff.args = list(label.sep = "\n", label.x.npc = 0.6, label.y = 1.2), palette = cancertype_color_vector) +#, label = "cell_line", repel = T) +
  xlab("% of NK cells in  cluster") +
  ylab("PRISM AUC") +
  facet_wrap(. ~ cluster, ncol = 1, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1.4)) +
  theme(legend.title = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_legend(nrow = 3))


ggsave("results/celllinepanel_integrated/cluster_fraction_prism_allclusters_manuscript_vertical.pdf", height = 10, width = 2.5)


# Pairwise correlations betweeen donors (Figure S1C)
cluster_freq_donors_mat <- cluster_freq_donors %>% 
  mutate(donor = gsub(".*-", "", Var1),
         cell_line = gsub("-.*", "", gsub("NK[0-9]", "NK", Var1))) %>% 
  filter(donor %in% c("NK1", "NK2", "NK3")) %>% 
  mutate(donor_cluster = paste(donor, gsub(" \\(.*", "", Var2), sep = "_")) %>% 
  tidyr::pivot_wider(id_cols = cell_line, names_from = donor_cluster, values_from = freq)

mat <- cluster_freq_donors_mat
rownames(mat) <- mat$cell_line
mat$cell_line <- NULL

cor_mat <- psych::corr.test(mat, method = "spearman", adjust = "none")

rownames(cor_mat$r) <- gsub("_", " ", rownames(cor_mat$r))
colnames(cor_mat$r) <- gsub("_", " ", colnames(cor_mat$r))

rownames(cor_mat$p) <- gsub("_", " ", rownames(cor_mat$p))
colnames(cor_mat$p) <- gsub("_", " ", colnames(cor_mat$p))


clusters <- gsub(".*_", "", colnames(mat))
df <- data.frame(clusters)
colnames(df) <- "Cluster"

scale_fill_manual(values = a_palette[c(1,2,6,8,5)])


ha <- HeatmapAnnotation(df = df,
                        col = list(Cluster = structure(a_palette[c(1,2,6,8,5)], 
                                                names = c("Resting", "Adaptive", "Type I IFN", "Activated", "Cytokine"))
                        ),
                        gap = unit(0.5, "mm"),
                        height = unit(0.3, "cm"),
                        simple_anno_size_adjust = T,
                        annotation_name_gp = gpar(fontsize = 6),
                        annotation_legend_param = list(
                          Cluster = list(title = "Cluster", title_gp = gpar(fontsize = 6), 
                                             labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(3, "mm"))))

row_ha <- HeatmapAnnotation(df = df,
                            which = "row",
                        col = list(Cluster = structure(a_palette[c(1,2,6,8,5)], 
                                                       names = c("Resting", "Adaptive", "Type I IFN", "Activated", "Cytokine"))
                        ),
                        gap = unit(0.5, "mm"),
                        width = unit(0.3, "cm"),
                        simple_anno_size_adjust = T,
                        show_legend = F,
                        annotation_name_gp = gpar(fontsize = 6),
                        annotation_legend_param = list(
                          Cluster = list(title = "Cluster", title_gp = gpar(fontsize = 6), 
                                         labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(4, "mm"))))

  
  

ht <- Heatmap(cor_mat$r,
        name = "Correlation",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(cor_mat$p[i, j] < 0.05)
            grid.text(sprintf("%.1f", cor_mat$r[i, j]), x, y, gp = gpar(fontsize = 6))
        },
        top_annotation = ha,
        left_annotation = row_ha,
        col = colorRamp2(seq(-1, 1, length.out = 9), pals::ocean.balance(13)[3:11]),
        rect_gp = gpar(col = "white", lwd = unit(0.4, "mm")),
        border_gp = gpar(col = "black", lty = 1, lwd = 1),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(title = "Correlation",
                                    title_gp = gpar(fontsize = 6),
                                    labels_gp = gpar(fontsize = 6),
                                    grid_height = unit(0.2, "cm"),
                                    grid_width = unit(2, "mm"),
                                    border = NA,
                                    legend_direction = "horizontal"))


pdf("results/celllinepanel_integrated/cluster_fraction_correlation_heatmap.pdf", height = 4.5, width = 4.5)
draw(ht, merge_legend = T, heatmap_legend_side = "bottom")
dev.off()

