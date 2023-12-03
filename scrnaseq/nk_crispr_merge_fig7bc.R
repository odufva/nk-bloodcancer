
# Merge data from two experiments with NK cells co-cultured with CRISPR KO cell lines and map to cell line panel NK cell reference

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

dir.create("results/crispr_combined")

# load K562/SUDHL4 co-culture experiment NK cell object
k562_sudhl4_seurat <- readRDS("results/k562_sudhl4_new//hash_seurat_singlet.rds")

# load NALM6 co-culture experiment NK cell object
nalm6_seurat <- readRDS("results/nalm6/hash_seurat_singlet.rds")

# load 26 cell line panel co-culture experiment NK cell object (only donor NK1)
celllinepanel_seurat <- readRDS("results/celllinepanel/hto_singlets/nk/batchcorrected/celllinepanel_nk1_batchcorrected_seurat.rds")

celllinepanel_seurat$cluster <- factor(celllinepanel_seurat$seurat_clusters,
                                 levels = c(0:4),
                                 labels = c("Resting (0)",
                                            "Adaptive (1)",
                                            "Activated (2)",
                                            "Type I IFN (3)",
                                            "Cytokine (4)"))# merge K562, SUDHL4, NALM6

k562_sudhl4_seurat$hash.ID.grouped <- gsub("\\.1|\\.2|\\.4|\\.5", "", k562_sudhl4_seurat$hash.ID)
k562_sudhl4_seurat$hash.ID.grouped <- gsub("-", " ", k562_sudhl4_seurat$hash.ID.grouped)
k562_sudhl4_seurat$hash.ID.grouped <- gsub("PMAIP", "PMAIP1", k562_sudhl4_seurat$hash.ID.grouped)


nalm6_seurat$hash.ID.grouped <- paste0("NALM6 ", nalm6_seurat$hash.ID)
nalm6_seurat$hash.ID.grouped <- gsub("-prestained|-old-method|.1|.2|.4|.5", "", nalm6_seurat$hash.ID.grouped)
nalm6_seurat$hash.ID.grouped <- gsub("-", " ", nalm6_seurat$hash.ID.grouped)
nalm6_seurat$hash.ID.grouped <- gsub("NALM6 NK only", "NK only", nalm6_seurat$hash.ID.grouped)

ksn_seurat <- merge(k562_sudhl4_seurat, nalm6_seurat)


celllinepanel_seurat <- ScaleData(celllinepanel_seurat, assay = 'RNA')

# Run UMAP on cell line panel data with returbn.model = T to be able to map other data to it
celllinepanel_seurat <- RunUMAP(celllinepanel_seurat, reduction = "pca", dims = 1:20, return.model = TRUE)


# Map K562, SUDHL4, NALM6 data to cell line panel reference
anchors <- FindTransferAnchors(
  reference = celllinepanel_seurat,
  query = ksn_seurat,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:20
)

ksn_seurat <- MapQuery(
  anchorset = anchors,
  query = ksn_seurat,
  reference = celllinepanel_seurat,
  refdata = list(
    cluster = "cluster"),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

saveRDS(ksn_seurat, "results/crispr_combined/k562_sudhl4_nalm6_seurat.rds")

DimPlot(ksn_seurat, group.by = "predicted.cluster", label = T, repel = T) &
  theme_bw(base_size = 12) &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

DimPlot(ksn_seurat, split.by = "hash.ID", ncol = 6) &
  theme_bw(base_size = 12) &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

DimPlot(ksn_seurat, split.by = "hash.ID.grouped", ncol = 6) &
  theme_bw(base_size = 12) &
  xlab("UMAP 1") &
  ylab("UMAP 2") &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

plotdata <- FetchData(ksn_seurat, vars = c("refUMAP_1", "refUMAP_2", "hash.ID.grouped"))
plotdata <- plotdata %>% mutate(hash.ID.grouped = factor(hash.ID.grouped, levels = c("NK only", "K562 sgCtrl", "K562 sgJAK1", "K562 sgB2M", "K562 sgKCNH2", "K562 sgGFI1B", "K562 sgNCR3LG1", "K562 sgCD58", 
                                                                                     "SUDHL4 sgCtrl", "SUDHL4 sgMETTL17", "SUDHL4 sgYTHDF2", "SUDHL4 sgFADD", "SUDHL4 sgPMAIP1", "SUDHL4 sgBID",
                                                                                     "NALM6 sgCtrl", "NALM6 sgSPN", "NALM6 sgCMIP", "NALM6 sgSPPL3", "NALM6 sgFADD")))

# PLot density UMAP
ggplot(plotdata, aes(x = refUMAP_1, y = refUMAP_2)) +
  theme_bw() +
  stat_density_2d(aes(x = refUMAP_1, y = refUMAP_2, color = ..level..), contour_var = "ndensity", h = 2) +
  scale_color_viridis_c(name = "Density", option = "magma", begin = 0.1, end = 1) +
  scale_x_continuous(limits = c(-5, 7)) +
  scale_y_continuous(limits = c(-7, 5.5)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  guides(color = "none") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  facet_wrap(. ~ hash.ID.grouped)
ggsave("results/crispr_combined/perturbation_density_umaps.pdf", height = 5, width = 6)


# main figure plot examples (Figure 7C)
selected <- c("NK only", "K562 sgCD58", "K562 sgCtrl", "K562 sgJAK1")

plotdata_selected <- plotdata %>%
  filter(hash.ID.grouped %in% selected) %>%
  mutate(hash.ID.grouped = factor(hash.ID.grouped, levels = selected))

plotdata_selected_background1 <- plotdata_selected %>% mutate(hash.ID.grouped == "NK only")
plotdata_selected_background2 <- plotdata_selected %>% mutate(hash.ID.grouped == "K562 sgCD58")
plotdata_selected_background3 <- plotdata_selected %>% mutate(hash.ID.grouped == "K562 sgCtrl")
plotdata_selected_background1 <- plotdata_selected %>% mutate(hash.ID.grouped == "K562_sgJAK1")

plotdata_background <- FetchData(celllinepanel_seurat, vars = c("UMAP_1", "UMAP_2"))

ggplot(plotdata_selected, aes(x = refUMAP_1, y = refUMAP_2)) +
  theme_bw() +
  stat_density_2d(data = plotdata_background, aes(x = UMAP_1, y = UMAP_2, fill = ..density..), contour = F, geom = "raster") +
  stat_density_2d(aes(x = refUMAP_1, y = refUMAP_2, color = ..level..), contour_var = "ndensity", h = 2) +
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
  facet_wrap(. ~ hash.ID.grouped, ncol = 1)
ggsave("results/crispr_combined/perturbation_density_umaps_selected.pdf", height = 6, width = 2)


# NK activation module score
nk_cluster_deg <- fread("results/celllinepanel/hto_singlets/nk/batchcorrected/clusters_deg.txt", data.table = F)

nkactivation <- nk_cluster_deg %>% filter(cluster == "Activated (2)") %>% top_n(50, desc(p_val)) %>% select(gene) %>% tibble::deframe()


ksn_seurat <- AddModuleScore(ksn_seurat, features = list(nkactivation), name = c("nk_activation_score"), assay = "RNA")

actscore <- FetchData(ksn_seurat, vars = c("hash.ID.grouped", "nk_activation_score1"))

actscore_k562 <- actscore %>% filter(grepl("K562", hash.ID.grouped))
actscore_sudhl4 <- actscore %>% filter(grepl("SUDHL4", hash.ID.grouped))
actscore_nalm6 <- actscore %>% filter(grepl("NALM6", hash.ID.grouped))

actscore_k562$nk_activation_score_normalized <- actscore_k562$nk_activation_score1 - median(actscore_k562$nk_activation_score1[actscore_k562$hash.ID.grouped=="K562 sgCtrl"])
actscore_sudhl4$nk_activation_score_normalized <- actscore_sudhl4$nk_activation_score1 - median(actscore_sudhl4$nk_activation_score1[actscore_sudhl4$hash.ID.grouped=="SUDHL4 sgCtrl"])
actscore_nalm6$nk_activation_score_normalized <- actscore_nalm6$nk_activation_score1 - median(actscore_nalm6$nk_activation_score1[actscore_nalm6$hash.ID.grouped=="NALM6 sgCtrl"])

# Compare to corresponding control for each cell line
actscore_normalized <- rbind(actscore_k562, actscore_sudhl4, actscore_nalm6)

actscore_order <- actscore_normalized %>% group_by(hash.ID.grouped) %>% summarize(median = median(nk_activation_score_normalized)) %>% arrange(median) %>% select(hash.ID.grouped) %>% tibble::deframe()

actscore_normalized <- actscore_normalized %>% mutate(hash.ID.grouped = factor(hash.ID.grouped, levels = rev(actscore_order))) %>% group_by(hash.ID.grouped) %>% mutate(median = median(nk_activation_score_normalized)) %>% ungroup()

# Test significance
test_wilcox_pairwise_k562 <- rstatix::wilcox_test(nk_activation_score1 ~ hash.ID.grouped, ref.group = "K562 sgCtrl", p.adjust.method = "BH", data = actscore_k562)
test_wilcox_pairwise_k562 <- test_wilcox_pairwise_k562 %>% mutate(y.position = 1) %>% mutate(p = gsub("ns", "", p), p.adj.signif = gsub("ns", "", p.adj.signif))

test_wilcox_pairwise_sudhl4 <- rstatix::wilcox_test(nk_activation_score1 ~ hash.ID.grouped, ref.group = "SUDHL4 sgCtrl", p.adjust.method = "BH", data = actscore_sudhl4)
test_wilcox_pairwise_sudhl4 <- test_wilcox_pairwise_sudhl4 %>% mutate(y.position = 1) %>% mutate(p = gsub("ns", "", p), p.adj.signif = gsub("ns", "", p.adj.signif))

test_wilcox_pairwise_nalm6 <- rstatix::wilcox_test(nk_activation_score1 ~ hash.ID.grouped, ref.group = "NALM6 sgCtrl", p.adjust.method = "BH", data = actscore_nalm6)
test_wilcox_pairwise_nalm6 <- test_wilcox_pairwise_nalm6 %>% mutate(y.position = 1) %>% mutate(p = gsub("ns", "", p), p.adj.signif = gsub("ns", "", p.adj.signif))

test_wilcox_pairwise <- rbind(test_wilcox_pairwise_k562, test_wilcox_pairwise_sudhl4, test_wilcox_pairwise_nalm6)
test_wilcox_pairwise$p.adj <- p.adjust(test_wilcox_pairwise$p, method = "BH")

# Plot figure 7B
ggplot(actscore_normalized, aes(y = nk_activation_score_normalized, x = hash.ID.grouped)) +
  geom_boxplot(aes(fill = median)) +
  scale_fill_gradientn("Median", colors = pals::ocean.curl(12)[3:10]) +
  coord_flip() +
  xlab("") +
  ylab("Normalized activation score") +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.3, "cm")) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  scale_y_continuous(limits = c(-0.6, 1.1))

ggsave("results/crispr_combined/activationscore_boxplot_separatecontrols.pdf", height = 5, width = 6)

