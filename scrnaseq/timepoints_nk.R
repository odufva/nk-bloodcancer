
# NK cell co-culture time course experiment, hashing scRNA-seq analysis
# NK cell analysis


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
library(slingshot)
library(awtools)
library(patchwork)
library(viridis)
library(fgsea)
library(circlize)
library(ComplexHeatmap)

dir.create("results/timepoints/hto_singlets/nk")

"%ni%" <- Negate("%in%")

hash_seurat_singlet <- readRDS("results/timepoints/timepoints_seurat.rds")

hash_seurat_nk <- subset(hash_seurat_singlet, seurat_clusters %in% c("0", "6", "12", "15"))

# Select the top 1000 most variable features
hash_seurat_nk <- FindVariableFeatures(hash_seurat_nk, selection.method = "mean.var.plot")

# Scale RNA data and regress out cell cycle
hash_seurat_nk <- CellCycleScoring(hash_seurat_nk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
hash_seurat_nk <- ScaleData(hash_seurat_nk, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(hash_seurat_nk))

# Run PCA
hash_seurat_nk <- RunPCA(hash_seurat_nk, features = VariableFeatures(hash_seurat_nk))

ElbowPlot(hash_seurat_nk, ndims = 50, reduction = "pca")


# Select the top 20 PCs for clustering and UMAP
hash_seurat_nk <- FindNeighbors(hash_seurat_nk, reduction = "pca", dims = 1:20)
hash_seurat_nk <- FindClusters(hash_seurat_nk, resolution = 0.4, verbose = FALSE)
hash_seurat_nk <- RunUMAP(hash_seurat_nk, reduction = "pca", dims = 1:20)


## Visualize
DimPlot(hash_seurat_nk, group.by = "seurat_clusters", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/nk/umap.png", width = 6, height = 4)

DimPlot(hash_seurat_nk, group.by = "Phase", label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(3)) + labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/nk/umap_cycle.png", width = 6, height = 4)

DimPlot(hash_seurat_nk, group.by = "seurat_clusters", split.by = "hash.ID", label = T, repel = T, ncol = 8) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/nk/umap_celllines.png", width = 14, height = 12)

DimPlot(hash_seurat_nk, group.by = "orig.ident", label = T, repel = T) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())
ggsave("results/timepoints/hto_singlets/nk/umap_batch.png", width = 6, height = 4)

# DEGs of clusters
Idents(hash_seurat_nk) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_nk, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/nk/clusters_deg.txt", sep = "\t", quote = F, row.names = F)

# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_nk,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/nk/cluster_deg_dotplot.pdf", height = 4, width = 18)

# Heatmap of top cluster DEGs
top25 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

p <- DoHeatmap(hash_seurat_nk,
               features = top25$gene,
               raster = F,
               angle = 0) + NoLegend()
ggsave("results/timepoints/hto_singlets/nk/cluster_deg_heatmap.pdf", p, height = 15, width = 15)

hash_seurat_nk$seurat_clusters <- factor(hash_seurat_nk$seurat_clusters, levels = c("0", "4", "3", "2", "1"),
                                         labels = c("Resting", "Cytokine", "Early activated", "Activated", "Type I IFN"))

# Save object
saveRDS(hash_seurat_nk, "results/timepoints/hto_singlets/nk/timepoints_seurat_nk.rds")

## ----------------------------

sampleorder <- c("K562", "JURKAT", "GDM1", "RI1")

timepoints <- c("1h", "3h", "6h", "12h", "24h")

grid <- expand.grid(sampleorder, timepoints) %>% arrange(Var1)
plotorder <- apply(grid, 1, paste, collapse="-")
plotorder <- c("NK-0h-1", "NK-0h-2", plotorder)

hash_seurat_nk_plot <- subset(hash_seurat_nk, hash.ID %in% plotorder)
hash_seurat_nk_plot$hash.ID <- factor(hash_seurat_nk$hash.ID, levels = plotorder)


DimPlot(hash_seurat_nk_plot, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 5) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = rev(getPalette(7))) +
  labs("UMAP 1", "UMAP 2")
ggsave("results/timepoints/hto_singlets/nk/umap_celllines_order.pdf", width = 10, height = 8)

## ----------------------------

# Plots for manuscript

# Load object
hash_seurat_nk <- readRDS("results/timepoints/hto_singlets/nk/timepoints_nk_seurat.rds")

## Visualize
DimPlot(hash_seurat_nk, group.by = "seurat_clusters", label = F, repel = T) +
  theme_bw(base_size = 12) &
  xlab("") &
  ylab("") &
  ggtitle("") &
  scale_color_manual(values = a_palette[c(1,5,3,6,8)]) &
  theme_void() &
  theme(legend.position = "bottom") &
  NoLegend()
  #guides(color = guide_legend(nrow = 3, byrow = T, override.aes = list(size = 2)))

ggsave("results/timepoints/hto_singlets/nk/umap_clusters.pdf", width = 3, height = 3)

p1 <- DimPlot(hash_seurat_nk, group.by = "seurat_clusters", label = F, repel = T) +
  theme_bw(base_size = 12) &
  xlab("") &
  ylab("") &
  ggtitle("") &
  scale_color_manual(values = a_palette[c(1,5,3,6,8)]) &
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave("results/timepoints/hto_singlets/nk/umap_clusters_rect.pdf", width = 6, height = 4)

# color by time point

hash_seurat_nk$timepoint <- gsub(".*-", "", gsub("-1$|-2$", "", hash_seurat_nk$hash.ID))

hash_seurat_nk$timepoint <- factor(hash_seurat_nk$timepoint, levels = c("0h", "1h", "3h", "6h", "12h", "24h"))

DimPlot(hash_seurat_nk, group.by = "timepoint", label = F, repel = T) +
  theme_bw(base_size = 12) &
  xlab("") &
  ylab("") &
  ggtitle("") &
  scale_color_manual(values = c("grey50", rev(pals::ocean.matter(6)[2:6]))) &
  theme_void() &
  theme(legend.position = "bottom")&
  guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 2)))

ggsave("results/timepoints/hto_singlets/nk/umap_clusters_timepoint.pdf", width = 3, height = 3.5)




# grid arrangement
sampleorder <- c("K562", "JURKAT", "RI1", "GDM1")

timepoints <- c("1h", "3h", "6h", "12h", "24h")

grid <- expand.grid(sampleorder, timepoints) %>% arrange(Var1)
plotorder <- apply(grid, 1, paste, collapse="-")

hash_seurat_nk_coculture <- subset(hash_seurat_nk, hash.ID %in% plotorder)
hash_seurat_nk_coculture$hash.ID <- factor(hash_seurat_nk_coculture$hash.ID, levels = plotorder)

DimPlot(hash_seurat_nk_coculture, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 5) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_color_manual(values = a_palette[c(1,5,3,6,8)]) &
  xlab("") &
  ylab("") &
  guides(color = "none")

ggsave("results/timepoints/hto_singlets/nk/umap_celllines_coculture_order.pdf", width = 10, height = 8)


DimPlot(hash_seurat_nk_coculture, group.by = "seurat_clusters", split.by = "hash.ID", label = F, repel = T, ncol = 5) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = a_palette[c(1,5,3,6,8)]) &
  xlab("") &
  ylab("") &
  guides(color = "none")

ggsave("results/timepoints/hto_singlets/nk/umap_celllines_coculture_order_labels.pdf", width = 10, height = 8)

p2 <- DimPlot(subset(hash_seurat_nk, hash.ID %in% c("NK-0h-1", "NK-0h-2")), group.by = "seurat_clusters", label = F, repel = T) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_color_manual(values = a_palette[c(1,5,3,6,8)]) &
  xlab("") &
  ylab("") &
  ggtitle("NK-only-0h") &
  guides(color = "none")

ggsave("results/timepoints/hto_singlets/nk/umap_nkonly.pdf", width = 4, height = 4)


# merge 0h samples for plotting

hash_seurat_nk_0hmerged <- hash_seurat_nk
hash_seurat_nk_0hmerged$hash.ID <- gsub("0h-1|0h-2", "0h", hash_seurat_nk$hash.ID)
saveRDS(hash_seurat_nk_0hmerged, "results/timepoints/hto_singlets/nk/imepoints_nk_0hmerged_seurat.rds")


# faceted plot with order

plotdata <- FetchData(hash_seurat_nk_coculture, vars = c("UMAP_1", "UMAP_2", "hash.ID", "seurat_clusters"))
plotdata_selected <- plotdata %>%
  mutate(timepoint = gsub(".*-", "", hash.ID),
         cell_line = gsub("-.*", "", hash.ID)) %>% 
  mutate(timepoint = factor(timepoint, levels = timepoints),
         cell_line = factor(cell_line, levels = sampleorder))


p3 <- ggplot(plotdata_selected, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual("", values = a_palette[c(1,5,3,6,8)]) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "top") +
  guides(fill = "none",
         colour = guide_legend(override.aes = list(size=2))) +
  xlab("") +
  ylab("") +
  facet_grid(rows = vars(cell_line), cols = vars(timepoint), switch = "y")
ggsave("results/timepoints/hto_singlets/nk/umap_celllines_coculture_order_labels_grid.pdf", p3, width = 6, height = 5)

layout <- "
ACCC
BCCC
"

p1 + p2 + p3 + plot_layout(design = layout, guides = "collect") & theme(legend.position = "top") 
ggsave("results/timepoints/hto_singlets/nk/umap_celllines__order_labels_grid_clusters.pdf", width = 7.5, height = 5)


# include 0h

plotdata <- FetchData(hash_seurat_nk_0hmerged, vars = c("UMAP_1", "UMAP_2", "hash.ID", "seurat_clusters"))
plotdata_selected <- plotdata %>%
  filter(!hash.ID %in% c("K562-1h", "K562-0h")) %>% 
  mutate(timepoint = gsub(".*-", "", hash.ID),
         cell_line = gsub("-.*", "", hash.ID)) %>%
  mutate(timepoint = ifelse(hash.ID == "NK-0h", "1h", timepoint),
         cell_line = ifelse(hash.ID == "NK-0h", "K562", cell_line)) %>% 
  filter(!is.na(timepoint)) %>% 
  mutate(timepoint = factor(timepoint, levels = timepoints),
         cell_line = factor(cell_line, levels = sampleorder))


p4 <- ggplot(plotdata_selected, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.25) +
  theme_bw() +
  scale_color_manual("", values = a_palette[c(1,5,3,6,8)]) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "top") +
  guides(fill = "none",
         colour = guide_legend(override.aes = list(size=2))) +
  xlab("") +
  ylab("") +
  facet_grid(rows = vars(cell_line), cols = vars(timepoint), switch = "y")
ggsave("results/timepoints/hto_singlets/nk/umap_celllines_coculture_order_labels_grid_with0h.pdf", p4, width = 6, height = 5)


# heatmap of top cluster DEGs

Idents(hash_seurat_nk) <- "seurat_clusters"
all_deg <- FindAllMarkers(hash_seurat_nk, test.use = "t", only.pos = T, return.thresh = 0.05)
fwrite(all_deg, "results/timepoints/hto_singlets/nk/clusters_deg_clusternames.txt", sep = "\t", quote = F, row.names = F)

all_deg <- fread("results/timepoints/hto_singlets/nk/clusters_deg_clusternames.txt")

top5 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

feats <- c("GZMA", "FCGR3A", "LTB",
           "IFNG", "TNF", "CCL3", "CCL4", "CSF2", "CD69", "FOS", "JUN", "EGR1", "NR4A2", "NFKBIA", "NFKBID", "NFKBIZ",
           "GZMB", "SERPINB9", "XCL1", "XCL2", "REL", "BIRC3",
           "TNFRSF9", "CRTAM", "TNFRSF18", "HAVRC2", "TIGIT",
           "ISG15", "ISG20", "MX1", "IFIT1", "IFIT3", "LAG3")

DoHeatmap(subset(hash_seurat_nk, downsample = 200), features = feats, label = T, draw.lines = F, group.bar = T, raster = F,
          group.colors = a_palette[c(1,5,3,6,8)],
          angle = 0,
          hjust = 0.5,
          size = 4) + NoLegend() +
  scale_fill_viridis_c(option = "viridis")
  
ggsave("results/timepoints/hto_singlets/nk/heatmap_clusters.pdf", width = 6, height = 4)


hash_seurat_nk_downsampled <- subset(hash_seurat_nk, downsample = 200)
mat_gexp <- t(FetchData(hash_seurat_nk_downsampled, vars = feats, slot = "scale.data"))
df_annot <- FetchData(hash_seurat_nk_downsampled, vars = c("seurat_clusters"), slot = "data")

cell_order = rownames(df_annot)[order(df_annot$seurat_clusters)]

mat_gexp <- mat_gexp[,cell_order]
df_annot <- data.frame(Cluster = df_annot[cell_order,])

rownames(mat_gexp) <- gsub("FCGR3A", "FCGR3A (CD16)",
                           gsub("CCL3", "CCL3 (MIP-1a)",
                                gsub("CCL4", "CCL4 (MIP-1b)",
                                     gsub("TNFRSF9", "TNFRSF9 (4-1BB)",
                                          gsub("TNFRSF18", "TNFRSF18 (GITR)",
                                               gsub("HAVRC2", "HAVCR2 (TIM-3)", 
                                                    gsub("CSF2", "CSF2 (GM-CSF)", rownames(mat_gexp))))))))

cols <- c("#2A363B", "#FF847C", "#99B898", "#E84A5F", "#96281B")
names(cols) <- c("Resting",
                 "Cytokine",
                 "Early activated",
                 "Activated",
                 "Type I IFN"
                 )

ha <- HeatmapAnnotation(df = df_annot,
                        col = list(#Pseudotime = colorRamp2(seq(min(df_annot$Pseudotime), max(df_annot$Pseudotime), length.out = 9), brewer.pal("RdPu", n = 9)),
                                   Cluster = cols),
                        show_legend = F,
                        annotation_name_side = "left",
                        annotation_legend_param = list(#Pseudotime = list(
                          # title_gp = gpar(fontsize = 10),
                          # labels_gp = gpar(fontsize = 10),
                          # grid_height = unit(0.2, "cm"),
                          # grid_width = unit(2, "mm"),
                          # title_position = "topcenter",
                          # legend_direction = "horizontal"),
                          Cluster = list(
                            title_gp = gpar(fontsize = 10),
                            labels_gp = gpar(fontsize = 10),
                            grid_height = unit(0.2, "cm"),
                            grid_width = unit(3, "mm"),
                            title_position = "topleft",
                            legend_direction = "horizontal",
                            nrow = 1)))

ht <- Heatmap(mat_gexp,
              col = colorRamp2(seq(quantile(mat_gexp, 0.05), quantile(mat_gexp, 0.95), length = 10), viridis::viridis_pal()(10)),
              #col = brewer.pal(name = "RdPu", n = 9),
              border_gp = gpar(col = "black", lty = 1),
              use_raster = T,
              top_annotation = ha,
              cluster_columns = F,
              cluster_rows = F,
              show_column_names = F,
              show_row_names = T,
              row_names_side = "left",
              row_names_gp = gpar(fontface = "italic"),
              #column_names_gp = gpar(fontsize = 5),
              #row_title_gp = gpar(fontsize = 5),
              show_heatmap_legend = T,
              raster_quality = 5,
              #raster_resize_mat = T,
              heatmap_legend_param = list(title = "Scaled expression",
                                          title_gp = gpar(fontsize = 10),
                                          labels_gp = gpar(fontsize = 10),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"),
                                          tick_length = unit(0, "mm"),
                                          border = "black",
                                          title_position = "topcenter",
                                          legend_direction = "horizontal"))

pdf("results/timepoints/hto_singlets/nk/complexheatmap_clusters.pdf", height = 5.5, width = 5.5)
draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE, annotation_legend_side = "bottom")
dev.off()


# Dot plot of top cluster DEGs
top10 <- all_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(hash_seurat_nk,
        features = unique(top10$gene),
        cols = "RdBu") + RotatedAxis()
ggsave("results/timepoints/hto_singlets/nk/cluster_deg_dotplot_labels.pdf", height = 4, width = 18)


# GSEA

# DEGs without thresholds
all_deg <- FindAllMarkers(hash_seurat_nk, test.use = "t", only.pos = F, logfc.threshold = 0, return.thresh = 1, min.pct = 0)
fwrite(all_deg, "results/timepoints/hto_singlets/nk/clusters_deg_clusternames_all.txt", sep = "\t", quote = F, row.names = F)


pathways_reactome <- gmtPathways("/Users/odufva/Documents/Labra/Scripts_programs/msigdb/msigdb_v7.0_GMTs/c2.cp.reactome.v7.0.symbols.gmt")
pathways_hallmark <- gmtPathways("/Users/odufva/Documents/Labra/Scripts_programs/msigdb/msigdb_v7.0_GMTs/h.all.v7.0.symbols.gmt")
pathways_biocarta <- gmtPathways("/Users/odufva/Documents/Labra/Scripts_programs/msigdb/msigdb_v7.0_GMTs/c2.cp.biocarta.v7.0.symbols.gmt")

dir.create("results/timepoints/hto_singlets/nk/gsea/")

# load DE results
data <- fread("results/timepoints/hto_singlets/nk/clusters_deg_clusternames_all.txt", data.table = F)


run_gsea <- function(DATA, NAME, CLUSTER, PATHWAYS){
  
  ranks <- DATA %>% 
    filter(cluster == CLUSTER) %>%
    dplyr::select(gene, avg_log2FC, p_val) %>%
    #mutate(p_val = ifelse(p_val == 0, 1e-305, p_val)) %>%
    #mutate(signed_p = -log10(p_val)*sign(avg_log2FC)) %>%
    #dplyr::select(gene, signed_p) %>%
    dplyr::select(gene, avg_log2FC) %>%
    tibble::deframe()
  
  # GSEA
  fgseaRes <- fgsea(PATHWAYS, ranks, maxSize=500, eps = 0)
  fgseaRes <- fgseaRes %>%
    mutate(pos_neg = ifelse(ES > 0, "pos", "neg")) %>%
    arrange(pos_neg, pval) %>%
    select(-pos_neg) %>% 
    mutate(cluster = CLUSTER)
  
  # table of top pathways
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  pdf(paste0("results/timepoints/hto_singlets/nk/gsea/", NAME, "_", CLUSTER, ".pdf"), height = 7, width = 20)
  plotGseaTable(PATHWAYS[topPathways], ranks, fgseaRes, gseaParam = 0.5)
  dev.off()
  
  return(fgseaRes)
}


clusters <- unique(data$cluster)

result_hallmark <- lapply(clusters, run_gsea, DATA = data, NAME = "HALLMARK", PATHWAYS = pathways_hallmark) %>% bind_rows()
fwrite(result_hallmark, "results/timepoints/hto_singlets/nk/gsea/hallmark.txt", quote = F, row.names = F, sep = "\t")

result_reactome <- lapply(clusters, run_gsea, DATA = data, NAME = "REACTOME", PATHWAYS = pathways_reactome) %>% bind_rows()
fwrite(result_reactome, "results/timepoints/hto_singlets/nk/gsea/reactome.txt", quote = F, row.names = F, sep = "\t")

result_biocarta <- lapply(clusters, run_gsea, DATA = data, NAME = "BIOCARTA", PATHWAYS = pathways_biocarta) %>% bind_rows()
fwrite(result_biocarta, "results/timepoints/hto_singlets/nk/gsea/biocarta.txt", quote = F, row.names = F, sep = "\t")



## Cluster lineplot over time

df_all <- as.data.frame(prop.table(table(hash_seurat_nk$hash.ID, hash_seurat_nk$seurat_clusters), 1))

df_coculture <- df_all %>% filter(!grepl("0h", Var1))
df_nk <- df_all %>% filter(Var1 %in% c("NK-0h-1", "NK-0h-2")) %>% group_by(Var2) %>% summarize(Freq = mean(Freq))

df_nk_k562 <- df_nk %>% mutate(Var1 = "K562-0h")
df_nk_jurkat <- df_nk %>% mutate(Var1 = "JURKAT-0h")
df_nk_gdm1 <- df_nk %>% mutate(Var1 = "GDM1-0h")
df_nk_ri1 <- df_nk %>% mutate(Var1 = "RI1-0h")

df <- rbind(df_coculture, df_nk_k562, df_nk_jurkat, df_nk_gdm1, df_nk_ri1)


df <- df %>% 
  mutate(cell_line = gsub("\\-.*", "", Var1),
         timepoint = gsub(".*-", "", Var1)) %>% 
  mutate(timepoint = factor(timepoint, levels = c("0h", "1h", "3h", "6h", "12h", "24h"))) %>% 
  mutate(cluster = factor(Var2)) %>% 
  mutate(cell_line = factor(cell_line, levels = c("K562", "JURKAT", "RI1", "GDM1"))) %>% 
  mutate(pct = Freq*100)

timepoints <- data.frame(
timepoint_continuous = c(0, 1, 3, 6, 12, 24),
timepoint = c("0h", "1h", "3h", "6h", "12h", "24h"))

df <- left_join(df, timepoints)

ggplot(df, aes(x = timepoint_continuous, y = pct, color = cluster, fill = cluster)) +
  geom_smooth(se = F, method = "loess", span = 0.6) +
  stat_smooth(se=FALSE, geom="area", method = 'loess', alpha = 0.25, span = 0.6) +
  scale_y_continuous(limits = c(0,107), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12, 24), expand = c(0,0)) +
  facet_wrap(. ~ cell_line, ncol = 1, switch = "both") +
  scale_color_manual("", values = a_palette[c(1,5,3,6,8)]) +
  scale_fill_manual("", values = a_palette[c(1,5,3,6,8)]) +
  ylab("%") +
  xlab("Time (h)") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "top",
        legend.justification = 0.5)+
  guides(fill = guide_legend(nrow = 2, byrow = T),
         color = guide_legend(nrow = 2, byrow = T))

ggsave("results/timepoints/hto_singlets/nk/lineplot_clusters.pdf", width = 6, height = 6.5)



