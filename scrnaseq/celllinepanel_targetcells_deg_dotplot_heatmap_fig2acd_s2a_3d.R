
# Plot target cell DEG dotplot and heatmap (Figures 2A, 2C, 2D, and S2A)


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
library(tidyr)
library(ComplexHeatmap)
library(circlize)


hash_seurat_targets_combined <- readRDS("results/celllinepanel_combined/targets/hash_seurat_targets_combined_mainclusters.rds")


# target cells

# FindMarkers t-test

# load DE results
data <- fread("results/celllinepanel_combined/targets/deg/deg_mainclusters_all.txt", data.table = F)

data %>% filter(p_val < 0.05) %>% group_by(cell_line) %>% summarize(count = n()) %>% summarize(mean_count = mean(count))
data %>% filter(p_val_adj < 0.05) %>% group_by(cell_line) %>% summarize(count = n()) %>% summarize(mean_count = mean(count))


# define core NK-induced genes
core_genes <- data %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC > 0.25) %>% 
  group_by(gene) %>% 
  summarize(cell_line_count = n()) %>% 
  arrange(desc(cell_line_count)) %>% 
  filter(cell_line_count > 12) %>%
  select(gene) %>%
  tibble::deframe()


# heatmap of most variable genes
topvargenes <- data %>%
  filter(avg_log2FC > 0) %>% 
  filter(p_val_adj < 0.05) %>% 
  group_by(gene) %>% summarize(variance = var(avg_log2FC)) %>% mutate(variance = ifelse(is.na(variance), 0, variance)) %>% arrange(desc(variance)) %>%
  top_n(200, variance) %>% dplyr::select(gene) %>% tibble::deframe()

df <- data %>% 
  filter(gene %in% topvargenes) %>% 
  mutate(gene = factor(gene, levels = topvargenes))

mat <- df %>% dplyr::select(gene, cell_line, avg_log2FC) %>% pivot_wider(names_from = cell_line, values_from = avg_log2FC) %>% as.data.frame()
rownames(mat) <- mat$gene
mat$gene <- NULL
mat <- as.matrix(mat)
mat[is.na(mat)] <- 0

# row clusters and annotation
split = data.frame(cluster = cutree(hclust(dist(mat, method = "euclidean"), method = "ward.D2"), k = 8))#,

# label clusters
split$title <- "IFNy 1"
split$title[split$cluster==2] <- "IFNy 4"
split$title[split$cluster==3] <- "IFNy 2"
split$title[split$cluster==4] <- "IFNy 3"
split$title[split$cluster==5] <- "IFNy/\nType I IFN"
split$title[split$cluster==6] <- "HLA II/\nMyeloid"
split$title[split$cluster==7] <- "Cell line-\nspecific"
split$title[split$cluster==8] <- "Type I IFN"

split$title <- factor(split$title, levels = c("IFNy 1", "IFNy 2",  "IFNy 3", "IFNy 4", "IFNy/\nType I IFN", "Type I IFN", "HLA II/\nMyeloid", "Cell line-\nspecific"))

row_ha <- rowAnnotation(Cluster = split$title,
                        col = list(Cluster = structure(brewer.pal(8, "BrBG"), names =  c("IFNy 1", "IFNy 2",  "IFNy 3", "IFNy 4", "IFNy/\nType I IFN", "Type I IFN", "HLA II/\nMyeloid", "Cell line-\nspecific"))),
                        show_annotation_name = F,
                        show_legend = F,
                        na_col = "white")


# column annotation

# read cancer types and colors
cols <- fread("data/nk_crispr_colors.txt", data.table = F) %>% dplyr::rename(cancer = cancer_type)
cols <- rbind(cols, data.frame(color = c("lightblue", "grey70"), cancer = c("CML", "NK cells only")))

cancertypes <- readxl::read_excel("cell_line_subtypes.xlsx")

cancertype_color <- merge(cancertypes, cols)

cancertype_color_df <- cancertype_color %>% 
  filter(cell_line %in% colnames(mat)) %>% 
  mutate(gene = "Cancer type") %>% 
  dplyr::select(gene, cell_line, cancer) %>% 
  mutate(cell_line = factor(cell_line, levels = colnames(mat))) %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  mutate(cancer_main = factor(cancer_main, levels = c("T-ALL", "AML/CML", "BCL", "B-ALL", "MM"))) %>% 
  mutate(cancer_color = gsub("AML/CML", "AML", gsub("BCL", "DLBCL", cancer_main)))

cancertype_color_vector <- c("#3C899E", "#E41A1C", "#9F5196", "#FFA60F", "#DF6F32")
names(cancertype_color_vector) <- c("AML/CML", "T-ALL", "BCL", "MM", "B-ALL") 

cancertype_color_df <- cancertype_color_df[match(colnames(mat), cancertype_color_df$cell_line),]
cancertype_color_df$cancer_main <- factor(cancertype_color_df$cancer_main, levels = names(cancertype_color_vector))

ha <- HeatmapAnnotation(`Cancer type` = cancertype_color_df$cancer_main,
                        col = list(`Cancer type` = cancertype_color_vector),
                        annotation_legend_param =  list(title_gp = gpar(fonttype = "regular"),
                                                        direction = "horizontal",
                                                        title_position = "topcenter"))

# genes to label

core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")
core_subset <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1")
hla2 <- c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1")
myeloid <- c("LYZ", "S100A8", "S100A9", "NCF1", "CTSZ")
ligands <- c("TNFRSF1B", "LGALS9", "CLEC2B", "ICAM1", "BTN3A2", "CALR")
cytokines <- c("TNFSF10", "CXCL8", "CXCL9", "CXCL10", "LTB")
type1ifn <- c("MX2", "OAS1", "OAS2", "OAS3", "IRF7", "IFIT1", "IFI44",
              "IFI27", "ISG15")
specific <- c("IDO1", "CCL2", "CCL3", "TGFB1")


plotgenes <- c(core_subset, hla2, myeloid, ligands, cytokines, type1ifn, specific)

genelabels <- rowAnnotation(genelabels = anno_mark(at = which(rownames(mat) %in% plotgenes), labels = rownames(mat)[rownames(mat) %in% plotgenes],
                                                   labels_gp = gpar(fontsize = 10, fontface = "italic"),
                                                   link_width = unit(20, "mm"))
)

ht <- Heatmap(mat,
              name = "Average\nlog2FC\nNK-treated vs.\nuntreated",
              row_split = split$title,
              row_title_rot = 0,
              cluster_row_slices = FALSE, 
              col = colorRamp2(seq(-quantile(mat, 0.975), quantile(mat, 0.975), length.out = 11), rev(brewer.pal(11, "RdBu"))),
              top_annotation = ha,
              left_annotation = row_ha,
              show_row_names = T,
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D2",
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "ward.D2",
              row_dend_width = unit(3, "cm"),
              heatmap_legend_param = list(title_gp = gpar(fonttype = "regular"),
                                          direction = "horizontal",
                                          title_position = "topcenter"))



pdf("results/celllinepanel_combined/targets/targets_deg_heatmap_top200variance.pdf", width = 10, height = 30)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)
dev.off()

# selected genes labeled (Figure S2A)
ht <- Heatmap(mat,
              name = "Average\nlog2FC\nNK-treated vs.\nuntreated",
              row_split = split$title,
              row_title_rot = 0,
              cluster_row_slices = FALSE, 
              col = colorRamp2(seq(-quantile(mat, 0.975), quantile(mat, 0.975), length.out = 11), rev(brewer.pal(11, "RdBu"))),
              top_annotation = ha,
              left_annotation = row_ha,
              right_annotation = genelabels,
              show_row_names = F,
              row_names_gp = gpar(fontsize = 4),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D2",
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "ward.D2",
              row_dend_width = unit(3, "cm"),
              heatmap_legend_param = list(title_gp = gpar(fonttype = "regular"),
                                          direction = "horizontal",
                                          title_position = "topcenter"))



pdf("results/celllinepanel_combined/targets/targets_deg_heatmap_top200variance_small.pdf", width = 9, height = 12)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = T)
dev.off()



## ------------------------------


# Dot plot (Figure 2A)
# genes filtered by p value

core_genes <- data %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC > 0.25) %>% 
  group_by(gene) %>% 
  summarize(cell_line_count = n()) %>% 
  arrange(desc(cell_line_count)) %>% 
  filter(cell_line_count > 13) %>%
  select(gene) %>%
  tibble::deframe()



# selected genes 
core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")
hla2 <- c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1")
myeloid <- c("LYZ", "NCF1", "CTSZ")
ligands <- c("LGALS9")
cytokines <- c("CXCL8", "CXCL10", "TNFSF10")
type1ifn <- c("MX2", "OAS1", "OAS2", "IRF7")

plotgenes <- c(core, hla2, myeloid, ligands, cytokines, type1ifn)

nkresponse <- fread("results/celllinepanel_combined/targets/targets_differential_modulescore_mainclusters.txt", data.table = F)

cancertypes <- readxl::read_excel("data/cell_line_subtypes.xlsx")

cellline_order <- nkresponse %>% 
  left_join(cancertypes) %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  group_by(cancer_main) %>% 
  mutate(cancer_mean = mean(avg_log2FC)) %>% 
  arrange(desc(cancer_mean), desc(avg_log2FC)) %>% 
  dplyr::select(cell_line) %>%
  tibble::deframe()

df_plot <- data %>% 
  filter(abs(avg_log2FC) > 0.1) %>%
  filter(gene %in% plotgenes) %>% 
  mutate(gene = factor(gene, levels = plotgenes)) %>% 
  mutate(cell_line = factor(cell_line, levels = cellline_order)) %>% 
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-310, p_val_adj),
         p_val = ifelse(p_val == 0, 1e-310, p_val))


# read GSEA results
gsea <- fread("results/celllinepanel_combined/gsea/hallmark.txt", data.table = F)

gsea_plot <- gsea %>% mutate(pathway_sentenceCase = gsub("_", " ",sub(".*?_", "", pathway))) %>%
  filter(padj < 0.05) %>%
  filter(pathway %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE")) %>% 
  mutate(log10_fdr = -log10(padj)) %>%
  mutate(celL_line = factor(cell_line, levels = cellline_order)) %>% 
  dplyr::rename(gene = pathway, p_val = pval, p_val_adj = padj, avg_log2FC = NES) %>% 
  mutate(pct.1 = NA, pct.2 = NA, p_adj = NA) %>% 
  dplyr::select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cell_line, p_adj)


df <- rbind(df_plot, gsea_plot) %>% 
  mutate(gene_prettyname = gsub("Il6 jak stat3", "IL6-JAK-STAT3",
                                gsub("Il2 stat5", "IL2-STAT5",
                                     gsub("Dna", "DNA",
                                          gsub("Myc", "MYC",
                                               gsub("E2f", "E2F",
                                                    gsub("Kras", "KRAS",
                                                         gsub("Tnfa", "TNFA",
                                                              gsub("G2m", "G2M",
                                                                   gsub("Nf-kb|nfkb", "NF-kB",
                                                                        gsub("Ifng", "IFNy response",
                                                                             gsub("Nkresponse", "NK cell response",
                                                                                  stringr::str_to_sentence(gsub("_", " ",
                                                                                                                sub(".*?_", "", 
                                                                                                                    gsub("nfkb", "NF-kB response", gene)))))))))))))))) %>%
  mutate(gene_prettyname = ifelse(grepl("HALLMARK", gene) | grepl("response", gene_prettyname), gene_prettyname, gene))


p_dotplot <- ggplot(df, aes(x = cell_line, y = gene, size = -log10(p_val_adj), fill = avg_log2FC)) +
  geom_point(data = df[!grepl("Cancer|HALLMARK", df$gene),], pch = 21, color = "white") +
  scale_fill_distiller("Fold\nchange\n(log2)", palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = max(abs(as.numeric(df$avg_log2FC[!grepl("HALLMARK", df$gene)]))) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top")) +
  geom_point(data = df[!grepl("Cancer|HALLMARK", df$gene),], pch = 1, color = ifelse(df[!grepl("Cancer|HALLMARK", df$gene),]$p_val_adj < 0.05, "grey50", "white")) +
  scale_size("FDR\n(-log10)", range = c(2, 4.5), guide = guide_legend(title.position = "top",
                                                                      override.aes = list(pch = 16))) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "bottom") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  ylab("") +
  xlab("")

p_dotplot_gsea <- ggplot(df[grepl("HALLMARK", df$gene),], aes(x = cell_line, y = gene_prettyname, size = -log10(p_val_adj), fill = avg_log2FC)) +
  geom_point(pch = 21, color = "white") +
  scale_fill_gradientn("NES", colors = pals::ocean.thermal(9)[2:8],
                       limits = max(abs(as.numeric(df$avg_log2FC[grepl("HALLMARK", df$gene)]))) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top")) +
  geom_point(pch = 1, color = ifelse(df[grepl("HALLMARK", df$gene),]$p_val_adj < 0.05, "grey50", "white")) +
  scale_size("FDR\n(-log10)", range = c(1, 4.5), guide = guide_legend(title.position = "top",
                                                                      override.aes = list(pch = 16))) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "bottom") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  ylab("") +
  xlab("")


# cancer type color dots
# read cancer types and colors
cols <- fread("data/nk_crispr_colors.txt", data.table = F) %>% dplyr::rename(cancer = cancer_type)
cols <- rbind(cols, data.frame(color = "lightblue", cancer = "CML"))

cancertypes <- readxl::read_excel("data/cell_line_subtypes.xlsx")

cancertype_color <- merge(cancertypes, cols)

cancertype_color_df <- cancertype_color %>% 
  filter(cell_line %in% cellline_order) %>% 
  mutate(gene = "Cancer type") %>% 
  dplyr::select(gene, cell_line, cancer) %>% 
  mutate(cell_line = factor(cell_line, levels = cellline_order))  %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  mutate(cancer_main = factor(cancer_main, levels = c("T-ALL", "AML/CML", "BCL", "B-ALL", "MM"))) %>% 
  mutate(cancer_color = gsub("BCL", "DLBCL", cancer_main))

cancertype_color_vector <- cancertype_color %>% filter(cell_line %in% cellline_order) %>% dplyr::select(cancer, color) %>% unique() %>% dplyr::select(color) %>% tibble::deframe()
names(cancertype_color_vector) <- cancertype_color %>% filter(cell_line %in% cellline_order) %>% dplyr::select(cancer, color) %>% unique() %>% dplyr::select(cancer) %>% tibble::deframe()

cancertype_color_vector <- cancertype_color_vector[c("T-ALL", "AML", "DLBCL", "MM", "B-ALL")]
names(cancertype_color_vector) <- c("T-ALL", "AML/CML", "DLBCL", "MM", "B-ALL")

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
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))



# differential core NK response barplot
nkresponse_cancertypes <- merge(cancertype_color_df, nkresponse, by = "cell_line")

nkresponse_cancertypes <- nkresponse_cancertypes %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  mutate(cancer_main = factor(cancer_main, levels = c("T-ALL", "AML/CML", "BCL", "B-ALL", "MM"))) %>% 
  mutate(cancer_color = gsub("BCL", "DLBCL", cancer_main))


p_nkresponse <- ggplot(nkresponse_cancertypes, aes(y = avg_log2FC, x = cell_line, fill = cancer_color)) +
  geom_col() +
  theme_cowplot() +
  scale_fill_manual("Cancer type", values = cancertype_color_vector,
                    guide = guide_legend(title.position = "top")) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1, 1.5)) +
  xlab("") +
  ylab("") +
  labs(fill = "Cancer type") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text = element_text(size = 7)) +
  guides(fill = "none")


# IFNB1 bar plot

ifnb1 <- FetchData(hash_seurat_targets_combined, vars = c("hash.ID", "IFNB1"))

ifnb1_freq <- ifnb1 %>%
  mutate(IFNB1_expressed = ifelse(IFNB1>0, 1, 0)) %>% 
  mutate(cell_line = gsub("-.*", "", hash.ID),
         donor = ifelse(grepl("NK", hash.ID), gsub("expanded", "NK1", gsub(".*-", "", hash.ID)), "Untreated"),
         treatment = ifelse(grepl("NK", hash.ID), "NK-treated", "Untreated")) %>%
  filter(!grepl("NK", cell_line), !is.na(cell_line)) %>%
  group_by(treatment, donor, cell_line) %>%
  summarize(sum = sum(IFNB1_expressed), count = n()) %>%
  mutate(freq = sum / count) %>% 
  mutate(pct = 100*freq) %>% 
  group_by(treatment, cell_line) %>%
  mutate(mean_pct = mean(pct))


ifnb1_freq_plot <- ifnb1_freq %>% mutate(cell_line = factor(cell_line, levels = cellline_order)) %>% filter(treatment == "NK-treated")

p_ifnb1 <- ggplot(ifnb1_freq_plot, aes(x = cell_line)) +
  geom_col(data = ifnb1_freq_plot[ifnb1_freq_plot$donor=="NK1",], aes(y = mean_pct), fill = "grey40") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2), breaks = c(0,0.5,1)) +
  ylab("") +
  xlab("")

# combined barplot + dotplot (Figure 2A)
(p_nkresponse / p_dotplot / p_dotplot_gsea / p_ifnb1 / p_cancertypedots & theme(plot.margin = unit(c(1,0,0,0), "mm"))) + plot_layout(heights = c(0.5,7,0.35,0.45,0.2), guides = "collect")
ggsave("results/celllinepanel_combined/targets/targets_deg_dotplot_manuscript_selectedgenes.pdf", width = 7, height = 7.75)



# Heatmap of cytokine log2FC (Figure 2A)

# top 10 cytokines correlating with core NK response
result_cor_cytokines <- fread("data/cytokines_correlation_corenkresponse.txt")

top10 <- result_cor_cytokines %>% slice_min(p, n = 10) %>% select(cytokine) %>% tibble::deframe()


# cytokine NK treatment lof2FC (with NK baseline added)
result_all_withnk <- fread("../NK_resistance Heme CollabPaper/Analysis/Cytokine_analysis/cytokines_differential_nk_added.txt")

mat_cytokine <- result_all_withnk %>% filter(cytokine %in% top10) %>% pivot_wider(names_from = cell_line, values_from = log2fc, id_cols = cytokine) %>% as.data.frame()
rownames(mat_cytokine) <- mat_cytokine$cytokine
mat_cytokine$cytokine <- NULL
mat_cytokine <- mat_cytokine[top10,cellline_order]

mat_cytokine_scaled <- t(apply(mat_cytokine, 1, scale))
rownames(mat_cytokine_scaled) <- rownames(mat_cytokine)
colnames(mat_cytokine_scaled) <- colnames(mat_cytokine)

ht <- Heatmap(mat_cytokine_scaled,
        cluster_rows = T,
        cluster_columns = F,
        border_gp = gpar(col = "black", lty = 1),
        col = pals::ocean.deep(9),
        show_column_names = F,
        show_row_names = T,
        row_names_side = "left",
        show_heatmap_legend = T,
        raster_quality = 3,
        heatmap_legend_param = list(title = "Scaled\nconcentration\nlog2 fold change",
                                    title_gp = gpar(fontsize = 10),
                                    labels_gp = gpar(fontsize = 10),
                                    grid_height = unit(0.2, "cm"),
                                    grid_width = unit(2, "mm"),
                                    tick_length = unit(0, "mm"),
                                    border = "black",
                                    title_position = "topleft",
                                    legend_direction = "vertical"))

pdf("results/celllinepanel_combined/targets/cytokines_heatmap_targetdotplot_manuscript.pdf", width = 9, height = 2.25)
ht
dev.off()


## -------------------------------------

## Boxplot of core NK response across cancer types (Figure)

nkresponse_maincancertypes <- nkresponse_cancertypes %>% 
  mutate(cancer_main = gsub("AML|CML", "AML/CML", gsub("MCL|DLBCL", "BCL", cancer))) %>% 
  mutate(cancer_main = factor(cancer_main, levels = c("T-ALL", "AML/CML", "BCL", "B-ALL", "MM")))

comparisons = list(c("BCL", "T-ALL"), c("BCL", "AML/CML"), c("BCL", "MM"), c("BCL", "B-ALL"),
                   c("T-ALL", "AML/CML"), c("T-ALL", "MM"), c("T-ALL", "B-ALL"),
                   c("AML/CML", "MM"),  c("AML/CML", "B-ALL"),
                   c("MM", "B-ALL"))

# test data normality
shapiro.test(nkresponse_maincancertypes$avg_log2FC)

test_wilcox_pairwise <- ggpubr::compare_means(avg_log2FC ~ cancer_main, comparisons = comparisons, p.adjust.method = "BH", method = 't.test', data = nkresponse_maincancertypes)
test_wilcox_pairwise <- test_wilcox_pairwise %>% arrange(p) %>% mutate(y.position = c(1.4, 1.55, rep(NA, 8)),
                                                                       p.signif = ifelse(p < 0.05, "*", "ns"))

colors_inducedgenes <- c("AML/CML" = "#3C899E", "T-ALL" = "#E41A1C", "BCL" = "#9F5196", "MM" = "#FFA60F", "B-ALL" = "#DF6F32")

ggplot(nkresponse_maincancertypes, aes(x = cancer_main, y = avg_log2FC)) +
  geom_boxplot(aes(fill = cancer_main), outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  theme_cowplot() +
  scale_fill_manual("Cancer type", values = colors_inducedgenes,
                    guide = guide_legend(title.position = "top")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.7), breaks = c(0,0.5,1,1.5)) +
  xlab("") +
  ylab("Core NK response signature\naverage log2 fold change") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.y = 1, label.x = "B-ALL", method = "kruskal.test") +
  ggpubr::stat_pvalue_manual(test_wilcox_pairwise, label = "p.signif")

ggsave("results/celllinepanel_combined/targets/targets_nkresponse_boxplot.pdf", width = 2.5, height = 3.75)



# Correlation of core NK response with PRISM AUC and NK cell cluster fraction (Figures 2D and 3D)

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

# cluster fraction correlation with PRISM AUC
nk_prism_auc <- fread("data/NK_PRISM_heme_pctviability_auc_norm.txt")
nk_prism_fm_auc <- nk_prism_auc %>%
  filter(Condition == "0") %>%
  dplyr::select(auc_norm, cell_line) %>% 
  as.data.frame()

cancertype_color_vector <- cancertype_color_vector[c("AML/CML", "T-ALL", "DLBCL", "MM", "B-ALL")]
names(cancertype_color_vector) <- gsub("DLBCL", "BCL", names(cancertype_color_vector))


cluster_freq_pert_prism <- cluster_freq_pert %>%
  dplyr::rename(cell_line = condition, cluster = Var2) %>%
  filter(cluster %in% c("Activated (2)", "Type I IFN (3)")) %>%
  group_by(cell_line) %>% 
  summarize(freq_sum = sum(mean_freq)) %>% 
  left_join(nk_prism_fm_auc) %>%
  left_join(cancertype_color_df)


ggpubr::ggscatter(cluster_freq_pert_prism, y = "auc_norm", x = "freq_sum", shape = 21, size = 3, fill = "cancer_main",
                  add = "reg.line", cor.coef = T, cor.coeff.args = list(label.sep = "\n"), cor.method = "spearman", palette = cancertype_color_vector,#, label = "cell_line", repel = T,
                  font.label = c(10, "plain"), cor.coef.coord = c(50, 1.1)) +
  xlab("% of NK cells in activated/\ntype I IFN clusters") +
  ylab("PRISM AUC") +
  theme(legend.title = element_blank(),
        plot.margin = unit(c(0,0.3,0,0), "cm")) +
  guides(fill="none")
ggsave("results/celllinepanel_integrated/cluster_fraction_prism.pdf", height = 3, width = 3.3)


# cluster fraction correlation with target core NK response (Figure 2D)
nkresponse_maincancertypes_cluster_freq <- nkresponse_maincancertypes %>% 
  left_join(cluster_freq_pert_prism[,c("cell_line", "freq_sum")], by = "cell_line")

colors_nkresponse <- c("AML/CML" = "#3C899E", "T-ALL" = "#E41A1C", "BCL" = "#9F5196", "MM" = "#FFA60F", "B-ALL" = "#DF6F32")

ggpubr::ggscatter(nkresponse_maincancertypes_cluster_freq, y = "avg_log2FC", x = "freq_sum", shape = 21, size = 3, fill = "cancer_main",
                  add = "reg.line", cor.coef = T, cor.method = "spearman", palette = colors_nkresponse,
                  font.label = c(10, "plain"), cor.coef.coord = c(50, 1.5), cor.coeff.args = list(label.sep = c("\n"))) +
  xlab("% of NK cells in activated/\ntype I IFN clusters") +
  ylab("Core NK-induced signature\naverage log2 fold change") +
  theme(legend.title = element_blank(),
        plot.margin = unit(c(0,0.3,0,0), "cm")) +
  guides(fill = "none")
ggsave("results/celllinepanel_combined/targets/cluster_fraction_corenk.pdf", height = 3, width = 3.3)


# Correlation of core NK response with PRISM AUC (Figure 3D)
nkresponse_maincancertypes <- nkresponse_maincancertypes %>%
  left_join(nk_prism_fm_auc)

ggpubr::ggscatter(nkresponse_maincancertypes, x = "auc_norm", y = "avg_log2FC", shape = 21, size = 3, fill = "cancer_main",
                  add = "reg.line", cor.coef = T, cor.method = "spearman", palette = colors_nkresponse, label = "cell_line", repel = T) +
  ylab("Core NK-induced signature\naverage log2 fold change") +
  xlab("PRISM AUC") +
  theme(legend.title = element_blank())
ggsave("results/celllinepanel_combined/targets/targets_prism_corenk_manuscript.pdf", height = 5, width = 5)
