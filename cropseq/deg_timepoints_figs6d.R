
# Compare DEG in MM1S timepoint CROP-seq

library(data.table)
library(dplyr)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(emdbook)

# load DEG lists
result_nonk_24h <- fread("results/mm1s/deg/singlet/pert_deg_all.txt") %>%
  filter(condition %in% c("no NK", "NK 1:4")) %>% 
  mutate(condition = gsub("NK 1:4", "NK 1:4 24h", condition))
result_nk_3h_6h <- fread("results/mm1s/timepoints/deg/singlet/pert_deg.txt")

data_all <- rbind(result_nonk_24h, result_nk_3h_6h)

data <- data_all %>% filter(p_val_adj < 0.05, avg_log2FC > 0.1)

count_deg <- function(PERTURBATION){
  
  data_pert <- data %>% filter(perturbation == PERTURBATION)
  
  genes_all <- unique(data_pert$gene)
  genes_nonk <- data_pert$gene[data_pert$condition == "no NK"]
  genes_3h <- data_pert$gene[data_pert$condition == "NK 1:4 3h"]
  genes_6h <- data_pert$gene[data_pert$condition == "NK 1:4 6h"]
  genes_24h <- data_pert$gene[data_pert$condition == "NK 1:4 24h"]
  
  genes_nonk_unique <- genes_nonk[!genes_nonk %in% c(genes_3h, genes_6h, genes_24h)]
  genes_3h_unique <- genes_3h[!genes_3h %in% c(genes_nonk, genes_6h, genes_24h)]
  genes_6h_unique <- genes_6h[!genes_6h %in% c(genes_nonk, genes_3h, genes_24h)]
  genes_24h_unique <- genes_24h[!genes_24h %in% c(genes_nonk, genes_3h, genes_6h)]
  
  # gene counts
  count_nonk <- length(genes_nonk)
  count_3h <- length(genes_3h)
  count_6h <- length(genes_6h)
  count_24h <- length(genes_24h)
  
  count_nonk_unique <- length(genes_nonk_unique)
  count_3h_unique <- length(genes_3h_unique)
  count_6h_unique <- length(genes_6h_unique)
  count_24h_unique <- length(genes_24h_unique)
  
  
  df <- data.frame(perturbation = PERTURBATION,
                   condition = c("no NK", "3h", "6h", "24h"),
                   all_count = c(count_nonk, count_3h, count_6h, count_24h),
                   unique_count = c(count_nonk_unique, count_3h_unique, count_6h_unique, count_24h_unique))
  df$common_count <- df$all_count-df$unique_count
  
  return(df)
  
}

perturbations <- unique(data_all$perturbation)
counts <- lapply(perturbations, count_deg) %>% bind_rows()

# list unique DEG

unique_deg <- function(PERTURBATION){
  
  data_pert <- data %>% filter(perturbation == PERTURBATION)
  
  genes_all <- unique(data_pert$gene)
  genes_nonk <- data_pert$gene[data_pert$condition == "no NK"]
  genes_3h <- data_pert$gene[data_pert$condition == "NK 1:4 3h"]
  genes_6h <- data_pert$gene[data_pert$condition == "NK 1:4 6h"]
  genes_24h <- data_pert$gene[data_pert$condition == "NK 1:4 24h"]
  
  genes_nonk_unique <- genes_nonk[!genes_nonk %in% c(genes_3h, genes_6h, genes_24h)]
  genes_3h_unique <- genes_3h[!genes_3h %in% c(genes_nonk, genes_6h, genes_24h)]
  genes_6h_unique <- genes_6h[!genes_6h %in% c(genes_nonk, genes_3h, genes_24h)]
  genes_24h_unique <- genes_24h[!genes_24h %in% c(genes_nonk, genes_3h, genes_6h)]
  genes_3h_6h_unique <- genes_3h[!genes_3h %in% c(genes_nonk, genes_24h)]
  
  list("no NK" = genes_nonk_unique, "3h" = genes_3h_unique, "6h" = genes_6h_unique, "24h" = genes_24h_unique, "3h_6h" = genes_3h_6h_unique)
}

unique_deg("JAK1")
unique_deg("JAK2")
unique_deg("STAT1")
unique_deg("IFNGR2")
unique_deg("PCGF5")


# Bar plot

# pivot longer
counts_df <- counts %>%
  tidyr::pivot_longer(cols = contains("count"), names_to = "category", values_to = "count") %>%
  filter(category != "all_count") %>% 
  mutate(condition = factor(condition, levels = c("no NK", "3h", "6h", "24h")),
         category = factor(category, levels = c("unique_count", "common_count"), labels = c("Unique to one timepoint", "Common across timepoints")))


ggplot(counts_df, aes(x = condition, y = count, fill = category)) +
  geom_col() +
  facet_grid(cols = vars(perturbation)) +
  ylab("Number of DEG") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0,0))


# matrix of all DEG counts
count_mat <- counts %>% select(-unique_count, -common_count) %>% tidyr::pivot_wider(names_from = condition, values_from = all_count) %>% as.data.frame()
rownames(count_mat) <- count_mat$perturbation
count_mat$perturbation <- NULL
count_mat <- as.matrix(count_mat)

# order by max number of DEGs
pert_order <- rownames(count_mat)[rev(order(rowMaxs(count_mat)))]
count_mat <- count_mat[pert_order,]


# matrix of common DEG counts
count_mat_common <- counts %>% select(-unique_count, -all_count) %>% tidyr::pivot_wider(names_from = condition, values_from = common_count) %>% as.data.frame()
rownames(count_mat_common) <- count_mat_common$perturbation
count_mat_common$perturbation <- NULL
count_mat_common <- as.matrix(count_mat_common)

# order by max number of DEGs
count_mat_common <- count_mat_common[pert_order,]

# fraction of common DEG
count_mat_common_fraction <- count_mat_common/count_mat

# remove NaN
is.nan.data.frame <- function(x){ do.call(cbind, lapply(x, is.nan)) }
count_mat_common_fraction[is.nan(count_mat_common_fraction)] <- 0

# Heatmap
ha <- HeatmapAnnotation(`Median` = anno_barplot(colMedians(count_mat)),
                        `Total` = anno_barplot(colSums(count_mat)),
                        annotation_name_side = "left",
                        annotation_name_rot = 0,
                        gap = unit(0.25, "cm"))


# top5 genes unique to 3-6h and 24 h
plotgenes <- c("IRF1", "PARP9", "GBP1", "AIM2", "PARP14", "HLA-A", "HLA-B", "HLA-F", "HLA-C" ,"PSME2")

mat_nonk <- data %>% mutate(gene = factor(gene, levels = plotgenes), perturbation = factor(perturbation, levels = pert_order)) %>% filter(gene %in% plotgenes, condition == "no NK") %>% tidyr::pivot_wider(names_from = gene, id_cols = perturbation, values_from = avg_log2FC, names_expand = T, id_expand = T) %>% as.data.frame()
rownames(mat_nonk) <- mat_nonk$perturbation
mat_nonk$perturbation <- NULL
mat_nonk[,plotgenes]

mat_3h <- data %>% mutate(gene = factor(gene, levels = plotgenes), perturbation = factor(perturbation, levels = pert_order)) %>% filter(gene %in% plotgenes, condition == "NK 1:4 3h") %>% tidyr::pivot_wider(names_from = gene, id_cols = perturbation, values_from = avg_log2FC, names_expand = T, id_expand = T) %>% as.data.frame()
rownames(mat_3h) <- mat_3h$perturbation
mat_3h$perturbation <- NULL
mat_3h[,plotgenes]

mat_6h <- data %>% mutate(gene = factor(gene, levels = plotgenes), perturbation = factor(perturbation, levels = pert_order)) %>% filter(gene %in% plotgenes, condition == "NK 1:4 6h") %>% tidyr::pivot_wider(names_from = gene, id_cols = perturbation, values_from = avg_log2FC, names_expand = T, id_expand = T) %>% as.data.frame()
rownames(mat_6h) <- mat_6h$perturbation
mat_6h$perturbation <- NULL
mat_6h[,plotgenes]

mat_24h <- data %>% mutate(gene = factor(gene, levels = plotgenes), perturbation = factor(perturbation, levels = pert_order)) %>% filter(gene %in% plotgenes, condition == "NK 1:4 24h") %>% tidyr::pivot_wider(names_from = gene, id_cols = perturbation, values_from = avg_log2FC, names_expand = T, id_expand = T) %>% as.data.frame()
rownames(mat_24h) <- mat_24h$perturbation
mat_24h$perturbation <- NULL
mat_24h[,plotgenes]



ht_counts <-  Heatmap(count_mat,
                      column_title = "Count of\nDEG",
                     top_annotation = ha,
                     cluster_columns = F,
                     cluster_rows = F,
                     row_names_side = "left",
                     row_names_gp = gpar(fontface = "italic"),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if(!is.na(count_mat[i, j]))
                         grid.text( count_mat[i, j], x, y, gp = gpar(fontsize = 10))},
                     na_col = "grey95",
                     col = colorRamp2(c(0, lseq(1, max(count_mat, na.rm = T), length.out = 7)), c("white", brewer.pal(name = "Reds", n = 9)[c(2:8)])),
                     show_heatmap_legend = T,
                     heatmap_legend_param = list(title = "Count of\nDEG",
                                                 title_gp = gpar(fontface = "plain"),
                                                 grid_height = unit(0.4, "cm"),
                                                 grid_width = unit(2, "mm"),
                                                 tick_length = unit(0, "mm"),
                                                 border = "black",
                                                 title_position = "topcenter",
                                                 legend_direction = "horizontal"),
                     border_gp = gpar(col = "black", lty = 1, lwd = 1))


ht_common_fraction <-Heatmap(count_mat_common_fraction,
                             column_title = "Fraction of\ncommon DEG",
                     cluster_columns = F,
                     cluster_rows = F,
                     row_names_side = "left",
                     row_names_gp = gpar(fontface = "italic"),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if(!is.na(count_mat_common_fraction[i, j]))
                         grid.text( signif(count_mat_common_fraction[i, j], digits = 2), x, y, gp = gpar(fontsize = 10))},
                     na_col = "grey95",
                     col = colorRamp2(c(0, seq(0, max(count_mat_common_fraction, na.rm = T), length.out = 7)), c("white", brewer.pal(name = "Oranges", n = 9)[c(2:8)])),
                     show_heatmap_legend = T,
                     heatmap_legend_param = list(title = "Fraction of\ncommon DEG",
                                                 title_gp = gpar(fontface = "plain"),
                                                 grid_height = unit(0.4, "cm"),
                                                 grid_width = unit(2, "mm"),
                                                 tick_length = unit(0, "mm"),
                                                 border = "black",
                                                 title_position = "topcenter",
                                                 legend_direction = "horizontal"),
                     border_gp = gpar(col = "black", lty = 1, lwd = 1),
                     width = 7)

color <- colorRamp2(seq(-0.75, 0.75, length.out = 8), rev(brewer.pal(name = "RdBu", n = 9)[c(1:8)]))

ht_nonk <- Heatmap(mat_nonk,
                  column_title = "no NK",
                  cluster_columns = F,
                  cluster_rows = F,
                  row_names_side = "left",
                  row_names_gp = gpar(fontface = "italic"),
                  na_col = "grey95",
                  col = color,
                  show_heatmap_legend = F,
                  heatmap_legend_param = list(title = "Average log2 fold change",
                                              title_gp = gpar(fontface = "plain"),
                                              grid_height = unit(0.4, "cm"),
                                              grid_width = unit(2, "mm"),
                                              tick_length = unit(0, "mm"),
                                              border = "black",
                                              title_position = "topcenter",
                                              legend_direction = "horizontal"),
                  border_gp = gpar(col = "black", lty = 1, lwd = 1))

ht_3h <- Heatmap(mat_3h,
                 column_title = "3h",
                   cluster_columns = F,
                   cluster_rows = F,
                   row_names_side = "left",
                   row_names_gp = gpar(fontface = "italic"),
                   na_col = "grey95",
                    col = color,
                   show_heatmap_legend = F,
                 heatmap_legend_param = list(title = "Average log2\nfold change",
                                             title_gp = gpar(fontface = "plain"),
                                             grid_height = unit(0.4, "cm"),
                                             grid_width = unit(2, "mm"),
                                             tick_length = unit(0, "mm"),
                                             border = "black",
                                             title_position = "topcenter",
                                             legend_direction = "horizontal"),
                   border_gp = gpar(col = "black", lty = 1, lwd = 1))

ht_6h <- Heatmap(mat_6h,
                 column_title = "6h",
                   cluster_columns = F,
                   cluster_rows = F,
                   row_names_side = "left",
                   row_names_gp = gpar(fontface = "italic"),
                   na_col = "grey95",
                 col = color,
                   show_heatmap_legend = F,
                 heatmap_legend_param = list(title = "Average log2 fold change",
                                             title_gp = gpar(fontface = "plain"),
                                             grid_height = unit(0.4, "cm"),
                                             grid_width = unit(2, "mm"),
                                             tick_length = unit(0, "mm"),
                                             border = "black",
                                             title_position = "topcenter",
                                             legend_direction = "horizontal"),
                   border_gp = gpar(col = "black", lty = 1, lwd = 1))

ht_24h <- Heatmap(mat_24h,
                  column_title = "24h",
                   cluster_columns = F,
                   cluster_rows = F,
                   row_names_side = "left",
                   row_names_gp = gpar(fontface = "italic"),
                   na_col = "grey95",
                  col = color,
                   show_heatmap_legend = T,
                  heatmap_legend_param = list(title = "Average log2\nfold change",
                                              title_gp = gpar(fontface = "plain"),
                                              grid_height = unit(0.4, "cm"),
                                              grid_width = unit(2, "mm"),
                                              tick_length = unit(0, "mm"),
                                              border = "black",
                                              title_position = "topcenter",
                                              legend_direction = "horizontal"),
                   border_gp = gpar(col = "black", lty = 1, lwd = 1))


ht_list <- ht_counts + ht_common_fraction + ht_nonk + ht_3h + ht_6h + ht_24h

pdf("results/mm1s/timepoints/deg/singlet/deg_count_heatmap.pdf", height = 6, width = 10)
draw(ht_list, heatmap_legend_side = "bottom") 
dev.off()




