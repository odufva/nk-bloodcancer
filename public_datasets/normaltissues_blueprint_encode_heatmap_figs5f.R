
# Plot heatmap of CRISPR screen hits with specific expression pattern in blood cancers vs solid tumors

library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# load data
load("data/blueprint_encode.rda")


# select genes based on differential analysis in CCLE cell line data
genelist <- c("CD48", "SPN", "RHOH", "MYB", "SELPLG", "TNFRSF1B", "PVR", "ULBP3")

# subset data to selected genes
gexp <- blueprint_encode$data[genelist,]

df <- data.frame(t(gexp),
                 main_type = blueprint_encode$main_types,
                 type = blueprint_encode$types,
                 sample = colnames(blueprint_encode$data),
                 check.names = F,
                 fix.empty.names = F)

# matrix of gene expression medians 
gexp_mat <- df %>%
  group_by(main_type) %>% 
  summarize(across(c(1:8), median, na.rm = T))

main_type <- gexp_mat$main_type
gexp_mat$main_type <- NULL

# scale data
gexp_mat_scaled <- t(apply(gexp_mat, 2, scale))
colnames(gexp_mat_scaled) <- main_type

heme <- c("DC", "Macrophages", "Monocytes", "Neutrophils", "Eosinophils", "Erythrocytes", "HSC", "CD4+ T-cells", "CD8+ T-cells", "NK cells", "B-cells")
solid <- colnames(gexp_mat_scaled)[!colnames(gexp_mat_scaled) %in% c(heme, "Neurons", "Skeletal muscle")]

gexp_mat_scaled <- gexp_mat_scaled[,c(heme, solid, "Neurons", "Skeletal muscle")]

# plot heatmap
ht_gexp <- Heatmap(gexp_mat_scaled,
                        name = "Expression",
                        col = colorRamp2(seq(-1.5, 2.5, length.out = 11), pals::ocean.deep(11)),
                        rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
                        column_names_side = "top",
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                        column_names_gp = gpar(fontsize = 7),
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        row_title_gp = gpar(fontsize = 5),
                        heatmap_legend_param = list(title = "Scaled\nmedian\nexpression",
                                                    title_gp = gpar(fontsize = 7),
                                                    labels_gp = gpar(fontsize = 7),
                                                    grid_height = unit(0.2, "cm"),
                                                    grid_width = unit(3, "mm"))
)

pdf("blueprint_encode_heatmap.pdf", height = 1.8, width = 4.25)
draw(ht_gexp)
decorate_heatmap_body("Expression", {
  
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.5))
  
})
dev.off()

