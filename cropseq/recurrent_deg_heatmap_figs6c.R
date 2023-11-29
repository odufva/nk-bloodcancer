
# Plot DEGs of all NK cell CROP-seq experiments

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)


dir.create("results/combine")

# load perturbation results
nalm6 <- fread("results/nalm6/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4 <- fread("results/sudhl4/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
k562 <- fread("results/k562/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "K562")
mm1s <- fread("results/mm1s/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S")
lp1 <- fread("results/lp1/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "LP1")
mm1s_a <- fread("results/mm1sa/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S CRISPRa")


# load gsea results
gsea_hallmark <- fread("results/combine/gsea/hallmark.txt", data.table = F) %>%
  rename(gene = pathway, p_val = pval, p_val_adj = padj, avg_log2FC = NES) %>% 
  mutate(pct.1 = NA, pct.2 = NA) %>% 
  select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, perturbation, cell_line, condition) %>% 
  mutate(condition = gsub("no_", "no ", gsub("NK_", "NK ", gsub("1_", "1:", condition))))

gsea_reactome <- fread("results/combine/gsea/reactome.txt", data.table = F) %>%
  rename(gene = pathway, p_val = pval, p_val_adj = padj, avg_log2FC = NES) %>% 
  mutate(pct.1 = NA, pct.2 = NA) %>% 
  select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, perturbation, cell_line, condition) %>% 
  mutate(condition = gsub("no_", "no ", gsub("NK_", "NK ", gsub("1_", "1:", condition))))

gsea_immunologic <- fread("results/combine/gsea/immunologic.txt", data.table = F) %>%
  rename(gene = pathway, p_val = pval, p_val_adj = padj, avg_log2FC = NES) %>% 
  mutate(pct.1 = NA, pct.2 = NA) %>% 
  select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, perturbation, cell_line, condition) %>% 
  mutate(condition = gsub("no_", "no ", gsub("NK_", "NK ", gsub("1_", "1:", condition))))


## ---------------------------------

# combine data
data <- rbind(k562, sudhl4, nalm6, mm1s, lp1)

data <- data %>% mutate(perturbation_cell_line = paste(perturbation, cell_line, sep = " "))

deg <- rbind(k562, sudhl4, nalm6, mm1s, lp1
             ) %>% mutate(cell_line_condition = paste(cell_line, condition, sep = " "))

deg_counts <- deg %>%
  mutate(cell_line_condition = factor(cell_line_condition)) %>% 
  group_by(cell_line_condition, perturbation) %>% 
  summarize(count = length(p_val[p_val_adj < 0.05])) %>% 
  mutate(signif = ifelse(count > 4, "At least\n5 DEG", "Less than\n5 DEG"))

# list perturbations with at least 5 DEGs
top_cell_line_condition_perturbation <- deg_counts %>% 
  mutate(cell_line = gsub("\\ NK.*|\\ no.*", "", cell_line_condition)) %>% 
  group_by(perturbation, cell_line) %>% 
  slice(which.max(count)) %>% 
  filter(signif == "At least\n5 DEG") %>% 
  ungroup() %>% 
  mutate(cell_line_condition_perturbation = paste(cell_line_condition, perturbation)) %>% 
  select(cell_line_condition_perturbation) %>% 
  tibble::deframe()

plotdata <- data %>%
  mutate(cell_line_condition_perturbation = paste(cell_line, condition, perturbation)) %>%
  filter(cell_line_condition_perturbation %in% top_cell_line_condition_perturbation | perturbation == "NK 1:16 vs no NK") %>%
  mutate(cell_line_condition_perturbation = factor(cell_line_condition_perturbation, levels = top_cell_line_condition_perturbation)) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>% 
  group_by(gene, direction) %>%
  mutate(count = n()) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(max_count = max(count)) %>% 
  filter(max_count > 5)

plotdata_wide <- dcast(plotdata, gene ~ perturbation_cell_line, value.var = "avg_log2FC", drop = F)
plotdata_wide[is.na(plotdata_wide)] <- 0
rownames(plotdata_wide) <- plotdata_wide$gene
plotdata_wide$gene <- NULL
plotdata_wide <- as.matrix(plotdata_wide)
plotdata_wide2 <- plotdata_wide

et_ratio_lookup_table <- plotdata %>%
  ungroup() %>%
  select(cell_line, perturbation_cell_line, condition, perturbation) %>%
  unique() %>%
  mutate(perturbation_cell_line = factor(perturbation_cell_line, levels = colnames(plotdata_wide))) %>% 
  arrange(perturbation_cell_line) %>% 
  mutate(perturbation_class = ifelse(perturbation %in% c("IFNGR2", "JAK1", "JAK2", "STAT1"), "IFNG signaling",
                                     ifelse(perturbation %in% c("RFXAP", "NLRC5"), "HLA regulator",
                                            ifelse(perturbation %in% c("TRAF2", "NFKBIA", "NFKBIB"), "NF-kB signaling",
                                                   ifelse(perturbation %in% c("FADD", "CASP8", "BID"), "Death receptor signaling",
                                            "TF/Other"))))) %>% 
  select(-perturbation_cell_line, -perturbation) %>% 
  mutate(cell_line = factor(cell_line, levels = c("K562", "SUDHL4", "MM1S", "MM1S CRISPRa", "LP1", "NALM6")),
         condition = factor(condition, levels = c("no NK", "NK 1:16", "NK 1:4"))) %>% 
  as.data.frame()

colnames(et_ratio_lookup_table) <- c("Cell line", "Condition", "Perturbation class")



ha <- HeatmapAnnotation(df = et_ratio_lookup_table,
                        col = list(`Cell line` = structure(LaCroixColoR::lacroix_palette("PeachPear", 6), 
                                                           names = c("K562", "SUDHL4", "MM1S", "MM1S CRISPRa", "LP1", "NALM6")),
                                   `Condition` = structure(pals::ocean.amp(9)[c(3,5,7)], 
                                                       names = c("no NK", "NK 1:16", "NK 1:4")),
                                   `Perturbation class` = structure(LaCroixColoR::lacroix_palette("PassionFruit", 6)[c(1:5)], 
                                                                    names = c("IFNG signaling", "HLA regulator", "NF-kB signaling", "Death receptor signaling", "TF/Other"))
                        ),
                        gap = unit(0.5, "mm"),
                        height = unit(0.75, "cm"),
                        simple_anno_size_adjust = T,
                        annotation_name_gp = gpar(fontsize = 6),
                        annotation_legend_param = list(`Cell line` = list(title = "Cell line", title_gp = gpar(fontsize = 6), 
                                                                      labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                                                       `Condition` = list(title = "Condition", title_gp = gpar(fontsize = 6), 
                                                                          labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                                                       `Perturbation class` = list(title = "Perturbation class", title_gp = gpar(fontsize = 6), 
                                                                          labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))))


ht_genes <-  Heatmap(plotdata_wide, 
                     col = colorRamp2(seq(-quantile(abs(plotdata_wide), 0.95), quantile(abs(plotdata_wide), 0.95), length = 9), rev(brewer.pal(name = "RdBu", n = 9))),
                     top_annotation = ha,
                     border_gp = gpar(col = "black", lty = 1, lwd = 1),
                     row_names_gp = gpar(fontsize = 6, fontface = "italic"),
                     column_names_gp = gpar(fontsize = 6),
                     show_heatmap_legend = T,
                     heatmap_legend_param = list(title = "Average\nfold change\n(log2)",
                                                 title_gp = gpar(fontsize = 6),
                                                 labels_gp = gpar(fontsize = 6),
                                                 grid_height = unit(0.2, "cm"),
                                                 grid_width = unit(2, "mm"),
                                                 legend_direction = "vertical"))

pdf("results/combine/dotplots/pert_heatmap_recurrent.pdf", height = 7, width = 7)
draw(ht_genes, merge_legend = T)
dev.off()

## ---------------------------------

## pathways

# Hallmark

data <- gsea_hallmark

data <- data %>% mutate(perturbation_cell_line = paste(perturbation, cell_line, sep = " "))


# top condition for each perturbation

plotdata <- data %>%
  mutate(cell_line_condition_perturbation = paste(cell_line, condition, perturbation)) %>%
  filter(cell_line_condition_perturbation %in% top_cell_line_condition_perturbation) %>%
  filter(perturbation_cell_line %in% colnames(plotdata_wide2)) %>% 
  mutate(perturbation_cell_line = factor(perturbation_cell_line, levels = colnames(plotdata_wide2))) %>% 
  filter(p_val_adj < 0.05) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>% 
  group_by(gene, direction) %>%
  mutate(count = n()) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(max_count = max(count)) %>% 
  filter(max_count > 5) 

plotdata_wide <- dcast(plotdata, gene ~ perturbation_cell_line, value.var = "avg_log2FC", drop = F)
plotdata_wide[is.na(plotdata_wide)] <- 0
rownames(plotdata_wide) <- plotdata_wide$gene
plotdata_wide$gene <- NULL
plotdata_wide <- as.matrix(plotdata_wide)

rownames(plotdata_wide) <- gsub("HALLMARK ", "", gsub("_", " ", rownames(plotdata_wide)))

ht_hallmark <- Heatmap(plotdata_wide, 
              col = colorRamp2(seq(-quantile(abs(plotdata_wide), 0.95), quantile(abs(plotdata_wide), 0.95), length = 7), pals::ocean.thermal(9)[2:8]),
              border_gp = gpar(col = "black", lty = 1, lwd = 1),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              show_heatmap_legend = T,
              heatmap_legend_param = list(title = "NES",
                                          title_gp = gpar(fontsize = 6),
                                          labels_gp = gpar(fontsize = 6),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"),
                                          legend_direction = "vertical"))

pdf("results/combine/dotplots/pert_heatmap_recurrent_pathways.pdf", height = 7, width = 16)
draw(ht_hallmark, heatmap_legend_side = "bottom", padding = unit(c(0.2, 0.2, 0.2, 8), "cm"))
dev.off()

# combined
pdf("results/combine/dotplots/pert_heatmap_recurrent_genes_hallmark.pdf", height = 7, width = 7)
draw(ht_genes %v% ht_hallmark, padding = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
dev.off()


## Reactome

data <- gsea_reactome

data <- data %>% mutate(perturbation_cell_line = paste(perturbation, cell_line, sep = " "))


# top condition for each perturbation

plotdata <- data %>%
  filter(perturbation_cell_line %in% c(paste0(c("KLF16", "RFXAP", "KIAA0922", "CMIP", "CHD7"), " NALM6"),
                                       paste0(c("IFNGR2", "FADD", "PCGF1", "RFXAP", "CASP8", "RUNX1", "STAT1", "YTHDF2", "METTL17", "CHD7", "JAK1", "JAK2", "STAG2", "PTEN", "BID"), " SUDHL4"),
                                       paste0(c("ARHGAP1", "NFKBIA", "NFKBIB", "JAK1", "JAK2", "PTEN", "GNA13", "PCGF5", "IFNGR2", "RFXAP", "NLRC5", "STAT1", "TRAF2"), " MM1S"),
                                       paste0(c("GFI1B", "JAK1", "JAK2", "PTPN2", "IFNGR2", "IRF1", "STAT1"), " K562"),
                                       paste0(c("IFNGR2", "JAK1", "JAK2","STAT1", "NLRC5", "RFXAP", "TRAF2", "PTEN", "GSK3B", "MYB", "MSI2", "CHD7"), " LP1"))) %>% 
  mutate(perturbation_cell_line = factor(perturbation_cell_line, levels = c(paste0(c("KLF16", "RFXAP", "KIAA0922", "CMIP", "CHD7"), " NALM6"),
                                                                            paste0(c("IFNGR2", "FADD", "PCGF1", "RFXAP", "CASP8", "RUNX1", "STAT1", "YTHDF2", "METTL17", "CHD7", "JAK1", "JAK2", "STAG2", "PTEN", "BID"), " SUDHL4"),
                                                                            paste0(c("ARHGAP1", "NFKBIA", "NFKBIB", "JAK1", "JAK2", "PTEN", "GNA13", "PCGF5", "IFNGR2", "RFXAP", "NLRC5", "STAT1", "TRAF2"), " MM1S"),
                                                                            paste0(c("GFI1B", "JAK1", "JAK2", "PTPN2", "IFNGR2", "IRF1", "STAT1"), " K562"),
                                                                            paste0(c("IFNGR2", "JAK1", "JAK2","STAT1", "NLRC5", "RFXAP", "TRAF2", "PTEN", "GSK3B", "MYB", "MSI2", "CHD7"), " LP1")))) %>% 
  filter(p_val_adj < 0.05) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>% 
  group_by(gene, direction) %>%
  mutate(count = n()) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(max_count = max(count)) %>% 
  filter(max_count > 5) %>% 
  filter(count > 19) %>% 
  group_by(gene, perturbation_cell_line) %>%
  top_n(1, desc(p_val))

plotdata_wide <- dcast(plotdata, gene ~ perturbation_cell_line, value.var = "avg_log2FC")
plotdata_wide[is.na(plotdata_wide)] <- 0
rownames(plotdata_wide) <- plotdata_wide$gene
plotdata_wide$gene <- NULL

rownames(plotdata_wide) <- gsub("REACTOME ", "", gsub("_", " ", rownames(plotdata_wide)))

df_annot <- data.frame(cell_line = gsub(".*\\ ", "", colnames(plotdata_wide)),
                       perturbation_class = ifelse(gsub(" .*", "", colnames(plotdata_wide)) %in% c("IFNGR2", "JAK1", "JAK2", "STAT1"), "IFNG signaling", "Other"))

ha <- HeatmapAnnotation(df = df_annot)

ht_reactome <- Heatmap(plotdata_wide, 
  col = rev(brewer.pal(name = "RdBu", n = 9)),
  top_annotation = ha,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = T,
  heatmap_legend_param = list(title = "NES",
                              title_gp = gpar(fontsize = 10),
                              labels_gp = gpar(fontsize = 10),
                              grid_height = unit(0.2, "cm"),
                              grid_width = unit(2, "mm"),
                              title_position = "topcenter",
                              legend_direction = "horizontal"))

pdf("results/combine/dotplots/pert_heatmap_recurrent_pathways_reactome.pdf", height = 8, width = 16)
draw(ht_reactome, heatmap_legend_side = "bottom", padding = unit(c(0.2, 0.2, 0.2, 8), "cm"))
dev.off()

## ------------------------------------

# top recurrent by perturbation (each perturbation only once, not in every cell line; Figure S6C)

recurrent_genes <- data %>%
  mutate(cell_line_condition_perturbation = paste(cell_line, condition, perturbation)) %>%
  filter(cell_line_condition_perturbation %in% top_cell_line_condition_perturbation | perturbation == "NK 1:16 vs no NK") %>%
  mutate(cell_line_condition_perturbation = factor(cell_line_condition_perturbation, levels = top_cell_line_condition_perturbation)) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(gene, perturbation) %>%
  top_n(1, p_val) %>% 
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>% 
  group_by(gene, direction) %>%
  mutate(count = n()) %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(max_count = max(count)) %>% 
  filter(max_count > 3) %>% 
  select(gene) %>% 
  tibble::deframe()

plotdata <- data %>% 
  mutate(cell_line_condition_perturbation = paste(cell_line, condition, perturbation)) %>%
  filter(cell_line_condition_perturbation %in% top_cell_line_condition_perturbation | perturbation == "NK 1:16 vs no NK") %>%
  mutate(cell_line_condition_perturbation = factor(cell_line_condition_perturbation, levels = top_cell_line_condition_perturbation)) %>%
  filter(p_val_adj < 0.05) %>%
  filter(gene %in% recurrent_genes)


plotdata_wide <- dcast(plotdata, gene ~ perturbation_cell_line, value.var = "avg_log2FC", drop = F)
plotdata_wide[is.na(plotdata_wide)] <- 0
rownames(plotdata_wide) <- plotdata_wide$gene
plotdata_wide$gene <- NULL
plotdata_wide <- as.matrix(plotdata_wide)
plotdata_wide2 <- plotdata_wide

et_ratio_lookup_table <- plotdata %>%
  ungroup() %>%
  select(cell_line, perturbation_cell_line, condition, perturbation) %>%
  unique() %>%
  mutate(perturbation_cell_line = factor(perturbation_cell_line, levels = colnames(plotdata_wide))) %>% 
  arrange(perturbation_cell_line) %>% 
  mutate(perturbation_class = ifelse(perturbation %in% c("IFNGR2", "JAK1", "JAK2", "STAT1"), "IFNG signaling",
                                     ifelse(perturbation %in% c("RFXAP", "NLRC5"), "HLA regulator",
                                            ifelse(perturbation %in% c("TRAF2", "NFKBIA", "NFKBIB"), "NF-kB signaling",
                                                   ifelse(perturbation %in% c("FADD", "CASP8", "BID"), "Death receptor signaling",
                                                          "TF/Other"))))) %>% 
  select(-perturbation_cell_line, -perturbation) %>% 
  mutate(cell_line = factor(cell_line, levels = c("K562", "SUDHL4", "MM1S", "MM1S CRISPRa", "LP1", "NALM6")),
         condition = factor(condition, levels = c("no NK", "NK 1:16", "NK 1:4"))) %>% 
  as.data.frame()

colnames(et_ratio_lookup_table) <- c("Cell line", "Condition", "Perturbation class")



ha <- HeatmapAnnotation(df = et_ratio_lookup_table,
                        col = list(`Cell line` = structure(LaCroixColoR::lacroix_palette("PeachPear", 6), 
                                                           names = c("K562", "SUDHL4", "MM1S", "MM1S CRISPRa", "LP1", "NALM6")),
                                   `Condition` = structure(pals::ocean.amp(9)[c(3,5,7)], 
                                                           names = c("no NK", "NK 1:16", "NK 1:4")),
                                   `Perturbation class` = structure(LaCroixColoR::lacroix_palette("PassionFruit", 6)[c(1:5)], 
                                                                    names = c("IFNG signaling", "HLA regulator", "NF-kB signaling", "Death receptor signaling", "TF/Other"))
                        ),
                        gap = unit(0.5, "mm"),
                        height = unit(0.75, "cm"),
                        simple_anno_size_adjust = T,
                        annotation_name_gp = gpar(fontsize = 6),
                        annotation_legend_param = list(`Cell line` = list(title = "Cell line", title_gp = gpar(fontsize = 6), 
                                                                          labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                                                       `Condition` = list(title = "Condition", title_gp = gpar(fontsize = 6), 
                                                                          labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                                                       `Perturbation class` = list(title = "Perturbation class", title_gp = gpar(fontsize = 6), 
                                                                                   labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))))


# plot Figure S6C
ht_genes <-  Heatmap(plotdata_wide, 
                     col = colorRamp2(seq(-quantile(abs(plotdata_wide), 0.95), quantile(abs(plotdata_wide), 0.95), length = 9), rev(brewer.pal(name = "RdBu", n = 9))),
                     top_annotation = ha,
                     border_gp = gpar(col = "black", lty = 1, lwd = 1),
                     row_names_gp = gpar(fontsize = 6, fontface = "italic"),
                     column_names_gp = gpar(fontsize = 6),
                     show_heatmap_legend = T,
                     heatmap_legend_param = list(title = "Average\nfold change\n(log2)",
                                                 title_gp = gpar(fontsize = 6),
                                                 labels_gp = gpar(fontsize = 6),
                                                 grid_height = unit(0.2, "cm"),
                                                 grid_width = unit(2, "mm"),
                                                 legend_direction = "vertical"))

pdf("results/combine/dotplots/pert_heatmap_recurrent.pdf", height = 7, width = 7)
draw(ht_genes, merge_legend = T)
dev.off()


