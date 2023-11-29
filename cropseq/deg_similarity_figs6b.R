
# Plot similarity of DEGs of all NK cell CROP-seq experiments (Figure S6B)

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(OrderedList)
library(ComplexHeatmap)
library(circlize)


dir.create("results/combine/similarity")

# load perturbation results
nalm6 <- fread("results/nalm6/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4 <- fread("results/sudhl4/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
mm1s <- fread("results/mm1s/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S")
k562 <- fread("results/k562/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "K562")
lp1 <- fread("results/lp1/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "LP1")
mm1sa <- fread("results/mm1sa/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S CRISPRa")


# calculate similarity scores for between all functional perturbations across cell lines
# for each perturbation, top condition with most DEGs

deg <- rbind(k562, sudhl4, nalm6, mm1s, lp1, mm1sa) %>% mutate(cell_line_condition = paste(cell_line, condition, sep = " "))

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


# filter common genes

data <- data %>%
  mutate(cell_line_condition_perturbation = paste(cell_line, condition, perturbation)) %>% 
  filter(cell_line_condition_perturbation %in% top_cell_line_condition_perturbation)

common_genes <- Reduce(intersect, split(data$gene, data$cell_line_condition_perturbation))

data <- data %>% filter(gene %in% common_genes)

calcSimilarityScore <- function(CELLLINE_CONDITION_PERTURBATION_1, CELLLINE_CONDITION_PERTURBATION_2, MAXRANK = 1000){
  
  list1 <- data %>% filter(cell_line_condition_perturbation == CELLLINE_CONDITION_PERTURBATION_1) %>% arrange(avg_log2FC) %>% select(gene) %>% tibble::deframe()
  list2 <- data %>% filter(cell_line_condition_perturbation == CELLLINE_CONDITION_PERTURBATION_2) %>% arrange(avg_log2FC) %>% select(gene) %>% tibble::deframe()
  
  y <- compareLists(list1, list2)
  
  score <- y$scores[y$nn == MAXRANK]
  
  df <- data.frame(cell_line_condition_perturbation_1 = CELLLINE_CONDITION_PERTURBATION_1, cell_line_condition_perturbation_2 = CELLLINE_CONDITION_PERTURBATION_2, similarity_score = score)
  return(df)
  
}

calcSimilarityScores <- function(CELLLINE_CONDITION_PERTURBATION, MAXRANK = 1000){
  
  result <- lapply(perturbations, calcSimilarityScore, CELLLINE_CONDITION_PERTURBATION_1 = CELLLINE_CONDITION_PERTURBATION, MAXRANK = MAXRANK) %>% bind_rows()
  return(result)
}

result <- lapply(top_cell_line_condition_perturbation, calcSimilarityScores) %>% bind_rows()

fwrite(result, "results/combine/similarity/similarity_scores_allpert.txt", row.names = F, quote = F, sep = "\t")


# visualize heatmap
result <- fread("results/combine/similarity/similarity_scores_allpert.txt", data.table = F)

mat <- result %>% tidyr::pivot_wider(names_from = cell_line_condition_perturbation_2, values_from = similarity_score) %>% as.data.frame()
rownames(mat) <- mat$cell_line_condition_perturbation_1
mat$cell_line_condition_perturbation_1 <- NULL
mat <- as.matrix(mat)

et_ratio_lookup_table <- data %>%
  ungroup() %>%
  select(cell_line_condition_perturbation, cell_line, condition, perturbation) %>%
  unique() %>%
  mutate(cell_line_condition_perturbation = factor(cell_line_condition_perturbation, levels = colnames(mat))) %>% 
  mutate(perturbation_cell_line = paste(perturbation, cell_line)) %>% 
  arrange(cell_line_condition_perturbation) %>% 
  mutate(perturbation_class = ifelse(perturbation %in% c("IFNGR2", "JAK1", "JAK2", "STAT1"), "IFNG signaling",
                                     ifelse(perturbation %in% c("RFXAP", "NLRC5"), "HLA regulator",
                                            ifelse(perturbation %in% c("TRAF2", "NFKBIA", "NFKBIB"), "NF-kB signaling",
                                                   ifelse(perturbation %in% c("FADD", "CASP8", "BID"), "Death receptor signaling",
                                                          "TF/Other"))))) %>% 
  select(-cell_line_condition_perturbation, -perturbation, -cell_line) %>% 
  mutate(condition = factor(condition, levels = c("no NK", "NK 1:16", "NK 1:4"))) %>% 
  as.data.frame()

colnames(mat) <- et_ratio_lookup_table$perturbation_cell_line
rownames(mat) <- et_ratio_lookup_table$perturbation_cell_line

et_ratio_lookup_table <- et_ratio_lookup_table %>% select(-perturbation_cell_line)

colnames(et_ratio_lookup_table) <- c("Condition", "Perturbation class")

ha <- HeatmapAnnotation(df = et_ratio_lookup_table,
                        col = list(`Condition` = structure(pals::ocean.amp(9)[c(3,5,7)], 
                                                           names = c("no NK", "NK 1:16", "NK 1:4")),
                                   `Perturbation class` = structure(LaCroixColoR::lacroix_palette("PassionFruit", 6)[c(1:5)], 
                                                                    names = c("IFNG signaling", "HLA regulator", "NF-kB signaling", "Death receptor signaling", "TF/Other"))
                        ),
                        gap = unit(0.5, "mm"),
                        height = unit(0.5, "cm"),
                        simple_anno_size_adjust = T,
                        annotation_name_gp = gpar(fontsize = 6),
                        annotation_legend_param = list(`Condition` = list(title = "Condition", title_gp = gpar(fontsize = 6), 
                                                                          labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                                                       `Perturbation class` = list(title = "Perturbation class", title_gp = gpar(fontsize = 6), 
                                                                                   labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))))



row_ha <- rowAnnotation(df = et_ratio_lookup_table,
                        col = list(`Condition` = structure(pals::ocean.amp(9)[c(3,5,7)], 
                                                           names = c("no NK", "NK 1:16", "NK 1:4")),
                                   `Perturbation class` = structure(LaCroixColoR::lacroix_palette("PassionFruit", 6)[c(1:5)], 
                                                                    names = c("IFNG signaling", "HLA regulator", "NF-kB signaling", "Death receptor signaling", "TF/Other"))
                        ),
                        gap = unit(0.5, "mm"),
                        width = unit(0.5, "cm"),
                        simple_anno_size_adjust = T,
                        annotation_name_gp = gpar(fontsize = 6),
                        annotation_legend_param = list(`Condition` = list(title = "Condition", title_gp = gpar(fontsize = 6), 
                                                                          labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                                                       `Perturbation class` = list(title = "Perturbation class", title_gp = gpar(fontsize = 6), 
                                                                                   labels_gp = gpar(fontsize = 6), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))))


ht_similarity <- Heatmap(mat,
                         name = "Similarity score",
                         top_annotation = ha,
                         right_annotation = row_ha,
                         col = colorRamp2(seq(quantile(abs(mat), 0.025), quantile(abs(mat), 0.975), length.out = 9), rev(pals::ocean.matter(9))),
                         rect_gp = gpar(col = "white", lwd = unit(0.4, "mm")),
                         border_gp = gpar(col = "black", lty = 1, lwd = 1),
                         row_names_gp = gpar(fontsize = 6),
                         column_names_gp = gpar(fontsize = 6),
                         heatmap_legend_param = list(title = "Similarity score",
                                                     title_gp = gpar(fontsize = 6),
                                                     labels_gp = gpar(fontsize = 6),
                                                     grid_height = unit(0.2, "cm"),
                                                     grid_width = unit(2, "mm"),
                                                     at = c(quantile(abs(mat), 0.025), quantile(abs(mat), 0.975)),
                                                     labels = c("Low", "High"),
                                                     border = NA,
                                                     legend_direction = "vertical"))



pdf("results/combine/similarity/allpert_heatmap.pdf", height = 7, width = 8)
draw(ht_similarity, merge_legend = T)
dev.off()

