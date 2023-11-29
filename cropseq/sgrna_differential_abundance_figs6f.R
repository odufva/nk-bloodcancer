
# Analyse sgRNA differential abundance from CROP-seq data

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
library(scales)
library(patchwork)
library(viridis)

theme_set(theme_classic(base_size = 12))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

# K562

# load object
crop_seurat <- readRDS("results/k562/k562_crop_seurat_singlet.rds")

# counts of perturbed cells
counts <- as.matrix(table(crop_seurat$sgrna_name, crop_seurat$orig.ident))
counts <- prop.table(counts, 2)

counts <- as.data.frame(counts)
colnames(counts) <- c("sgrna", "condition", "norm_count")

counts <- counts %>%
  mutate(gene = gsub("sg|\\..*", "", sgrna))

test_gene <- function(GENE, DATA, CONDITION1, CONDITION2){
  
  data <- DATA
  counts1 <- data[data$condition == CONDITION1 & data$gene == GENE, "norm_count"]
  counts2 <- data[data$condition == CONDITION2 & data$gene == GENE, "norm_count"]
  res <- t.test(counts1, counts2, paired = T)
  log2fc <- log2(mean(counts2)/mean(counts1))
  df <- data.frame(gene = GENE, condition1 = CONDITION1, condition2 = CONDITION2,
                   log2fc = log2fc, 
                   log2fc_1 = log2(counts2[1]/counts1[1]),
                   log2fc_2 = log2(counts2[2]/counts1[2]),
                   log2fc_3 = log2(counts2[3]/counts1[3]),
                   log2fc_4 = log2(counts2[4]/counts1[4]),
                   log2fc_4 = log2(counts2[5]/counts1[5]),
                   log2fc_4 = log2(counts2[6]/counts1[6]),
                   p = res$p.value)
  return(df)
  
}

result1 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:16") %>%
  bind_rows() %>% arrange(log2fc)
result1$FDR <- p.adjust(result1$p, method = "BH")
result1

dir.create("results/k562/differential_abundance")
dir.create("results/k562/differential_abundance/singlet")
fwrite(result1, "results/k562/differential_abundance/singlet/differential_abundance.txt", sep = "\t")

# plots

# 1:16 
data_plot <- result1 %>%
  mutate(gene = factor(gene, levels = result1$gene[order(result1$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave("results/k562/differential_abundance/singlet/1_16_vs_nonk_barplot.pdf", height = 4, width = 6)

## --------------------------------

# SUDHL4 

# load object
crop_seurat <- readRDS("results/sudhl4/sudhl4_crop_seurat_singlet.rds")

# counts of perturbed cells
counts <- as.matrix(table(crop_seurat$sgrna_name, crop_seurat$orig.ident))
counts <- prop.table(counts, 2)

counts <- as.data.frame(counts)
colnames(counts) <- c("sgrna", "condition", "norm_count")

counts <- counts %>%
  mutate(gene = gsub("sg|\\..*", "", sgrna))

test_gene <- function(GENE, DATA, CONDITION1, CONDITION2){
  
  data <- DATA
  counts1 <- data[data$condition == CONDITION1 & data$gene == GENE, "norm_count"]
  counts2 <- data[data$condition == CONDITION2 & data$gene == GENE, "norm_count"]
  res <- t.test(counts1, counts2, paired = T)
  log2fc <- log2(mean(counts2)/mean(counts1))
  df <- data.frame(gene = GENE, condition1 = CONDITION1, condition2 = CONDITION2,
                   log2fc = log2fc, 
                   log2fc_1 = log2(counts2[1]/counts1[1]),
                   log2fc_2 = log2(counts2[2]/counts1[2]),
                   log2fc_3 = log2(counts2[3]/counts1[3]),
                   log2fc_4 = log2(counts2[4]/counts1[4]),
                   log2fc_4 = log2(counts2[5]/counts1[5]),
                   log2fc_4 = log2(counts2[6]/counts1[6]),
                   p = res$p.value)
  return(df)
  
}

result1 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:16") %>%
  bind_rows() %>% arrange(log2fc)
result1$FDR <- p.adjust(result1$p, method = "BH")
result1

result2 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:4") %>%
  bind_rows() %>% arrange(log2fc)
result2$FDR <- p.adjust(result2$p, method = "BH")
result2

result <- rbind(result1, result2)
dir.create("results/sudhl4/differential_abundance/singlet")
fwrite(result, "results/sudhl4/differential_abundance/singlet/differential_abundance.txt", sep = "\t")

# plots

# 1:16 
data_plot <- result1 %>%
  mutate(gene = factor(gene, levels = result1$gene[order(result1$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave("results/sudhl4/differential_abundance/singlet/1_16_vs_nonk_barplot.pdf", height = 4, width = 6)


# 1:4
data_plot <- result2 %>%
  mutate(gene = factor(gene, levels = result2$gene[order(result2$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")
ggsave("results/sudhl4/differential_abundance/singlet/1_4_vs_nonk_barplot.pdf", height = 4, width = 6)


## --------------------------------

# NALM6

# load object
crop_seurat <- readRDS("results/nalm6/nalm6_crop_seurat_singlet.rds")

# counts of perturbed cells
counts <- as.matrix(table(crop_seurat$sgrna_name, crop_seurat$orig.ident))
counts <- prop.table(counts, 2)

counts <- as.data.frame(counts)
colnames(counts) <- c("sgrna", "condition", "norm_count")

counts <- counts %>%
  mutate(gene = gsub("sg|\\..*", "", sgrna))

test_gene <- function(GENE, DATA, CONDITION1, CONDITION2){
  
  data <- DATA
  counts1 <- data[data$condition == CONDITION1 & data$gene == GENE, "norm_count"]
  counts2 <- data[data$condition == CONDITION2 & data$gene == GENE, "norm_count"]
  res <- t.test(counts1, counts2, paired = T)
  log2fc <- log2(mean(counts2)/mean(counts1))
  df <- data.frame(gene = GENE, condition1 = CONDITION1, condition2 = CONDITION2,
                   log2fc = log2fc, 
                   log2fc_1 = log2(counts2[1]/counts1[1]),
                   log2fc_2 = log2(counts2[2]/counts1[2]),
                   log2fc_3 = log2(counts2[3]/counts1[3]),
                   log2fc_4 = log2(counts2[4]/counts1[4]),
                   log2fc_4 = log2(counts2[5]/counts1[5]),
                   log2fc_4 = log2(counts2[6]/counts1[6]),
                   p = res$p.value)
  return(df)
  
}

result1 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:16") %>%
  bind_rows() %>% arrange(log2fc)
result1$FDR <- p.adjust(result1$p, method = "BH")
result1

result2 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:4") %>%
  bind_rows() %>% arrange(log2fc)
result2$FDR <- p.adjust(result2$p, method = "BH")
result2

result <- rbind(result1, result2)
dir.create("results/nalm6/differential_abundance/singlet")
fwrite(result, "results/nalm6/differential_abundance/singlet/differential_abundance.txt", sep = "\t")

# plots

# 1:16 
data_plot <- result1 %>%
  mutate(gene = factor(gene, levels = result1$gene[order(result1$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave("results/nalm6/differential_abundance/singlet/1_16_vs_nonk_barplot.pdf", height = 4, width = 6)


# 1:4
data_plot <- result2 %>%
  mutate(gene = factor(gene, levels = result2$gene[order(result2$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")
ggsave("results/nalm6/differential_abundance/singlet/1_4_vs_nonk_barplot.pdf", height = 4, width = 6)


## --------------------------------

# MM1S

# load object
crop_seurat <- readRDS("results/mm1s/mm1s_crop_seurat_singlet_2.rds")


# counts of perturbed cells
counts <- as.matrix(table(crop_seurat$sgrna_name, crop_seurat$orig.ident))
counts <- prop.table(counts, 2)

counts <- as.data.frame(counts)
colnames(counts) <- c("sgrna", "condition", "norm_count")

counts <- counts %>%
  mutate(gene = gsub("sg|\\..*", "", sgrna))

test_gene <- function(GENE, DATA, CONDITION1, CONDITION2){
  
  data <- DATA
  counts1 <- data[data$condition == CONDITION1 & data$gene == GENE, "norm_count"]
  counts2 <- data[data$condition == CONDITION2 & data$gene == GENE, "norm_count"]
  res <- t.test(counts1, counts2, paired = T)
  log2fc <- log2(mean(counts2)/mean(counts1))
  df <- data.frame(gene = GENE, condition1 = CONDITION1, condition2 = CONDITION2,
                   log2fc = log2fc, 
                   log2fc_1 = log2(counts2[1]/counts1[1]),
                   log2fc_2 = log2(counts2[2]/counts1[2]),
                   log2fc_3 = log2(counts2[3]/counts1[3]),
                   log2fc_4 = log2(counts2[4]/counts1[4]),
                   log2fc_4 = log2(counts2[5]/counts1[5]),
                   log2fc_4 = log2(counts2[6]/counts1[6]),
                   p = res$p.value)
  return(df)
  
}

result1 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:16") %>%
  bind_rows() %>% arrange(log2fc)
result1$FDR <- p.adjust(result1$p, method = "BH")
result1

result2 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:4") %>%
  bind_rows() %>% arrange(log2fc)
result2$FDR <- p.adjust(result2$p, method = "BH")
result2

result <- rbind(result1, result2)
dir.create("results/mm1s/differential_abundance")
dir.create("results/mm1s/differential_abundance/singlet")
fwrite(result, "results/mm1s/differential_abundance/singlet/differential_abundance.txt", sep = "\t")

# plots

# 1:16 
data_plot <- result1 %>%
  mutate(gene = factor(gene, levels = result1$gene[order(result1$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave("results/mm1s/differential_abundance/singlet/1_16_vs_nonk_barplot_qc3.pdf", height = 4, width = 6)


# 1:4
data_plot <- result2 %>%
  mutate(gene = factor(gene, levels = result2$gene[order(result2$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")
ggsave("results/mm1s/differential_abundance/singlet/1_4_vs_nonk_barplot_qc3.pdf", height = 4, width = 6)


## --------------------------------

# LP1

# load object
crop_seurat <- readRDS("results/lp1/lp1_crop_seurat_singlet_2.rds")

# counts of perturbed cells
counts <- as.matrix(table(crop_seurat$sgrna_name, crop_seurat$orig.ident))
counts <- prop.table(counts, 2)

counts <- as.data.frame(counts)
colnames(counts) <- c("sgrna", "condition", "norm_count")

counts <- counts %>%
  mutate(gene = gsub("sg|\\..*", "", sgrna))

test_gene <- function(GENE, DATA, CONDITION1, CONDITION2){
  
  data <- DATA
  counts1 <- data[data$condition == CONDITION1 & data$gene == GENE, "norm_count"]
  counts2 <- data[data$condition == CONDITION2 & data$gene == GENE, "norm_count"]
  res <- t.test(counts1, counts2, paired = T)
  log2fc <- log2(mean(counts2)/mean(counts1))
  df <- data.frame(gene = GENE, condition1 = CONDITION1, condition2 = CONDITION2,
                   log2fc = log2fc, 
                   log2fc_1 = log2(counts2[1]/counts1[1]),
                   log2fc_2 = log2(counts2[2]/counts1[2]),
                   log2fc_3 = log2(counts2[3]/counts1[3]),
                   log2fc_4 = log2(counts2[4]/counts1[4]),
                   log2fc_4 = log2(counts2[5]/counts1[5]),
                   log2fc_4 = log2(counts2[6]/counts1[6]),
                   p = res$p.value)
  return(df)
  
}

result1 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:16") %>%
  bind_rows() %>% arrange(log2fc)
result1$FDR <- p.adjust(result1$p, method = "BH")
result1

result2 <- lapply(unique(counts$gene), test_gene, DATA = counts, CONDITION1 = "no NK", CONDITION2 = "NK 1:4") %>%
  bind_rows() %>% arrange(log2fc)
result2$FDR <- p.adjust(result2$p, method = "BH")
result2

result <- rbind(result1, result2)
dir.create("results/lp1/differential_abundance")
dir.create("results/lp1/differential_abundance/singlet")
fwrite(result, "results/lp1/differential_abundance/singlet/differential_abundance_qc3.txt", sep = "\t")

# plots

# 1:16 
data_plot <- result1 %>%
  mutate(gene = factor(gene, levels = result1$gene[order(result1$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave("results/lp1/differential_abundance/singlet/1_16_vs_nonk_barplot_qc3.pdf", height = 4, width = 6)


# 1:4
data_plot <- result2 %>%
  mutate(gene = factor(gene, levels = result2$gene[order(result2$log2fc)])) %>%
  tidyr::pivot_longer(cols = log2fc:FDR, names_to = "variable", values_to = "value")

ggplot(data_plot[data_plot$variable=="log2fc",], aes(x = value, y = gene, fill = value)) +
  geom_col() +
  geom_point(data = data_plot[grepl("log2fc_", data_plot$variable),]) +
  scale_fill_viridis(option = "magma", begin = 0.2, end = 0.8) +
  xlab("Fold change (log2)") +
  ylab("") +
  theme_cowplot() +
  theme(legend.position = "none")
ggsave("results/lp1/differential_abundance/singlet/1_4_vs_nonk_barplot_qc3.pdf", height = 4, width = 6)



## --------------------------------

# Analyze all differential abundances together and compare to genome-wide screens

# load differential abundance results
nalm6_counts <- fread("results/nalm6/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4_counts <- fread("results/sudhl4/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
k562_counts <- fread("results/k562/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "K562")
mm1s_counts <- fread("results/mm1s/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "MM1S")
lp1_counts <- fread("results/lp1/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "LP1")

counts <- rbind(nalm6_counts, sudhl4_counts, k562_counts, mm1s_counts, lp1_counts)


# bar plot of differential abundance

# get activating vs inhibitory classification from genome-wide screen data

# load screen results
crispr <- fread("crispr_mageck_combined.txt", data.table = F)

# counts of cells

plot_counts <- counts %>%
  filter(condition2 == "NK 1:4") %>% 
  rename(perturbation = gene, condition = condition2, ) %>% 
  mutate(perturbation_cell_line = paste(perturbation, cell_line)) %>% 
  mutate(signif = ifelse(p < 0.05, "p < 0.05", "ns")) %>% 
  left_join(crispr, by = c("perturbation" = "gene", "cell_line" = "cell_line"), suffix = c("_cropseq", "_crispr")) %>% 
  mutate(effect = ifelse(lfc > 0, "Enriched", "Depleted")) %>% 
  mutate(effect = ifelse(p_crispr > 0.05, "ns", effect)) %>% 
  mutate(cell_line = factor(cell_line, levels = c("SUDHL4", "MM1S", "LP1", "NALM6")))

plot_counts <- plot_counts %>% 
  tidyr::pivot_longer(cols = contains("log2fc_"), names_to = "replicate", values_to = "log2fc_replicate")

library(tidytext)


barplot_counts <- ggplot(plot_counts, aes(x = reorder_within(perturbation, log2fc, cell_line), y = log2fc, fill = effect)) +
  geom_col(data = plot_counts[plot_counts$replicate=="log2fc_1",]) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_manual("Effect in\ngenome-scale\nCRISPR screen", values = c(brewer.pal(11, "RdBu")[c(8,3)], "grey70"),
                    breaks = c("Depleted", "Enriched", "ns")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_reordered() +
  xlab("") +
  ylab("Fold change (log2)\naverage sgRNA abundance\nNK-treated vs. untreated") +
  facet_grid(cols = vars(cell_line), scales = "free_x")

# plot Figure S6F
barplot_counts
dir.create("results/combine/differential_abundance")
ggsave("results/combine/differential_abundance/pert_counts_barplot_allconditions.pdf", height = 3.5, width = 20)


