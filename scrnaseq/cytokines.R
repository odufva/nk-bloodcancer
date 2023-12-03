
# Analyze cytokine data from 26 cell lines co-cultured with NK cells

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(data.table)
library(gridExtra)
library(patchwork)
library(Seurat)
library(psych)
library(tidyr)

# load data 
data_1_raw <- read_excel("data/Plate 1.xls")
data_2_raw <- read_excel("data/230331 Mustjoki cytokine_results.xlsx", sheet = 5)
data_3_raw <- read_excel("data/230331 Mustjoki cytokine_results.xlsx", sheet = 7)

# load annotations
annot <-fread("data/NK_PRISM_heme_pctviability.txt")

annot <- annot %>%
  mutate(Subtype_simple = ifelse(grepl("AML", Subtype), "AML",
                                 ifelse(grepl("ALL), B-cell", Subtype), "B-ALL",
                                        ifelse(grepl("ALL), T-cell", Subtype), "T-ALL",
                                               ifelse(grepl("Multiple", Subtype), "MM",
                                                      ifelse(grepl("CLL", Subtype), "CLL",
                                                             ifelse(grepl("Mantle", Subtype), "MCL",
                                                                    ifelse(grepl("B-cell, Hodgkins", Subtype), "CHL",
                                                                           ifelse(grepl("Burkitt", Subtype), "BL",
                                                                                  ifelse(grepl("ALCL", Subtype), "ALCL",
                                                                                         ifelse(grepl("DLBCL", Subtype), "DLBCL", 
                                                                                                ifelse(grepl("B-cell, Non", Subtype), "BCL other",
                                                                                                       ifelse(grepl("Cutaneous", Subtype), "CTCL",
                                                                                                              ifelse(grepl("^T-cell", Subtype), "TCL other", Subtype)))))))))))))) %>% 
  mutate(Subtype_main = ifelse(Subtype_simple %in% c("DLBCL", "BCL other", "CHL", "BL", "MCL", "CLL"), "BCL",
                               ifelse(Subtype_simple %in% c("CTCL", "ALCL", "ATL", "TCL other"), "TCL", Subtype_simple))) %>% 
  select(cell_line, Subtype_simple, Subtype_main) %>% 
  unique()

# add missing cell lines
annot <- rbind(annot,
               data.frame(cell_line = c("K562", "MOLM14", "SUDHL4"),
                          Subtype_simple = c("CML", "AML", "BCL"),
                          Subtype_main = c("CML", "AML", "BCL")
               ))

# load colors
cols <- fread("../CCLE_featurematrix_NK_PRISM/nk_crispr_colors.txt", data.table = F) %>% dplyr::rename(cancer = cancer_type)
cols <- rbind(cols, data.frame(color = c("lightblue", "grey70"), cancer = c("CML", "NK cells only")))

cols_heatmap <- cols
cols_heatmap$cancer <- gsub("DLBCL", "BCL", cols_heatmap$cancer)
cols_heatmap <- cols_heatmap[cols_heatmap$cancer %in% c("CML", "AML", "T-ALL", "BCL", "MM", "B-ALL", "NK cells only"),]

cols_heatmap_vector <- as.character(cols_heatmap$color)
names(cols_heatmap_vector) <- cols_heatmap$cancer
cols_heatmap_vector <- cols_heatmap_vector[c("CML", "AML", "T-ALL", "BCL", "MM", "B-ALL", "NK cells only")]


# select data of interest
data1 <- data_1_raw[30:86,-c(1:2)]
data2 <- data_2_raw[19:93,-c(1:2)]
data3 <- data_3_raw[19:85,-c(1:2)]

colnames(data1) <- gsub("Hu | \\([0-9]*\\)", "", data_1_raw[8,-c(1,2)])
colnames(data2) <- gsub("Hu | \\([0-9]*\\)", "", data_2_raw[8,-c(1,2)])
colnames(data3) <- gsub("Hu | \\([0-9]*\\)", "", data_3_raw[8,-c(1,2)])

data <- rbind(data1, data2, data3)

# remove asterisks and "OOR <"

clean_data <- function(COLUMN){
  
  vector <- data[,COLUMN][[1]]
  
  vector <- as.character(gsub("\\*", "", vector)) # remove asterisks
  
  min_vector <- min(as.numeric(vector[!grepl("OOR", vector)])) # get minimum value
  
  max_vector <- max(as.numeric(vector[!grepl("OOR", vector)])) # get maximum value
  
  vector[vector == "OOR <"] <- min_vector # replace "OOR <" with minimum value
  
  vector[vector == "OOR >"] <- max_vector # replace "OOR >" with maximum value
  
  vector <- as.numeric(vector)
  
}

data_clean <- lapply(2:28, clean_data) %>% bind_cols() %>% as.data.frame()

colnames(data_clean) <- colnames(data)[2:28]

cytokine_names_complete <- c("IL-1b (IL1B)", "IL-1ra (IL1RN)", "IL-2 (IL2)", "IL-4 (IL4)", "IL-5 (IL5)", "IL-6 (IL6)", "IL-7 (IL7)", "IL-8 (IL8)", "IL-9 (IL9)", "IL-10 (IL10)", "IL-12 (IL12A)", "IL-13 (IL13)",
                             "IL-15 (IL15)", "IL-17 (IL17A)", "Eotaxin", "FGF basic (FGF2)", "G-CSF (CSF3)", "GM-CSF (CSF2)", "IFN-g (IFNG)", "IP-10 (CXCL1)0", "MCP-1 (CCL2)", "MIP-1a (CCL3)", "MIP-1b (CCL4)", "RANTES (CCL5)", "TNF-a (TNF)", "VEGF")

cytokine_names <- c("IL-1b", "IL-1ra", "IL-2", "IL-4", "IL-5", "IL-6", "IL-7", "IL-8", "IL-9", "IL-10", "IL-12", "IL-13",
                    "IL-15", "IL-17", "Eotaxin", "FGF2", "G-CSF", "GM-CSF", "IFNg", "IP-10 (CXCL10)", "MCP-1 (CCL2)", "MIP-1a (CCL3)", "PDGF-BB", "MIP-1b (CCL4)", "RANTES (CCL5)", "TNF-a", "VEGF")


colnames(data_clean) <- cytokine_names

# remove batch controls
data_clean <- data_clean[c(1:131, 133:197),]


row_names <- data[c(1:131, 133:197),1][[1]]
row_names[53:57] <- c("NK1_1", "NK1_2", "NK1_3", "NK1_4", "NK1_5")
row_names <- gsub("RIVA", "RI1", row_names)


# prepare metadata

sample_annot <- read_excel("230331 Mustjoki cytokine_results_annotations.xlsx", sheet = 8)


cell_line = gsub("_.*", "", row_names)
cell_line[grepl("NK", cell_line)] = "NK"


df_annot <- data.frame(sample_id = row_names,
                       cell_line = cell_line
                       )

df_annot <- df_annot %>% left_join(sample_annot, by = c("sample_id" = "Sample"))

row_names <- gsub("-", "_", gsub("NK$|NK-ctrl|NK2587", "NK1", gsub("NK2770", "NK2", gsub("NK2763", "NK3", gsub("NK2768", "NK4", gsub("NK2771", "NK5", gsub("NK2931", "NK6", row_names)))))))
rownames(data_clean) <- row_names

df_annot$sample_id_final <- row_names


# data prepared

# data for supplement

df_annot$sample_id_final == rownames(data_clean)

data_supplement <- cbind(df_annot[,c("sample_id_final", "cell_line", "Batch_main", "Batch")], data_clean)

fwrite(data_supplement, "cytokines_samples.txt", sep = "\t")


# PCA

pca <- prcomp(data_clean, center = TRUE, scale. = TRUE)

pca_df <- data.frame(pca$x)
pca_df <- cbind(pca_df, df_annot)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  geom_text_repel(aes(label = sample_id))


p2 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = Batch)) +
  geom_point() +
  geom_text_repel(aes(label = sample_id))

p1 + p2

ggsave("cytokines_pca.pdf", height = 6, width = 16)


data_clean_targets <- data_clean[!grepl("NK", row_names),]

pca_targets <- prcomp(data_clean_targets, center = TRUE, scale. = TRUE)

pca_df_targets <- data.frame(pca_targets$x)
pca_df_targets$sample <- rownames(data_clean_targets)

ggplot(pca_df_targets, aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(aes(label = sample), max.overlaps = 20)

ggsave("cytokines_pca_targets.pdf", height = 6, width = 10)


# Prepare annotations

# Transpose data for differential tests
data <- t(data_clean)[,1:57]


# Vectors of sample groups for differential tests
untreated <- colnames(data)[!grepl("NK", colnames(data))]
nk_treated <- colnames(data)[grepl("_NK", colnames(data))]

untreated == gsub("_NK1", "", nk_treated) # check equality

# Cancer type annotation
cancer_type <- annot$Subtype_main[match(untreated, annot$cell_line)]
cancer_type[c(1:6)] <- c("CML", "AML", "AML", "AML", "BCL", "BCL")
cancer_type <- c(cancer_type, cancer_type, rep("NK cells only", 5))

# Treatment annotation
treatment <- c(rep("Untreated", 26),  rep("NK-treated", 26),  rep("NK cells only", 5))

# Data frame for heatmap annotations
heatmap_annot <- cbind(cancer_type, treatment)
rownames(heatmap_annot) <- colnames(data)
heatmap_annot <- as.data.frame(heatmap_annot)

heatmap_annot$treatment <- factor(heatmap_annot$treatment, levels = c("Untreated", "NK-treated", "NK cells only"))
heatmap_annot$cancer_type <- factor(heatmap_annot$cancer_type, levels = c("CML", "AML", "T-ALL", "BCL", "MM", "B-ALL", "NK cells only"))

# Sample order for heatmap and barplots
column_order <- rownames(heatmap_annot[order(heatmap_annot$treatment, heatmap_annot$cancer_type, gsub("-NK", "", rownames(heatmap_annot))),])

## -----------------------------------------------------------

# Prepare data for bar plots

data_clean_sample <- data_clean
data_clean_sample$sample <- rownames(data_clean)
data_clean_sample$batch <- df_annot$Batch

df <- data_clean_sample %>% tidyr::pivot_longer(cols = !contains(c("batch", "sample")), names_to = "cytokine", values_to = "concentration")

# remove time point data
df <- df %>%
  filter(batch %in% c("Main_batch1", "Batch5", "Batch6", "Batch7")) %>% 
  mutate(cell_line = gsub("_.*", "", gsub("NK.*", "NK cells only", sample)),
                    condition = ifelse(grepl("NK", sample), "NK-treated", "Untreated")) %>% 
mutate(cell_line = factor(cell_line, levels = c(column_order[1:26], "NK cells only")),
       condition = factor(condition, levels = c("Untreated", "NK-treated")))



df$cancer_type <- annot$Subtype_main[match(df$cell_line, annot$cell_line)]
df$cancer_type[df$cell_line=="NK"] <- "NK cells only"


# calculate mean and sd

df <- df %>% 
  group_by(cytokine, cell_line, condition) %>% 
  mutate(mean = mean(concentration), sd = sd(concentration))



ggplot(df, aes(x = sample, y = concentration)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ cytokine, scales = "free_y", ncol = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fwrite(df, "cytokine_celllinepanel_processed.txt", sep = "\t")

# supplement table with means

df_supplement_mean <- df %>%
  select(cell_line, condition, cytokine, mean) %>% unique() %>%
  pivot_wider(id_cols = c("cell_line", "condition"), names_from = cytokine, values_from = mean) 

fwrite(df_supplement_mean, "cytokines_mean.txt", sep = "\t")


# Bar plots 

df <- fread("cytokine_celllinepanel_processed.txt")


cytokine_barplot <- function(CYTOKINE){

  df_plot <- df %>% filter(cytokine == CYTOKINE)
  
  ggplot(df_plot, aes(x = cell_line, y = concentration, fill = cancer_type)) +
    geom_bar(stat = "summary", fun = "mean") +
    geom_errorbar(data = df_plot[df_plot$batch=="Main_batch1",], aes(ymin=mean,ymax=mean+sd),
                  width=0.3,
                  size=0.3) +
    geom_point() +
    facet_wrap(. ~ condition, ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(color = "black"),
          panel.grid.minor.y = element_blank() ,
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line( size=.1, color="grey50" ),
          strip.background = element_blank()) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(values = cols_heatmap_vector) +
    ylab("Concentration (pg/ml)") +
    xlab("") +
    labs(fill = "") +
    guides(fill = guide_legend(override.aes = list(shape = NA))) +
    ggtitle(CYTOKINE)
  
}


p1 <- lapply(cytokine_names, cytokine_barplot)
m1 <- marrangeGrob(p1, nrow = 2, ncol = 3, top = NULL)
ggsave("cytokine_barplots_alldonors.pdf", m1, height = 10, width = 20)


## ----------------------------------------


# Test differential abundance


test_differential <- function(CYTOKINE, DATA, GROUP1, GROUP2, PAIRED = F){
  
  group1 <- as.numeric(DATA[DATA$cytokine==CYTOKINE&DATA$condition==GROUP1,]$concentration)
  group2 <- as.numeric(DATA[DATA$cytokine==CYTOKINE&DATA$condition==GROUP2,]$concentration)
  
  
  if (sd(group1) != 0 & sd(group2) != 0){
  
    res <- t.test(group1, group2, paired = PAIRED)
    
    result <- data.frame(cytokine = CYTOKINE,
                       log2fc = log2(mean(group1)/mean(group2)),
                       mean1 = mean(group1),
                       mean2 = mean(group2),
                       p = res$p.value)
  
  }
  
  else {
    
    result <- data.frame(cytokine = CYTOKINE,
                         log2fc = log2(mean(group1)/mean(group2)),
                         mean1 = mean(group1),
                         mean2 = mean(group2),
                         p = NA)
  }
  
  return(result)
  
}


test_cellline <- function(CELLLINE, DATA = df){
  
  # DATA = df
  # CELLLINE = "K562"
  
  # subset to one cell line
  df_cellline  <- DATA[DATA$cell_line == CELLLINE,]
  
  # apply differential test function over all cytokines
  res_cellline <- lapply(cytokines, test_differential, DATA = df_cellline, GROUP1 = "NK-treated", GROUP2 = "Untreated") %>% bind_rows()
  res_cellline$p_adj <- p.adjust(res_cellline$p, method = "BH")
  res_cellline$cell_line <- CELLLINE
  res_cellline <- res_cellline %>% arrange(p)
  
  return(res_cellline)
  
}


cytokines <- unique(df$cytokine)
celllines <- unique(df$cell_line)
celllines <- celllines[celllines!="NK cells only"]


# test with original values
result_all <- lapply(celllines, test_cellline) %>% bind_rows

fwrite(result_all, "cytokines_differential.txt", sep = "\t")


# test with values where NK cell median has been added to unreated
nk_median <- df %>%
  filter(cell_line == "NK cells only") %>% 
  group_by(cytokine) %>%
  summarize(nk_median_concentration = median(concentration))
  

df_nk_median_untreated <- df %>%
  filter(condition == "Untreated") %>% 
  left_join(nk_median) %>% 
  mutate(concentration_withnk = concentration + nk_median_concentration)

df_nktreated <- df %>% 
  filter(condition == "NK-treated")

df_nk_median <- rbind(df_nktreated, df_nk_median_untreated) %>%
  mutate(concentration = ifelse(condition == "Untreated", concentration_withnk, concentration))


result_all_withnk <- lapply(celllines, test_cellline, DATA = df_nk_median) %>% bind_rows

fwrite(result_all_withnk, "cytokines_differential_nk_added.txt", sep = "\t")


# Dot plot

df_dotplot <- result_all_withnk %>% 
  mutate(log10p = -log10(p)) %>% 
  filter(p < 0.1)

ggplot(df_dotplot, aes(x = cell_line, y = cytokine, color = log2fc)) +
  geom_point() +
  scale_color_distiller("Fold\nchange\n(log2)", palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                       type = "div", limits = max(abs(as.numeric(df_dotplot$log2fc))) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top")) +
  scale_size("P value\n(-log10)", range = c(2, 4.5), guide = guide_legend(title.position = "top",
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
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  ylab("") +
  xlab("")


## ----------------------------------------

# Test differential abundance with paired tests using means of replicates

df_mean <- df %>% 
  filter(cell_line != "NK cells only") %>% 
  ungroup() %>% 
  group_by(cytokine, cell_line, condition) %>% 
  summarize(concentration_mean = mean(concentration))

test_differential_wilcox <- function(CYTOKINE, DATA, GROUP1, GROUP2, PAIRED = F){
  
  group1 <- as.numeric(DATA[DATA$cytokine==CYTOKINE&DATA$condition==GROUP1,]$concentration_mean)
  group2 <- as.numeric(DATA[DATA$cytokine==CYTOKINE&DATA$condition==GROUP2,]$concentration_mean)
  
  
    res <- wilcox.test(group1, group2, paired = PAIRED)
    
    result <- data.frame(cytokine = CYTOKINE,
                         log2fc = log2(mean(group1)/mean(group2)),
                         mean1 = mean(group1),
                         mean2 = mean(group2),
                         p = res$p.value)
    
  return(result)
  
}


result_treated_untreated <- lapply(cytokines, test_differential_wilcox, DATA = df_mean, GROUP1 = "NK-treated", GROUP2 = "Untreated", PAIRED = T) %>% bind_rows
result_treated_untreated$p_adj <- p.adjust(result_treated_untreated$p, method = "BH")
result_treated_untreated <- result_treated_untreated %>% arrange(p)


# Test differential abundance with paired tests using means of replicates with NK median added to untreated

df_mean_nk_median_untreated <- df_mean %>%
  filter(condition == "Untreated") %>% 
  left_join(nk_median) %>% 
  mutate(concentration_mean_withnk = concentration_mean + nk_median_concentration)

df_mean_nktreated <- df_mean %>% 
  filter(condition == "NK-treated")

df_mean_nk_median <- rbind(df_mean_nktreated, df_mean_nk_median_untreated) %>%
  mutate(concentration_mean = ifelse(condition == "Untreated", concentration_mean_withnk, concentration_mean))

result_treated_untreated_withnk <- lapply(cytokines, test_differential_wilcox, DATA = df_mean_nk_median, GROUP1 = "NK-treated", GROUP2 = "Untreated", PAIRED = T) %>% bind_rows
result_treated_untreated_withnk$p_adj <- p.adjust(result_treated_untreated_withnk$p, method = "BH")
result_treated_untreated_withnk <- result_treated_untreated_withnk %>% arrange(p)



## ----------------------------------------------

# Correlate cytokines with scRNA-seq clusters

combined_seurat <- readRDS("results/celllinepanel_integrated/celllinepanel_nk.rds")


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


# 3 donors merged

min(table(combined_seurat_3donors$hash.ID))
Idents(combined_seurat_3donors) <- "hash.ID"

combined_seurat_3donors$hash.ID.grouped <- gsub("[0-9]$", "", combined_seurat_3donors$hash.ID)


# cluster fractions

Idents(combined_seurat_3donors) <- "hash.ID"

cluster_freq <- data.frame(table(Idents(combined_seurat_3donors), combined_seurat_3donors$predicted.cluster))
cluster_freq <- cluster_freq %>%
  filter(Var1 != "697") %>% 
  group_by(Var1, Var2) %>%
  summarize(count = sum(Freq)) %>%
  mutate(freq = count / sum(count))

cluster_freq_mean <- cluster_freq %>% 
  mutate(condition = gsub("^NK.*", "NK cells only", gsub("-NK.", "", Var1))) %>% 
  group_by(condition, Var2) %>% 
  summarize(mean_freq = mean(freq)) %>% 
  rename(cluster = Var2)

sample_order <- cluster_freq_mean %>%
  filter(cluster %in% c("Resting (0)", "Adaptive (1)")) %>%
  group_by(condition) %>% 
  summarize(freq_sum = sum(mean_freq)) %>%
  arrange(freq_sum) %>%
  dplyr::select(condition) %>%
  tibble::deframe()

cluster_freq_pert <- cluster_freq_mean %>%
  mutate(condition = factor(condition, levels = rev(sample_order))) %>% 
  mutate(mean_freq = 100*mean_freq)

cluster_freq_pert_resting_adaptive <- cluster_freq_mean %>%
  filter(cluster %in% c("Resting (0)", "Adaptive (1)")) %>%
  group_by(condition) %>% 
  summarize(mean_freq = sum(mean_freq)) %>% 
  mutate(cluster = "Resting + Adaptive") %>% 
  select(condition, cluster, mean_freq)

cluster_freq_pert <- rbind(cluster_freq_pert, cluster_freq_pert_resting_adaptive)

# merge with cytokine differential abundance

result_all_withnk_scrna <- result_all_withnk %>%
  left_join(cluster_freq_pert, by = c("cell_line" = "condition"))

result_all_withnk_scrna$cancer_type <- annot$Subtype_main[match(result_all_withnk_scrna$cell_line, annot$cell_line)]
result_all_withnk_scrna$cancer_type[result_all_withnk_scrna$cell_line=="NK cells only"] <- "NK cells only"

fwrite(result_all_withnk_scrna, "cytokines_log2fc_nkclusters.txt", sep = "\t")


correlate_cytokines <- function(CYTOKINE, DATA){
  
  vector1 <- DATA$mean_freq[DATA$cytokine == CYTOKINE]
  vector2 <- DATA$log2fc[DATA$cytokine == CYTOKINE]

  res <- corr.test(vector1, vector2, adjust = "BH", method = "pearson")
  
  result <- data.frame(cytokine = CYTOKINE,
                       cor = res$r,
                       p = res$p)
  
  return(result)
  
}

loop_clusters <- function(CLUSTER){
  
  df_cytokine <- result_all_withnk_scrna[result_all_withnk_scrna$cluster == CLUSTER,]
  
  res_cytokine <- lapply(cytokines, correlate_cytokines, DATA = df_cytokine) %>% bind_rows()
  res_cytokine$p_adj  <- p.adjust(res_cytokine$p, method = "BH")
  res_cytokine$cluster  <- CLUSTER
  res_cytokine <- res_cytokine %>% arrange(p)
  return(res_cytokine)
  
}

clusters <- unique(cluster_freq_pert$cluster)

result_cor_cytokines <- lapply(clusters, loop_clusters) %>% bind_rows()

fwrite(result_cor_cytokines, "cytokines_correlation_nkclusters.txt", sep = "\t")

  
  
# Plot scatter plot

cytokine_scatterplot <- function(CYTOKINE){
  
  df_plot_scatter <- result_all_withnk_scrna %>% filter(cluster != "Resting + Adaptive") %>% filter(cytokine == CYTOKINE)
  
  ggpubr::ggscatter(df_plot_scatter, x = "mean_freq", y = "log2fc", shape = 21, size = 3, fill = "cancer_type",
                    add = "reg.line", cor.coef = T, cor.method = "pearson", label = "cell_line", repel = T, cor.coef.coord = c(0, max(df_plot_scatter$log2fc)+((max(df_plot_scatter$log2fc)-min(df_plot_scatter$log2fc))/10))) +
    ylab("Cytokine log2 fold change\nNK-treated vs. untreated") +
    xlab("% of NK cells in cluster") +
    facet_wrap(. ~ cluster, scales = "free_x", ncol = 5) +
    scale_fill_manual(values = cols_heatmap_vector) +
    theme_bw() +
    theme(legend.title = element_text(face = "plain"),
          axis.text = element_text(color = "black"),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggtitle(CYTOKINE)
  
}
  

p1 <- lapply(cytokine_names, cytokine_scatterplot)
m1 <- marrangeGrob(p1, nrow = 3, ncol = 1, top = NULL)
ggsave("cytokine_scatterplots_nkclusters.pdf", m1, height = 15, width = 20)



## -----------------------------------------

# Heatmap of cytokine log2FC ordered by correlation with NK activation

result_all_withnk <- fread("cytokines_differential_nk_added.txt")
result_cor_cytokines <- fread("cytokines_correlation_nkclusters.txt")


# top 10 cytokines correlating with NK activation
top5_resting_adaptive <- result_cor_cytokines %>% filter(cluster == "Resting + Adaptive") %>% slice_min(p, n = 5) %>% select(cytokine) %>% tibble::deframe()

# cytokine NK treatment lof2FC (with NK baseline added)
mat_cytokine <- result_all_withnk %>% filter(cytokine %in% top5_resting_adaptive) %>% pivot_wider(names_from = cell_line, values_from = log2fc, id_cols = cytokine) %>% as.data.frame()
rownames(mat_cytokine) <- mat_cytokine$cytokine
mat_cytokine$cytokine <- NULL
mat_cytokine <- mat_cytokine[top5_resting_adaptive,rev(sample_order[sample_order != "NK cells only"])]

mat_cytokine_scaled <- t(apply(mat_cytokine, 1, scale))
rownames(mat_cytokine_scaled) <- rownames(mat_cytokine)
colnames(mat_cytokine_scaled) <- colnames(mat_cytokine)

ht <- Heatmap(mat_cytokine_scaled,
              cluster_rows = F,
              cluster_columns = F,
              border_gp = gpar(col = "black", lty = 1),
              col = pals::ocean.deep(9),
              use_raster = F,
              show_column_names = F,
              show_row_names = T,
              row_names_side = "left",
              show_heatmap_legend = T,
              raster_quality = 3,
              heatmap_legend_param = list(title = "Scaled concentration\nlog2 fold change",
                                          title_gp = gpar(fontsize = 10),
                                          labels_gp = gpar(fontsize = 10),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"),
                                          tick_length = unit(0, "mm"),
                                          border = "black",
                                          title_position = "topcenter",
                                          legend_direction = "horizontal"))

pdf("cytokines_heatmap_nkactivation_manuscript.pdf", width = 6, height = 1)
draw(ht, heatmap_legend_side = "left")
dev.off()


# Plot unscaled values
ht <- Heatmap(mat_cytokine,
              cluster_rows = F,
              cluster_columns = F,
              border_gp = gpar(col = "black", lty = 1),
              col = pals::ocean.deep(9),
              use_raster = F,
              show_column_names = F,
              show_row_names = T,
              row_names_side = "left",
              show_heatmap_legend = T,
              raster_quality = 3,
              heatmap_legend_param = list(title = "Concentration\nlog2 fold change",
                                          title_gp = gpar(fontsize = 10),
                                          labels_gp = gpar(fontsize = 10),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"),
                                          tick_length = unit(0, "mm"),
                                          border = "black",
                                          title_position = "topcenter",
                                          legend_direction = "horizontal"))

pdf("cytokines_heatmap_nkactivation_manuscript_unscaled.pdf", width = 6, height = 1)
draw(ht, heatmap_legend_side = "left")
dev.off()


# Correlate cytokines with scRNA-seq target cells core NK response score FC
nkresponse <- fread("results/celllinepanel_combined/targets/targets_differential_modulescore_mainclusters.txt", data.table = F)


# merge with cytokine differential abundance

result_all_withnk_scrna_nkresponse <- result_all_withnk_scrna %>%
  left_join(nkresponse, by = ("cell_line"))

correlate_cytokines <- function(CYTOKINE, DATA){

  vector1 <- DATA$avg_log2FC[DATA$cytokine == CYTOKINE]
  vector2 <- DATA$log2fc[DATA$cytokine == CYTOKINE]
  
  res <- corr.test(vector1, vector2, adjust = "BH", method = "pearson")
  
  result <- data.frame(cytokine = CYTOKINE,
                       cor = res$r,
                       p = res$p)
  
  return(result)
  
}

result_cor_cytokines <- lapply(cytokines, correlate_cytokines, DATA = result_all_withnk_scrna) %>% bind_rows()
result_cor_cytokines$p_adj  <- p.adjust(result_cor_cytokines$p, method = "BH")
result_cor_cytokines <- result_cor_cytokines %>% arrange(p)

fwrite(result_cor_cytokines, "cytokines_correlation_corenkresponse.txt", sep = "\t")

## ----------------------------------------------------------

# Time point data bar plots

data_clean_sample <- data_clean[grepl("h$|B4", rownames(data_clean)),]
data_clean_sample$sample <- rownames(data_clean_sample)

df <- data_clean_sample %>% tidyr::pivot_longer(cols = !contains("sample"), names_to = "cytokine", values_to = "concentration")

df <- df %>% mutate(cell_line = gsub("_.*", "", gsub("NK1_B4", "NK only", sample)),
                    timepoint = gsub(".*NK1_", "", gsub(".*B4", "0h", sample)),
                    sample = factor(sample, levels = rownames(data_clean_sample))) %>% 
  mutate(cell_line = factor(cell_line, levels = c("K562", "JURKAT", "RI1", "GDM1", "NK only")),
         timepoint = factor(timepoint, levels = c("0h", "1h", "3h", "6h", "12h", "24h")))



df$cancer_type <- annot$Subtype_main[match(df$cell_line, annot$cell_line)]
df$cancer_type[df$cell_line=="NK cells only"] <- "NK cells only"

# supplement table
df_timepoint_supplement <- df %>% pivot_wider(id_cols = c("cell_line", "timepoint"), values_from = concentration, names_from = cytokine)

fwrite(df_timepoint_supplement, "cytokines_timecourse.txt", sep = "\t")

ggplot(df, aes(x = cell_line, y = concentration, fill = timepoint)) +
  geom_bar(position = position_dodge2(preserve = "single"), stat = "identity") +
  facet_wrap(. ~ cytokine, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.y = element_blank() ,
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line( size=.1, color="grey50" ),
        strip.background = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = rev(pals::ocean.matter(7)[c(1:6)])) +
  ylab("Concentration (pg/ml)") +
  xlab("") +
  labs(fill = "")
  
ggsave("cytokine_timepoint_barplots.pdf", height = 10, width = 10)


# line plot (Figure 1I)

df_plot <- df %>%
  filter(cytokine == "IFNg") %>% 
  mutate(timepoint = as.numeric(gsub("h", "", timepoint)))

ggplot(df_plot, aes(x = timepoint, y = concentration, color = cell_line, group = cell_line)) +
  geom_smooth(se = F, method = "loess", span = 0.6) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12, 24), expand = c(0,0)) +
  scale_color_manual("", values = rev(cartography::carto.pal("wine.pal", 6)[2:5])) +
  ylab("IFN-y (pg/ml)") +
  xlab("Time (h)") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.justification = 0.5)+
  guides(color = "none") +
  geom_text_repel(data = df_plot[df_plot$timepoint=="24",], aes(label = cell_line), nudge_y = 2)

ggsave("cytokine_ifng_timepoint_lineplot.pdf", width = 6, height = 1.97)



