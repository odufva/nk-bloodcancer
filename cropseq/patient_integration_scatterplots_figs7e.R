
# Match DEGs of all NK cell CROP-seq experiments to DEGs from patients with vs without mutations in perturbed genes (Figure S7E)

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(viridis)
library(ggrastr)


dir.create("results/combine/patient_integration")

# load results
nalm6 <- fread("results/nalm6/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4 <- fread("results/sudhl4/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
k562 <- fread("results/k562/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "K562")
mm1s <- fread("results/mm1s/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S")
lp1 <- fread("results/lp1/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "LP1")

data <- rbind(k562, sudhl4, nalm6, mm1s, lp1)

data <- data %>% mutate(perturbation_cell_line = paste(perturbation, cell_line, sep = " "))

# load patient DEG lists (mut vs WT)

files <- list.files("../NK_resistance Heme CollabPaper/Analysis", recursive = T, pattern = "DEG.txt", full.names = T)
files <- files[!grepl("vs", files)]
read_files <- function(FILE){
  mutation  <- gsub("_.*", "", gsub(".*results_DE.", "", FILE))
  dataset  <- sub(".*?_", "", gsub("_DEG.txt*", "", gsub(".*results_DE.", "", FILE)))
  fread(FILE, data.table = F) %>% mutate(perturbation = mutation, dataset = dataset) %>% rename(gene = Gene)
}

deg <- lapply(files, read_files) %>% bind_rows()

data_deg <- merge(data, deg, by = c("gene", "perturbation"))

integrated <- data_deg %>%
  filter(sign(avg_log2FC) == sign(logFC)) %>%
  filter(p_val < 0.05) %>%
  filter(P.Value < 0.05) %>%
  arrange(p_val)

fwrite(integrated, "results/combine/patient_integration/cropseq_patient_integration.txt", sep = "\t")

integrated_supplement <- integrated
colnames(integrated_supplement)[2:10] <- paste(colnames(integrated_supplement)[2:10], "CROPseq", sep = "_")
colnames(integrated_supplement)[11:17] <- paste(colnames(integrated_supplement)[11:17], "patient", sep = "_")
fwrite(integrated_supplement, "results/combine/patient_integration/cropseq_patient_integration_supplement.txt", sep = "\t")

integrated_top <- data_deg %>%
  filter(sign(avg_log2FC) == sign(logFC)) %>%
  filter(p_val < 0.05) %>%
  filter(P.Value < 0.05) %>%
  group_by(gene, perturbation_cell_line, dataset) %>%
  top_n(1, desc(p_val)) %>%
  arrange(p_val)

fwrite(integrated_top, "results/combine/patient_integration/cropseq_patient_integration_topcondition.txt", sep = "\t")


## scatter plots for Figure S7E

# only positive FC

cutoff_cropseq <- 0.25
cutoff_patient <- 1

data_plot <- integrated_top %>%
  filter(!grepl("Chapuy", dataset)) %>%
  filter(avg_log2FC > 0) %>% 
  mutate(dataset = gsub("CoMMpass", "MM", gsub("_Reddy|TCGA_", "", dataset))) %>%
  mutate(label = ifelse(avg_log2FC > cutoff_cropseq & logFC > cutoff_patient, "show_label", "hide_label"))

ggplot(data_plot, aes(x = avg_log2FC, y = logFC, size = -log10(p_val), color = -log10(P.Value))) +
  geom_point_rast(aes(alpha = ifelse(label == "show_label", 1, 1))) +
  geom_point(data = data_plot[data_plot$label=="show_label",], aes(fill = -log10(P.Value)), color = "black", shape = 21) +
  scale_color_viridis(limits = c(1, 12), option = "magma", direction = -1) +
  scale_fill_viridis(limits = c(1, 12), option = "magma", direction = -1) +
  geom_text_repel(data = data_plot[data_plot$label=="show_label",], aes(label = ifelse(label == "show_label",
                                                                                       paste(gene, "\n", perturbation_cell_line, "\n", dataset), "")),
                  size = 2,
                  color = "black",
                  max.overlaps = 100) +
  theme_cowplot() +
  xlab("Fold change (log2)\nCRISPR perturbation vs control") +
  ylab("Fold change (log2)\nPatients with vs without mutation") +
  labs(size = "P value\n(-log10)\nCROP-seq",
       color = "P value\n(-log10)\nPatients") +
  guides(alpha = F, fill = F) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))
  

ggsave("results/combine/patient_integration/scatter_selected_positive.pdf", height = 5, width = 8)

# only neagtive FC

cutoff_cropseq <- 0.2
cutoff_patient <- 0.8

data_plot <- integrated_top %>%
  filter(!grepl("Chapuy", dataset)) %>%
  filter(avg_log2FC < 0) %>% 
  mutate(dataset = gsub("CoMMpass", "MM", gsub("_Reddy|TCGA_", "", dataset))) %>%
mutate(label = ifelse(avg_log2FC < -cutoff_cropseq & logFC < -cutoff_patient, "show_label", "hide_label"))

ggplot(data_plot, aes(x = avg_log2FC, y = logFC, size = -log10(p_val), color = -log10(P.Value))) +
  geom_point_rast(aes(alpha = ifelse(label == "show_label", 1, 1))) +
  geom_point(data = data_plot[data_plot$label=="show_label",], aes(fill = -log10(P.Value)), color = "black", shape = 21) +
  scale_color_viridis(limits = c(1, 12), option = "magma", direction = -1) +
  scale_fill_viridis(limits = c(1, 12), option = "magma", direction = -1) +
  geom_text_repel(data = data_plot[data_plot$label=="show_label",], aes(label = ifelse(label == "show_label",
                                                                                       paste(gene, "\n", perturbation_cell_line, "\n", dataset), "")),
                  size = 2,
                  color = "black",
                  max.overlaps = 100) +
  theme_cowplot() +
  xlab("Fold change (log2)\nCRISPR perturbation vs control") +
  ylab("Fold change (log2)\nPatients with vs without mutation") +
  labs(size = "P value\n(-log10)\nCROP-seq",
       color = "P value\n(-log10)\nPatients") +
  guides(alpha = F, fill = F) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01))


ggsave("results/combine/patient_integration/scatter_selected_negative.pdf", height = 5, width = 8)


