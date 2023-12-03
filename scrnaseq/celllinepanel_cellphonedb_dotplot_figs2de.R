
# NK cell co-culture with panel of 26 cell lines
# Plot cellphonedb results (Figure S2D-E)

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

dir.create("results/celllinepanel/cellphonedb/")

# load significant means
data <- fread("data/significant_means.txt", data.table = F)

data_long_supplement <- data %>% 
  select(1:12, contains("|effector")) %>% 
  tidyr::pivot_longer(cols = contains("effector"), names_to = "interacting_cells", values_to = "mean") %>% 
  filter(!is.na(mean)) %>% 
  mutate(interacting_cells = gsub("\\ |_$|-like|-producing", "", interacting_cells))
fwrite(data_long_supplement, "results/celllinepanel/cellphonedb/interactions_nkexpanded.txt", row.names = F, quote = F, sep = "\t")

data_long <- data %>% 
  filter(annotation_strategy == "curated") %>% # only curated interactions
  select(1:12, contains("|effector")) %>% 
  tidyr::pivot_longer(cols = contains("effector"), names_to = "interacting_cells", values_to = "mean") %>% 
  mutate(interacting_cells = gsub("_$|-like|-producing", "", interacting_cells))


# Dot plot of interactions
data_plot <- data_long %>% 
  filter(!is.na(mean)) %>% 
  mutate(treatment = ifelse(grepl("NK_expanded", interacting_cells), "NK", "No NK")) %>% 
  mutate(celltype_a = gsub("\\|.*", "", interacting_cells),
         celltype_b = gsub(".*\\|", "", interacting_cells)) %>% 
  filter(!grepl("^effector|target_NA$", celltype_a)) %>% 
  mutate(interacting_cells = gsub("effector_|target_|_NK_expanded| _[0-9]_", "", interacting_cells))
  
ggplot(data_plot, aes(x = interacting_cells, y = interacting_pair, size = mean)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ treatment + celltype_b, scales = "free_x") +
  xlab("") +
  ylab("")

ggsave("results/celllinepanel/cellphonedb/dotplot_allinteractions.pdf", height = 40, width = 40)



## Differential interactions clusters 0 and 1 + untreated targets vs clusters 2-5 + NK-treated targets
data_diff_interacting <- data_plot %>% 
  filter(treatment == "NK" & grepl("Activated|IFN|Cytokine", celltype_b)) %>% 
  group_by(interacting_pair, celltype_b) %>% 
  summarize(count = n()) %>% 
  summarize(mean_count = mean(count))

data_diff_noninteracting <- data_plot %>% 
  filter(treatment == "No NK" & !grepl("Activated|IFN|Cytokine", celltype_b)) %>% 
  group_by(interacting_pair, celltype_b) %>% 
  summarize(count = n()) %>% 
  summarize(mean_count = mean(count))

data_diff <- data_diff_interacting %>% left_join(data_diff_noninteracting, by = c("interacting_pair"), suffix = c("_int", "_noint")) %>% 
  mutate(mean_count_noint = ifelse(is.na(mean_count_noint), 0, mean_count_noint)) %>% 
  mutate(diff_count = mean_count_int - mean_count_noint) %>% 
  arrange(desc(diff_count))

data_diff_plot <- data_diff %>% 
  filter(diff_count > 4)

data_diff_plot$interacting_pair <- factor(data_diff_plot$interacting_pair, levels = unique(data_diff_plot$interacting_pair))


data_plot_diffinteractions <- data_plot %>% 
  filter((treatment == "NK" & grepl("Activated|IFN|Cytokine", celltype_b)) | (treatment == "No NK" & !grepl("Activated|IFN|Cytokine", celltype_b))) %>% 
  filter(interacting_pair %in% unique(data_diff_plot$interacting_pair))

ggplot(data_plot_diffinteractions, aes(x = interacting_cells, y = interacting_pair, size = mean)) +
  geom_point() +
  scale_x_discrete(drop = F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ treatment + celltype_b, scales = "free_x") +
  xlab("") +
  ylab("")


# Differential interactions between clusters 2-4 and 0-1          
data_plot_act <- data_plot %>% 
  filter(grepl("Activated|IFN|Cytokine", celltype_b)) %>% 
  group_by(interacting_pair, celltype_b, treatment) %>% 
  summarize(count = n()) %>% 
  group_by(interacting_pair) %>% 
  summarize(mean_count = mean(count))

data_plot_noact <- data_plot %>% 
  filter(!grepl("Activated|IFN|Cytokine", celltype_b)) %>% 
  group_by(interacting_pair, celltype_b, treatment) %>% 
  summarize(count = n()) %>% 
  group_by(interacting_pair) %>% 
  summarize(mean_count = mean(count))

data_diff_act_clusters <- data_plot_act %>% left_join(data_plot_noact, by = c("interacting_pair"), suffix = c("_act", "_noact")) %>% 
  mutate(mean_count_noact = ifelse(is.na(mean_count_noact), 0, mean_count_noact)) %>% 
  mutate(diff_count = mean_count_act - mean_count_noact) %>% 
  arrange(desc(diff_count)) 

data_diff_act_clusters_plot <- data_diff_act_clusters %>% 
  filter(diff_count > 4) %>% 
  select(interacting_pair) %>% 
  tibble::deframe()

# Test significance
test_fisher <- function(interaction){

  df <- as.data.frame(data_diff_act_clusters)
  rownames(df) <- data$interacting_pair
  
  counts <- matrix(c(df[interaction, "mean_count_act"],
                     df[interaction, "mean_count_noact"],
                     26-df[interaction, "mean_count_act"],
                     26-df[interaction, "mean_count_noact"]),
                   nrow = 2,
                        dimnames =
                          list(c("Activated clusters", "Non-activated clusters"),
                               c("Interaction", "No interaction")))
  
  res <- fisher.test(counts)
  result <- data.frame(interacting_pair = interaction, p = res$p.value)
  return(result)
  
}

result_clusters <- lapply(unique(data_diff_act_clusters$interacting_pair), test_fisher) %>% bind_rows()


# Differential interactions between NK-treated vs untreated
data_plot_act <- data_plot %>% 
  filter(treatment=="NK") %>% 
  group_by(interacting_pair, celltype_b, treatment) %>% 
  summarize(count = n()) %>% 
  group_by(interacting_pair) %>% 
  summarize(mean_count = mean(count))

data_plot_noact <- data_plot %>% 
  filter(treatment=="No NK") %>% 
  group_by(interacting_pair, celltype_b, treatment) %>% 
  summarize(count = n()) %>% 
  group_by(interacting_pair) %>% 
  summarize(mean_count = mean(count))

data_diff_act_treatment <- data_plot_act %>% left_join(data_plot_noact, by = c("interacting_pair"), suffix = c("_act", "_noact")) %>% 
  mutate(mean_count_noact = ifelse(is.na(mean_count_noact), 0, mean_count_noact)) %>% 
  mutate(diff_count = mean_count_act - mean_count_noact) %>% 
  arrange(desc(diff_count)) 

data_diff_act_treatment_plot <- data_diff_act_treatment %>% 
  filter(diff_count > 4) %>% 
  select(interacting_pair) %>% 
  tibble::deframe()

plotgenes <- c(intersect(data_diff_act_clusters_plot, data_diff_act_treatment_plot),
setdiff(data_diff_act_clusters_plot, data_diff_act_treatment_plot),
setdiff(data_diff_act_treatment_plot, data_diff_act_clusters_plot))

data_plot_diffinteractions <- data_plot %>% 
  group_by(treatment, celltype_b, interacting_pair) %>% 
  summarize(count = n()) %>% 
  filter(interacting_pair %in% plotgenes) %>% 
  mutate(interacting_pair = gsub("_", " | ", interacting_pair)) %>% 
  mutate(interacting_pair = factor(interacting_pair, levels = rev(gsub("_", " | ", plotgenes)))) %>% 
  mutate(celltype_b = gsub("\\ _.*", "", gsub("effector_", "", celltype_b))) %>% 
  mutate(celltype_b = factor(celltype_b, levels = c("Activated", "Type I IFN", "Cytokine", "Adaptive", "Resting"),
                             labels = c("Activated", "Type I IFN", "Cytokine", "Adaptive", "Resting"))) %>% 
mutate(interacting_pair_a = gsub("\\|.*", "", interacting_pair),
       interacting_pair_b = gsub(".*\\|", "", interacting_pair))


ggplot(data_plot_diffinteractions, aes(x = treatment, y = interacting_pair, size = count, color = count)) +
  geom_point() +
  facet_grid(. ~ celltype_b, scales = "free_x") +
  xlab("") +
  ylab("") +
  scale_color_distiller("Number of cell lines\nwith interaction", palette = "Reds", direction = 1) +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        legend.margin = margin(-1,-1,-1,-1),
        legend.key.height = unit(c(0.3), "cm"),
        legend.box = "vertical",
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = "black", size = 11),
        strip.text = element_text(angle = 90, hjust = 0, size = 11),
        plot.margin = unit(c(0,1,0,0), "cm")) +
  guides(size = guide_legend(title = ""))
  
ggsave("results/celllinepanel/cellphonedb/dotplot_top_differential_interactions.pdf", height = 11, width = 5)


# separate plots for NK clusters and target treatment

# NK (Figure S2D)
pair_order <- c("TNFSF9 | TNFRSF9", "TNFSF4 | TNFRSF4",
                "TNFRSF10A | TNFSF10", "TNFRSF10B | TNFSF10", "FAS | FASLG",
                "LGALS9 | HAVCR2", "PVR | TIGIT", "CD47 | SIRPG",
                "TGFB1 | TGFBR2", "TGFB3 | TGFBR2",
                "EPHB4 | EFNB1", "EPHB2 | EFNB1",  "EPHB2 | EFNA5", "EFNB2 | EPHA4", "EFNA3 | EPHA4", "EFNA4 | EPHA4", "EFNA5 | EPHA4", 
                "ICAM1 | aMb2 complex", "ICAM1 | aXb2 complex", "JAM3 | aMb2 complex",
                "FCER2 | aMb2 complex", "FCER2 | aXb2 complex", "COL9A2 | a1b1 complex", "COL27A1 | a1b1 complex", "C3 | aMb2 complex",
                "NRG2 | ERBB3", "CD55 | ADGRE5", "TNFRSF17 | TNFSF13B", "CD28 | CD80",
                "NOTCH2 | JAG2", "NOTCH4 | JAG2", 
                "HLA-C | KIR2DL1")

data_plot_diffinteractions_nk <- data_plot %>% 
  group_by(treatment, celltype_b, interacting_pair) %>% 
  summarize(count = n()) %>% 
  filter(interacting_pair %in% data_diff_act_clusters_plot) %>% 
  mutate(interacting_pair = gsub("_", " | ", gsub("beta receptor2", "BR2", interacting_pair))) %>% 
  mutate(interacting_pair = factor(interacting_pair, levels = rev(gsub("_", " | ", gsub("beta receptor2", "BR2", data_diff_act_clusters_plot))))) %>% 
  mutate(interacting_pair = factor(interacting_pair, levels = rev(pair_order))) %>% 
  mutate(celltype_b = gsub("\\ _.*|_.*", "", gsub("effector_", "", celltype_b))) %>% 
  mutate(celltype_b = factor(celltype_b, levels = c("Activated", "Type I IFN", "Cytokine", "Adaptive", "Resting"),
                             labels = c("Activated", "Type I IFN", "Cytokine", "Adaptive", "Resting"))) %>% 
  mutate(interacting_pair_a = gsub("\\|.*", "", interacting_pair),
         interacting_pair_b = gsub(".*\\|", "", interacting_pair))


ggplot(data_plot_diffinteractions_nk, aes(x = treatment, y = interacting_pair, size = count, color = count)) +
  geom_point() +
  facet_grid(. ~ celltype_b, scales = "free_x") +
  xlab("") +
  ylab("") +
  scale_color_distiller("Number of cell lines\nwith interaction", palette = "Reds", direction = 1) +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.margin = margin(-1,-1,-1,-1),
        legend.key.width = unit(c(0.3), "cm"),
        legend.box = "vertical",
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = "black", size = 11),
        strip.text = element_text(angle = 90, hjust = 0, size = 11),
        plot.margin = unit(c(0,1,0,0), "cm")) +
  guides(size = guide_legend(title = "Number of cell lines\nwith interaction"),
         color = guide_legend(title = ""))

ggsave("results/celllinepanel/cellphonedb/dotplot_top_differential_interactions_nk.pdf", height = 7, width = 7)



# target (Figure S2E)
pair_order <- c("HLA-E | KLRC1", "HLA-E | KLRC2", "HLA-E | KLRK1", "HLA-B | KIR3DL2", "HLA-C | KIR2DL3",
                "HLA-F | KIR3DL2", "HLA-F | LILRB1", "FAS | FASLG", "CXCL10 | CXCR3",
                "TGFB1 | TGFBR1", "TGFB3 | TGFBR2", "TGFB3 | TGFBR3", "TGFB1 | TGFBR3", "TGFB1 | TGFBR2",
                "TNFSF9 | TNFRSF9", "CLEC2B | KLRF1", "ICOSLG | ICOS", "LGALS9 | HAVCR2", "NECTIN2 | NECTIN3",
                "ICAM1 | aMb2 complex", "ICAM3 | aLb2 complex", "EFNA3 | EPHA4", "EPHB4 | EFNB1")

data_plot_diffinteractions_target <- data_plot %>% 
  group_by(treatment, celltype_b, interacting_pair) %>% 
  summarize(count = n()) %>% 
  filter(interacting_pair %in% data_diff_act_treatment_plot) %>% 
  mutate(interacting_pair = gsub("_", " | ", gsub("beta receptor", "BR", interacting_pair))) %>% 
  mutate(interacting_pair = factor(interacting_pair, levels = rev(gsub("_", " | ", pair_order)))) %>% 
  mutate(celltype_b = gsub("\\ _.*", "", gsub("effector_", "", celltype_b))) %>% 
  mutate(celltype_b = factor(celltype_b, levels = c("Activated", "Type I IFN", "Cytokine", "Adaptive", "Resting"),
                             labels = c("Activated", "Type I IFN", "Cytokine", "Adaptive", "Resting"))) %>% 
  mutate(interacting_pair_a = gsub("\\|.*", "", interacting_pair),
         interacting_pair_b = gsub(".*\\|", "", interacting_pair))


ggplot(data_plot_diffinteractions_target, aes(x = celltype_b, y = interacting_pair, size = count, color = count)) +
  geom_point() +
  facet_grid(. ~ treatment, scales = "free_x") +
  xlab("") +
  ylab("") +
  scale_color_distiller("Number of cell lines\nwith interaction", palette = "Reds", direction = 1) +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.margin = margin(-1,-1,-1,-1),
        legend.key.width = unit(c(0.3), "cm"),
        legend.box = "vertical",
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color = "black", size = 11),
        strip.text = element_text(angle = 0, hjust = 0.5, size = 11),
        plot.margin = unit(c(0,1,0,0), "cm")) +
  guides(size = guide_legend(title = "Number of cell lines\nwith interaction"),
         color = guide_legend(title = ""))

ggsave("results/celllinepanel/cellphonedb/dotplot_top_differential_interactions_target.pdf", height = 6, width = 7)


