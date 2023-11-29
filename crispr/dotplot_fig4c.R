
# Dot plot of NK cell CRISPR screen logFC (Figure 4C)

# load libraries
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(data.table)
library(ggsci)
library(ggrepel)
library(viridis)
library(fgsea)
library(dplyr)
library(tibble)
library(ggpointdensity)
library(RColorBrewer)
library(ggnewscale)
library(wesanderson)

# load screen results
crispr <- fread("crispr_mageck_combined.txt", data.table = F)

# gene lists based on categories
apc <- c("HLA-E", "HLA-A", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5", "RFXAP", "RFXANK")
ifny <- c("JAK1", "STAT1")
activating_receptors <- c("NCR3LG1", "CD58", "ULBP1", "ULBP2", "ULBP3", "CD48", "TNFSF9", "PVR", "MICA")
death_receptor_apoptosis <- c("TNFRSF1B", "TNFRSF10A", "TNFRSF10B", "TNFRSF10D", "FAS", "FADD", "BID", "BAX", "PMAIP1", "CASP8", "CFLAR", "TRAF2", "RNF31", "XIAP", "BIRC3")
adhesion <- c("ICAM1", "SPN", "SELPLG", "CD44")
fucosylation <- c("GMDS", "FUT8", "SLC35C1")
mucin <- c("MUC1", "MUC21")
tf_chromatin <- c("ARID1A", "GFI1B", "TAL1", "YTHDF2", "RUNX1", "CHD7", "CMIP", "PCGF5", "RBBP4", "FOXA1", "IRF4", "MYB", "MSI2", "STAG2", "KLF16")
other <- c("PTEN", "NFKBIA", "GNA13", "ZC3HAV1", "GSK3B", "SPPL3", "KIAA0922", "ARHGAP1", "MAPK14", "HMBG1")

genelist <- c(apc, ifny, activating_receptors, death_receptor_apoptosis,adhesion, fucosylation, mucin, tf_chromatin, other)

crispr <- crispr %>% filter(gene %in% genelist)

p_threshold <- 0.001
lfc_threshold_pos <- 0.75
lfc_threshold_neg <- -0.75
lfc_threshold <- 0.75

activating_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc>lfc_threshold_pos)

activating_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc<lfc_threshold_neg)

inhibitory_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc<lfc_threshold_neg)

inhibitory_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc>lfc_threshold_pos)

activating <- rbind(activating_lof, activating_gof)
inhibitory <- rbind(inhibitory_lof, inhibitory_gof)

activating_top <- activating %>%
  filter(p < p_threshold, abs(lfc)>lfc_threshold) %>%
  group_by(cell_line) %>%
  top_n(20, desc(p)) %>%
  ungroup() %>%
  dplyr::select(gene) %>%
  unique() %>%
  tibble::deframe()

inhibitory_top <- inhibitory %>%
  filter(p < p_threshold, abs(lfc)>lfc_threshold) %>%
  group_by(cell_line) %>%
  top_n(20, desc(p)) %>%
  ungroup() %>%
  dplyr::select(gene) %>%
  unique() %>%
  tibble::deframe()

genelist_top <- unique(c(activating_top, inhibitory_top))


plotdata <- crispr %>%
  filter(gene %in% genelist) %>%
  filter(p < 0.05) %>%
  mutate(gene = factor(gene, levels = genelist)) %>%
  mutate(lof_gof = ifelse(grepl("GOF", cell_line), "GOF", "LOF")) %>% 
  mutate(lof_gof = factor(lof_gof, levels = c("LOF", "GOF"))) %>% 
  mutate(cell_line = gsub("_", " ", gsub("_GOF", "", cell_line))) %>% 
  mutate(cell_line = factor(cell_line, levels = c("K562", "MOLM14", "SUDHL4", "NALM6", "MM1S", "KMS11" ,"LP1", "LP1 KHYG1"),
                            labels = c("K562", "MOLM14", "SUDHL4", "NALM6", "MM1S", "KMS11" ,"LP1", "LP1 (KHYG1)")))

ggplot(plotdata, aes(x = cell_line, y = gene, size = -log10(p))) +
  geom_point(aes(fill = lfc), pch = 21) +
  scale_fill_distiller("Fold change log2)\nLOF screens", palette = "RdBu", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                        type = "div", limits = max(abs(plotdata$lfc[!plotdata$lof_gof == "GOF"])) * c(-1, 1), guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +#, limits = LIMITS) +
  ggnewscale::new_scale_fill() +
  geom_point(data = plotdata[plotdata$lof_gof == "GOF",], aes(fill = lfc), pch = 21) +
  scale_fill_distiller("Fold change (log2)\nGOF screens", palette = "PiYG", values = c(0, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1),
                        type = "div", limits = max(abs(plotdata$lfc[plotdata$lof_gof == "GOF"])) * c(-1, 1), direction = 1, guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +#, limits = LIMITS) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        legend.key.width = unit(0.3, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        plot.title = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed")) +
  labs(size = "P value (-log10)") +
  guides(size = guide_legend(override.aes = list(pch = 16), title.position = "top", title.hjust = 0.5)) +
  ylab("") +
  xlab("") +
  facet_grid(. ~ lof_gof, scales = "free_x", space = "free_x")

ggsave("nk_crispr_dotplot.pdf", height = 20, width = 3.75)



