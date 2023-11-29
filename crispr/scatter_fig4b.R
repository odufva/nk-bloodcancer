
# Scatter plot of NK cell CRISPR screen logFC (Figure 4B)

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
library(RColorBrewer)
library(ggnewscale)
library(biomaRt)

# load screen results
crispr <- fread("crispr_mageck_combined.txt", data.table = F)

# remove controls
crispr <- crispr[!grepl("Control|ORGENES|-Mar", crispr$gene),]

# add random number for genes for plotting
set.seed(24)
position <- data.frame(gene = unique(crispr$gene), position_rank = sample(c(1:length(unique(crispr$gene)))))
crispr <- merge(crispr, position)


plotdata <- crispr %>%
  mutate(signed_p = -log10(p) * sign(lfc)) %>%
  mutate(cell_line = factor(cell_line, levels = c("K562", "MOLM14", "SUDHL4", "NALM6", "MM1S", "LP1", "KMS11", "MM1S_GOF", "KMS11_GOF", "KMS11_KHYG1_GOF", "LP1_KHYG1_GOF", "KMS11_KHYG1")))

activating_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(lfc>0)

activating_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(lfc<0)

inhibitory_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(lfc<0)

inhibitory_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(lfc>0)

activating <- rbind(activating_lof, activating_gof)
inhibitory <- rbind(inhibitory_lof, inhibitory_gof)

plotdata_activating <- activating %>%
  group_by(gene) %>%
  top_n(1, desc(p)) %>%
  mutate(signed_p = -log10(p) * sign(lfc)) %>%
  filter(p < 0.0001, abs(lfc) > 0.75) %>%
  mutate(cell_line = factor(cell_line, levels = c("K562", "MOLM14", "SUDHL4", "NALM6", "MM1S", "LP1", "KMS11", "MM1S_GOF", "KMS11_GOF", "KMS11_KHYG1_GOF", "LP1_KHYG1_GOF", "KMS11_KHYG1")))

plotdata_inhibitory <- inhibitory %>%
  group_by(gene) %>%
  top_n(1, desc(p)) %>%
  mutate(signed_p = -log10(p) * sign(lfc)) %>%
  filter(p < 0.0001, abs(lfc) > 0.75) %>%
  mutate(cell_line = factor(cell_line, levels = c("K562", "MOLM14", "SUDHL4", "NALM6", "MM1S", "LP1", "KMS11", "MM1S_GOF", "KMS11_GOF", "KMS11_KHYG1_GOF", "LP1_KHYG1_GOF", "KMS11_KHYG1")))

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

colors <- c("#F54F29", "#FF974F", "#FFD393", "#9C9B7A", "#AD5472", "#F07C6C", "#C0CA55", "#6C9380", "grey80")

colors <- c(LaCroixColoR::lacroix_palette("PeachPear", 6), LaCroixColoR::lacroix_palette("PassionFruit", 6)[c(1,3)], "grey70")
colors[3] <- brewer.pal(9, "YlOrRd")[3]

names(colors) <- c("Activating receptor ligand", "Death receptor apoptosis", "TF/epigenetic", "Adhesion", "Antigen presentation", "IFNy",  "Mucin", "Fucosylation", "Other")

dat <- rbind(plotdata_activating, plotdata_inhibitory) %>%
  mutate(signed_p = ifelse(grepl("GOF", cell_line), signed_p * -1, signed_p)) %>% 
  mutate(category = ifelse(gene %in% apc, "Antigen presentation",
                           ifelse(gene %in% ifny, "IFNy",
                                  ifelse(gene %in% activating_receptors, "Activating receptor ligand",
                                         ifelse(gene %in% death_receptor_apoptosis, "Death receptor apoptosis",
                                                ifelse(gene %in% adhesion, "Adhesion",
                                                       ifelse(gene %in% fucosylation, "Fucosylation",
                                                              ifelse(gene %in% mucin, "Mucin",
                                                                     ifelse(gene %in% tf_chromatin, "TF/epigenetic", "Other")))))))))

ggplot(dat, aes(x = position_rank, y = signed_p, size = abs(lfc))) +
  ggrastr::geom_point_rast(data = plotdata[sample(nrow(plotdata), 20000),], size = 0.1, color = "grey70") +
  geom_point(data = dat[dat$category == "Other",], aes(color = category), pch = 16) +
  geom_point(data = dat[dat$category != "Other",], aes(color = category), pch = 16) +
  geom_text_repel(aes(label = gene),
                  size = 4,
                  color = "black",
                  max.overlaps = 30,
                  box.padding = 0.5,
                  point.padding = 1,
                  fontface = "italic") +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(-7, 7)) +
  scale_color_manual("", values = colors) +
  scale_fill_manual("", values = colors) +
  scale_size_area(max_size = 8) +
  ylab("Signed P value (-log10)") +
  xlab("Gene") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 3, title.position = "top", nrow = 3),
         size = guide_legend(override.aes = list(pch = 16), title.position = "top",
                             title = "Fold change (log2)"))

ggsave("nk_crispr_scatter.pdf", height = 10, width = 10)

