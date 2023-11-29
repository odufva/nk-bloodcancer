
# Plot dot plot combining univariable survival data from different datasets (Figure S5G)

library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
library(data.table)
library(patchwork)

# load univariable survival results
aml_hemap <- fread("result_Hemap_AML.tsv") %>% mutate(Cancer_type = "AML")
aml_beat <- fread("result_beatAML.tsv") %>% mutate(Cancer_type = "AML")
aml_tcga <- fread("result_TCGA_AML.tsv") %>% mutate(Cancer_type = "AML")
mm_hemap <- fread("result_Hemap_MM.tsv") %>% mutate(Cancer_type = "MM")
mm_commpass <- fread("result_CoMMpass.tsv") %>% mutate(Cancer_type = "MM")
dlbcl_hemap_chop <- fread("result_Hemap_DLBCL_RCHOP.tsv") %>% mutate(Cancer_type = "DLBCL")
dlbcl_chapuy <- fread("result_GSE98588.tsv") %>% mutate(Cancer_type = "DLBCL")
dlbcl_reddy <- fread("result_Reddy_DLBCL.tsv") %>% mutate(Cancer_type = "DLBCL")

# load CRISPR screen data
crispr <- fread("crispr_mageck_combined.txt", data.table = F)


data <- rbind(aml_hemap, mm_hemap, aml_beat, aml_tcga, mm_commpass, dlbcl_hemap_chop, dlbcl_chapuy, dlbcl_reddy)

# select significant features
signif_feats <- data %>%
  filter(Adj.P < 0.05) %>% 
  dplyr::select(Feature) %>% 
  unique() %>% 
  tibble::deframe()

# prepare data frame for plotting
df_plot <- data %>%
  filter(Feature %in% signif_feats, Feature != "Age-1") %>% 
  mutate(log10p = -log10(P),
         coef = `exp(coef)`) %>% 
  mutate(log2hr = log2(coef)) %>% 
  mutate(Name = factor(Name, levels = c("CoMMpass", "CoMMpass_bortezomib_IMID", "CoMMpass_bortezomib",  "Hemap_MM", "GSE16716,GSE24080_Hemap_MM", "GSE19784_Hemap_MM", "BeatAML", "TCGA_AML", "Hemap_AML", "GSE10358 AML", "GSE10358,GSE12662 AML", "GSE6891,GSE14468 AML", "GSE12417 AML", "Hemap_DLBCL_RCHOP", "DLBCL_GSE98588", "Reddy_DLBCL"),
                       labels = c("MM (CoMMpass)", "MM (CoMMpass bortezomib+IMID)", "MM (CoMMpass bortezomib)",  "MM (Hemap)", "MM (Hemap GSE16716,GSE24080)", "MM (Hemap GSE19784)", "AML (BeatAML)", "AML (TCGA)", "AML (Hemap)", "AML (Hemap GSE10358)", "AML (Hemap GSE10358,GSE12662)", "AML (Hemap GSE6891,GSE14468", "AML (Hemap GSE12417)", "DLBCL (Hemap RCHOP)", "DLBCL (Chapuy et al.)", "DLBCL (Reddy et al.)"))) %>% 
  mutate(Type = ifelse(grepl("HLA", Feature), "CRISPRhit",
                       ifelse(grepl("GCB|ABC", Feature), "Subtype", Type)))

# ordger genes based on median log2 hazard ratio
genes_order_top_combined <- df_plot %>%
  group_by(Feature) %>% 
  dplyr::summarize(median_log2hr = median(log2hr, na.rm = T)) %>% 
  arrange(median_log2hr) %>% 
  select(Feature) %>% 
  tibble::deframe()

df_plot_genes <- df_plot %>%
  filter(Type == "CRISPRhit") %>% 
  mutate(Feature = factor(Feature, levels = genes_order_top_combined))


# plot gene expression panel
p_genes <- ggplot(df_plot_genes, aes(x = Name, y = Feature, size = log10p, color = log2hr)) +
  geom_point() +
  geom_point(data = df_plot_genes[df_plot_genes$Adj.P < 0.1,], pch = 21, aes(fill = log2hr), color = "black") +
  scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                        type = "div", limits = quantile(abs(df_plot$log2hr), 0.95) * c(-1, 1), oob = squish) +
  scale_fill_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                       type = "div", limits = quantile(abs(df_plot$log2hr), 0.95) * c(-1, 1), oob = squish) +
  scale_size(limits = c(0, max(df_plot$log10p))) +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "right", 
        legend.key.width = unit(0.3, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        plot.title = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(color = "black", fill = "white")) +
  labs(color = "Hazard ratio (log2)", size = "P (-log10)") +
  guides(fill = "none")


# plot subtype panel
df_plot_subtype <- df_plot %>%
  filter(Type == "Subtype") %>% 
  mutate(Feature = factor(Feature, levels = rev(c("WHSC1_FGFR3_Ig", "Hyperdiploid_gain1q", "MDS-like", "ABC", "GCB")),
                          labels = rev(c("WHSC1_FGFR3_Ig", "Hyperdiploid_gain1q", "MDS-like", "ABC", "GCB"))))

p_subtype <- ggplot(df_plot_subtype, aes(x = Name, y = Feature, size = log10p, color = log2hr)) +
  geom_point() +
  geom_point(data = df_plot_subtype[df_plot_subtype$Adj.P < 0.1,], pch = 21, aes(fill = log2hr), color = "black") +
  scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                        type = "div", limits = quantile(abs(df_plot$log2hr), 0.95) * c(-1, 1), oob = squish) +
  scale_fill_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                       type = "div", limits = quantile(abs(df_plot$log2hr), 0.95) * c(-1, 1), oob = squish) +
  scale_size(limits = c(0, max(df_plot$log10p))) +
  scale_x_discrete(drop = F) +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "right", 
        legend.key.width = unit(0.3, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        plot.title = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "plain"),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(color = "black", fill = "white")) +
  labs(color = "Hazard ratio (log2)", size = "P (-log10)") +
  guides(fill = "none")


# plot clinical risk classification panel
df_plot_clinical <- df_plot %>%
  filter(Type == "Clinical") %>% 
  mutate(Feature = factor(Feature, levels = rev(c("ISS3", "ISS1", "ELN2017_Adverse", "ELN2017_Favorable", "isTransformed", "isDenovo", "BM_Transplant", "IPI_4to5", "IPI_0to1", "Age", "Age-1")),
                          labels = rev(c("ISS3", "ISS1", "ELN2017 Adverse", "ELN2017 Favorable", "Socondary AML", "De novo AML", "BM transplant", "IPI 4-5", "IPI 0-1", "Age", "Age-1"))))

p_clinical <- ggplot(df_plot_clinical, aes(x = Name, y = Feature, size = log10p, color = log2hr)) +
  geom_point() +
  geom_point(data = df_plot_clinical[df_plot_clinical$Adj.P < 0.1,], pch = 21, aes(fill = log2hr), color = "black") +
  scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                        type = "div", limits = quantile(abs(df_plot$log2hr), 0.95) * c(-1, 1), oob = squish) +
  scale_fill_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                       type = "div", limits = quantile(abs(df_plot$log2hr), 0.95) * c(-1, 1), oob = squish) +
  scale_size(limits = c(0, max(df_plot$log10p))) +
  scale_x_discrete(drop = F) +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "right", 
        legend.key.width = unit(0.3, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        plot.title = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(face = "plain"),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(color = "black", fill = "white")) +
  labs(color = "Hazard ratio (log2)", size = "P (-log10)") +
  guides(fill = "none")


# get activating vs inhibitory classification for target genes from genome-wide screen data
df_crispr_target <- crispr %>%
  filter(gene %in% unique(df_plot_genes$Feature)) %>% 
  mutate(effect = ifelse(lfc > 0, "Enriched", "Depleted")) %>% 
  select(gene, effect, p) %>%
  mutate(effect = ifelse(p > 0.01, "ns", effect)) %>%
  group_by(gene) %>% 
  arrange(p) %>% 
  top_n(1, desc(p)) %>% 
  unique() %>% 
  filter(!is.na(effect)) %>% 
  mutate(CRISPR = "CRISPR effect") %>%
  mutate(gene = factor(gene, levels = genes_order_top_combined))


gw_colors <- c("Enriched" = brewer.pal(11, "RdBu")[3], "Depleted" = brewer.pal(11, "RdBu")[9], "ns" = "grey70")

# plot colored dots of CRISPR screen direction
p_gw_target <- ggplot(df_crispr_target, aes(x = CRISPR, y = gene, color = effect)) +
  geom_point() +
  scale_color_manual("", values = gw_colors, na.value = "white") +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black", face = "italic"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5))


layout <- "
ABBBB
#CCCC
#DDDD
"

p_gw_target + p_genes + p_subtype + p_clinical + plot_layout(guides = "collect", design = layout, heights = c(6,0.75,1.5))

ggsave("univariable_dotplot_cancertypes_genes_subtypes_clinical.pdf", height = 13, width = 8.5)

## -------------------------------------

# Save supplement table (Table S5D)

df_supplement <- df_plot %>% 
  select(-log10p, -coef, -log2hr) %>% 
  rename(Cohort = Name, Hazard_ratio = `exp(coef)`)

fwrite(df_supplement, "univariable_supplement.txt", sep = "\t")
    



