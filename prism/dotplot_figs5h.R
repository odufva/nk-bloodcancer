
# Dot plot of recurrent NK PRISM multiomic and GSEA correlations (Figure 5H)

library(data.table)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)

# load correlations
all <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_all.tsv", data.table = F) %>% mutate(cancer_type = "All")
mm <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_MM.tsv", data.table = F) %>% mutate(cancer_type = "MM")
bcl <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_BCL.tsv", data.table = F) %>% mutate(cancer_type = "BCL")
tall <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_T-ALL.tsv", data.table = F) %>% mutate(cancer_type = "T-ALL")
ball <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_B-ALL.tsv", data.table = F) %>% mutate(cancer_type = "B-ALL")
aml <- fread("results_20Q4_complete/NK_PRISM_heme_correlations_AML.tsv", data.table = F) %>% mutate(cancer_type = "AML")

data_main <- rbind(all, mm, bcl, tall, ball, aml)

# RPPA

data_rppa <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:RPPA") %>% 
  filter(p < 0.05) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("\\.", "", gsub("_", " ", gsub("_Caution", "", gsub("^N:....:", "", featureB))))) %>%
  group_by(featureB_plot) %>%
  filter(n()>1)

# miRNA

data_mirna <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:MIRN") %>% 
  filter(p < 0.05) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("^N:....:", "", featureB)) %>%
  group_by(featureB_plot) %>%
  filter(n()>2)

# LCMS

data_lcms <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:LCMS") %>% 
  filter(p < 0.05) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("hilic", "HILIC", gsub("tag", "TAG", gsub("lpc", "LPC", stringr::str_to_sentence(gsub("\\.|\\.\\.", " ", gsub("^N:....:", "", featureB))))))) %>%
  group_by(featureB_plot) %>%
  filter(n()>2)

# GDSC

data_gdsc <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:GDSC") %>% 
  filter(p < 0.05) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("^N:....:", "", featureB)) %>%
  group_by(featureB_plot) %>%
  filter(n()>1)


# GEXP

data_gexp <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:GEXP") %>% 
  filter(p < 0.05) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("^N:....:", "", featureB)) %>%
  group_by(featureB_plot) %>%
  filter(n()>2)


# GNAB

gnab_feats <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:GNAB") %>% 
  filter(p < 0.05) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("^B:....:", "", featureB)) %>%
  filter(cancer_type != "T-ALL") %>% # exclude T-ALL because of hypermutated cell lines -> large number of correlations
  group_by(featureB_plot) %>%
  filter(n()>1) %>% 
  select(featureB_plot) %>% 
  tibble::deframe()

data_gnab <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:GNAB") %>% 
  filter(p < 0.05) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("^B:....:", "", featureB)) %>%
  filter(featureB_plot %in% gnab_feats)


# CNVR

data_cnvr_neg <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:CNVR") %>% 
  filter(p < 0.05, cor < 0) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("^N:....:", "", featureB)) %>%
  group_by(featureB_plot) %>%
  filter(n()>2)

data_cnvr_pos <- data_main %>% filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:CNVR") %>% 
  filter(p < 0.05, cor > 0) %>%
  mutate(log10_p = -log10(p)) %>%
  mutate(featureB_plot = gsub("^N:....:", "", featureB)) %>%
  group_by(featureB_plot) %>%
  filter(n()>2)

# GSEA
gsea_all <- fread("results_gsea/NK_PRISM_all_HALLMARK_GSEA.txt", data.table = F) %>% mutate(cancer_type = "All")
gsea_mm <- fread("results_gsea/NK_PRISM_MM_HALLMARK_GSEA.txt", data.table = F) %>% mutate(cancer_type = "MM")
gsea_bcl <- fread("results_gsea/NK_PRISM_BCL_HALLMARK_GSEA.txt", data.table = F) %>% mutate(cancer_type = "BCL")
gsea_tall <- fread("results_gsea/NK_PRISM_TALL_HALLMARK_GSEA.txt", data.table = F) %>% mutate(cancer_type = "T-ALL")
gsea_ball <- fread("results_gsea/NK_PRISM_BALL_HALLMARK_GSEA.txt", data.table = F) %>% mutate(cancer_type = "B-ALL")
gsea_aml <- fread("results_gsea/NK_PRISM_AML_HALLMARK_GSEA.txt", data.table = F) %>% mutate(cancer_type = "AML")

gsea <- rbind(gsea_all, gsea_mm, gsea_bcl, gsea_tall, gsea_ball, gsea_aml)

gsea <- gsea %>%
  mutate(pathway_sentenceCase = gsub("Il6 jak stat3", "IL6-JAK-STAT3",
                                     gsub("Il2 stat5", "IL2-STAT5",
                                          gsub("Dna", "DNA",
                                               gsub("Myc", "MYC",
                                                    gsub("E2f", "E2F",
                                                         gsub("Kras", "KRAS",
                                                              gsub("Tnfa", "TNFA",
                                                                   gsub("G2m", "G2M",
                                                                        gsub("nfkb", "NF-kB",
                                                                             gsub("Mtorc1", "MTORC1",
                                                                                  gsub("Uv", "UV",
                                                                                       gsub("Pi3k akt mtor", "PI3K-AKT-MTOR",
                                                                                            gsub("Tgf", "TGF",
                                                                                                 stringr::str_to_sentence(gsub("_", " ",
                                                                                                                               sub(".*?_", "", pathway)))))))))))))))))




signif <- gsea %>% filter(padj < 0.05) %>% select(pathway) %>% unique()

signif_sentenceCase <- signif %>%
  mutate(pathway_sentenceCase = gsub("Il6 jak stat3", "IL6-JAK-STAT3",
                                     gsub("Il2 stat5", "IL2-STAT5",
                                          gsub("Dna", "DNA",
                                               gsub("Myc", "MYC",
                                                    gsub("E2f", "E2F",
                                                         gsub("Kras", "KRAS",
                                                              gsub("Tnfa", "TNFA",
                                                                   gsub("G2m", "G2M",
                                                                        gsub("nfkb", "NF-kB",
                                                                             gsub("Mtorc1", "MTORC1",
                                                                                  gsub("Uv", "UV",
                                                                                       gsub("Pi3k akt mtor", "PI3K-AKT-MTOR",
                                                                                            gsub("Tgf", "TGF",
                                                                                                 stringr::str_to_sentence(gsub("_", " ",
                                                                                                                               sub(".*?_", "", pathway)))))))))))))))))



data_gsea <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(log10_fdr = -log10(padj))

nes_order <- data_gsea %>% 
  group_by(pathway_sentenceCase) %>%
  filter(all(NES > 0) | all(NES < 0)) %>% 
  summarize(NES_mean = mean(NES), n_signif = n(), sum_NES = sum(NES)) %>% 
  filter(n_signif > 1) %>% 
  arrange(desc(sum_NES)) %>% 
  select(pathway_sentenceCase) %>% 
  tibble::deframe()

data_gsea <- data_gsea %>%
  filter(pathway_sentenceCase %in% nes_order) %>%
  mutate(featureA = NA, featureB = NA, signifCode = NA, test.method = NA, test.group = NA, datapairs = "PRSM:GSEA", log10_p = -log10(pval)) %>% 
  rename(cor = NES, p = pval, adj.p = padj, featureB_plot = pathway_sentenceCase) %>% 
  select(colnames(data_rppa)) %>% 
  mutate(featureB_plot = factor(featureB_plot, levels = nes_order))

# order features
feature_order <- c(levels(data_gsea$featureB_plot),
                   unique(data_gnab$featureB_plot),
                   unique(data_mirna$featureB_plot),
                   unique(data_rppa$featureB_plot), 
                   unique(data_lcms$featureB_plot))#,
feature_order <- feature_order[!feature_order %in% c("AMPK_pT172", "hsa-miR-548d-5p")] # filter out discordant hits

data_plot <- rbind(data_gsea, data_mirna, data_gnab, data_rppa, data_lcms) %>%
  filter(!featureB_plot %in% c("AMPK_pT172", "hsa-miR-548d-5p")) %>% # filter out discordant hits
  mutate(featureB_plot = factor(featureB_plot, levels = rev(feature_order)))

ggplot(data_plot, aes(x = cancer_type, y = featureB_plot, color = cor, size = log10_p)) +
  geom_point(aes(color = cor)) +
  scale_color_distiller("R", palette = "RdBu", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(data_plot$cor[data_plot$datapairs!="PRSM:GSEA"])) * c(-1, 1),
                        guide = guide_colorbar(title.position = "top")) +
  ggnewscale::new_scale_color() +
  geom_point(data = data_plot[data_plot$datapairs=="PRSM:GSEA",], aes(color = cor)) +
  scale_color_gradientn(name = "NES", colors = pals::ocean.thermal(9)[2:8],
                        guide = guide_colorbar(title.position = "top")) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.box = "vertical",
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  guides(size = guide_legend(title.position = "top")) +
  labs(size = "P (-log10)")

ggsave("CCLE_NK_PRISM_multiomic_dotplot.pdf", height = 11, width = 5)


# only GSEA

feature_order <- c(levels(data_gsea$featureB_plot))


data_plot <- rbind(data_gsea) %>%
  mutate(featureB_plot = factor(featureB_plot, levels = rev(feature_order)))

ggplot(data_plot, aes(x = cancer_type, y = featureB_plot, color = cor, size = log10_p)) +
  geom_point(aes(color = cor)) +
  scale_color_distiller("R", palette = "RdBu", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(data_plot$cor[data_plot$datapairs!="PRSM:GSEA"])) * c(-1, 1),
                        guide = guide_colorbar(title.position = "top")) +
  ggnewscale::new_scale_color() +
  geom_point(data = data_plot[data_plot$datapairs=="PRSM:GSEA",], aes(color = cor)) +
  scale_color_gradientn(name = "NES", colors = pals::ocean.thermal(9)[2:8],
                        guide = guide_colorbar(title.position = "top")) +
  scale_y_discrete(drop = F) +
  guides(x = guide_axis(angle = 45)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text = element_text(color = "black", size = 12),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.grid.major = element_line(size=0.25, colour="grey80", linetype = "dashed"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.box = "vertical",
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,-5,0,-5)) +
  guides(size = guide_legend(title.position = "top")) +
  labs(size = "P (-log10)")

ggsave("CCLE_NK_PRISM_GSEA_dotplot.pdf", height = 5, width = 5)

ggsave("CCLE_NK_PRISM_GSEA_dotplot_small.pdf", height = 5, width = 4.25)

