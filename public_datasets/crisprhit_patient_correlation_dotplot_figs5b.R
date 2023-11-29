
## Plots of NK CRISPR screen hit expression vs genetic alterations across cancer types

# load libraries
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(ggsci)
library(parallel)
library(grid)
library(ggpubr)
library(data.table)
library(tibble)
library(parallel)
library(readxl)
library(RColorBrewer)
library(gridExtra)
library(patchwork)
library(scales)

# laod plotting scripts
source("../functions/plotting_functions_may2020.R")
source("../CCLE_featurematrix_NK_PRISM/compute.pairwise.R")
source("../CCLE_featurematrix_NK_PRISM/regulon.feats.R")
source("../functions/functions_statistics.R")

# ------------------------------------------------------------------

# load screen results
data <- fread("crispr_mageck_combined.txt", data.table = F)

genelist <- data %>%
  filter(!grepl("NonTargeting|hsa-mir", Gene)) %>%
  filter(p < 0.0001, abs(lfc) > 0.75) %>%
  select(Gene) %>%
  unique() %>%
  deframe()

# load correlation data and clean annotations
data_mm <- read.table("CoMMpass_NK_CRISPR_correlations.tsv", header = TRUE, sep = "\t")

df_mm <- data_mm %>%
  filter(featureA %in% paste0("N:GEXP:", genelist)) %>%
  mutate(log10.adj.p = -log10(adj.p),
         featureA_short = gsub("^.......", "", featureA),
         featureB_short = gsub("B:GNAB:", "B:",
                               gsub("B:SAMP:cancermap_subtypes_", "A:",
                                    gsub("N:CNVR:", "C:", featureB))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("finalFusion", featureB), "Structural variation (SV)",
                                         ifelse(grepl("CLIN", featureB), "Clinical",
                                                ifelse(grepl("B:GNAB", featureB), "Nonsynonymous mutation (MUT)", 
                                                       ifelse(grepl("subtype", featureB), "Subtype",
                                                              ifelse(grepl("CNVR", featureB), "CNV", ""))))))
  ) %>%
  mutate(featureAB_short = paste(featureA_short, featureB_short, sep = "_")) %>%
  mutate(data = "CoMMpass") %>%
  filter(grepl("@", featureB)) %>%
  filter(!grepl("Score|_and_|_vs_|cancermap_cluster|CGA|RNASeq|CNVR:[0-9]|MUT:FGFR3|MUT:CCND1|MUT:IG", featureB_short)) %>%
  mutate(featureB_notype = gsub("A:|B:|C:|D:", "", featureB_short)) %>%
  mutate(featureB_notype = factor(featureB_notype, levels = unique(featureB_notype)[order(unique(featureB_short))])) %>%
  mutate(featureA_short = factor(featureA_short, levels = genelist)) %>%
  filter(adj.p < 0.05)


data_dlbcl_reddy <- read.table("DLBCL_Reddy_NK_CRISPR_correlations.tsv", header = TRUE, sep = "\t")

df_dlbcl_reddy <- data_dlbcl_reddy %>%
  filter(featureA %in% paste0("N:GEXP:", genelist)) %>%
  mutate(log10.adj.p = -log10(adj.p),
         featureA_short = gsub("^.......", "", featureA),
         featureB_short = gsub("B:GNAB:", "B:",
                               gsub("B:CLIN:ABC_GCB_RNAseq_", "A:",
                                    gsub("N:CNVR:", "C:", featureB))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("finalFusion", featureB), "Structural variation (SV)",
                                         ifelse(grepl("CLIN", featureB), "Clinical",
                                                ifelse(grepl("B:GNAB", featureB), "Nonsynonymous mutation (MUT)", 
                                                       ifelse(grepl("subtype", featureB), "Subtype",
                                                              ifelse(grepl("CNVR", featureB), "CNV", ""))))))
  ) %>%
  mutate(featureAB_short = paste(featureA_short, featureB_short, sep = "_")) %>%
  mutate(data = "DLBCL Reddy") %>%
  filter(grepl("@", featureB)) %>%
  filter(!grepl("Score|_and_|_vs_|Cells", featureB_short)) %>%
  mutate(featureB_notype = gsub("A:|B:|C:|D:", "", featureB_short)) %>%
  mutate(featureB_notype = factor(featureB_notype, levels = unique(featureB_notype)[order(unique(featureB_short))])) %>%
  mutate(featureA_short = factor(featureA_short, levels = genelist)) %>%
  filter(adj.p < 0.05)


data_dlbcl_chapuy <- read.table("DLBCL_Chapuy_NK_CRISPR_correlations.tsv", header = TRUE, sep = "\t")

df_dlbcl_chapuy <- data_dlbcl_chapuy %>%
  filter(featureA %in% paste0("N:GEXP:", genelist)) %>%
  mutate(log10.adj.p = -log10(adj.p),
         featureA_short = gsub("^.......", "", featureA),
         featureB_short = gsub("B:GNAB:", "B:",
                               gsub("B:SAMP:COO_byGEP_", "A:",
                                    gsub(".:CNVR:", "C:", featureB))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("finalFusion", featureB), "Structural variation (SV)",
                                         ifelse(grepl("CLIN", featureB), "Clinical",
                                                ifelse(grepl("B:GNAB", featureB), "Nonsynonymous mutation (MUT)", 
                                                       ifelse(grepl("subtype", featureB), "Subtype",
                                                              ifelse(grepl("CNVR", featureB), "CNV", ""))))))
  ) %>%
  mutate(featureAB_short = paste(featureA_short, featureB_short, sep = "_")) %>%
  mutate(data = "DLBCL Chapuy") %>%
  filter(grepl("|@", featureB)) %>%
  filter(!grepl("Score|_and_|_vs_|Cells", featureB_short)) %>%
  mutate(featureB_notype = gsub("A:|B:|C:|D:", "", featureB_short)) %>%
  mutate(featureB_notype = factor(featureB_notype, levels = unique(featureB_notype)[order(unique(featureB_short))])) %>%
  mutate(featureA_short = factor(featureA_short, levels = genelist)) %>%
  filter(adj.p < 0.05)


data_dlbcl_tcga <- read.table("TCGA_DLBCL_NK_CRISPR_correlations.tsv", header = TRUE, sep = "\t")

df_dlbcl_tcga <- data_dlbcl_tcga %>%
  filter(grepl(paste(paste0("N:GEXP:", genelist), collapse = "|"), featureA)) %>%
  mutate(log10.adj.p = -log10(adj.p),
         featureA_short = gsub("^.......|:chr.*", "", featureA),
         featureB_short = gsub(":chr.*", "",
                               gsub("B:GNAB", "B:",
                                    gsub("B:CLIN:ABC_GCB_RNAseq_", "A:",
                                         gsub(".:CNVR:", "C:", 
                                              gsub(".*cg........_", "", featureB))))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("finalFusion", featureB), "Structural variation (SV)",
                                         ifelse(grepl("CLIN", featureB), "Clinical",
                                                ifelse(grepl("B:GNAB", featureB), "Nonsynonymous mutation (MUT)", 
                                                       ifelse(grepl("subtype", featureB), "Subtype",
                                                              ifelse(grepl("CNVR", featureB), "CNV", ""))))))
  ) %>%
  mutate(featureAB_short = paste(featureA_short, featureB_short, sep = "_")) %>%
  mutate(data = "DLBCL TCGA") %>%
  filter(test.group %in% c("regulatory.features")) %>%
  filter(!grepl("\\,|mean|1stExon|UTR|Body", featureB_short)) %>%
  mutate(featureB_notype = gsub("A:|B:|C:|D:", "", featureB_short)) %>%
  mutate(featureB_notype = factor(featureB_notype, levels = unique(featureB_notype)[order(unique(featureB_short))])) %>%
  mutate(featureA_short = factor(featureA_short, levels = genelist)) %>%
  filter(adj.p < 0.05) %>% 
  filter(cor < 0) # only negatice correlations relevant


data_aml_tcga <- read.table("TCGA_AML_NK_CRISPR_correlations.tsv", header = TRUE, sep = "\t")

df_aml_tcga <- data_aml_tcga %>%
  filter(grepl(paste(paste0("N:GEXP:", genelist), collapse = "|"), featureA)) %>%
  mutate(log10.adj.p = -log10(adj.p),
         featureA_short = gsub("^.......|:chr.*", "", featureA),
         featureB_short = gsub(":chr.*||\\|GENETICS.*|:::::", "",
                               gsub("B:GNAB:", "B:",
                                    gsub("B:SAMP:cancermap_", "A:",
                                         gsub("B:SAMP:I\\(", "C:",
                                              gsub(".*cg........_", "", featureB))))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("finalFusion", featureB), "Structural variation (SV)",
                                         ifelse(grepl("CLIN", featureB), "Clinical",
                                                ifelse(grepl("B:GNAB", featureB), "Nonsynonymous mutation (MUT)", 
                                                       ifelse(grepl("subtype", featureB), "Subtype",
                                                              ifelse(grepl("CNVR", featureB), "CNV", ""))))))
  ) %>%
  mutate(featureAB_short = paste(featureA_short, featureB_short, sep = "_")) %>%
  mutate(data = "AML TCGA") %>%
  filter(test.group %in% c("regulatory.features")) %>%
  filter(!grepl("\\,|mean|_and|_vs|1stExon|UTR|Body", featureB_short)) %>%
  mutate(featureB_notype = gsub("A:|B:|C:|D:", "", featureB_short)) %>%
  mutate(featureB_notype = factor(featureB_notype, levels = unique(featureB_notype)[order(unique(featureB_short))])) %>%
  mutate(featureA_short = factor(featureA_short, levels = genelist)) %>%
  filter(adj.p < 0.05) %>% 
  filter(cor < 0) # only negatice correlations relevant


data_aml_beat <- read.table("BeatAML_NK_CRISPR_correlations.tsv", header = TRUE, sep = "\t")

df_aml_beat <- data_aml_beat %>%
  filter(featureA %in% paste0("N:GEXP:", genelist)) %>%
  mutate(log10.adj.p = -log10(adj.p),
         featureA_short = gsub("^.......|:chr.*", "", featureA),
         featureB_short = gsub(":chr.*||\\|GENETICS.*", "",
                               gsub("B:GNAB:", "C:",
                                    gsub("B:SAMP:cancermap_cluster_TCGA_", "A:",
                                         gsub("B:CLIN:finalFusion_", "D:",
                                              gsub("B:CLIN:priorMDS_TRUE", "B:priorMDS", 
                                                   gsub("B:CLIN:FAB_Blast_Morphology_", "B:FAB_", featureB)))))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("finalFusion", featureB), "Structural variation (SV)",
                                         ifelse(grepl("CLIN", featureB), "Clinical",
                                                ifelse(grepl("B:GNAB", featureB), "Nonsynonymous mutation (MUT)", 
                                                       ifelse(grepl("subtype", featureB), "Subtype",
                                                              ifelse(grepl("CNVR", featureB), "CNV", ""))))))
  ) %>%
  mutate(featureAB_short = paste(featureA_short, featureB_short, sep = "_")) %>%
  mutate(data = "BeatAML") %>%
  filter(test.group %in% c("regulatory.features")) %>%
  filter(!grepl("\\,|mean|and|vs", featureB_short)) %>%
  mutate(featureB_notype = gsub("A:|B:|C:|D:", "", featureB_short)) %>%
  mutate(featureB_notype = factor(featureB_notype, levels = unique(featureB_notype)[order(unique(featureB_short))])) %>%
  mutate(featureA_short = factor(featureA_short, levels = genelist)) %>%
  filter(adj.p < 0.05)


# ------------------------------------------------------------------


# faceted plot
df_plot <- rbind(df_aml_tcga, df_dlbcl_tcga, df_dlbcl_reddy, df_mm) %>% 
  mutate(featureB_notype = gsub("_", " ", gsub("\\@.*", "", featureB_notype))) %>% 
  mutate(featureB_notype = factor(featureB_notype, levels = c("TSS200 NShore", "TSS200 Island", "TSS200 SShore", "TSS1500 NShelf", "TSS1500 NShore", "TSS1500 Island", "TSS1500 SShore", "TRAF2", "ARID1A", "PTEN", "CD58")),
         data = factor(data, levels = c("AML TCGA", "DLBCL TCGA", "DLBCL Reddy", "CoMMpass"), labels = c("AML (TCGA)", "DLBCL (TCGA)", "DLBCL (Reddy et al.)", "MM (CoMMpass)")))

  ggplot(df_plot, aes(x = featureB_notype, y = featureA_short, fill = cor, size = log10.adj.p)) +
  geom_point(pch = 21) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), values = seq(0, 1, length.out = 11), limits = c(-0.75, 0.75), oob = squish) +
  scale_size(limits = c(0.1,60)) +
  scale_y_discrete(drop = T) +
  guides(x = guide_axis(angle = 45),
         size = guide_legend(override.aes = list(pch = 16))) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major=element_line(size=0.25, colour="grey80", linetype = "dashed"),
        plot.title = element_text(face = "plain", hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text = element_text(color = "black", size = 12)) +
  labs(fill = "R (Spearman)", size = "FDR (-log10)") +
  facet_grid(. ~ data, space = "free", scales = "free_x")

ggsave("NK_CRISPR_pancancer_correlations_meth_cnvr.pdf", height = 8, width = 6.5)

# write table for supplement (Table S5B)
df_plot_supplement <- df_plot %>% 
  select(-featureB_short, -featureAB_short, -test.group, -log10.adj.p) %>% 
  rename(featureB_short = featureB_notype) %>% 
  relocate(featureA_short, .before = featureB_short)

fwrite(df_plot_supplement, "NK_CRISPR_pancancer_correlations_meth_cnvr.txt", row.names = F, quote = F, sep = "\t")

