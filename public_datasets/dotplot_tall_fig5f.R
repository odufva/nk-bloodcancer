
# Compute pairwise correlations for NK PRISM sensitivity correlates in Liu et al T-ALL data and plot dot plot (Figure 5F)

# load libraries
library(data.table)
library(readxl)
library(dplyr)
library(GSVA)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(cowplot)


# load data
load("data/TALL_Liu_fm.Rdata")

# read correlation results
cor <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_T-ALL.tsv", data.table = F)

# GSVA gene set scores with PRISM correlating genes
genelist_sens = cor %>%
  filter(datapairs == "PRSM:GEXP" & featureA == "N:PRSM:AUC") %>%
  arrange(cor) %>%
  slice_min(cor, n = 50) %>%
  select(featureB) %>%
  tibble::deframe()

genelist_res = cor %>%
  filter(datapairs == "PRSM:GEXP" & featureA == "N:PRSM:AUC") %>%
  arrange(cor) %>%
  slice_max(cor, n = 50) %>%
  select(featureB) %>%
  tibble::deframe()

genesets = list(`N:SAMP:PRISM_NK_SENSITIVE_TALL` = genelist_sens,
                `N:SAMP:PRISM_NK_RESISTANT_TALL` = genelist_res)

fm_gsva <- fm[,!is.na(colSums(fm[grepl("^N:GEXP", rownames(fm)),]))] # remove samples with all gexp NA

# run GSVA
gsva_results <- gsva(as.matrix(fm_gsva), genesets)

# join GSVA results to fm
fm_genesets <- do.call(plyr::rbind.fill,list(fm, as.data.frame(gsva_results)))
rownames(fm_genesets) <- c(rownames(fm), as.character(rownames(gsva_results)))

genelist = c(genelist_sens, genelist_res, as.character(rownames(gsva_results)))
genelist = c(as.character(rownames(gsva_results)), "N:GEXP:FAS", "N:GEXP:PVR", "N:GEXP:ULBP1", "N:GEXP:CD44", "N:GEXP:DACH1", "N:GEXP:PSMG1") # add genes found in Fig5A analysis

# test significance of differential expression between subtypes
test_wilcox <- function(FEATURE1, FEATURE2){

  res <- wilcox.test(as.numeric(fm_genesets[FEATURE1, fm_genesets[FEATURE2,]=="1"]),
                     as.numeric(fm_genesets[FEATURE1, fm_genesets[FEATURE2,]=="0"]))

  exp <- median(as.numeric(fm_genesets[FEATURE1, fm_genesets[FEATURE2,]=="1"]))

  result <- data.frame(feature1 = FEATURE1, feature2 = FEATURE2, p = res$p.value, exp = exp)
  return(result)
}

test_wilcox_genes <- function(GENE, CLASS_LIST){
  result <- lapply(CLASS_LIST, test_wilcox, FEATURE1 = GENE) %>% bind_rows()
  result$FDR <- p.adjust(result$p, method = "BH")
  return(result)
}

# select genes based on analysis in Fig5A
genes <- c("N:GEXP:FAS", "N:GEXP:PVR", "N:GEXP:ULBP1", "N:SAMP:PRISM_NK_SENSITIVE_TALL")
test_classes <- rownames(fm)[grepl("B:CLIN:group", rownames(fm))]

df <- lapply(genes, test_wilcox_genes, CLASS_LIST = test_classes) %>% bind_rows()

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


df_scaled <- df %>%
  filter(!grepl("Unknown", feature2)) %>%
  mutate(feature1 = gsub("^.:....:", "", feature1),
         feature2 = gsub("^.:....:group_", "", feature2),
         log10_fdr = -log10(FDR)) %>%
  mutate(feature2 = gsub("NKX2_1", "NKX2.1", gsub("LMO2_LYL1", "LMO2/LYL1", feature2))) %>%
  mutate(feature1 = gsub("PRISM_NK_SENSITIVE_TALL", "NK sens", feature1)) %>%
  mutate(feature1 = factor(feature1, levels = rev(c("NK sens", "ULBP1", "PVR", "FAS"))),
         feature2 = factor(feature2, levels = c("TAL1", "TAL2", "LMO1/2", "NKX2.1", "TLX1", "TLX3", "HOXA", "LMO2/LYL1"))) %>%
  group_by(feature1) %>%
  mutate(exp_scaled = scale_this(exp))

ggplot(df_scaled, aes(y = feature2, x = feature1, color = exp_scaled, size = log10_fdr)) +
  geom_point() +
  scale_color_distiller(palette = "RdBu", values = seq(0, 1, length.out = 11),
                        type = "div", limits = max(abs(df_scaled$exp_scaled)) * c(-1, 1)) +
  scale_y_discrete(drop = F, position = "right") +
  guides(x = guide_axis(angle = 90)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(panel.border = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(face = "italic"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_blank(),
        plot.margin = unit(c(3,1,1,1), "cm"),
  ) +
  labs(color = "Median expression\n(Z-score)", size = "FDR (-log10)")

ggsave("TALL_Liu_subtype_dotplot_small.pdf", height = 4, width = 4)





