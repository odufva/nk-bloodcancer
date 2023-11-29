
# Dot plot of NK PRISM sensitivity correlates in MMRF CoMMpass dataset

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

# load feature matrix
fm = get(load("data/MM_COMPASS_FM.Rdata"))

# load PRISM correlation results
cor = fread("results_20Q4_complete/NK_PRISM_heme_correlations_MM.tsv", data.table = F)

genes= cor %>%
  filter(datapairs == "PRSM:GEXP" & featureA == "N:PRSM:AUC") %>%
  arrange(p) %>%
  slice_min(p, n = 50) %>%
  mutate(gene = gsub("^.:....:", "", featureB)) %>%
  select(gene) %>%
  tibble::deframe()

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

genesets = list(`N:SAMP:PRISM_NK_SENSITIVE_MM` = genelist_sens,
                `N:SAMP:PRISM_NK_RESISTANT_MM` = genelist_res)

fm_gsva <- fm[,!is.na(colSums(fm[grepl("^N:GEXP", rownames(fm)),]))] # remove samples with all gexp NA

# run GSVA
gsva_results <- gsva(as.matrix(fm_gsva), genesets)

# join GSVA results to fm
fm_genesets <- do.call(plyr::rbind.fill,list(fm, as.data.frame(gsva_results)))
rownames(fm_genesets) <- c(rownames(fm), as.character(rownames(gsva_results)))

# test significance of differential expression between subtypes
test_wilcox <- function(FEATURE1, FEATURE2){
  
  res <- wilcox.test(as.numeric(fm_genesets[FEATURE1, fm_genesets[FEATURE2,]%in%"1"]),
                     as.numeric(fm_genesets[FEATURE1, fm_genesets[FEATURE2,]%in%"0"]))
  
  exp <- median(as.numeric(fm_genesets[FEATURE1, fm_genesets[FEATURE2,]%in%"1"]))
  
  result <- data.frame(feature1 = FEATURE1, feature2 = FEATURE2, p = res$p.value, exp = exp)
  return(result)
}

test_wilcox_genes <- function(GENE, CLASS_LIST){
  result <- lapply(CLASS_LIST, test_wilcox, FEATURE1 = GENE) %>% bind_rows()
  result$FDR <- p.adjust(result$p, method = "BH")
  return(result)
}

# eelct genes/signatures to plot
genes <- c("N:GEXP:CFLAR", "N:SAMP:HLAIScore", "N:SAMP:PRISM_NK_SENSITIVE_MM")
test_classes <- rownames(fm)[grepl("subtypes_", rownames(fm))]
test_classes <- test_classes[!grepl("and|vs", test_classes)]

df <- lapply(genes, test_wilcox_genes, CLASS_LIST = test_classes) %>% bind_rows()

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# prepare data frame for plotting
df_scaled <- df %>%
  mutate(feature1 = gsub("^.:....:", "", feature1),
         feature2 = gsub("^.:....:cancermap_subtypes_", "", feature2),
         log10_fdr = -log10(FDR)) %>%
  mutate(feature2 = gsub("WHSC1_FGFR3", "WHSC1/FGFR3", gsub("_Ig", "", gsub("_amp1q", " (amp 1q)", gsub("_Aberrated", "-aberrated", feature2))))) %>%
  mutate(feature1 = gsub("PRISM_NK_SENSITIVE_MM", "NK sensitivity signature", feature1)) %>%
  mutate(feature1 = factor(feature1, levels = rev(c("NK sensitivity signature", "CFLAR", "HLA I"))),
         feature2 = factor(feature2, levels = c("CCND1", "Hyperdiploid", "Hyperdiploid (amp 1q)", "MAF", "WHSC1/FGFR3", "TRAF3-aberrated"))) %>%
  group_by(feature1) %>%
  mutate(exp_scaled = scale_this(exp))

# dot plot
df_scaled <- df %>%
  mutate(feature1 = gsub("^.:....:", "", feature1),
         feature2 = gsub("^.:....:cancermap_subtypes_", "", feature2),
         log10_fdr = -log10(FDR)) %>%
  mutate(feature2 = gsub("WHSC1_FGFR3", "WHSC1", gsub("_Ig|_Aberrated", "", gsub("_amp1q", " (amp 1q)", gsub("Hyperdiploid", "HD", feature2))))) %>%
  mutate(feature1 = gsub("PRISM_NK_SENSITIVE_MM", "NK sens", gsub("HLAIScore", "HLA I", feature1))) %>%
  mutate(feature1 = factor(feature1, levels = rev(c("NK sens", "HLA I", "CFLAR"))),
         feature2 = factor(feature2, levels = c("CCND1", "HD", "HD (amp 1q)", "MAF", "WHSC1", "TRAF3"))) %>%
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

ggsave("CoMMpass_NK_PRISM_subtype_dotplot_small.pdf", height = 4, width = 4)



