
# Differential expression of CRISPR screen hits between heme and solid tumor cell lines

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
library(tidyr)
library(tibble)
library(ggpointdensity)
library(RColorBrewer)
library(ggnewscale)
library(wesanderson)

# load screen results
data <- fread("data/crispr_mageck_combined.txt", data.table = F)

# load feature matrix
fm <- readRDS("data/CCLE.fm_20Q4_complete.rds")

genelist <- data %>%
  filter(p < 0.001, abs(lfc)>0.75) %>%
  filter(Gene %in% gsub("N:GEXP:", "", rownames(fm))) %>%
  select(Gene) %>%
  unique() %>%
  deframe()


ccle_annot <- fread("../NK_PRISM/sample_info.csv")

extrafeatures=NULL
l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 3, cores = 2, adjust.method = "BH")

dir.create("results_CCLE_all")
fwrite(results, "results_CCLE_all/NK_CRISPR_CCLE_correlations.tsv", sep ="\t")


mat_gexp <- as.data.frame(t(fm[paste0("N:GEXP:", genelist),]))
mat_gexp$cell_line <- rownames(mat_gexp)
mat_gexp <- mat_gexp %>%
  pivot_longer(!cell_line, names_to = "gene", values_to = "gexp") %>%
  mutate(gene = gsub("N:GEXP:", "", gene)) %>%
  filter(!is.na(gexp)) %>%
  filter(gene %in% genelist)

mat_gexp <- merge(mat_gexp, ccle_annot, by.x = "cell_line", by.y = "DepMap_ID")

# Test differential expression heme vs solid
mat_gexp <- mat_gexp %>%
  mutate(heme_solid = ifelse(primary_disease %in% c("Leukemia", "Lymphoma", "Myeloma"), "Hematological", "Solid"))
  

de_wilcox <- function(GENE){
  
  heme <- as.numeric(mat_gexp[mat_gexp$gene %in% GENE & mat_gexp$heme_solid %in% "Hematological", "gexp"])
  solid <- as.numeric(mat_gexp[mat_gexp$gene %in% GENE & mat_gexp$heme_solid %in% "Solid", "gexp"])
  
  res <- wilcox.test(heme, solid)
  
  median_hematological = median(heme, na.rm = T) + 0.001
  median_solid = median(solid, na.rm = T) + 0.001
  
  result <- data.frame(gene = GENE,
                       median_hematological = median_hematological,
                       median_solid = median_solid,
                       log2fc = log2(median_hematological/median_solid),
                       p = res$p.value)
  return(result)
}


result <- lapply(genelist, de_wilcox) %>% bind_rows()

result$FDR <- p.adjust(result$p, method = "bonferroni")

result <- result %>%
  arrange(p) %>%
  mutate(log10p = -log10(p),
         log10fdr = -log10(FDR))

# Save supplementary table (Table S5C)
result_supplement <- result %>% select(-log10p, -log10fdr)
fwrite(result_supplement, "ccle_heme_solid_supplement.txt", sep = "\t")

# Volcano plot
cutoff <- result$FDR < 0.05
cutoff2 <- result$FDR < 10^-20 & abs(result$log2fc) > 2
selected <- F

ggplot(result, aes(x = log2fc, y = log10fdr)) +
  geom_point(color = ifelse(cutoff2, "grey50", "grey70")) +
  geom_point(data = result[selected|cutoff2,], size = 4, pch = 21, aes(fill = log2fc)) +
  scale_fill_gradientn(colors = c(rev(brewer.pal(9, "YlGnBu")), brewer.pal(9, "YlOrRd"))) +
  geom_text_repel(aes(label = ifelse(selected|cutoff2, as.character(gene), "")), point.padding = 6) +
  ylab("Adjusted p value (-log10)") +
  xlab("Fold change (log2)") +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.justification = "center", legend.key.height = unit(0.2, "cm")) +
  scale_x_continuous(limits = c(-8,8)) +
  guides(fill = guide_colorbar(title = ""))

ggsave("NK_CRISPR_CCLE_heme_solid_volcanoplot.pdf", height = 4, width = 6)

