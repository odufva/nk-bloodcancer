
# Plot DEGs of all NK cell CROP-seq experiments

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(patchwork)


dir.create("results/combine")

# load perturbation results
nalm6_deg <- fread("results/nalm6/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4_deg <- fread("results/sudhl4/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
k562_deg <- fread("results/k562/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "K562")
mm1s_deg <- fread("results/mm1s/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S")
lp1_deg <- fread("results/lp1/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "LP1")
mm1s_a_deg <- fread("results/mm1sa/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S CRISPRa")

# load score results
nalm6_score <- fread("results/nalm6/deg/modulescore/revision/pert_scores.txt", data.table = F) %>% filter(gene == "nkresponse") %>% mutate(cell_line = "NALM6")
sudhl4_score <- fread("results/sudhl4/deg/modulescore/revision/pert_scores.txt", data.table = F) %>% filter(gene == "nkresponse") %>% mutate(cell_line = "SUDHL4")
k562_score <- fread("results/k562/deg/modulescore/revision/pert_scores.txt", data.table = F) %>% filter(gene == "nkresponse") %>% mutate(cell_line = "K562")
mm1s_score <- fread("results/mm1s/deg/modulescore/revision/pert_scores.txt", data.table = F) %>% filter(gene == "nkresponse") %>% mutate(cell_line = "MM1S")
lp1_score <- fread("results/lp1/deg/modulescore/revision/pert_scores.txt", data.table = F) %>% filter(gene == "nkresponse") %>% mutate(cell_line = "LP1")

# load nk vs no nk results
nalm6_1_16 <- fread("results/nalm6/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "NALM6", condition = "NK 1:16") %>% dplyr::rename(perturbation = comparison)
sudhl4_1_16 <- fread("results/sudhl4/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "SUDHL4", condition = "NK 1:16") %>% dplyr::rename(perturbation = comparison)
k562_1_16 <- fread("results/k562/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% mutate(cell_line = "K562", condition = "NK 1:16") %>% dplyr::rename(perturbation = comparison)
mm1s_1_16 <- fread("results/mm1s/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% mutate(cell_line = "MM1S", condition = "NK 1:16", comparison = "NK 1:16 vs no NK") %>% dplyr::rename(perturbation = comparison)
lp1_1_16 <- fread("results/lp1/deg/singlet/nk_nonk_deg.txt", data.table = F) %>% filter(comparison == "NK 1:16 vs no NK") %>% select(-comparison) %>% mutate(cell_line = "LP1", comparison = "NK 1:16 vs no NK", condition = "NK 1:16") %>% dplyr::rename(perturbation = comparison)

# load nk vs no nk score results
nalm6_1_16_score <- fread("results/nalm6/deg/modulescore/revision/condition_scores.txt", data.table = F) %>% filter(condition == "NK 1:16", gene == "nkresponse") %>% mutate(cell_line = "NALM6", perturbation = "NK 1:16 vs no NK")
sudhl4_1_16_score <- fread("results/sudhl4/deg/modulescore/revision/condition_scores.txt", data.table = F) %>% filter(condition == "NK 1:16", gene == "nkresponse") %>% mutate(cell_line = "SUDHL4", perturbation = "NK 1:16 vs no NK")
mm1s_1_16_score <- fread("results/mm1s/deg/modulescore/revision/condition_scores.txt", data.table = F) %>% filter(condition == "NK 1:16", gene == "nkresponse") %>% mutate(cell_line = "MM1S", perturbation = "NK 1:16 vs no NK",
                                                                                                                                                                  p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj),
                                                                                                                                                                  p_val = ifelse(p_val == 0, 1e-300, p_val))
lp1_1_16_score <- fread("results/lp1/deg/modulescore/revision/condition_scores.txt", data.table = F) %>% filter(condition == "NK 1:16", gene == "nkresponse") %>% mutate(cell_line = "LP1", perturbation = "NK 1:16 vs no NK",
                                                                                                                                                                p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj),
                                                                                                                                                                p_val = ifelse(p_val == 0, 1e-300, p_val))
k562_1_16_score <- fread("results/k562/deg/modulescore/revision/condition_scores.txt", data.table = F) %>% filter(condition == "NK 1:16", gene == "nkresponse") %>% mutate(cell_line = "K562", perturbation = "NK 1:16 vs no NK")

# load gsea results
gsea_hallmark <- fread("results/combine/gsea/hallmark.txt", data.table = F) %>%
  dplyr::rename(gene = pathway, p_val = pval, p_val_adj = padj, avg_log2FC = NES) %>% 
  mutate(pct.1 = NA, pct.2 = NA) %>% 
  select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, perturbation, cell_line, condition) %>% 
  mutate(condition = gsub("no_", "no ", gsub("NK_", "NK ", gsub("1_", "1:", condition))))

# load nk vs no nk gsea results
gsea_hallmark_nk <- fread("results/combine/gsea/hallmark_nkvsnonk.txt", data.table = F) %>%
  dplyr::rename(gene = pathway, p_val = pval, p_val_adj = padj, avg_log2FC = NES) %>% 
  mutate(pct.1 = NA, pct.2 = NA, perturbation = "NK 1:16 vs no NK", condition = "NK 1:16") %>% 
  select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, perturbation, cell_line, condition) %>% 
  mutate(condition = gsub("no_", "no ", gsub("NK_", "NK ", gsub("1_", "1:", condition))))


# load cell cycle phase results
nalm6_phase <- fread("results/nalm6/phase_analysis/phase_analysis.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4_phase <- fread("results/sudhl4/phase_analysis/phase_analysis.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
k562_phase <- fread("results/k562/phase_analysis/phase_analysis.txt", data.table = F) %>% mutate(cell_line = "K562")
mm1s_phase <- fread("results/mm1s/phase_analysis/phase_analysis.txt", data.table = F) %>% mutate(cell_line = "MM1S")
lp1_phase <- fread("results/lp1/phase_analysis/phase_analysis.txt", data.table = F) %>% mutate(cell_line = "LP1")

phase <- rbind(nalm6_phase, sudhl4_phase, k562_phase, mm1s_phase, lp1_phase)

# load differential abundance results
nalm6_counts <- fread("results/nalm6/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4_counts <- fread("results/sudhl4/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
k562_counts <- fread("results/k562/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "K562")
mm1s_counts <- fread("results/mm1s/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "MM1S")
lp1_counts <- fread("results/lp1/differential_abundance/singlet/differential_abundance.txt", data.table = F) %>% mutate(cell_line = "LP1")

counts <- rbind(nalm6_counts, sudhl4_counts, k562_counts, mm1s_counts, lp1_counts)

# get activating vs inhibitory classification from genome-wide screen data

# load screen results
crispr <- fread("crispr_mageck_combined.txt", data.table = F)

# combine data
data <- rbind(k562_deg, sudhl4_deg, nalm6_deg, mm1s_deg, lp1_deg, mm1s_a_deg,
              k562_score, sudhl4_score, nalm6_score, mm1s_score, lp1_score,
              nalm6_1_16, sudhl4_1_16, k562_1_16, mm1s_1_16, lp1_1_16,
              nalm6_1_16_score, sudhl4_1_16_score, k562_1_16_score, mm1s_1_16_score, lp1_1_16_score,
              gsea_hallmark, gsea_hallmark_nk)

data <- data %>% mutate(perturbation_cell_line = paste(perturbation, cell_line, sep = " "))


# numbers of DEG

dir.create("results/combine/deg_counts/")

deg <- rbind(k562_deg, sudhl4_deg, nalm6_deg, mm1s_deg, lp1_deg, mm1s_a_deg) %>% mutate(cell_line_condition = paste(cell_line, condition, sep = " "))

deg_counts <- deg %>%
  mutate(cell_line_condition = factor(cell_line_condition)) %>% 
  group_by(cell_line, condition, cell_line_condition, perturbation) %>% 
  summarize(count = length(p_val[p_val_adj < 0.05])) %>% 
  mutate(signif = ifelse(count > 4, "At least\n5 DEG", "Less than\n5 DEG"))

# Percent DEG change with vs without NK

percent_deg_change <- deg_counts %>% 
  filter(perturbation %in% unique(deg_counts$perturbation[deg_counts$signif=="At least\n5 DEG"])) %>% # only peturbations with a phenotype
  mutate(condition_binary = ifelse(condition %in% c("NK 1:16", "NK 1:4"), "With NK", "No NK")) %>% 
  group_by(perturbation, condition_binary) %>% 
  summarize(mean_count = mean(count)) %>% 
  tidyr::pivot_wider(values_from = mean_count, names_from = condition_binary) %>% 
  mutate(percent_change = 100*((`With NK` - `No NK`)/`With NK`)) %>% 
  arrange(desc(percent_change)) %>% 
  mutate(change = ifelse(percent_change > 50, "Yes", "No"))


pert_order <- deg_counts %>% 
  group_by(perturbation) %>%
  summarize(max_count = max(count)) %>% 
  arrange(desc(max_count)) %>% 
  select(perturbation) %>% 
  tibble::deframe() %>%
  as.character()

cell_lines <- c("K562", "SUDHL4", "MM1S", "LP1", "NALM6", "MM1S CRISPRa")
row_order <- c(paste("no NK", cell_lines),
               paste("NK 1:16", cell_lines),
               paste("NK 1:4", cell_lines[!grepl("K562", cell_lines)]))

conditions <- c("no NK", "NK 1:16", "NK 1:4")
row_order <- c(paste(cell_lines[1], c("no NK", "NK 1:16")),
               paste(cell_lines[2], conditions),
               paste(cell_lines[3], conditions),
               paste(cell_lines[4], conditions),
               paste(cell_lines[5], conditions),
               paste(cell_lines[6], conditions))

plotdata_wide <- dcast(deg_counts, cell_line_condition ~ perturbation, value.var = "count", drop = F)
rownames(plotdata_wide) <- plotdata_wide$cell_line_condition
plotdata_wide$cell_line_condition <- NULL

plotdata_wide <- plotdata_wide[row_order, pert_order]

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

change <- as.data.frame(percent_deg_change[match(colnames(plotdata_wide), percent_deg_change$perturbation), "change"])
colnames(change) <- "Majority of DEG NK-induced"

row_ha <- rowAnnotation(df = change,
                        col = list(`Majority of DEG NK-induced` = structure(c("black", "white"), names = c("Yes", "No"))),
                        na_col = "white",
                        show_legend = F)

# vertical heatmap (Figure S6A)
ht_genes <-  Heatmap(t(plotdata_wide),
                     cluster_columns = F,
                     cluster_rows = F,
                     left_annotation = row_ha,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if(!is.na(t(plotdata_wide)[i, j]))
                         grid.text( t(plotdata_wide)[i, j], x, y, gp = gpar(fontsize = 10))
                     },
                     na_col = "grey95",
                     col = colorRamp2(c(0, lseq(1, max(plotdata_wide, na.rm = T), length.out = 6)), c("white", brewer.pal(name = "Reds", n = 9)[c(2:7)])),
                     show_heatmap_legend = F,
                     border_gp = gpar(col = "black", lty = 1, lwd = 1),
                     row_names_gp = gpar(fontface = "italic"),
                     row_names_side = "left")

pdf("results/combine/deg_counts/deg_count_heatmap_vertical.pdf", height = 15, width = 7)
ht_genes
dev.off()


## --------------------------

## Dot plots

dir.create("results/combine/dotplots/revision")

# list genes/features to plot

nfkb <- c("BIRC3", "FAS", "NFKBIA", "LTB", "CD70", "CXCL10")

b_cell_development <- c("CD9", "DNTT", "MYB", "VPREB1", "VPREB3")

crispr_hits <- c("CD48", "CFLAR")

scores <- c("nkresponse")

core <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "TAP1", "TAPBP", "STAT1", "IRF1", "PSMB8", "PSMB9", "PSME1", "PSME2", "UBE2L6", "MT2A", "BST2", "GNLY")
hla2 <- c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1")

nk_induced_genes <- unique(c(core, hla2, nfkb, b_cell_development, crispr_hits, scores))

pathways = unique(gsea_hallmark$gene)


# -------------------------------------

# top condition for each perturbation, ordered by cell line (Figure 6E)

# list perturbations with at least 5 DEGs
top_cell_line_condition_perturbation <- deg_counts %>% 
  mutate(cell_line = gsub("\\ NK.*|\\ no NK.*", "", cell_line_condition)) %>% 
  group_by(perturbation, cell_line) %>% 
  slice(which.max(count)) %>% 
  filter(signif == "At least\n5 DEG") %>% 
  ungroup() %>% 
  mutate(cell_line_condition_perturbation = paste(cell_line_condition, perturbation)) %>% 
  select(cell_line_condition_perturbation) %>% 
  tibble::deframe()

pathways_selected <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")#, "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS")

plotdata <- data %>%
  mutate(cell_line_condition_perturbation = paste(cell_line, condition, perturbation)) %>%
  filter(gene %in% c(nk_induced_genes, pathways_selected) & (cell_line_condition_perturbation %in% top_cell_line_condition_perturbation | perturbation == "NK 1:16 vs no NK")) %>%
  filter(p_val < 0.05,
         perturbation_cell_line %in% c(paste0(c("KLF16", "RFXAP", "KIAA0922", "CMIP", "CHD7"), " NALM6"),
                                       paste0(c("IFNGR2", "JAK1", "JAK2", "STAT1", "FADD", "PCGF1", "RFXAP", "CASP8", "RUNX1", "YTHDF2", "METTL17", "CHD7", "STAG2", "PTEN", "BID"), " SUDHL4"),
                                       paste0(c("ARHGAP1", "NFKBIA", "NFKBIB", "JAK1", "JAK2", "PTEN", "GNA13", "PCGF5", "RFXAP", "NLRC5", "STAT1", "TRAF2"), " MM1S"),
                                       paste0(c("NLRC5", "SASH3"), " MM1S CRISPRa"),
                                       paste0(c("GFI1B", "JAK1", "JAK2", "PTPN2", "IFNGR2", "STAT1"), " K562"),
                                       paste0(c("IFNGR2", "JAK1", "JAK2","STAT1", "NLRC5", "RFXAP", "TRAF2", "CFLAR", "PTEN", "GSK3B", "MYB", "MSI2", "HMGB1", "CHD7", "RBBP4", "ARID1A"), " LP1"),
                                       paste0("NK 1:16 vs no NK ", c("K562", "SUDHL4", "NALM6", "MM1S", "LP1"))))


perturbation_cell_line_order <- c(paste0("NK 1:16 vs no NK ", c("K562", "SUDHL4", "MM1S", "LP1", "NALM6")),
                                  paste0(c("IFNGR2", "JAK1", "JAK2", "STAT1", "GFI1B"), " K562"),
                                  paste0(c("IFNGR2", "STAT1", "JAK1", "JAK2", "RFXAP", "STAG2", "PTEN", "METTL17", "FADD", "CASP8", "BID", "PCGF1", "YTHDF2", "RUNX1", "CHD7"), " SUDHL4"),
                                  paste0(c("JAK1", "JAK2", "STAT1", "NLRC5", "RFXAP", "TRAF2", "NFKBIA", "NFKBIB", "GNA13", "PTEN", "PCGF5", "ARHGAP1"), " MM1S"),
                                  paste0(c("NLRC5"), " MM1S CRISPRa"),
                                  paste0(c("IFNGR2", "JAK1", "JAK2","STAT1", "NLRC5", "RFXAP", "TRAF2", "CFLAR", "PTEN", "GSK3B", "MYB", "MSI2", "HMGB1", "CHD7", "RBBP4", "ARID1A"), " LP1"),
                                  paste0(c("RFXAP", "KLF16", "KIAA0922", "CMIP", "CHD7"), " NALM6")
)

nk_induced_genes_prettynames <- gsub("ifng", "IFNy response",
                                     gsub("nfkb", "NF-kB response",
                                          gsub("nkresponse", "NK cell response", nk_induced_genes)))

pathways_selected_prettynames <- gsub("Il6 jak stat3", "IL6-JAK-STAT3",
                                      gsub("Il2 stat5", "IL2-STAT5",
                                           gsub("Dna", "DNA",
                                                gsub("Myc", "MYC",
                                                     gsub("E2f", "E2F",
                                                          gsub("Kras", "KRAS",
                                                               gsub("Tnfa", "TNFA",
                                                                    gsub("G2m", "G2M",
                                                                         gsub("nfkb|Nfkb", "NF-kB",
                                                                              gsub("Ifng", "IFNy",
                                                                                   gsub("Nkresponse", "NK cell response",
                                                                                        stringr::str_to_sentence(gsub("_", " ",
                                                                                                                      sub(".*?_", "", pathways_selected))))))))))))))

df <- plotdata %>%
  mutate(gene_prettyname = gsub("Il6 jak stat3", "IL6-JAK-STAT3",
                                gsub("Il2 stat5", "IL2-STAT5",
                                     gsub("Dna", "DNA",
                                          gsub("Myc", "MYC",
                                               gsub("E2f", "E2F",
                                                    gsub("Kras", "KRAS",
                                                         gsub("Tnfa", "TNFA",
                                                              gsub("G2m", "G2M",
                                                                   gsub("Nf-kb|nfkb", "NF-kB",
                                                                        gsub("Ifng", "IFNy response",
                                                                             gsub("Nkresponse", "NK cell response",
                                                                                  stringr::str_to_sentence(gsub("_", " ",
                                                                                                                sub(".*?_", "", 
                                                                                                                    gsub("nfkb", "NF-kB response", gene)))))))))))))))) %>%
  mutate(gene_prettyname = ifelse(grepl("HALLMARK", gene) | grepl("response", gene_prettyname), gene_prettyname, gene)) %>% 
  ungroup() %>%
  filter(perturbation_cell_line %in% perturbation_cell_line_order) %>% 
  mutate(perturbation_cell_line = factor(perturbation_cell_line, levels = perturbation_cell_line_order)) %>%
  mutate(gene_prettyname = factor(gene_prettyname, levels = c(nk_induced_genes_prettynames, pathways_selected_prettynames))) %>% 
  mutate(cell_line = ifelse(perturbation == "NK 1:16 vs no NK", "NK vs no NK", cell_line)) %>% 
  mutate(cell_line = factor(cell_line, levels = c("NK vs no NK", "K562", "SUDHL4", "MM1S", "MM1S CRISPRa", "LP1", "NALM6"))) %>% 
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj)) %>% 
  mutate(p_val = ifelse(p_val == 0, 1e-300, p_val))

# get activating vs inhibitory classification for perturbations from genome-wide screen data

df_crispr <- df %>% left_join(crispr, by = c("perturbation" = "gene", "cell_line" = "cell_line")) %>% 
  mutate(effect = ifelse(lfc > 0, "Enriched", "Depleted")) %>% 
  select(perturbation, condition, cell_line, perturbation_cell_line, cell_line_condition_perturbation, effect) %>% 
  mutate(CRISPR = "CRISPR effect",
         ETratio = "          E:T ratio") %>% 
  mutate(cell_line = factor(cell_line, levels = c("NK vs no NK", "K562", "SUDHL4", "MM1S", "MM1S CRISPRa", "LP1", "NALM6")))


# get activating vs inhibitory classification for target genes from genome-wide screen data

df_crispr_target <- df %>% left_join(crispr, by = c("gene" = "gene", "cell_line" = "cell_line")) %>% 
  mutate(effect = ifelse(lfc > 0, "Enriched", "Depleted")) %>% 
  select(gene, effect, p) %>%
  mutate(effect = ifelse(p > 0.01, "ns", effect)) %>%
  group_by(gene) %>% 
  arrange(p) %>% 
  top_n(1, desc(p)) %>% 
  unique() %>% 
  filter(!is.na(effect)) %>% 
  mutate(CRISPR = "CRISPR effect") 

df_crispr_target <- rbind(df_crispr_target, data.frame(gene = c("NK cell response", pathways_selected_prettynames), effect = "ns", p = 0, CRISPR = "CRISPR effect")) %>% 
  mutate(gene = factor(gene, levels = rev(c(nk_induced_genes_prettynames, pathways_selected_prettynames))))



p_main <- ggplot(df, aes(x = perturbation_cell_line, y = gene_prettyname, size = -log10(p_val), fill = avg_log2FC)) +
  geom_point(pch = 21, color = "white") +
  scale_fill_distiller("Fold change (log2)", palette = "RdBu",
                       type = "div", limits = quantile(abs(df$avg_log2FC), 0.75) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top"),
                       oob = squish) +
  ggnewscale::new_scale_fill() +
  geom_point(data = df[df$gene %in% c("nkresponse"),], aes(fill = avg_log2FC), pch = 21, color = "white") +
  scale_fill_distiller("Fold change (log2)", palette = "BrBG", 
                       type = "div", limits = quantile(abs(df$avg_log2FC[df$gene %in% c("nkresponse")]), 0.85) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top"),
                       oob = squish) +
  geom_point(pch = 1, color = ifelse(df$p_val_adj < 0.05, "black", "white"), alpha = 0.5) +
  scale_size("FDR (-log10)", range = c(2.5, 4.5), guide = guide_legend(title.position = "top", override.aes = list(pch = 16, color = "black"))) +
  ggnewscale::new_scale_fill() +
  geom_point(data = df[grepl("HALLMARK", df$gene),], aes(fill = avg_log2FC), pch = 21, color = "white") +
  scale_fill_gradientn("NES", colors = pals::ocean.thermal(9)[2:8],
                       limits = max(abs(df$avg_log2FC[grepl("HALLMARK", df$gene)])) * c(-1, 1),
                       guide = guide_colorbar(title.position = "top"),
                       oob = squish) +
  geom_point(pch = 1, color = ifelse(df$p_val_adj < 0.05, "black", "white")) +
  scale_size("FDR (-log10)", range = c(2.5, 4.5), guide = guide_legend(title.position = "top", override.aes = list(pch = 16, color = "black"))) +
  scale_y_discrete(limits = rev, drop = F) +
  scale_x_discrete(position = "top", breaks = perturbation_cell_line_order, labels = c(gsub("NK 1:16 vs no NK ", "", perturbation_cell_line_order[1:5]), gsub(" K562| SUDHL4| NALM6| MM1S| LP1", "", perturbation_cell_line_order[-c(1:5)]))) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
    panel.border = element_rect(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom", 
    plot.title = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_text(color = "black"),
    plot.margin = unit(c(0,1,0,0), "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,-5,0,-5),
    legend.title.align = 0.5) +
  ylab("") +
  xlab("") +
  facet_grid(. ~ cell_line, scales = "free_x", space = "free", switch = "x")

gw_colors <- c("Enriched" = brewer.pal(11, "RdBu")[3], "Depleted" = brewer.pal(11, "RdBu")[9], "ns" = "grey70")

p_gw <- ggplot(df_crispr, aes(x = perturbation_cell_line, y = CRISPR, color = effect)) +
  geom_point() +
  scale_color_manual("", values = gw_colors, na.value = "white") +
  scale_x_discrete(position = "top", breaks = perturbation_cell_line_order, labels = c(gsub("NK 1:16 vs no NK ", "", perturbation_cell_line_order[1:5]), gsub(" K562| SUDHL4| NALM6| MM1S| LP1", "", perturbation_cell_line_order[-c(1:5)]))) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.placement = "outside") +
  facet_grid(. ~ cell_line, scales = "free_x", space = "free", switch = "x")


gw_colors_target <- c("Enriched" = brewer.pal(11, "RdBu")[3], "Depleted" = brewer.pal(11, "RdBu")[9], "ns" = "white")

p_gw_target <- ggplot(df_crispr_target, aes(x = CRISPR, y = gene, color = effect)) +
  geom_point() +
  scale_color_manual("", values = gw_colors_target, na.value = "white") +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black", face = "italic", size = 11),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5))

etratio_colors <- structure(pals::ocean.amp(9)[c(2,5,8)], 
                            names = c("no NK", "NK 1:16", "NK 1:4"))

p_etratio <- ggplot(df_crispr, aes(x = perturbation_cell_line, y = ETratio, color = condition)) +
  geom_point() +
  scale_color_manual("", values = etratio_colors, na.value = "white") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.placement = "outside") +
  facet_grid(. ~ cell_line, scales = "free_x", space = "free", switch = "x")

# plot Figure 6E
(plot_spacer() + p_gw + plot_layout(widths = c(2,100))) / (plot_spacer() + p_etratio + plot_layout(widths = c(2,100))) / (p_gw_target + p_main + plot_layout(widths = c(1,100))) + plot_layout(heights = c(1,1,100))

ggsave("results/combine/dotplots/pert_dotplot_celllineorder_facet.pdf", height = 9, width = 13)



