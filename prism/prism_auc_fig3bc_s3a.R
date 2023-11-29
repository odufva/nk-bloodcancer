
# Calculate normalized AUC for NK PRISM data
# Percent viability

# load libraries
library(drc)
library(DescTools)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(data.table)
library(dplyr)
library(cowplot)
library(readxl)
library(tidyr)

# load processed data
nk_prism_reps <- fread("data/NK_PRISM_heme_pctviability_biologicalreplicates.txt")
nk_prism <- fread("data/NK_PRISM_heme_pctviability.txt") %>%
  filter(!is.na(pctviability), !cell_line %in% "HUT78") %>% # remove NA and non-heme cell line
  mutate(Condition = as.numeric(gsub("CTRL", "0", gsub("D1_NK_|D1_", "", Condition))))

## -----------------------

# function for AUC
calculate_auc <- function(CELLLINE){
  
  d <- data.frame(nk_prism[nk_prism$cell_line == CELLLINE,])
  auc = DescTools::AUC(d[,"Condition"],d[,"pctviability"])
  return(data.frame(cell_line = CELLLINE, auc = auc))
}

cell_lines <- unique(nk_prism$cell_line)

# apply AUC function across cell lines
auc <- do.call(rbind, lapply(cell_lines, calculate_auc))

auc <- auc %>% arrange(auc)

pctviability_auc <- merge(nk_prism, auc)

# normalize to 1 = 100 % viability in all E:T ratios
pctviability_auc <- pctviability_auc %>%
  mutate(auc_norm = auc/500)

write.table(pctviability_auc, "NK_PRISM_heme_pctviability_auc_norm.txt", quote = F, row.names = F, sep = "\t")

pctviability_auc %>% filter(Condition == 0) %>% 
  ggplot(aes(x = reorder(cell_line, auc_norm), y = auc_norm)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -----------------------

# plot line plots arranged by AUC

nk_prism_reps <- nk_prism_reps %>%
  pivot_longer(!c(Condition, Biological_replicate), names_to = "cell_line", values_to = "pctviability") %>%
  filter(grepl("HAEMATOPOIETIC", cell_line)) %>%
  mutate(cell_line = gsub("_.*", "", cell_line)) %>%
  filter(!is.na(pctviability), !cell_line %in% "HUT78") %>% # remove NA
  mutate(Condition = as.character(gsub("CTRL", "0", gsub("D1_NK_|D1_", "", Condition))))


data_plot_subtype <- pctviability_auc %>%
  mutate(Subtype_simple = ifelse(grepl("AML", Subtype), "AML",
                                 ifelse(grepl("ALL), B-cell", Subtype), "B-ALL",
                                        ifelse(grepl("ALL), T-cell", Subtype), "T-ALL",
                                               ifelse(grepl("Multiple", Subtype), "MM",
                                                      ifelse(grepl("CLL", Subtype), "CLL",
                                                             ifelse(grepl("Mantle", Subtype), "MCL",
                                                                    ifelse(grepl("B-cell, Hodgkins", Subtype), "CHL",
                                                                           ifelse(grepl("Burkitt", Subtype), "BL",
                                                                                  ifelse(grepl("ALCL", Subtype), "ALCL",
                                                                                         ifelse(grepl("DLBCL", Subtype), "DLBCL", 
                                                                                                ifelse(grepl("B-cell, Non", Subtype), "BCL other",
                                                                                                       ifelse(grepl("Cutaneous", Subtype), "CTCL",
                                                                                                              ifelse(grepl("^T-cell", Subtype), "TCL other", Subtype)))))))))))))) %>% 
  mutate(Subtype_main = ifelse(Subtype_simple %in% c("DLBCL", "BCL other", "CHL", "BL", "MCL", "CLL"), "BCL",
                               ifelse(Subtype_simple %in% c("CTCL", "ALCL", "ATL", "TCL other"), "TCL", Subtype_simple)))


# faceted line plot ordered by cancer type (Figure S3A)

data_plot_subtype <- data_plot_subtype %>% mutate(Condition = as.character(Condition))

complete <- nk_prism %>% filter(Condition == 5) %>% dplyr::select(cell_line) %>% tibble::deframe()

cell_lines <- as.character(auc$cell_line)

nk_prism_reps_plot <- nk_prism_reps %>%
  left_join(data_plot_subtype, by = c("cell_line", "Condition"), suffix = c("_rep", "_mean")) %>% 
  mutate(cell_line = factor(cell_line, levels = cell_lines)) %>% 
  filter(cell_line %in% complete) %>% # include only cell lines with complete data from all condition
  mutate(Subtype_main = factor(Subtype_main, levels = c("AML", "TCL", "T-ALL", "BCL", "MM", "B-ALL")))

cell_lines <- nk_prism_reps_plot %>% arrange(Subtype_main, auc) %>% dplyr::select(cell_line) %>% unique() %>% tibble::deframe()

nk_prism_reps_plot <- nk_prism_reps_plot %>% mutate(cell_line = factor(cell_line, levels = cell_lines))

cols <- fread("../CCLE_featurematrix_NK_PRISM/nk_crispr_colors.txt", data.table = F)
cols$cancer_main <- gsub("ALCL", "TCL", gsub("DLBCL", "BCL", as.character(cols$cancer)))

cols_vector <- cols$color
names(cols_vector) <- cols$cancer_main

subtype_simple_df <- nk_prism_reps_plot %>% dplyr::select(cell_line, Subtype_simple, auc_norm) %>% unique()
subtype_simple_df <- subtype_simple_df[match(cell_lines, subtype_simple_df$cell_line),]

label_df_type <- data.frame(x = 0.6, y = 300,
                            cell_line = cell_lines,
                            subtype_simple = subtype_simple_df$Subtype_simple)

label_df_auc <- data.frame(x = 0.6, y = 250,
                           cell_line = cell_lines,
                           auc_norm = paste("AUC =", signif(subtype_simple_df$auc_norm, 2)))

ggline(nk_prism_reps_plot, x = "Condition", y = "pctviability_rep", add = c("mean_sd", "jitter"),
       point.size = 0.25, add.params = list(size = 0.75, width = 0.4), facet.by = "cell_line", ncol = 9, color = "Subtype_main") +
  ylab("% viability\n(normalized to without NK cells)") +
  xlab("Effector:target ratio") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 310), breaks = c(0, 50, 100, 150, 200)) +
  scale_color_manual(values = cols_vector[c("AML", "TCL", "T-ALL", "BCL", "MM", "B-ALL")]) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey70") +
  guides(color = guide_legend(nrow = 1)) +
  geom_text(aes(x, y, label = subtype_simple), data = label_df_type, color = "grey50", hjust = 0) +
  geom_text(aes(x, y, label = auc_norm), data = label_df_auc, color = "grey50", hjust = 0)

ggsave("NK_PRISM_heme_pctviability_lineplots_facet_bycancertype.pdf", height = 12, width = 14)


## ------------------------------------------

# plot example line plots (Figure 3B)
colors <- fread("../CCLE_featurematrix_NK_PRISM/nk_crispr_colors.txt", data.table = F)
cols_example <- colors$color[colors$cancer_type %in% c("B-ALL", "AML")]
names(cols_example) <- c("THP1 (AML)", "697 (B-ALL)")

nk_prism_reps_plot %>% 
  filter(cell_line %in% c("THP1", "697")) %>%
  mutate(cell_line = gsub("THP1", "THP1 (AML)", gsub("697", "697 (B-ALL)", cell_line))) %>% 
  
  ggline(x = "Condition", y = "pctviability_rep", add = c("mean_sd", "jitter"),
         point.size = 0.25, point.color = "black", add.params = list(size = 0.75, width = 0.2), color = "cell_line") +
  ylab("% viability") +
  xlab("Effector:target ratio") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.25, "cm"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        legend.text = element_text(color = "black")) +#,
  #axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 150), breaks = c(0, 50, 100, 150)) +
  scale_color_manual(values = cols_example) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey70") +
  guides(color = guide_legend(nrow = 1))
ggsave("NK_PRISM_heme_pctviability_lineplots_example.pdf", height = 3.5, width = 3.5)



# boxplot of main cancer types (Figure 3C)

data_plot_subtype$Subtype_main <- data_plot_subtype$Subtype_simple
data_plot_subtype$Subtype_main[data_plot_subtype$Subtype_simple %in% c("DLBCL", "BCL other", "CHL", "BL", "MCL", "CLL")] <- "BCL"
data_plot_subtype$Subtype_main[data_plot_subtype$Subtype_simple %in% c("CTCL", "ALCL", "ATL", "TCL other")] <- "TCL"

data_plot_subtype_order <- data_plot_subtype %>%
  filter(Condition == "5") %>%
  group_by(Subtype_main) %>%
  summarize(auc_median = median(auc_norm)) %>%
  arrange(desc(auc_median))

names(cols) <- gsub("ALCL", "TCL", gsub("DLBCL", "BCL", as.character(colors$cancer_type)))


data_plot_final <- data_plot_subtype %>%
  filter(Condition == "5") %>%
  mutate(Subtype_main = factor(Subtype_main, levels = rev(data_plot_subtype_order$Subtype_main)))

comparisons = list(c("B-ALL", "MM"), c("B-ALL", "BCL"), c("B-ALL", "T-ALL"), c("B-ALL", "TCL"), c("B-ALL", "AML"),
                   c("MM", "BCL"), c("MM", "T-ALL"), c("MM", "TCL"), c("MM", "AML"),
                   c("BCL", "T-ALL"), c("BCL", "TCL"), c("BCL", "AML"),
                   c("T-ALL", "TCL"), c("T-ALL", "AML"),
                   c("TCL", "AML"))

test_wilcox_pairwise <- ggpubr::compare_means(auc_norm ~ Subtype_main, comparisons = comparisons, p.adjust.method = "BH", method = 'wilcox.test', data = data_plot_final)
test_wilcox_pairwise <- test_wilcox_pairwise %>% mutate(y.position = c(NA, NA, 1.4, 1.55, 1.7, NA, NA,  NA, NA, NA,NA, NA, NA, NA, NA),
                                                        p.signif = ifelse(p.adj < 0.001, "***",
                                                                          ifelse(p.adj < 0.01, "**", 
                                                                                 ifelse(p.adj < 0.05, "*", "ns"))))

ggplot(data_plot_final, aes(x = Subtype_main, y = auc_norm)) +
  geom_boxplot(aes(fill = Subtype_main), outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cols) +
  guides(fill = F) +
  ylab("PRISM AUC") +
  xlab("") +
  ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.y = 1, label.x = "AML") +
  ggpubr::stat_pvalue_manual(test_wilcox_pairwise, label = "p.signif", size = 6)


ggsave("NK_PRISM_heme_boxplot_subtype_main_AUC_normalized_pctviability.pdf", height = 3.5, width = 3.5)  



## ----------------------------

# supplement tables

# S5A (annotations, barcodes, AUC)

barcodes <- readxl::read_excel("PRISM_barcodes.xlsx")
barcodes <- barcodes %>% 
  mutate(cell_line = toupper(gsub("\\-|\\ |\\.|\\/", "", Name))) %>% 
  dplyr::select(cell_line, Barcode = Sequence)

data_plot_subtypes_supplement <- data_plot_subtypes %>%
  filter(Condition == 0) %>% 
  left_join(barcodes, by = "cell_line") %>% 
  dplyr::select(Cell_line = cell_line, AUC_normalized = auc_norm, DepMap_ID,
                CCLE_Name, COSMICID, RRID, Barcode, Cancer_type = Subtype, Cancer_type_simple = Subtype_simple, Cancer_type_main = Subtype_main,
                Subtype = subtype, Subtype_2 = subtype2, Subtype_source = source.y)

fwrite(data_plot_subtypes_supplement, "prism_norm_auc_annotations_supplement.txt", quote = F, row.names = F, sep = "\t")

# S5B (percent viability)
nk_prism_mean <- nk_prism[,c("cell_line", "Condition", "pctviability")]
class(nk_prism_mean$Condition) <- "character"

nk_prism_reps_supplement <- nk_prism_reps %>% 
  left_join(nk_prism_mean, by = c("cell_line", "Condition")) %>% 
  dplyr::select(Cell_line = cell_line, Effector_target_ratio = Condition, Replicate = Biological_replicate,
                Percent_viability_replicate = pctviability.x,
                Percent_viability_mean = pctviability.y) %>% 
  arrange(Cell_line, Effector_target_ratio, Replicate)

fwrite(nk_prism_reps_supplement, "prism_pctviability_supplement.txt", quote = F, row.names = F, sep = "\t")

