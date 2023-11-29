
# Plot stacked bar plot of NK CRISPR screen hit genetic alterations across cancer types (Figures 4I and S5A)

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
library(matrixStats)
library(ggplot2)
library(cowplot)

# load screen results
crispr <- fread("crispr_mageck_combined.txt", data.table = F)

genelist <- crispr %>%
  filter(p < 0.00005 & FDR < 0.2) %>%
  select(Gene) %>%
  unique() %>%
  filter(!grepl("hsa|Non|ORGENES", Gene)) %>% # remove miRNAs and control sgRNAs
  deframe()

# save gene list for cBioPortal
write.table(genelist, "NK_CRISPR_genelist.txt", quote = F, row.names = F, sep = "\t")

# use the above produced gene list to query mutations in cBioPortal

# load sample matrix from cBioPortal after query
mut <- fread("sample_matrix_20220416.txt", data.table = F)

mut <- mut %>%
  mutate(study = gsub("\\:.*", "", `studyID:sampleId`),
         sample = gsub(".*\\:", "", `studyID:sampleId`),
         cancer_type = toupper(gsub("laml", "aml", gsub("_.*", "", `studyID:sampleId`))))

samplecounts <- as.numeric(table(mut$cancer_type))

mut_pct <- mut %>%
  group_by(cancer_type) %>%
  summarize_at(vars(UBE2QL1:BCR), sum)

cancer_type <- mut_pct$cancer_type
mut_pct <- mut_pct[,-1]
mut_pct <- 100*mut_pct/samplecounts
mut_pct <- as.data.frame(t(mut_pct))
colnames(mut_pct) <- cancer_type
mut_pct$gene <- rownames(mut_pct)


geneorder <- mut_pct$gene[order(rowSums(mut_pct[,1:5], na.rm = T))]

mut_pct_long <- mut_pct %>%
  pivot_longer(ALL:MM, names_to = "cancer_type", values_to = "percentage") %>%
  mutate(gene = factor(gene, levels = geneorder))


colors <- fread("data/nk_crispr_colors.txt", data.table = F)
cols <- as.character(colors$color)
names(cols) <- gsub("B-ALL", "ALL", colors$cancer_type)
cols["CLL"] <- brewer.pal(11, "PuOr")[8]


topgenes <- mut_pct_long %>%
  group_by(gene) %>%
  summarize(percentage_sum = sum(percentage)) %>%
  filter(percentage_sum > 1.5) %>%
  select(gene) %>%
  tibble::deframe()
  
mut_pct_long %>% filter(gene %in% topgenes) %>% 

ggplot(aes(x = percentage, y = gene, fill = cancer_type)) +
  geom_col() +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0)) +
  ylab("") +
  xlab("% patients with mutation") +
  labs(fill = "Cancer type") +
  scale_fill_manual(values = cols[unique(mut_pct_long$cancer_type)]) +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("NK_CRISPR_patient_ns_mut_manuscript.pdf", height = 3.5, width = 5)


## --------------------------------------

# mutation types

# load sample matrix from cBioPortal

mut <- fread("mutations_20220416.txt", data.table = F)

mut <- mut %>%
  mutate(cancer_type = toupper(gsub("laml", "aml", gsub("_.*", "", `STUDY_ID`))))

samplecounts <- as.numeric(table(mut$cancer_type))

mut_long <- mut %>%
  pivot_longer(cols = UBE2QL1:BCR, names_to = "gene", values_to = "mutation") %>% 
  filter(mutation != "WT")

write.table(mut_long, "mutations_long_20220416.txt", quote = F, row.names = F, sep = "\t")


# mutation types barplot

mut_long_simple <- mut_long %>%
  mutate(mutation_simple = ifelse(grepl("fs", mutation), "Frameshift",
                                         ifelse(grepl("\\*$", mutation), "Stopgain",
                                                      ifelse(grepl("WT", mutation), "WT", "Missense_other"))))
mut_count <- mut_long_simple %>%
  group_by(gene, mutation_simple) %>%
  summarize(count = n())


mut_count %>% filter(gene %in% topgenes) %>% mutate(gene = factor(gene, levels = topgenes), mutation_simple = factor(mutation_simple, levels = c("Stopgain", "Frameshift", "Missense/other"))) %>% 
  
  ggplot(aes(x = count, y = gene, fill = mutation_simple)) +
  geom_col() +
  theme_cowplot() +
  scale_x_continuous(expand = c(0,0)) +
  ylab("") +
  xlab("Mutations") +
  labs(fill = "Mutation type") +
  scale_fill_manual(values = brewer.pal(11, "RdGy")[c(11, 3, 8)]) +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("NK_CRISPR_patient_ns_mut_types_apr2022_manuscript.pdf", height = 3.5, width = 5)





