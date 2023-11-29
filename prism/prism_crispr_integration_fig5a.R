  
# Plot cancer type-specific CCLE NK PRISM gene expression correlations (Figure 5A)

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(viridis)

# load cancer type-specific correlation results
aml <- fread("results_20Q4_complete_pctviability//NK_PRISM_heme_correlations_AML.tsv", data.table = F) %>% mutate(cancer_type = "AML")
bcl <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_BCL.tsv", data.table = F) %>% mutate(cancer_type = "BCL")
mm <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_MM.tsv", data.table = F) %>% mutate(cancer_type = "MM")
ball <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_B-ALL.tsv", data.table = F) %>% mutate(cancer_type = "B-ALL")
tall <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_T-ALL.tsv", data.table = F) %>% mutate(cancer_type = "T-ALL")
all <- fread("results_20Q4_complete_pctviability/NK_PRISM_heme_correlations_all.tsv", data.table = F) %>% mutate(cancer_type = "All")

data <- rbind(aml, bcl, mm, ball, tall, all)

data_gexp <- data %>%
  filter(featureA == "N:PRSM:AUC", datapairs == "PRSM:GEXP") %>%
  mutate(signed_p = sign(cor) * -log10(p)) %>%
  mutate(featureB = gsub("N:GEXP:", "", featureB))

# CRISPR hits 

crispr <- fread("crispr_mageck_combined.txt", data.table = F)

p_threshold <- 0.0001
lfc_threshold_pos <- 0.75
lfc_threshold_neg <- -0.75

activating_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc>lfc_threshold_pos) %>%
  select(Gene) %>%
  unique() %>%
  tibble::deframe()

activating_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc<lfc_threshold_neg) %>%
  select(Gene) %>%
  unique() %>%
  tibble::deframe()

inhibitory_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc<lfc_threshold_neg) %>%
  select(Gene) %>%
  unique() %>%
  tibble::deframe()

inhibitory_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc>lfc_threshold_pos) %>%
  select(Gene) %>%
  unique() %>%
  tibble::deframe()

activating <- c(activating_lof, activating_gof)
inhibitory <- c(inhibitory_lof, inhibitory_gof)

# plot (Figure 5A)

data_gexp_act <- data_gexp %>%
  filter(featureB %in% activating, p < 0.05, cor < 0)

data_gexp_inh <- data_gexp %>%
  filter(featureB %in% inhibitory, p < 0.05, cor > 0)

data_gexp_act_inh <- rbind(data_gexp_act, data_gexp_inh)

data_gexp$cancer_type <- factor(data_gexp$cancer_type, levels = rev(c("All", "AML", "T-ALL", "BCL", "MM", "B-ALL")))

set.seed(40)

p <- ggplot(data_gexp[sample(nrow(data_gexp), 50000),], aes(x = signed_p, y = cancer_type)) +
  ggrastr::geom_point_rast(position = position_jitter(seed = 42, height = 0.4), size = 0.1, color = "grey70") +
  geom_point(position = position_jitter(seed = 40, height = 0.4), data = data_gexp_act, pch = 21, size = 4, color = "black", aes(fill = signed_p)) +#fill = brewer.pal(11, "RdBu")[3]) +
  scale_fill_viridis(option = "magma", begin = 0.4, end = 1, direction = 1) +
  ggnewscale::new_scale_fill() +
  geom_point(position = position_jitter(seed = 40, height = 0.4), data = data_gexp_inh, pch = 21, size = 4, color = "black", aes(fill = signed_p)) +# fill = brewer.pal(11, "RdBu")[9]) +
  scale_fill_viridis(option = "viridis", begin = 0.4, end = 1, direction = -1) +
  geom_text_repel(position = position_jitter(seed = 40, height = 0.5), data = data_gexp_act,
                  aes(x = signed_p, y = cancer_type, label = featureB),
                  max.overlaps = 30,
                  box.padding = 0.5,
                  point.padding = 0.1,
                  fontface = "italic",
                  size = 5) +
  geom_text_repel(position = position_jitter(seed = 40, height = 0.5), data = data_gexp_inh,
                  aes(x = signed_p, y = cancer_type, label = featureB),
                  max.overlaps = 30,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  fontface = "italic",
                  size = 5) +
  theme_bw() +
  theme(axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom", 
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 16),
        plot.margin = unit(c(0,1,0,0), "cm")) +
  xlab("Signed P value (-log10)") +
  ylab("") +
  scale_x_continuous(limits = c(-4.5, 4.5)) +
  theme(legend.position = "none")

ggsave("CCLE_NK_PRISM_CRISPR_jitterplot.pdf", p, height = 5.25, width = 15)


## ---------------------

# Supplementary table with combined CRISRP and PRISM hits (Table S6F)

# combine PRISM and CRISPR data

prism_crispr <- data_gexp %>% 
left_join(crispr, by = c("featureB" = "Gene"), suffix = c("_prism", "_crispr"))

# add CRISPR effect annotations
prism_crispr <- prism_crispr %>%
mutate(CRISPR_effect = ifelse(!grepl("GOF", cell_line) & p_crispr < p_threshold & lfc>lfc_threshold_pos, "LOF confers resistance", "n.s.")) %>%
mutate(CRISPR_effect = ifelse(grepl("GOF", cell_line) & p_crispr < p_threshold & lfc<lfc_threshold_neg, "GOF sensitizes", CRISPR_effect)) %>%
mutate(CRISPR_effect = ifelse(!grepl("GOF", cell_line) & p_crispr < p_threshold & lfc<lfc_threshold_neg, "LOF sensitizes", CRISPR_effect)) %>%
mutate(CRISPR_effect = ifelse(grepl("GOF", cell_line) & p_crispr < p_threshold & lfc>lfc_threshold_pos, "GOF confers resistance", CRISPR_effect))

# keep only single instances per gene with CRISPR effect annotations
prism_crispr_supplement <- prism_crispr %>% 
select(-signed_p, -lfc, -p_crispr, -FDR, -cell_line, -lfc_z) %>% 
unique() %>% 
group_by(cancer_type, featureB) %>% 
filter(!(CRISPR_effect == "n.s." & n() > 1)) %>% 
rename(p = p_prism) %>% 
mutate(signed_p = sign(cor)*-log10(p),
       cancer_type = factor(cancer_type, levels = c("All", "AML", "T-ALL", "BCL", "MM", "B-ALL"))) %>% 
arrange(cancer_type, signed_p) %>% 
dplyr::select(-test.group, -signed_p, -test.method, -signifCode)

prism_crispr_supplement_signif <- prism_crispr_supplement %>%
filter(p < 0.05, CRISPR_effect != "n.s.") %>% 
filter(cor < 0 & CRISPR_effect %in% c("LOF confers resistance", "GOF sensitizes") | cor > 0 & CRISPR_effect %in% c("GOF confers resistance", "LOF sensitizes") )
fwrite(prism_crispr_supplement_signif, "prism_crispr_supplement.txt")

prism_crispr_supplement_all <- prism_crispr_supplement %>% filter(cancer_type == "All")
prism_crispr_supplement_ball <- prism_crispr_supplement %>% filter(cancer_type == "B-ALL")
prism_crispr_supplement_tall <- prism_crispr_supplement %>% filter(cancer_type == "T-ALL")
prism_crispr_supplement_mm <- prism_crispr_supplement %>% filter(cancer_type == "MM")
prism_crispr_supplement_aml <- prism_crispr_supplement %>% filter(cancer_type == "AML")
prism_crispr_supplement_bcl <- prism_crispr_supplement %>% filter(cancer_type == "BCL")

fwrite(prism_crispr_supplement_all, "prism_crispr_supplement_all.txt")
fwrite(prism_crispr_supplement_ball, "prism_crispr_supplement_ball.txt")
fwrite(prism_crispr_supplement_tall, "prism_crispr_supplement_tall.txt")
fwrite(prism_crispr_supplement_mm, "prism_crispr_supplement_mm.txt")
fwrite(prism_crispr_supplement_aml, "prism_crispr_supplement_aml.txt")
fwrite(prism_crispr_supplement_bcl, "prism_crispr_supplement_bcl.txt")


