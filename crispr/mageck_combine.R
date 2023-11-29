
# Generate combined file with all CRISPR screen results

library(data.table)
library(dplyr)

# load screen results
k562 <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/results_OD/screen_result_tables/k562_3reps_paired.gene_summary.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = neg.lfc, p = ifelse(neg.lfc>0, pos.p.value, neg.p.value), FDR = ifelse(neg.lfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc, p, FDR) %>%
  mutate(cell_line = "K562") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

molm14 <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/results_OD/screen_result_tables/molm14_2reps_paired.gene_summary.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = neg.lfc, p = ifelse(neg.lfc>0, pos.p.value, neg.p.value), FDR = ifelse(neg.lfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc, p, FDR) %>%
  mutate(cell_line = "MOLM14") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

sudhl4 <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/results_OD/screen_result_tables/sudhl4_all_nk_vs_ctrl_nod0_apr20.gene_summary.txt", data.table = FALSE, check.names = T) %>% 
  mutate(lfc = neg.lfc, p = ifelse(neg.lfc>0, pos.p.value, neg.p.value), FDR = ifelse(neg.lfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc, p, FDR) %>%
  mutate(cell_line = "SUDHL4") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

nalm6 <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/results_OD/screen_result_tables/nk2378_nk2542_nalm6_nk_vs_ctrl.gene_summary.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = neg.lfc, p = ifelse(neg.lfc>0, pos.p.value, neg.p.value), FDR = ifelse(neg.lfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc, p, FDR) %>%
  mutate(cell_line = "NALM6") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

mm1s <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/BRUNELLO_NK_mageck/MM1.S_LOF/NK180823_BRUN_TRT-HJTYKBGX7_MM1S_BRUN_pNKw5_vs_CTRL-HJTYKBGX7_MM1S_BRUN_CTRLw5.mageck.genes.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = avgfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "MM1S") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

mm1s_a <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/CALABRESE_NK_mageck_v190523/MITSIADES_CALA_TRT-MM1S_CALA_pNKw5_vs_CTRL-MM1S_CALA_CTRL3w5.norm.median.mageck.genes.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = avgfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "MM1S_GOF") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

kms11_a <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/CALABRESE_NK_mageck_v190523/KMS11_pNK_donor9/CALA_m5.9TRT-KMS11_CALA_NKw5_D9_vs_CTRL-KMS11_CALA_CTRLw5.norm.median.mageck.genes.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = avgfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "KMS11_GOF") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

kms11_khyg1_a <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/CALABRESE_NK_mageck_v190523/KMS11_KHYG1/_CALA_TRT-KMS11_CALA_KHYG1_vs_CTRL-KMS11_CALA_CTRLw5.norm.median.mageck.genes.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = avgfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "KMS11_KHYG1_GOF") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

lp1_a <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/CALABRESE_NK_mageck_v190523/LP1_KHYG1/LP1_CALA_TRT-LP1_CALA_coNKw7_KHYG1_vs_CTRL-LP1_CALA_CTRLw7.norm.median.mageck.genes.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = avgfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "LP1_KHYG1_GOF") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

lp1 <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/BRUNELLO_NK_mageck/LP1_LOF/202001007_LP1_BRUN_NKdonor9w2_vs_CTRL-LP1_BRUN_CTRLw2.norm.median.mageck.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = neg.lfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "LP1") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

kms11_khyg1 <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/BRUNELLO_NK_mageck/KMS11_LOF/HHFCTBGXH_BRUN_NK_v0.5.9_TRT-KHYG1_E2T_6to1_KMS11_BrunelloGW_vs_CTRL-KHYG1_cntr_KMS11_BrunelloGW.norm.total.mageck.genes.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = neg.lfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "KMS11_KHYG1") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

kms11 <- fread("/Users/odufva/Dropbox/NK_resistance Heme CollabPaper/Results_SG/Heme_CRISPR_data/BRUNELLO_NK_mageck/KMS11_LOF/HHFCTBGXH_BRUN_NK_v0.5.9_TRT-pNK_E2T_1to1_KMS11_BrunelloGW_vs_CTRL-pNK_cntr_KMS11_BrunelloGW.norm.median.mageck.genes.txt", data.table = FALSE, check.names = T) %>%
  mutate(lfc = neg.lfc, p = ifelse(avgfc>0, pos.p.value, neg.p.value), FDR = ifelse(avgfc>0, pos.fdr, neg.fdr)) %>%
  select(Gene = id, lfc = avgfc, p, FDR) %>%
  mutate(cell_line = "KMS11") %>%
  mutate(lfc_z = (lfc - mean(lfc))/sd(lfc))

# combine data
crispr <- rbind(k562, molm14, sudhl4, nalm6, mm1s, lp1, kms11, kms11_khyg1, mm1s_a, kms11_a, kms11_khyg1_a, lp1_a)

fwrite(crispr, "crispr_mageck_combined.txt")

