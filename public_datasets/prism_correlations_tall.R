
# Compute pairwise correlations for NK PRISM sensitivity correlates in Liu et al T-ALL data

# load scripts
source("../CCLE_featurematrix_NK_PRISM/regulon.feats.R")
source("../CCLE_featurematrix_NK_PRISM/compute.pairwise.R")
source("../CCLE_featurematrix_NK_PRISM/functions_statistics.R")

# load libraries
library(data.table)
library(parallel)
library(readxl)
library(dplyr)
library(GSVA)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


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
genelist = c(as.character(rownames(gsva_results)), "N:GEXP:FAS", "N:GEXP:PVR", "N:GEXP:ULBP1", "N:GEXP:CD44", "N:GEXP:DACH1", "N:GEXP:PSMG1")  # add genes found in Fig5A analysis

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^B:CNVR|^N:GEXP", rownames(fm_genesets), value=T))

l.regulon.gene=regulon.feats(fm_genesets, genelist)
results=pairwise.correlation(l.regulon.gene, fm_genesets, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

FDR = 0.1
filter.table = data.frame("pair"=unique(results$datapairs), "FDR"=rep(FDR, length(unique(results$datapairs))), stringsAsFactors = F)

type1=gsub(":.*.", "", filter.table[,1])
type2=gsub(".*.:", "", filter.table[,1])
filter.table$FDR[type1==type2]=0.001
filter.table$FDR[type1=="SAMP"&type2=="SAMP"]=0.1 # these can be very different
filter.table$FDR[type1=="GNAB"|type2=="GNAB"]=0.2 # mutations can be interesting, high P-value kept
filter.table$FDR[type1=="METH"|type2=="METH"]=5e-2 # needs to be strong association
filter.table$FDR[type1=="CNVR"&type2=="CNVR"]=0 # remove, highly correlated
filter.table$FDR[type1=="METH"&type2=="METH"]=0 # remove, highly correlated

results2=filter.pairwise.res(results, filter.table = filter.table)

# write result tables
fwrite(results, "Liu_TALL_NK_PRISM_TALL_correlations_pctviability.tsv", sep = "\t")
