
## Correlations of NK CRISPR screen hit expression vs genetic alterations across cancer types

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

# laod plotting scripts
source("../functions/plotting_functions_may2020.R")
source("../CCLE_featurematrix_NK_PRISM/compute.pairwise.R")
source("../CCLE_featurematrix_NK_PRISM/regulon.feats.R")
source("../functions/functions_statistics.R")


# load screen results
data <- fread("crispr_mageck_combined.txt", data.table = F)

genelist <- data %>%
  filter(!grepl("NonTargeting|hsa-mir", Gene)) %>%
  filter(p < 0.0001) %>%
  select(Gene) %>%
  unique() %>%
  deframe()

# CoMMpass

# load feature matrix
fm=get(load("data/MM_COMPASS_FM.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results2, "CoMMpass_NK_CRISPR_correlations.tsv", sep = "\t")

# ------------------------------------------------------------------

# Reddy DLBCL

# load feature matrix
fm=get(load("data/REDDY_DLBCL_fm.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results2, "DLBCL_Reddy_NK_CRISPR_correlations.tsv", sep = "\t")

# ------------------------------------------------------------------

# Chapuy DLBCL

# load feature matrix
fm=get(load("data/GSE98588_fm.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)#, cnv_annot)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results2, "DLBCL_Chapuy_NK_CRISPR_correlations.tsv", sep = "\t")

# ------------------------------------------------------------------

# TCGA DLBCL

# load feature matrix
fm=get(load("../TCGA_AML_DLBCL/DUFVA_TCGA_DLBCL_FM_meth.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)#, cnv_annot)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results2, "TCGA_DLBCL_NK_CRISPR_correlations.tsv", sep = "\t")

# ------------------------------------------------------------------

# TCGA AML

# load feature matrix
fm=get(load("data//DUFVA_TCGA_AML_FM_meth.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)#, cnv_annot)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results2, "TCGA_AML_NK_CRISPR_correlations.tsv", sep = "\t")

# ------------------------------------------------------------------

# TCGA AML

# load feature matrix
fm=get(load("data/DUFVA_TCGA_AML_FM_meth.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results2, "TCGA_AML_NK_CRISPR_correlations.tsv", sep = "\t")

# ------------------------------------------------------------------

# BeatAML

# load feature matrix
fm=get(load("data/BeatAML_fm.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results2, "BeatAML_NK_CRISPR_correlations.tsv", sep = "\t")

# ------------------------------------------------------------------

# load feature matrix
fm=get(load("data/TALL_Liu_fm.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results, "TALL_Liu_NK_CRISPR_correlations.tsv", sep = "\t")

## ---------------------------

# load feature matrix
fm=get(load("data/BALL_Gu_fm.Rdata"))

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:SAMP:|^N:SAMP:|^N:CNVR|^B:CNVR", rownames(fm), value=T))

l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores = 2, adjust.method = "BH")

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
fwrite(results, "BALL_Gu_NK_CRISPR_correlations.tsv", sep = "\t")


