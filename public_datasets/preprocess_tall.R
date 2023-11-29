
# Generate T-ALL Liu et al. feature matrix

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)

source("../prism/plotting_functions.R")
source("../prism/functions_generate_fm.R")

# load supplementary tables from Liu et al. Nature Genetics 2017
clin <- read_excel("41588_2017_BFng3909_MOESM2_ESM.xlsx", sheet = 1)
gexp <- read_excel("41588_2017_BFng3909_MOESM2_ESM.xlsx", sheet = 5)
mut <- read_excel("41588_2017_BFng3909_MOESM2_ESM.xlsx", sheet = 8)

# prepare gexp matrix
gexp_mat <- as.matrix(gexp[,-1])
rownames(gexp_mat) <- gexp$Gene
gexp_mat <- log2(gexp_mat+0.01)

clin_match <- clin[clin$RNAseq_id_D %in% colnames(gexp_mat),]
gexp_mat <- gexp_mat[,colnames(gexp_mat) %in% clin_match$RNAseq_id_D]
colnames(gexp_mat) <- clin_match$USI[match(colnames(gexp_mat), clin_match$RNAseq_id_D)]
gexp_mat <- gexp_mat[,!is.na(colnames(gexp_mat))]

# prepare mutation matrix
mut_mat <- mut %>%
  select(gene, sample) %>%
  mutate(mutation = 1) %>%
  pivot_wider(names_from = sample, values_from = mutation, values_fn = mean) %>%
  as.data.frame()
  
rownames(mut_mat) <- mut_mat$gene
mut_mat$gene <- NULL
mut_mat <- as.data.frame(t(mut_mat))
mut_mat[is.na(mut_mat)] <- 0


samples <- intersect(colnames(gexp_mat), rownames(mut_mat))
clin_match=clin_match[match(samples,clin_match$USI),]
mut_mat=mut_mat[samples,]
gexp_mat=gexp_mat[,samples]

clindatfm=make.features(df = as.data.frame(clin_match), datatype="CLIN", prefix="", make.pairwise = F)
colnames(clindatfm)=samples

mutfm=make.features(mut_mat, datatype="GNAB", prefix="", make.pairwise = F)

gexpfm=data.frame(gexp_mat)
rownames(gexpfm)=paste0("N:GEXP:", rownames(gexp_mat))

l.fm=list(clindatfm, mutfm, gexpfm)

fm=rbindlist(l.fm, use.names=F, fill=F)

fm=data.frame(fm, stringsAsFactors=F)
rownames(fm)=unlist(lapply(l.fm, rownames))
annot=data.frame(clin_match, stringsAsFactors = F)

save(fm, file="TALL_Liu_fm.Rdata")
save(annot, file="TALL_Liu_annot.Rdata")
