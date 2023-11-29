
## DE analysis of samples with vs without mutations in NK CRISPR screen hits in CoMMpass data

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
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(limma)
library(edgeR)


# load screen data

# CRISPR hits 

data <- fread("crispr_mageck_combined.txt", data.table = F)

genelist <- data %>%
  filter(p < 0.0001) %>%
  select(Gene) %>%
  unique() %>%
  tibble::deframe()

# load feature matrix
fm=get(load("MM_COMPASS_FM.Rdata"))

# subset to patients with gexp data
fm <- fm[, !is.na(fm["N:GEXP:KRAS",])]

# find genes with at least 5 alterations in CoMMpass data
genelist_gnab <- paste0("B:GNAB:", genelist)

fm_selected <- fm[rownames(fm) %in% c(genelist_gnab),]
fm_selected <- fm_selected[rowSums(fm_selected, na.rm = T)>4,]

# get gexp data
gexp = fm[grepl("N:GEXP", rownames(fm)),]


# function for DE test, volcanoplot, result table

run_de <- function(ALT, CUTOFF_FDR=1, CUTOFF_PVAL=1, GENES = NULL, TOPGENES=F, NTOPGENES=30){
  
  class = c(rep("wt", dim(gexp)[2]))
  class[fm[ALT,] %in% "1"] = "mut"
  
  design <- model.matrix(~0+class)
  colnames(design) <- c("wt","mut")
  
  cont=makeContrasts(wt - mut, levels = design)
  fit <- lmFit(gexp, design)
  fit <- contrasts.fit(fit, cont)
  
  fit2 <- eBayes(fit, trend = T)
  
  res=topTable(fit2, p.value = 0.05, number = 999999)
  res_all=topTable(fit2, p.value = 1, number = 999999)
  res_all$Gene <- gsub("N:GEXP:", "", rownames(res_all))
  res_all <- res_all[,c(7,1:6)]
  
  # save result table
  fwrite(res_all, paste0("results_DE/", gsub("B:GNAB:", "", ALT), "_CoMMpass_DEG.txt"), sep = "\t")
  
  ## prepare data for volcano plot
  res_all$log10.P.Val <- -log10(res_all$P.Value)
  
  # plot volcano plot
  plotdata <- res_all
  
  if (TOPGENES == T) {
    labelgenes <- plotdata %>%
      top_n(-NTOPGENES, wt = P.Value) %>%
      dplyr::select(Gene) %>%
      tibble::deframe()
  }
  
  else
    
  { labelgenes <- plotdata$Gene }
  
    ggplot(plotdata, aes(x=logFC, y=log10.P.Val, label = Gene)) +
      geom_point(aes(color = ifelse(logFC > 0, "up", "down"))) +
      geom_point(data = plotdata[plotdata$adj.P.Val>0.05,], color = "grey80") +
      geom_text_repel(aes(label=ifelse(adj.P.Val <= CUTOFF_FDR & P.Value <= CUTOFF_PVAL & Gene %in% labelgenes, as.character(Gene),'')),
                      point.padding = 0.2,
                      max.overlaps = 50) +
      ylab("p value (-log10)") +
      xlab("Fold change (log2)") +
      scale_color_manual(values = c(brewer.pal(11, "RdBu")[10], brewer.pal(11, "RdBu")[2])) +
      ggtitle(gsub("B:GNAB:", "", ALT)) +
      guides(color = F) +
      theme_cowplot() +
      theme(plot.title = element_text(hjust = 0, face = "plain"))
  }
  

# CRISPR hit DE analysis
dir.create("results_DE")
p1 <- lapply(rownames(fm_selected), run_de, TOPGENES = T, NTOPGENES=50)
m1 <- marrangeGrob(p1, nrow = 2, ncol = 2, top = NULL)
ggsave("results_DE/CoMMpass_NK_CRISPR_DE_volcanoplots.pdf", m1, height = 20, width = 20)

