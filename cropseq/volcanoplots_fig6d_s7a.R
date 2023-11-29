
# Plot volcano plots of KO vs Ctrl in CROP-seq data (Figures 6D and S7A)

# load libraries
library(patchwork)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(gridExtra)

# load perturbation results
nalm6 <- fread("results/nalm6/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "NALM6")
sudhl4 <- fread("results/sudhl4/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "SUDHL4")
k562 <- fread("results/k562/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "K562")
mm1s <- fread("results/mm1s/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S")
lp1 <- fread("results/lp1/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "LP1")
mm1s_a <- fread("results/mm1sa/deg/singlet/pert_deg_all.txt", data.table = F) %>% mutate(cell_line = "MM1S CRISPRa")

# combine data
data <- rbind(k562, sudhl4, nalm6, mm1s, lp1, mm1s_a)
              
# list cell lines and conditions with most DEGs

dir.create("results/combine/volcanoplots")

data <- rbind(k562, sudhl4, nalm6, mm1s, lp1, mm1s_a) %>%  mutate(cell_line_condition_perturbation = paste(cell_line, condition, perturbation, sep = " "))


## Volcano plots

# volcano plot function
volcanoplot <- function(CELLLINE_COND_PERT, CUTOFF_FDR=1, CUTOFF_PVAL=1, GENES = NULL, TOPGENES=T, NTOPGENES=30){
 
   plotdata <- data %>%
    filter(cell_line_condition_perturbation == CELLLINE_COND_PERT)
    
    CELLLINE = unique(plotdata$cell_line)
    PERTURBATION = unique(plotdata$perturbation)
    CONDITION = unique(plotdata$condition)

  if (TOPGENES == T) {
    labelgenes <- plotdata %>%
      top_n(-NTOPGENES, wt = p_val) %>%
      dplyr::select(gene) %>%
      tibble::deframe()
  }
  
  else
    
  { labelgenes <- plotdata$gene }
  
  
  if (is.null(GENES)) {
    
    ggplot(plotdata, aes(x=avg_log2FC, y=-log10(p_val), label = gene)) +
      geom_point(aes(color = ifelse(avg_log2FC > 0, "up", "down"))) +
      geom_point(data = plotdata[plotdata$p_val_adj>0.05,], color = "grey80") +
      geom_text_repel(aes(label=ifelse(p_val_adj <= CUTOFF_FDR & p_val <= CUTOFF_PVAL & gene %in% labelgenes, as.character(gene),'')),
                      point.padding = 0.2,
                      max.overlaps = 20) +
      xlab("") +
      ylab("") +
      scale_color_manual(values = c(brewer.pal(11, "RdBu")[10], brewer.pal(11, "RdBu")[2])) +
      ggtitle(paste0(PERTURBATION, " (", CELLLINE, " ", CONDITION, ")")) +
      guides(color = F) +
      theme_cowplot() +
      theme(plot.title = element_text(hjust = 0.5, face = "plain")) +
      scale_x_continuous(limits = c(max(abs(plotdata$avg_log2FC))*(-1), max(abs(plotdata$avg_log2FC))))
  }
  
  else
    
  {  ggplot(data[data$gene_to_test == COMPARISON,], aes(x=avg_log2FC, y=-log10(p_val), label = gene)) +
      geom_point(aes(color = ifelse(avg_log2FC > 0, "up", "down"))) +
      geom_point(data = plotdata[plotdata$p_val_adj>0.05,], color = "grey80") +
      geom_point(data = data[data$gene_to_test == COMPARISON & data$gene %in% GENES,], aes(fill = ifelse(avg_logFC > 0, "up", "down")),
                 color = "black", size = 3, pch = 21) +
      geom_text_repel(aes(label=ifelse(gene %in% GENES, as.character(gene), '')),
                      point.padding = 0.2,
                      max.overlaps = 20) +
      ylab("p value (-log10)") +
      xlab("Average fold change (log2)") +
      scale_color_manual(values = c(brewer.pal(11, "RdBu")[10], brewer.pal(11, "RdBu")[2])) +
      scale_fill_manual(values = c(brewer.pal(11, "RdBu")[2], brewer.pal(11, "RdBu")[10])) +
      ggtitle(COMPARISON) +
      guides(color = F, fill = F) +
      theme_cowplot() +
      theme(plot.title = element_text(hjust = 0, face = "plain")) +
      scale_x_continuous(limits = c(-0.6, 0.6))
  }
  
}

# volcanoplots for Figure S7A
selected_conditions <- c("SUDHL4 NK 1:16 JAK1",
                         "MM1S NK 1:16 NLRC5", "MM1S CRISPRa NK 1:16 NLRC5", "MM1S NK 1:16 RFXAP",
                         "K562 NK 1:16 GFI1B", "SUDHL4 NK 1:16 YTHDF2",  "MM1S no NK PCGF5", "MM1S NK 1:16 GNA13", "NALM6 NK 1:16 KIAA0922",
                         "MM1S no NK TRAF2",  "MM1S no NK NFKBIB",
                          "NALM6 NK 1:16 CMIP")

p <- lapply(selected_conditions, volcanoplot)
m <- marrangeGrob(p, top = NULL, ncol = 4, nrow = 3, layout_matrix = matrix(1:12, 3, 4, TRUE))

ggsave("results/combine/volcanoplots/volcanoplots_supplement_selected.pdf", m, height = 10, width = 20)



# volcano plot function for main figure (Figure 6D)
volcanoplot_mainfigure <- function(CELLLINE_COND_PERT, CUTOFF_FDR=0.05, CUTOFF_PVAL=0.05, GENES = NULL, TOPGENES=T, NTOPGENES=30){
  
  
  plotdata <- data %>%
    filter(cell_line_condition_perturbation == CELLLINE_COND_PERT)
  
  CELLLINE = unique(plotdata$cell_line)
  PERTURBATION = unique(plotdata$perturbation)
  CONDITION = unique(plotdata$condition)
  
  if (TOPGENES == T) {
    labelgenes <- plotdata %>%
      top_n(-NTOPGENES, wt = p_val) %>%
      dplyr::select(gene) %>%
      tibble::deframe()
  }
  
  else
    
  { labelgenes <- plotdata$gene }
  
  
  if (is.null(GENES)) {
    
    ggplot(plotdata, aes(x=avg_log2FC, y=-log10(p_val), label = gene)) +
      geom_point(aes(color = ifelse(avg_log2FC > 0, "up", "down"))) +
      geom_point(data = plotdata[plotdata$p_val_adj>0.05,], color = "grey80") +
      geom_text_repel(aes(label=ifelse(p_val_adj <= CUTOFF_FDR & p_val <= CUTOFF_PVAL & gene %in% labelgenes, as.character(gene),'')),
                      point.padding = 0.2,
                      label.padding = 0.5,
                      force = 2,
                      max.overlaps = 100) +
      ylab("p value (-log10)") +
      xlab("Average fold change (log2)") +
      scale_color_manual(values = c(brewer.pal(11, "RdBu")[10], brewer.pal(11, "RdBu")[2])) +
      ggtitle(paste0(PERTURBATION, " (", CELLLINE, " ", CONDITION, ")")) +
      guides(color = F) +
      theme_cowplot() +
      theme(plot.title = element_text(hjust = 0.5, face = "plain"),
            plot.margin = unit(c(0,0,0,0), "cm")) +
      scale_x_continuous(limits = c(max(abs(plotdata$avg_log2FC))*(-1), max(abs(plotdata$avg_log2FC))))
  }
  
  else
    
  {   ggplot(plotdata, aes(x=avg_log2FC, y=-log10(p_val), label = gene)) +
      geom_point(aes(color = ifelse(avg_log2FC > 0, "up", "down"))) +
      geom_point(data = plotdata[plotdata$p_val_adj>0.05,], color = "grey80") +
      geom_text_repel(aes(label=ifelse(p_val_adj <= CUTOFF_FDR & p_val <= CUTOFF_PVAL & gene %in% c(labelgenes, GENES), as.character(gene),'')),
                      point.padding = 0.2,
                      max.overlaps = 100) +
      ylab("p value (-log10)") +
      xlab("Average fold change (log2)") +
      scale_color_manual(values = c(brewer.pal(11, "RdBu")[10], brewer.pal(11, "RdBu")[2])) +
      ggtitle(paste0(PERTURBATION, " (", CELLLINE, " ", CONDITION, ")")) +
      guides(color = F) +
      theme_cowplot() +
      theme(plot.title = element_text(hjust = 0.5, face = "plain"),
            plot.margin = unit(c(0,0,0,0), "cm")) +
      scale_x_continuous(limits = c(max(abs(plotdata$avg_log2FC))*(-1), max(abs(plotdata$avg_log2FC))))
    
  }
  
}

selected_conditions <- c("SUDHL4 NK 1:16 IFNGR2", "MM1S no NK NFKBIA",  "LP1 no NK MYB", "LP1 no NK GSK3B")

p1 <- volcanoplot_mainfigure("MM1S no NK NFKBIA", NTOPGENES = 20) + xlab("")
p2 <- volcanoplot_mainfigure("LP1 no NK MYB", NTOPGENES = 10, GENES = c("HLA-DRA", "HLA-DRB1", "HES6", "CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DMA", "CD38", "HLA-E", "HLA-B", "CIITA")) + xlab("")
p3 <- volcanoplot_mainfigure("LP1 no NK GSK3B", NTOPGENES = 12)

p1 / p2 / p3

ggsave("results/combine/volcanoplots/volcanoplots_mainfigure_selected_3.pdf", height = 9, width = 5.5)

