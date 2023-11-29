
# Ridge plots of individual genes in CROP-seq and patient data using mixscape filtered data (Figure 6G)

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(viridis)
library(patchwork)
library(ggridges)
library(ggpubr)


dir.create("results/combine/patient_integration")

# load Seurat objects
mm1s_nonk <- readRDS("results/mm1s/mixscape/singlet/nonk_ko_nt.rds")
mm1s_1_16 <- readRDS("results/mm1s/mixscape/singlet/nk_1_16_ko_nt.rds")
mm1s_1_4 <- readRDS("results/mm1s/mixscape/singlet/nk_1_4_ko_nt.rds")
sudhl4_nonk <- readRDS("results/sudhl4/mixscape/singlet/nonk_ko_nt.rds")
sudhl4_1_16 <- readRDS("results/sudhl4/mixscape/singlet/nk_1_16_ko_nt.rds")
sudhl4_1_4 <- readRDS("results/sudhl4/mixscape/singlet/nk_1_4_ko_nt.rds")

# load patient data feature matrices
commpass_fm <- get(load("data/MM_COMPASS_FM.Rdata"))
commpass_df <- as.data.frame(t(commpass_fm)) 

reddy_fm=get(load("data/REDDY_DLBCL_fm.Rdata"))
reddy_df <- as.data.frame(t(reddy_fm))

aml_fm=get(load("data/DUFVA_TCGA_AML_FM_meth.Rdata"))
aml_df <- as.data.frame(t(aml_fm))

# load CROP-seq vs patient data analysis result
result <- fread("results/combine/patient_integration/cropseq_patient_integration_topcondition.txt", data.table = F)

result %>% 
  filter((gene == "HLA-E" & perturbation == "NLRC5") |
           (gene == "BIRC3" & perturbation == "NFKBIA") |
           (gene == "NFKB2" & perturbation == "NFKBIA") |
           (gene == "CD74" & perturbation == "NFKBIA") |
           (gene == "HLA-E" & perturbation == "NFKBIA") |
           (gene == "BIRC3" & perturbation == "TRAF2") |
           (gene == "NFKB2" & perturbation == "TRAF2") |
           (gene == "CXCL10" & perturbation == "TRAF2") |
           (gene == "CCND2" & perturbation == "TRAF2") |
           (gene == "RGS1" & perturbation == "TRAF2") |
           (gene == "ISG20" & perturbation == "PTEN") |
           (gene == "CFLAR" & perturbation == "PTEN") |
           (gene == "CD44" & perturbation == "PTEN") |
           (gene == "IRF4" & perturbation == "PTEN") |
           (gene == "CFLAR" & perturbation == "GNA13") |
           (gene == "ISG20" & perturbation == "RUNX1"))


# ridge plots of selected DEG in CROP-seq data

plot_ridgeplot <- function(data, perturbation, targetgene, cols){
  
  seurat <- subset(data, gene %in% c(perturbation, "NT"))
  
  plotdata <- FetchData(seurat, vars = c(targetgene, "gene"), slot = "data")
  plotdata$gene <- factor(plotdata$gene, levels = c("NT", perturbation))
  colnames(plotdata)[colnames(plotdata)!="gene"] <- "target_gene"
  
  ggplot(plotdata, aes(x = target_gene, y = gene, fill = gene)) +
    geom_density_ridges(alpha = 0.75, scale = 10000) +
    ylab("") +
    ggtitle("CROP-seq") +
    xlab(substitute(italic(x)~" expression", list(x = targetgene))) +
    theme_cowplot() +
    theme(axis.title.y = element_text(vjust = 0.5, angle = 0, face = "italic"),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "plain"),
          axis.line.y = element_blank()) +#,
    scale_fill_manual(values = cols) +
    guides(fill = "none")
  
}

p1 <- plot_ridgeplot(mm1s_nonk, "TRAF2", "BIRC3", brewer.pal(12, "Paired")[c(8,7)])
p2 <- plot_ridgeplot(mm1s_nonk, "TRAF2", "CXCL10", brewer.pal(12, "Paired")[c(8,7)])
p3 <- plot_ridgeplot(mm1s_nonk, "NFKBIA", "BIRC3", brewer.pal(12, "Paired")[c(8,7)])
p4 <- plot_ridgeplot(mm1s_1_4, "NFKBIA", "HLA-E", brewer.pal(12, "Paired")[c(8,7)])
p5 <- plot_ridgeplot(mm1s_1_16, "NLRC5", "HLA-E", brewer.pal(12, "Paired")[c(8,7)])
p6 <- plot_ridgeplot(sudhl4_nonk, "PTEN", "ISG20", brewer.pal(12, "Paired")[c(10,9)])
p7 <- plot_ridgeplot(sudhl4_1_16, "RUNX1", "ISG20", brewer.pal(12, "Paired")[c(10,9)])


# box plots of selected DEG in patient data

# function for boxplot of altered vs unaltered cases
alt_bplot <- function(DATA, GEXP, ALT, ALTTOPLOT=NULL, COLORS){

  df <- DATA[,c(GEXP, ALT)] 
  colnames(df) <- c("gexp", "alt") 
  
  # discretize CNVs
  if (grepl("CNVR", ALT)) {
    df$alt_discrete[df$alt < -1.5] <- "del"   
    df$alt_discrete[df$alt < -0.5 & df$alt >= -1.5] <- "loss"   
    df$alt_discrete[df$alt > 0.25  & df$alt <= 0.75] <- "gain"   
    df$alt_discrete[df$alt > 0.75] <- "amp"   
    df$alt_discrete[df$alt <= 0.25 & df$alt >= -0.5] <- "wt"   
    df$alt <- df$alt_discrete
    df$alt[df$alt%in%c(ALTTOPLOT)] <- "Mut"
    df$alt[!(df$alt%in%c("+"))] <- "WT"
  }
  
  xlabel <- gsub("N:GEXP:", "", gsub(":::::|:chr.*", "", GEXP))
  
  alt_short <- gsub("([0-9])P([0-9])", "\\1p\\2", gsub("([0-9])Q([0-9])", "\\1q\\2", gsub("Q$", "q", gsub("P$", "p", gsub("B:.....|N:.....", "", gsub("_", ".", gsub("SV_", "", gsub("_nonsynonymous:::::", "", gsub(":AMP", " amp ", gsub(":LOSS", " loss", gsub(":DEL", " del ", gsub(":SV", " SV  ", ALT))))))))))))
  alt_short <- ifelse(grepl("GNAB", ALT), paste0(alt_short, " mut "), alt_short)
  
  df <- df %>%
    mutate(gexp = as.numeric(as.character(gexp)),
           alt = gsub("^1", "Mut", 
                      gsub("^0", "WT", alt))) %>%
    mutate(alt = factor(alt, levels = sort(unique(alt)))) %>%
    na.omit()
  
  ggplot(df, aes(x = alt, y = gexp, fill = alt)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, color = "grey20", size = 1) +
    ylim(c(min(df$gexp)-(max(df$gexp)-min(df$gexp))*0.1, max(df$gexp)+(max(df$gexp)-min(df$gexp))*0.1)) +
    scale_fill_manual(values = COLORS) +
    ylab(substitute(italic(x)~" expression (log2)", list(x = xlabel))) +
    ggtitle("MM patients") +
    xlab("") +
    theme_cowplot() +
    theme(legend.position="none",
          axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "plain")) +
    theme(plot.margin = unit(c(0,2,0,1), "cm"))
    }

# plot selected alterations
commpass1 <- alt_bplot(commpass_df, "N:GEXP:BIRC3", "B:GNAB:TRAF2", COLORS = brewer.pal(12, "Paired")[c(8,7)])
commpass2 <- alt_bplot(commpass_df, "N:GEXP:CXCL10", "B:GNAB:TRAF2", COLORS = brewer.pal(12, "Paired")[c(8,7)])
commpass3 <- alt_bplot(commpass_df, "N:GEXP:BIRC3", "B:GNAB:NFKBIA", COLORS = brewer.pal(12, "Paired")[c(8,7)])
commpass4 <- alt_bplot(commpass_df, "N:GEXP:HLA-E", "B:GNAB:NFKBIA", COLORS = brewer.pal(12, "Paired")[c(8,7)])
commpass5 <- alt_bplot(commpass_df, "N:GEXP:HLA-E", "B:GNAB:NLRC5", COLORS = brewer.pal(12, "Paired")[c(8,7)])

reddy1 <- alt_bplot(reddy_df, "N:GEXP:ISG20", "B:GNAB:PTEN", COLORS = brewer.pal(12, "Paired")[c(10,9)])
aml1 <- alt_bplot(aml_df, "N:GEXP:ISG20:chr15:89179384:89199714:+:3669", "B:GNAB:RUNX1:chr21:36160098:37357047:-:y_n_somatic", COLORS = brewer.pal(12, "Paired")[c(2,1)])


# combined plots

(p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(ncol = 1)) |
  (commpass1 + commpass2 + commpass3 + commpass4 + commpass5 + reddy1 + aml1 + plot_layout(ncol = 1))
ggsave("results/combine/patient_integration/ridgeplots_mixscape.pdf", height = 10, width = 8)

# Figure 6G
p1 + commpass1 + p2 + commpass2 + p3 + commpass3 + p5 + commpass5 + plot_layout(nrow = 1, widths = c(1,0.5,1,0.5,1,0.5,1,0.5))
ggsave("results/combine/patient_integration/ridgeplots_mixscape_horizontal.pdf", height = 2.5, width = 18)

