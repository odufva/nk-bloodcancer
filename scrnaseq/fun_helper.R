
require(clusterProfiler)
require(org.Hs.eg.db)

if(me == "hru"){
  
  hallmark   <- read.gmt("/Users/hru/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")
  
  
}


if(me == "janihuuh"){
  
  hallmark   <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")
  
}


plotHypergeometric <- function(genes_df, universe_df, term_df){
  
  require(clusterProfiler)
  
  if(nrow(genes_df) == 0) return(NULL)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
  
  # out: df with enrichment results
  
  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  
  if(table(enrich@result$p.adjust < 0.05) %>% length() > 1){
    heatplot(enrich)
  }
  
  else(NULL)
  
}

getHypergeometric <- function(genes_df, universe_df, term_df){
  
  require(clusterProfiler)
  
  if(nrow(genes_df) == 0) return(NULL)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
  
  # out: df with enrichment results
  
  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  enrich <- do.call(rbind, enrich@result) %>% t %>% as.data.frame()
  enrich[,c(5:7, 9)] <- sapply(enrich[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
  return(enrich)
  
}



getDEGsgRNA <- function(seurat_object, sgrna_genes){
  
  Idents(seurat_object) <- "gene_functional"
  lapply(sgrna_genes, function(x) FindMarkers(seurat_object, ident.1 = x, ident.2 = "Control", test.use = "t", max.cells.per.ident = 100, logfc.threshold = 0.05, min.pct = 0.01) %>% 
           add_rownames(var = "gene") %>% mutate(gene_to_test = paste(x, "vs NT"))) %>%
    rbindlist() %>% mutate(p_val_adj = p.adjust(p_val, method = "BH") ) %>% arrange(p_val_adj) %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down"))
  
}

getGLMtest <- function(seurat_object, gene_to_test, sgrna_to_test){
  
  ## Build models for each gene
  
  cells.to.keep   <- seurat_object@meta.data %>% filter(gene_functional %in% c(sgrna_to_test, "Control")) %>% pull(barcode)
  seurat_temp     <- subset(crop_seurat, cells = cells.to.keep)
  genes_expressed <- rownames(seurat_temp)[rowSums(seurat_temp) > 0]
  
  if(gene_to_test %in% genes_expressed){
    
    y      <- seurat_temp@assays$RNA@counts[rownames(crop_seurat) == gene_to_test, ]
    group1 <- droplevels(as.factor(seurat_temp$gene_functional))
    group2 <- droplevels(as.factor(seurat_temp$orig.ident))
    group3 <- droplevels(as.factor(seurat_temp$cluster))
    
    tapply(y, group1, mean)
    tapply(y, group1, sd)
    
    tapply(y, group2, mean)
    tapply(y, group2, sd)
    
    fit0 <- glm(y ~ 1, family = poisson())
    fit1 <- glm(y ~ group1, family = poisson())
    fit2 <- glm(y ~ group1 + group2, family = poisson())
    fit3 <- glm(y ~ group1 + group2 + group3, family = poisson())
    
    # p.df <- anova(fit0, fit1, test="Chisq")
    # plot(p.df)
    # plot(fit0)
    
    p.df <- anova(fit0, fit1, test="Chisq") %>% broom::tidy() %>% as.data.frame() # %>% mutate(sgrna = sgrna_to_test, gene = gene_to_test)
    p.df <- p.df[2,]
    
    p.df2 <- anova(fit0, fit2, test="Chisq") %>% broom::tidy() %>% as.data.frame() # %>% mutate(sgrna = sgrna_to_test, gene = gene_to_test)
    p.df2 <- p.df2[2,]
    
    p.df3 <- anova(fit0, fit3, test="Chisq") %>% broom::tidy() %>% as.data.frame() # %>% mutate(sgrna = sgrna_to_test, gene = gene_to_test)
    p.df3 <- p.df3[2,]
    
    p.tot <- data.frame(p.sgrna = p.df$p.value, 
                        p.sgrna.ident = p.df2$p.value, 
                        p.sgrna.ident.cluster = p.df3$p.value, 
                        
                        # dev.sgrna = p.df$Deviance,
                        # dev.sgrna.ident = p.df2$Deviance, 
                        # dev.sgrna.ident.cluster = p.df3$Deviance, 
                        
                        sgrna = sgrna_to_test, gene = gene_to_test, 
                        mean_x = tapply(y, group1, mean)[1], 
                        mean_y = tapply(y, group1, mean)[2]) %>% mutate(log2fc = log2(mean_y/mean_x))
    return(p.tot)
    
  }
  
}


getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){
  
  hpca.se   <- celldex::HumanPrimaryCellAtlasData()
  blueprint <- celldex::BlueprintEncodeData()
  ## @ params
  ## cluster = possible cluster vec, if not provided, tries to find in meta.data$cluster
  ## method = if "cluster", then performs preds based on clusters, not cells
  ## sample = to subsample or not 
  
  if(!is.null(sample)){
    
    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])
    
  }
  
  sce       <- as.SingleCellExperiment(seurat_object)
  
  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)
    
    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$labels
      seurat_object$singler_blueprint_pred <- pred.blu$labels
      return(seurat_object)  
    }
    
    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }
    
  }
  
  
  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}




getDoublets <- function(seurat_object){
  
  require(scds)
  
  # Annotate doublet using co-expression based doublet scoring:
  sce_object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = seurat_object@assays$RNA@counts), colData = seurat_object@meta.data)
  sce_object <- cxds(sce_object)
  sce_object <- bcds(sce_object)
  sce_object <- cxds_bcds_hybrid(sce_object)
  
  ## Add into Seurat
  seurat_object$cxds_doublet_score   <- SingleCellExperiment::colData(sce_object)$cxds_score
  seurat_object$bcds_doublet_score   <- SingleCellExperiment::colData(sce_object)$bcds_score
  
  seurat_object$hybrid_doublet_score <- SingleCellExperiment::colData(sce_object)$hybrid_score
  seurat_object$cxds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$cxds_score - min(SingleCellExperiment::colData(sce_object)$cxds_score)) / max(SingleCellExperiment::colData(sce_object)$cxds_score)
  seurat_object$bcds_doublet_score_norm <- c(SingleCellExperiment::colData(sce_object)$bcds_score - min(SingleCellExperiment::colData(sce_object)$bcds_score)) / max(SingleCellExperiment::colData(sce_object)$bcds_score)
  return(seurat_object)
  
}






plotQC <- function(seurat_object, folder, min_mito = 0, max_mito = 10, min_ribo = 5, max_ribo = 50, min_features = 300, max_features = 10e3, min_counts = 3e3, max_counts = 50e4){
  
  qc_df <- seurat_object@meta.data %>% as.data.frame()
  qc_df$cluster = Idents(seurat_object)
  
  plotQcViolin(qc_df, var_to_plot = "nFeature_RNA", grouping = "cluster", min = min_features, max = max_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_nFeature_RNA.png"), width = 12, height = 6)
  
  plotQcViolin(qc_df, var_to_plot = "nCount_RNA", grouping = "cluster", min = min_counts, max = max_counts) + scale_y_log10() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_nCount_RNA.png"), width = 12, height = 6)
  
  plotQcViolin(qc_df, var_to_plot = "percent.mt", grouping = "cluster", min = min_mito, max = max_mito) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "violin_percent_mt.png"), width = 12, height = 6)
  
  plotQcViolin(qc_df, var_to_plot = "percent.ribo", grouping = "cluster", min = min_ribo, max = max_ribo) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(folder, "/violin_percent_ribo.png"), width = 12, height = 6)
  
}

getQC <- function(seurat_object, min_mito = 0, max_mito = 10, min_ribo = 5, max_ribo = 50,
                  min_features = 300, max_features = 10e3, min_counts = 3e3, max_counts = 50e4){
  
  ###################
  
  # min_mito     <- 0
  # max_mito     <- 10
  # 
  # min_ribo     <- 5
  # max_ribo     <- 50
  # 
  # min_features <- 300
  # max_features <- 10e3
  # 
  # min_counts   <- 3e3
  # max_counts   <- 50e4
  
  ###################
  
  seurat_object@meta.data$barcode <- colnames(seurat_object)
  
  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()
  
  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()
  
  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)
  
  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))
  
  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))
  
  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)
  
}



getQCexpansionNK <- function(seurat_object){
  
  ###################
  
  min_mito     <- 0
  max_mito     <- 15
  
  min_ribo     <- 10
  max_ribo     <- 50
  
  min_features <- 250
  max_features <- 6000
  
  min_counts   <- 1000
  max_counts   <- 30e3
  
  cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                    "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                    "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                    "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                    "TUBB", "TYMS", "UBE2C")
  
  ###################
  
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^RP", col.name = "percent.ribo")
  seurat_object@meta.data$barcode <- colnames(seurat_object)
  
  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()
  
  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()
  
  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)
  
  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))
  
  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))
  
  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)
  
}


extractClusterNumber <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }
  
  return(p)
  
}


getClusterPhenotypesXXX <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    "0"  = "0 ",
    "1"  = "1 " ,
    "2"  = "2 " ,
    "3"  = "3 " ,
    "4"  = "4 " ,
    "5"  = "5 " ,
    "6"  = "6 " ,
    "7"  = "7 " ,
    "8"  = "8 " ,
    "9"  = "9 " ,
    "10" = "10 " ,
    "11" = "11 " ,
    "12" = "12 " ,
    "13" = "13 " ,
    "14" = "14 " ,
    "15" = "15 " ,
    "16" = "16 " ,
    "17" = "17 " ,
    "18" = "18 " ,
    "19" = "19 " ,
    "20" = "20 " ,
    "21" = "21 " ,
    "22" = "22 ",
    "23" = "23 ",
    "24" = "24 ",
    "25" = "25 "))
  
  return(clusters)
  
}


plotClustering <- function(seurat_object){
  
  res           <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]
  
  q <- NULL; i <- 1
  
  for(clustering_column in clustering_columns){
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
  }
  
  data.frame(resolution = res, nClusters = do.call(q, what="c")) %>%
    ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
  
}

getClustering <- function(seurat_object){
  
  nPCs <- sum(seurat_object@reductions$pca@stdev > 2)
  
  ## Clustering
  res           <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:nPCs, force.recalc = T)
  seurat_object <- FindClusters(seurat_object, resolution = res, verbose = F, )
  
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]
  
  q <- NULL; i <- 1
  
  for(clustering_column in clustering_columns){
    
    message(clustering_column)
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
    
  }
  
  p <- data.frame(resolution = res, nClusters = do.call(q, what = "c")) %>%
    ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
  print(p)
  
  return(seurat_object)
  
}




preprocessSeurat <- function(orig_object, cells.to.use, vars.to.regress = NULL){
  
  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)
  
  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta
  
  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)
  
  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]
  
  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))
  
  ## Scale data
  object <- ScaleData(object, features = hvg, vars.to.regress = vars.to.regress)
  
  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  message(paste("nPCs:", nPCs))
  
  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)
  return(object)
  
}


facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractSeuratName <- function(str1){
  
  str1 <- substr(str1, 1, nchar(str1) - 27)
  extractFileName(str1)
}

extractTimepoint <- function(strs){
  
  strs2 <-NULL
  i <- 1
  for(str1 in strs){
    strs2[[i]] <- strsplit(str1, "[_]")[[1]][2]
    i <- i + 1
  }
  
  return(strs2)
  
}


plotQcViolin <- function(viz_df, var_to_plot, grouping, min, max){
  
  ## Plot univariate violin plots with filter thresholds
  
  # @ params:
  # viz_df = df that contains qc-analysis results and covariates of interest
  # var_to_plot = char, a column name that contains the variable to plot
  # grouping = char, a column name that contains the x-axis grouping
  # min = num, min value for variable
  # max = num, max value for variable
  
  viz_df_temp <- viz_df %>% dplyr::select(var_to_plot)
  
  label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
  label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table
  
  ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
    geom_violin(alpha = 0.5) +
    # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
    
    geom_hline(yintercept = min, linetype = "dotted") +
    geom_hline(yintercept = max, linetype = "dotted") +
    
    annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
    annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +
    
    labs(x = "", title = var_to_plot) + theme(legend.position = "none")
  
}






fixSeurat <- function(seurat_object){
  
  ## Fix meta data if it brokes
  
  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)
  
  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]
  
  rownames(meta.data) <- meta.data$barcode
  
  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)
  
  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg
  
  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)
  
}





extractCoarsePhenotype <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }
  
  return(p)
  
}




















