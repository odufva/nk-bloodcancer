
# Run cellphoneDB on interacting NK and cancer cells

dir.create("results/cellphonedb/", showWarnings = F)
dir.create("results/cellphonedb/input_files/", showWarnings = F)

pan_nk_seurat <- readRDS("results/celllinepanel/hto_singlets/nk/batchcorrected/hash_seurat_nk_batchcorrected.rds")

file_names        <- list.files("results/celllinepanel/hto_singlets/targets/mainclusters/", pattern = "rds", full.names = T)
scrnaseq_files    <- lapply(file_names, FUN = function(x){message(x); readRDS(x)})
pan_target_seurat <- merge(scrnaseq_files[[1]], scrnaseq_files[-1])

saveRDS(pan_nk_seurat, "results/cellphonedb/pan_nk_seurat.rds")
saveRDS(pan_target_seurat, "results/cellphonedb/pan_target_seurat.rds")

## Only expanded NK cells
cells.to.keep         <- pan_target_seurat@meta.data %>% filter(!grepl("PBMC", pan_target_seurat$hash.ID)) %>% pull(barcode)
pan_target_exp_seurat <- subset(pan_target_seurat, cells = cells.to.keep)

pan_nk_seurat$cellphonedb     <- gsub("\\(", "\\_", pan_nk_seurat$cluster)
pan_nk_seurat$cellphonedb     <- gsub("\\)", "\\_", pan_nk_seurat$cellphonedb)
pan_nk_seurat$cellphonedb     <- paste0("effector_", pan_nk_seurat$cellphonedb)

pan_target_seurat$cellphonedb <- gsub("\\-", "\\_", pan_target_exp_seurat$hash.ID)
pan_target_seurat$cellphonedb <- paste0("target_", pan_target_seurat$cellphonedb)

pan_nk_target_exp_seurat <- merge(pan_nk_seurat, pan_target_seurat)
pan_nk_target_exp_seurat <- pan_nk_target_exp_seurat %>% preprocessSeurat(cells.to.use = colnames(pan_nk_target_exp_seurat))
saveRDS(pan_nk_target_exp_seurat, "results/cellphonedb/pan_nk_target_exp_seurat.rds")

data.frame(pan_nk_target_exp_seurat@meta.data) %>% 
  mutate(cluster = cellphonedb) %>% 
  initCellphonedb(seurat_object = pan_nk_target_exp_seurat, 
                  name = "pan_expanded_NK_pan_target", 
                  sample_n = 10, 
                  folder = "results/cellphonedb/input_files/")


