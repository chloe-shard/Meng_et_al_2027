setwd("I:/Data/Aijia/TAOF_R/data") 

# Change GBO to 

GBO <- readRDS("scRNA_seq/GBO/GBO_malignant_qc_cellstate_hypsig_mergealias_region.rds")

# Log-Norm, Scale, FindFeature, UMAP, PCA etc before DietSeurat. (orig.ident-level) 
## (same procedues from 0_1, but analyse again just to ensure the consistent processing)
# split the RNA assay by patient → this creates one “layer” per patient 
GBO[["RNA"]] <- split(GBO[["RNA"]], f = GBO$orig.ident)
table(GBO$orig.ident)
# run standard RNA workflow on the GBOect 
GBO <- NormalizeData(GBO)
GBO <- FindVariableFeatures(GBO)
GBO <- ScaleData(GBO)
GBO <- RunPCA(GBO)
GBO <- FindNeighbors(GBO, dims = 1:30, reduction = "pca")
GBO <- FindClusters(GBO, resolution = 2) # 0.3 for broad, 2 for one population only   
GBO <- RunUMAP(GBO, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(GBO, label = TRUE) # malignant cluster regardless orig.ident 
DimPlot(GBO, reduction = "umap.unintegrated", group.by = c("orig.ident")) # malignant cluster by orig.ident separately
saveRDS(GBO, "scRNA_seq/GBO/GBO_malignant_qc_cellstate_hypsig_mergealias_region.rds")

# Extract the Seurat object of gene of interest ---------------------------------------------
alias_OGDD <- read_excel("markers/uniprotkb_GO_0016706_AND_reviewed.xlsx") 
alias_map_OGDD <- alias_OGDD %>% 
  filter(!`Entry Name` %in% c("PDIA1_HUMAN", "HBAP1_HUMAN")) %>%  # remove non-OGDD
  separate_rows(`Gene Names`, sep = " |,|;|/") %>%  
  mutate(`alias` = toupper(trimws(`Gene Names`))) %>%  
  group_by(Entry) %>%
  mutate(genes = first(alias)) %>% # Use first alias as the main gene 
  ungroup() %>% 
  select(alias, genes) %>%
  distinct()  # remove duplicates, 175 

alias_NOX <- read_excel("markers/uniprotkb_NADPH_Oxidase_AND_reviewed_tr_2026_01_10.xlsx")
alias_map_NOX <- alias_NOX %>% 
  filter(
    grepl("^NOX[1-5]", `Gene Names`) |
      grepl("^DUOX[12]", `Gene Names`) |
      grepl("^CYBB", `Gene Names`)) %>% 
  separate_rows(`Gene Names`, sep = " |,|;|/") %>%  
  mutate(`alias` = toupper(trimws(`Gene Names`))) %>%  
  group_by(Entry) %>%
  mutate(genes = first(alias)) %>% # Use first alias as the main gene 
  ungroup() %>% 
  select(alias, genes) %>%
  distinct()  # remove duplicates, 17

gene_of_interest <- unique(c(alias_map_OGDD$genes, alias_map_NOX$genes, "ACTB")) # 66 

GBO_slim <- DietSeurat(GBO, features = gene_of_interest)
table(GBO_slim$cellstate) 
table(GBO_slim$orig.ident)
summary(GBO_slim$hyp_sig)
table(GBO_slim$tumour_region)

saveRDS(GBO_slim, "scRNA_seq/slim/GBO_slim.rds") 

# after making all slim dataset, filter  remove recurrent gbm, orig.ident with low cells, add hyp_bin \(0_0)/
