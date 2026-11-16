setwd("I:/Data/Aijia/TAOF_R/data") 
library(Seurat)
library(readxl)
library(tidyr)
library(dplyr)

alias <- read_excel("markers/uniprotkb_GO_0016706_AND_reviewed.xlsx")
head(alias)

alias_map <- alias %>% 
  separate_rows(`Gene Names`, sep = " |,|;|/") %>%  
  mutate(`alias` = toupper(trimws(`Gene Names`))) %>%  
  group_by(Entry) %>%
  mutate(genes = first(alias)) %>% # Use first alias as the main gene 
  ungroup() %>% 
  select(alias, genes) %>%
  distinct()  # remove duplicates
gene_of_interest <- alias_map$alias # 182
gene <- unique(alias_map$genes) # 60

# Merge alias gene in counts  ---------------------------------------
Neftel <- JoinLayers(Neftel, assay = "RNA") # Join if previously processed for orig.ident-based normalisation
counts <- GetAssayData(Neftel, assay = "RNA", layer = "counts") # extract count layer 
valid_map <- alias_map %>% filter(alias %in% rownames(counts)) # aliases that are present in the seurat Neftelect 

for (i in seq_len(nrow(valid_map))) {
  alias <- valid_map$alias[i]
  gene  <- valid_map$genes[i]
  
  if (alias == gene) next # skip if alias is same as gene
  
  if (alias %in% rownames(counts) && gene %in% rownames(counts)) {
    counts[gene, ] <- counts[gene, ] + counts[alias, ] # If alias and gene both exist, sum them. 
  }
}


alias_to_remove <- setdiff(valid_map$alias, valid_map$genes)
counts <- counts[!rownames(counts) %in% alias_to_remove, ]

# Pad any missing genes with 0s (from original)
original_genes <- rownames(GetAssayData(Neftel, assay = "RNA", layer = "counts"))

missing_genes <- setdiff(original_genes, rownames(counts)) 

if (length(missing_genes) > 0) {
  zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(counts),
                     dimnames = list(missing_genes, colnames(counts)))
  counts <- rbind(counts, zero_mat)
}

counts <- counts[original_genes, ] # reorder to match Seurat’s expectations 
Neftel <- SetAssayData(Neftel, assay = "RNA", layer = "counts", new.data = counts) # replace the count layer with the updated counts

rm(alias_map, counts, valid_map, zero_mat)

# Normalize / clutering / annotation / filter for only malignant / save ------------
## Run this option if you have multiple sample/patients (helps to control batch effect)
# split the RNA assay by patient → this creates one “layer” per patient
Neftel[["RNA"]] <- split(Neftel[["RNA"]], f = Neftel$orig.ident)

# run standard RNA workflow on the Neftelect 
Neftel <- NormalizeData(Neftel)
Neftel <- FindVariableFeatures(Neftel)
Neftel <- ScaleData(Neftel)
Neftel <- RunPCA(Neftel)
Neftel <- FindNeighbors(Neftel, dims = 1:30, reduction = "pca")
Neftel <- FindClusters(Neftel, resolution = 0.3) # 0.3 for broad, 2 for one population only 
Neftel <- RunUMAP(Neftel, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated") 

DimPlot(Neftel, label = TRUE) #  cluster regardless orig.ident 
DimPlot(Neftel, reduction = "umap.unintegrated", group.by = c("orig.ident")) 

saveRDS(Neftel, "scRNA_seq/webtool/Neftel_qc_mergealis.rds") 



