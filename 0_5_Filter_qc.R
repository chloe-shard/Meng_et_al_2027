setwd("I:/Data/Aijia/TAOF_R/data") 
# filter orig.ident with low cells, merge Abdelfattah dataset, remove recurrent gbm, finally add hyp_bin \(0_0)/
library(purrr)
library(dplyr)

# Process select orig.ident >100 cells 
# Load datasets 
Abdelfattah <- readRDS("scRNA_seq/slim/Abdelfattah_slim.rds")
Chen <- readRDS("scRNA_seq/slim/Chen_slim.rds")
Couturier <- readRDS("scRNA_seq/slim/Couturier_slim.rds")
Ebert <- readRDS("scRNA_seq/slim/Ebert_slim.rds")
LeBlanc <- readRDS("scRNA_seq/slim/LeBlanc_slim.rds")
Neftel <- readRDS("scRNA_seq/slim/Neftel_slim.rds")
TS <- readRDS("scRNA_seq/slim/TS_slim.rds")
GBO <- readRDS("scRNA_seq/slim/GBO_slim.rds")

# Function to remove orig.ident with total cell count < 100
filter_seurat_low_cells <- function(seurat_obj) {
  cell_counts <- table(seurat_obj$orig.ident)
  orig_qc <- names(cell_counts[cell_counts >= 100])
  subset(seurat_obj, subset = orig.ident %in% orig_qc & cellstate != "None"
  )
} 

# Apply to all objects
seurat_list <- list(
  Abdelfattah = filter_seurat_low_cells(Abdelfattah),
  Chen = filter_seurat_low_cells(Chen),
  Couturier = filter_seurat_low_cells(Couturier),
  Ebert = filter_seurat_low_cells(Ebert),
  LeBlanc = filter_seurat_low_cells(LeBlanc),
  Neftel = filter_seurat_low_cells(Neftel),
  GBO = filter_seurat_low_cells(GBO),
  TS = filter_seurat_low_cells(TS)
) 

table(seurat_list$Abdelfattah$orig.ident) # Check no CNSTM-081 with only 18 cells 
table(seurat_list$Neftel$orig.ident) # Check no BT1187 with only 8 cells 

# only keep primary glioblastoma --------------------
table(seurat_list$Abdelfattah$Type)
# GBM Recurrent GBM 
# 16950          9545

table(seurat_list$Neftel$GBMType)
# Adult Pediatric 
# 4734      1830 


# Abdelfattah: keep only non-recurrent GBM
seurat_list$Abdelfattah <- subset(
  seurat_list$Abdelfattah,
  subset = Type != "Recurrent GBM"
)

# Neftel: keep only Adult
seurat_list$Neftel <- subset(
  seurat_list$Neftel,
  subset = GBMType != "Pediatric"
) 

# make 40 hypoxia bins based on hyp_sig per dataset ----------- 
add_hyp_bin <- function(seurat_obj, n_bins = 40) {
  seurat_obj$hyp_bin <- ntile(seurat_obj$hyp_sig, n_bins)
  seurat_obj$hyp_bin <- factor(
    seurat_obj$hyp_bin,
    levels = sort(unique(seurat_obj$hyp_bin))
  )
  seurat_obj
}

seurat_list <- map(seurat_list, add_hyp_bin)
dataset_order <- c("Abdelfattah", "Chen", "Couturier", "Ebert", "LeBlanc", "Neftel", "TS", "GBO")
seurat_list <- seurat_list[order(factor(names(seurat_list), levels = dataset_order))]
saveRDS(seurat_list, "scRNA_seq/slim/seurat_list.rds")

# Join layer() 
# seurat_list <- lapply(seurat_list, function(seu) {
#   seu <- JoinLayers(
#     object = seu,
#     assay  = "RNA",
#     layers = "data"
#   )
#   seu
# })
# 
# # extract for expr_mat
# expr_list <- lapply(seurat_list, function(seu) {
#   GetAssayData(
#     object = seu,
#     assay  = "RNA",
#     layer  = "data"
#   )
# })  

