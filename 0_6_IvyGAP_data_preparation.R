# Load libraries
library(Seurat)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(UCell)

# Create Seurat Object: Ivy Gap GBM RNA --------------------------
# Set working directory
setwd("I:/Data/Aijia/TAOF_R/data")

# Read metadata and expression data
meta <- read_excel("bulk/all/GSE107559_ivygap_columns-samples.xlsx")
expr <- read.csv("bulk/all/GSE107559_ivygap_fpkm_table.csv.gz",
                 row.names = 1,
                 check.names = FALSE)
genes <- read.csv("bulk/all/GSE107559_ivygap_rows-genes.csv.gz",
                  stringsAsFactors = FALSE)

# Assign gene symbols as rownames
rownames(expr) <- genes$gene_symbol

# Reorder metadata to match expression matrix
meta <- meta[match(colnames(expr), meta$rna_well_id), ]

# Create Seurat object
obj <- CreateSeuratObject(
  counts = expr,
  meta.data = meta
)

# Create FPKM assay
obj[["FPKM"]] <- CreateAssayObject(counts = expr)
DefaultAssay(obj) <- "FPKM"

# Store log1p(FPKM) in data slot, because FPKM data is normalised but not log-transformed
log_expr <- log1p(as.matrix(expr))
obj <- SetAssayData(
  obj,
  assay = "FPKM",
  layer = "data",
  new.data = log_expr
)

# Clean metadata
obj$nCount_RNA <- NULL
obj$nFeature_RNA <- NULL

# Orig.ident is the patient id
obj$orig.ident <- obj$tumor_name
obj$tumour_region <- sub("-.*", "", obj$structure_abbreviation)
table(obj$tumour_region, obj$structure_color)
obj <- subset(obj, subset = !(orig.ident %in% c("W20-2-1","W22-2-1","W27-2-1")))
table(obj$orig.ident)
obj$orig.ident <- obj$rna_well_id

DefaultAssay(obj) <- "FPKM"
obj[["RNA"]] <- NULL


# Save object
saveRDS(obj, "bulk/all/GSE107559_RNA.rds")


# score for hyp_sig and region ------------------------
setwd("H:/Data/Aijia/TAOF_R_OGDD/data")
library(UCell)
obj <- readRDS("bulk/all/GSE107559_RNA.rds")
Chloe <- read.csv("markers/Hypoxia_signatures_from_Chloe.csv")
Puchlaski <- read_excel("markers/Puchlaski2018.xlsx")
Ishika <- read.csv("markers/Ishika_IT.csv")

hyp_sig <- Chloe$HYPOXIA
region_marker <- Puchlaski %>%
  group_by(ROI) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()
IT <- Ishika %>%
  pull(Infiltrating_region) %>%
  na.omit() %>%
  unique() %>%
  as.character()
markers <- c(region_marker, list(IT = IT), list(hyp_sig = hyp_sig))


# Run UCell
obj <- AddModuleScore_UCell(
  obj,
  features = markers,
  name = "",
  maxRank = max(lengths(markers))
)

# saveRDS(obj, "bulk/all/GSE107559_RNA.rds")
