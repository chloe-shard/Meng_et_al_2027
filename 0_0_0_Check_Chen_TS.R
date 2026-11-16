# Check the data integrity 
library(Seurat)
library(ggplot2)
setwd("I:/Data/Aijia/TAOF_R/data")

table(GBO$cellstate)
# AC MES1 MES2 None NPC1 NPC2  OPC 
# 1338 1308 1378 1067 1161 1313 1199 

table(GBO$cellstate) 
# AC MES1 MES2 None NPC1 NPC2  OPC 
# 1388 1304  621  293  908 1043 1306 

summary(GBO$NPC2)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002565 0.604451 0.883844 1.051249 1.281871 6.248034 

summary(GBO$NPC2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.2578  0.7625  1.1200  1.5211 13.4800

summary(GBO$NPC1)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004545 0.683876 1.017300 1.084655 1.420540 3.838741 

summary(GBO$NPC1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3061  0.8141  1.1463  1.5867 11.1214 


summary(GBO$OPC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01445 0.65634 0.98836 1.08523 1.40181 4.88869 

summary(GBO$OPC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3564  0.7583  1.1590  1.4540 13.6042 

GBO$dataset <- "GBO"
GBO$dataset <- "GBO"
GBO$dataset <- "GBO"
MDAG$dataset <- "MDAG"
obj <- merge(GBO, y = list(GBO, GBO, MDAG), add.cell.ids = c("GBO", "GBO", "GBO", "MDAG"), project = "Check")


# 1) Normalize
obj <- NormalizeData(obj)

# 2) Find variable genes
obj <- FindVariableFeatures(obj)

# 3) Scale
obj <- ScaleData(obj)

# 4) PCA
obj <- RunPCA(obj, npcs = 30)

# 5) UMAP
obj <- RunUMAP(obj, dims = 1:30)
saveRDS(obj, "scRNA_seq/webtool/check_obj.rds")


# QC plots
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "dataset", pt.size = 0.1, ncol = 3) 

markers <- list(
  NPC2  = c("NDRG4", "SOX11"), 
  MES1 = c("CHI3L1","ANXA2"),
  MES2 = c("ATF3","PGK1"),
  AC   = c("GFAP","AQP4"),
  OPC  = c("APOD","OLIG1"),
  NPC1  = c("DLL3", "SOX4")
) 
for(state in names(markers)){
  genes <- markers[[state]]   
  genes <- genes[genes %in% rownames(obj)]
  
  p <- FeaturePlot(obj, features = genes, ncol = 3, pt.size = 0.1, split.by = "dataset") + 
    ggtitle(genes)
  print(p)
}


# all markers -------- 
library(readxl) 
library(tibble)

# state_markers <- read.table("markers/neftel_markers.txt", header = TRUE, row.names = NULL, sep = "\t")
# cellstate <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2")
# state_markers <- as.data.frame(state_markers) %>%
#   select(all_of(cellstate)) %>%
#   lapply(function(x) {
#     x <- x[ nzchar(x) & !is.na(x) ]  # remove empty + NA
#     unique(x)                       # remove duplicates
#   })
# 
# Chloe <- read.csv("markers/Hypoxia_signatures_from_Chloe.csv") 
# Puchlaski <- read_excel("markers/Puchlaski2018.xlsx")
# Ishika <- read.csv("markers/Ishika_IT.csv") 
# 
# hyp_sig <- Chloe$HYPOXIA
# region_marker <- Puchlaski %>%
#   group_by(ROI) %>% 
#   summarise(genes = list(unique(gene)), .groups = "drop") %>% 
#   deframe() 
# IT <- Ishika %>%
#   pull(Infiltrating_region) %>%
#   na.omit() %>%
#   unique() %>%
#   as.character() 
# markers <- c(state_markers, region_marker, list(IT = IT), list(hyp_sig = hyp_sig))
# saveRDS(markers, "markers/all_markers.rds")

markers <- readRDS("markers/all_markers.rds")
# test the existence of signatures expressing in each dataset ------------- 
GBO <- readRDS("scRNA_seq/GBO/GBO_malignant_qc_cellstate_hypsig_mergealias_region.rds") 



library(Matrix)
# GBO <- JoinLayers(GBO)
counts <- GetAssayData(GBO[["RNA"]], layer = "counts")
total_cells <- ncol(counts)
result_list <- lapply(names(markers), function(sig_name){
  
  genes <- markers[[sig_name]]
  genes_use <- intersect(genes, rownames(counts))
  
  mat <- counts[genes_use, , drop = FALSE]
  
  cells_express <- Matrix::rowSums(mat > 0)
  cells_not_express <- total_cells - cells_express
  
  data.frame(
    signature = sig_name,
    gene = genes_use,
    cells_express = cells_express,
    cells_not_express = cells_not_express,
    frac_express = cells_express / total_cells
  )
})
result_df <- do.call(rbind, result_list)
write.csv(result_df,
          file = "markers/analysis/GBO_markers_expression_summary.csv",
          row.names = FALSE)

# Visualise the markers not expressing ---------
library(dplyr)
path <- "markers/analysis"
files <- list.files(path, pattern = "*.csv", full.names = TRUE)
all_data <- lapply(files, read.csv) %>%
  bind_rows(.id = "dataset_id")

dataset_names <- gsub("_markers_expression_summary.csv", "",
                      basename(files))

all_data$dataset <- dataset_names[as.numeric(all_data$dataset_id)]
zero_expr <- all_data %>%
  filter(frac_express == 0)
head(zero_expr)


zero_summary <- zero_expr %>%
  group_by(dataset) %>%
  summarise(n_failed_genes = n())

zero_summary 

library(ggplot2) 

ggplot(zero_summary,
       aes(x = dataset, y = n_failed_genes)) +
  geom_col() +
  theme_classic() +
  labs(title = "Signatures with ZERO expression per dataset",
       x = "Dataset",
       y = "Number of genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

zero_summary_sig <- zero_expr %>%
  group_by(dataset, signature) %>%
  summarise(n_failed_genes = n(), .groups = "drop")

head(zero_summary_sig)

sig_colors <- c(
  "AC"="#E41A1C",
  "CT"="#377EB8",
  "CTmvp"="#4DAF4A",
  "CTpan"="#984EA3",
  "hyp_sig"="#FF7F00",
  "IT"="#A65628",
  "LE"="#F781BF",
  "MES2"="#999999",
  "NPC1"="#66C2A5",
  "NPC2"="#FC8D62",
  "OPC"="#8DA0CB"
)

ggplot(zero_summary_sig,
       aes(x = dataset,
           y = n_failed_genes,
           fill = signature)) +
  geom_col() +
  scale_fill_manual(values = sig_colors) +
  theme_classic() +
  labs(title = "Zero-expression genes per dataset (stacked by signature)",
       x = "Dataset",
       y = "Number of genes",
       fill = "Signature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



