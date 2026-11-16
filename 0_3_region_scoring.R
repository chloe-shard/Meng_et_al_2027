# region scoring is annotated on malignant populations only for detecting invasion pattern. 
# malignant cells were quality controled and annotated with Neftel state and hypoxia signature scores on merged alias based on raw counts. 

# Change GNS1 to 

setwd("I:/Data/Aijia/TAOF_R") 
library(Seurat)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(UCell)

# do region scoring per cell 
Puchlaski <- read_excel("data/markers/Puchlaski2018.xlsx")
Ishika <- read.csv("data/markers/Ishika_IT.csv")
region_marker <- Puchlaski %>%
  group_by(ROI) %>% 
  summarise(genes = list(unique(gene)), .groups = "drop") %>% 
  deframe() 
IT <- Ishika %>%
  pull(Infiltrating_region) %>%
  na.omit() %>%
  unique() %>%
  as.character()
region_marker <- c(region_marker, list(IT = IT))
names(region_marker)[names(region_marker) == "CTpan"] <- "PAN" 
names(region_marker)[names(region_marker) == "CTmvp"] <- "MVP" 


# Score cells -------------------------------------------------------------
# Define function to score cells
score_region_marker <- function(input.seurat){
  seurat.counts <- GetAssayData(input.seurat, assay = "RNA", layer = "counts")
  
  scores <- ScoreSignatures_UCell(
    seurat.counts,
    features = region_marker,
    name = "",
    maxRank = max(lengths(region_marker))
  )
  
  input.seurat <- AddMetaData(input.seurat, scores)
  
  # Add dominant putative region 
  input.seurat$tumour_region <- names(region_marker)[apply(scores[, names(region_marker)], 1, which.max)] 
  return(input.seurat)
}

# Apply the function to the Seurat objects
GNS1 <- score_region_marker(GNS1) 
table(GNS1$tumour_region) 

GNS1 <- saveRDS(GNS1, "data/scRNA_seq/GNS1/GNS1_malignant_qc_cellstate_hypsig_mergealias_region.rds")


