setwd("I:/Data/Aijia/TAOF_R/data") 

# Change LeBlanc to 

# Calculate LeBlanc  State score for datasets

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(readr)
library(dplyr)
library(UCell)
library(reshape2) 

# Check orig.ident are samples, and celltype only malignant for cellstate annotation. 
table(LeBlanc$orig.ident)
table(LeBlanc$celltype)

LeBlanc <- JoinLayers(LeBlanc, assay = "RNA")
LeBlanc.mal.df <- GetAssayData(object = LeBlanc[["RNA"]], layer = "data")

# 2. Read LeBlanc  subgroup marker genes ---------------------------------------
gbm.genes.tab <- read.table("markers/neftel_markers.txt", header = TRUE, row.names = NULL, sep = "\t")
gbm.genes.tab <- gbm.genes.tab[,!(colnames(gbm.genes.tab) %in% c("G1_S", "G2_M"))]
gbm.genes.lst <- melt(gbm.genes.tab, measure.vars = colnames(gbm.genes.tab))
gbm.genes.lst <- gbm.genes.lst[gbm.genes.lst$value != "",]
gbm.genes.uniq <- unique(gbm.genes.lst$value)
table(gbm.genes.lst$variable)

# 3. after set LeBlanc.mal.df of seurat LeBlancect then come back to continue----------------------------------
genes.aggregate.expr <- sort(rowMeans(LeBlanc.mal.df))
genes.aggregate.expr <- genes.aggregate.expr[genes.aggregate.expr > 0]
LeBlanc.mal.df <- LeBlanc.mal.df[row.names(LeBlanc.mal.df) %in% names(genes.aggregate.expr), ]
gbm_states <- unique(gbm.genes.lst$variable)
state_scores <- matrix(NA, nrow = ncol(LeBlanc.mal.df), ncol = length(gbm_states),
                       dimnames = list(colnames(LeBlanc.mal.df), gbm_states))
for (g in gbm_states) {
  cat(g, "\n")
  
  genes <- gbm.genes.lst[gbm.genes.lst$variable == g, ]$value
  
  LeBlanc.mal.state <- LeBlanc.mal.df[row.names(LeBlanc.mal.df) %in% genes, , drop = FALSE]
  LeBlanc.mal.nonstate <- LeBlanc.mal.df[!row.names(LeBlanc.mal.df) %in% genes, , drop = FALSE]
  
  mean_exp_tumor_state <- rowMeans(LeBlanc.mal.state)
  LeBlanc.mal.state.rel <- sweep(LeBlanc.mal.state, 1, mean_exp_tumor_state, "/")
  mean_rel_tumor_state <- colMeans(LeBlanc.mal.state.rel, na.rm = TRUE)
  
  mean_exp_tumor_nonstate <- rowMeans(LeBlanc.mal.nonstate)
  LeBlanc.mal.nonstate.rel <- sweep(LeBlanc.mal.nonstate, 1, mean_exp_tumor_nonstate, "/")
  mean_rel_tumor_nonstate <- colMeans(LeBlanc.mal.nonstate.rel, na.rm = TRUE)
  
  state.score <- mean_rel_tumor_state / mean_rel_tumor_nonstate
  state_scores[, g] <- state.score
}
state_scores_df <- as.data.frame(state_scores) # Add all state scores to the LeBlancect as metadata 


# 4. Add all state scores to the LeBlancect as metadata   -----------------------------------
LeBlanc  <- AddMetaData(LeBlanc , metadata = state_scores_df)

# 5 add dominant cell state based on max state score (if max > 1)   -----------------------------------
cellstate <- apply(state_scores_df, 1, function(x) ifelse(max(x, na.rm = TRUE) > 1, names(which.max(x)), "None"))
LeBlanc  <- AddMetaData(LeBlanc , metadata = cellstate, col.name = "cellstate")

# 6. Create a new metadata column 'cellstate_merge'   -----------------------------------
LeBlanc$cellstate_merge <- LeBlanc$cellstate
LeBlanc$cellstate_merge[LeBlanc$cellstate_merge %in% c("NPC1", "NPC2")] <- "NPC"
LeBlanc$cellstate_merge[LeBlanc$cellstate_merge %in% c("MES1", "MES2")] <- "MES"
table(LeBlanc$cellstate) 

# !!! Note: This needs to process on merged alias count. 
# run hyp_scores before normalisation based on new counts ------------------------
setwd("H:/Data/Aijia/TAOF_R/data")
hyp_df <- read.csv("markers/Hypoxia_signatures_from_Chloe.csv", stringsAsFactors = FALSE)
hyp_vec <- hyp_df$HYPOXIA
hyp_vec <- hyp_vec[!is.na(hyp_vec) & hyp_vec != ""]
hyp_sig <- list(hyp_sig = hyp_vec)


hyp_scores <- function(input.seurat){
  seurat.counts <- GetAssayData(input.seurat, assay="RNA", layer="counts")
  scores <- ScoreSignatures_UCell(seurat.counts, features= hyp_sig , name="")
  input.seurat <- AddMetaData(input.seurat, scores)
  return(input.seurat) 
}

LeBlanc <- hyp_scores(LeBlanc) 
summary(LeBlanc$hyp_sig) 


