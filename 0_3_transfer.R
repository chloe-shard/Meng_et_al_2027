setwd("I:/Data/Aijia/TAOF_R/data") 

# merge the neftel state scores, state annotation, hyp_sig from malignant cell dataset to the total cel data
## from malignant_cellstate_hypsig to all_qc_mergealias

Abdelfattah_all <- readRDS("scRNA_seq/webtool/Abdelfattah_qc_mergealias.rds") 
Abdelfattah_mal <- readRDS("scRNA_seq/slim/Abdelfattah_slim.rds")

colnames(Abdelfattah_all@meta.data)
colnames(Abdelfattah_mal@meta.data)
sig_cols <- c("cellstate", "MES2", "MES1", "AC", "OPC", "NPC1", "NPC2", "cellstate_merge", "hyp_sig", 
              "LE", "IT", "CT", "PAN", "MVP", "tumour_region")


Abdelfattah_all@meta.data[, sig_cols] <- NA # create empry columns in all dataset 
common_cells <- intersect(
  rownames(Abdelfattah_all@meta.data),
  rownames(Abdelfattah_mal@meta.data)
)

Abdelfattah_all@meta.data[common_cells, sig_cols] <-
Abdelfattah_mal@meta.data[common_cells, sig_cols]
# Should be NA for non-malignant cells
table(is.na(Abdelfattah_all$cellstate))

# Spot check
head(Abdelfattah_all@meta.data[common_cells, sig_cols])
head(Abdelfattah_mal@meta.data[common_cells, sig_cols])

saveRDS(Abdelfattah_all, "scRNA_seq/webtool/Abdelfattah_qc_mergealias_cellstate_hypsig.rds")

