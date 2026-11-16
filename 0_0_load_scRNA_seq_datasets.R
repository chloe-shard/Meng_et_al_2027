
setwd("I:/Data/Aijia/TAOF_R/data")
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(readr)
library(dplyr)
library(UCell)
library(reshape2) 
library(readxl)
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
gene_of_interest <- unique(alias_map$genes) # 60
gene_of_interest <- unique(c(alias_map$genes, "ACTB")) # 61, ACTB as housekeeping gene

# Abdelfattah ------------------------------------------------------------------
Abdelfattah [["percent.mt"]] <- PercentageFeatureSet(Abdelfattah , pattern = "^MT-") # run this if percent.mt not found. 
Abdelfattah  <- subset(Abdelfattah , subset = nFeature_RNA > 200 & percent.mt < 10) 

Idents(Abdelfattah) <- "Assignment"
Abdelfattah  <- subset(Abdelfattah , subset = Assignment == "Glioma") 
Idents(Abdelfattah) <- "orig.ident"
table(Abdelfattah$orig.ident)

CNSTM_070 <- subset(Abdelfattah, orig.ident %in% c("CNSTM-070"))
CNSTM_390 <- subset(Abdelfattah, orig.ident %in% c("CNSTM-390"))

# Abdelfattah_CNSTM ------------------------ done -----------------------------
CNSTM <- readRDS("data/scRNA_seq/webtool/Abdelfattah_CNSTM_malignant_qc_cellstate_hypsig_mergealias.rds")
CNSTM <- subset(CNSTM, orig.ident %in% c("CNSTM-068", "CNSTM-081", "CNSTM-096", "CNSTM-375", "CNSTM-379", "CNSTM-397"))
table(CNSTM$orig.ident) 
CNSTM <- subset(CNSTM, subset = orig.ident != "CNSTM-081") # only has 18 malignant cells 
table(CNSTM$Assignment) # Glioma 16595 -> 16577 


# Abdelfattah_MDAG ------------------------------- done ----------------------
MDAG <- readRDS("data/scRNA_seq/webtool/Abdelfattah_MDAG_malignant_qc_cellstate_hypsig_mergealias.rds")
MDAG <- subset(Abdelfattah, orig.ident %in% c("MDAG-1", "MDAG-10", "MDAG-12", "MDAG-4", "MDAG-6", "MDAG-7", "MDAG-9"))
table(MDAG$orig.ident) 
MDAG <- subset(MDAG, orig.ident != "MDAG-1")
table(MDAG$Assignment) # Glioma 12690 


# Chen ------------------------------------------------ done ------------------
Chen <- readRDS("data/scRNA_seq/webtool/Chen_malignant_qc_cellstate_hypsig_mergealias.rds")

# Chen$sample <- Chen$orig.ident
table(Chen$sample)
# Chen$orig.ident <- sub("-.*", "", Chen$sample)
table(Chen$orig.ident) 
table(Chen$cellkb) # Tumor 22171


# Couturier ------------------------------------------------------------------ done -----------------------------
Couturier <- readRDS("data/scRNA_seq/webtool/Couturier_malignant_qc_cellstate_hypsig_mergealias.rds")

# Couturier$sample <- Couturier$orig.ident
# Couturier$orig.ident <- sub("^(BT\\d+).*", "\\1", Couturier$orig.ident)
table(Couturier$orig.ident)
table(Couturier$cellkb) # Tumor 29874


# Ebert_malignant  ----------------------------- done  -----------------
Ebert <- readRDS("data/scRNA_seq/webtool/Ebert_malignant_qc_cellstate_hypsig_mergealias.rds")
table(Ebert$orig.ident)
table(Ebert$pooled_annot) # Tumor cell 731, Tumor cell/OPC 640, Tumor cell/APC 416
table(Ebert$cellstate)

# LeBlanc --------------------------------------- done -----------------
LeBlanc <- readRDS("data/scRNA_seq/webtool/LeBlanc_malignant_qc_cellstate_hypsig_mergealias.rds")
summary(LeBlanc$percent.mt)
# LeBlanc$orig.ident <- sub("^((JK\\d{3})).*", "\\1", LeBlanc$orig.ident)
table(LeBlanc$orig.ident) 
table(LeBlanc$orig.ident) 
table(LeBlanc$cell_type) # malignant 23545


# Neftel --------------------------- done  -----------------
Neftel <- readRDS("scRNA_seq/webtool/Neftel_qc_mergealias_cellstate_hypsig.rds")
Neftel [["percent.mt"]] <- PercentageFeatureSet(Neftel , pattern = "^MT-") # run this if percent.mt not found. 
Neftel  <- subset(Neftel , subset = nFeature_RNA > 200 & percent.mt < 10) 
Neftel$orig.ident <- Neftel$Sample
Neftel$Sample <- NULL
table(Neftel$orig.ident)
table(Neftel$CellAssignment) # Malignant 6863 
Neftel  <- subset(Neftel , subset = CellAssignment == "Malignant") 

Neftel  <- subset(Neftel , subset = GBMType == "Adult") 
saveRDS(Neftel, "scRNA_seq/webtool/Neftel_qc_mergealias_cellstate_hypsig.rds")


# GBO and TS are in the same dataset --------------------------- done -----------------
GBO <- readRDS("data/scRNA_seq/GBO/GBO_TS.rds")

GBO$orig.ident <- gsub("^Ana_", "", GBO$sample)
GBO$orig.ident <- ifelse(
  grepl("_[0-9]+$", GBO$orig.ident) & lengths(regmatches(GBO$orig.ident, gregexpr("_", GBO$orig.ident))) > 1,
  gsub("_[0-9]+$", "", GBO$orig.ident),
  GBO$orig.ident
)
table(GBO$orig.ident)
table(GBO$broad) # Tumor
summary(GBO$percent.mt)

# more GBO, more TS 

# GBO ------------------------------------------------------------- 
GBO <- readRDS("data/scRNA_seq/GBO/GBO_07_12_41_malignant_qc_cellstate_hypsig_mergealias.rds")
GBO  <- subset(GBO , subset = nFeature_RNA > 200 & percent.mt < 10) 
GBO <- subset(GBO, broad == "Tumor" & orig.ident %in% c(" GBO_03", "GBO_07", 
                                                        "GBO_12", "GBO_27", 
                                                        "GBO_32_R3", "GBO_36",
                                                        "GBO_37_W0_1", 
                                                        "GBO_37_W1", "GBO_38", "GBO_41", "GBO_43")) 
table(GBO$orig.ident)
table(GBO$broad)  # Tumor 40301 

GBO$sample <- GBO$orig.ident
GBO$orig.ident[GBO$orig.ident %in% c("Ana_GBO_41_1", "Ana_GBO_41_2")] <- "Ana_GBO_41"
table(GBO$orig.ident)
saveRDS(GBO, "data/GBO_slim.rds")

# TS ------------------------------------------------------------- 
TS <- readRDS("data/scRNA_seq/TS/TS_12_41_malignant_qc_cellstate_hypsig_mergealias.rds")
TS <- subset(GBO, broad == "Tumor" & orig.ident %in% c("Nga_BT102_SANTB00706", "Nga_BT103_SANTB00707", 
                                                       "Nga_BT107_SANTB00729", "Nga_BT108_SANTB00486", 
                                                       "Nga_BT48_SANTB00469", "Nga_BT64_SANTB00554", 
                                                       "Nga_BT90_SANTB00682", "Nga_BT95_SANTB00691"))
TS  <- subset(TS, subset = nFeature_RNA > 200 & percent.mt < 10) 
table(TS$orig.ident)
table(TS$broad) # Tumor 8764
saveRDS(TS, "data/scRNA_seq/TS/TS_malignant_qc.rds")

TS <- readRDS("data/scRNA_seq/TS/TS_12_41_malignant_qc_cellstate_hypsig_mergealias.rds")

# Cell line ----------------------------------------------------
obj <- readRDS("G:/TAOF/Raw_seq_data/Processed_rds_files/Single_cell/Parse_scRNAseq_2024/all_samples_tumour_GBO_cell_lines.rds")
# metadata <- read_excel("G:/TAOF/Raw_seq_data/Processed_rds_files/Single_cell/Parse_scRNAseq_2024/Parse_metadata.xlsx", sheet = 2)
# 
# samples <- metadata %>%
#   filter(Diagnosis == "Glioblastoma", 
#          Tumour_occurrence == "Primary", 
#          Sample_type = "GNS_Line") %>%
#   pull(Sample_ID) # 21 
#  [1] "Sakthi_G6T P8"                "Sakthi_G12T P8"               "Sakthi_G15T"                 
#  [4] "Sakthi_G4T RFP P22"           "Sakthi_G4T RFP CD47 KO P25"   "Sakthi_G4T RFP NTC P24"      
#  [7] "Sakthi_G5C RFP#1 P23"         "Sakthi_G5C RFP#1 NTC P35"     "Sakthi_G5C RFP#1 CD47 KO P36"
# [10] "Sakthi_G20T P8"               "Sakthi_G23T P9"               "Sakthi_G27T P12"             
# [13] "Sakthi_G33T P9"               "Sakthi_G35T P8"               "Sakthi_G40T P6"              
# [16] "Sakthi_G47T p9"               "Sakthi_G52T P7"               "Sakthi_G53T P8"              
# [19] "Sakthi_G54T"                  "Sakthi_G55T P9"               "Sakthi_G57T P8"
# take only primary glioblastoma cell line data   
obj <- subset(obj, sample %in% c("Sakthi_G12T_P8", "Sakthi_G15T", "Sakthi_G20T_P8", "Sakthi_G23T_P9", "Sakthi_G27T_P12",
                                     "Sakthi_G33T_P9", "Sakthi_G35T_P8", "Sakthi_G40T_P6", "Sakthi_G47T_p9", 
                                 "Sakthi_G4T_RFP_NTC_P24", "Sakthi_G4T_RFP_P22", 
                                     "Sakthi_G52T_P7", "Sakthi_G53T_P8", "Sakthi_G54T_", "Sakthi_G55T_P9", "Sakthi_G57T_P8", 
                                 "Sakthi_G5C_RFP_1_NTC_P35", "Sakthi_G5C_RFP_1_P23", "Sakthi_G6T_P8")) 
obj  <- subset(obj , subset = nFeature_RNA > 200 & percent.mt < 10) 
obj$orig.ident <- sub("^Sakthi_([^_]+).*", "\\1", obj$sample)
sort(table(obj$orig.ident), decreasing = TRUE)

saveRDS(obj, "data/scRNA_seq/cell_line/all_cell_lines_qc.rds") 
obj <- readRDS("data/scRNA_seq/cell_line/all_cell_lines_qc.rds")

GNS1 <- subset(obj, orig.ident %in% c("G27T", "G35T", "G55T")) # three cell line growing 
saveRDS(GNS1, "data/scRNA_seq/cell_line/GNS1_qc.rds") # this is merged alias. 
rm(obj)
GNS2 <- subset(obj, orig.ident %in% c()) # cell lines matching with analysed TS samples 




