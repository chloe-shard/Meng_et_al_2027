
save_dir_fig_4 <- "I:/Writing/02-Meng et al_2027/Figures/Fig 4"

Averaged_patient_proteomics_data <- read_csv("I:/Data/Chloe/Other Published datasets/Proteomics_Lam_et_al_2022/Parul_analysis/Averaged_patient_proteomics_data.csv")

df <- Averaged_patient_proteomics_data
rm(Averaged_patient_proteomics_data)

# Step 1: Extract gene names from the 'Fasta headers' column using a clear extraction method
df$Gene_Name <- sapply(df$`Fasta headers`, function(x) {
  # Extract the gene name using regular expression, which matches 'GN=' followed by the gene name
  match <- regmatches(x, regexpr("GN=[A-Za-z0-9]+", x))
  # Remove 'GN=' to keep only the gene name
  if(length(match) > 0) {
    gsub("GN=", "", match)  # Replace 'GN=' with an empty string to keep only the gene name
  } else {
    NA  # Return NA if no match is found
  }
})


# Step 2: Move the Gene_Name column to the first position to make it easier to work with
df <- df[, c("Gene_Name", setdiff(names(df), "Gene_Name"))]

# Step 3: Remove unnecessary columns, such as 'Protein IDs' and 'Fasta headers'
df <- df[, !colnames(df) %in% c("Protein IDs", "Fasta headers")]

# Step 4: Remove rows with missing or malformed Gene_Name values (if any)
df <- df[!is.na(df$Gene_Name) & df$Gene_Name != "", ]


# Step 5: Remove any existing row names to ensure a clean dataframe
rownames(df) <- NULL

# Step 1: Remove any existing names attribute from the 'Gene_Name' column
df$Gene_Name <- unname(df$Gene_Name)



############### Subset based on genes of interest ################ note threshold score >20 (anything less is removed from from analysis)
# Step 2: Subset the dataframe based on the genes of interest
# Check for different names 
# # protein_of_interest <- Purine_and_Pyrimidine_metabolism[["Gene"]]
# protein_of_interest <- c( # KDM7A is named as JHDM1D in Ivy Gap GBM bulk RNA seq data. 
#   "ACTB", # house-keeping 
#   "KDM4C", "JHDM1D", "JMJD6", "EGLN1", "EGLN3", 
#   "P4HA2", "ALKBH5", "PLOD2", "KDM3A", "P4HA1",  
#   "PLOD1", "PHYH", "ASPH", "P4HTM", "BBOX1", 
#   "PLOD3", "KDM5D", "ADO", "KDM8", "KDM2B", 
#   "PHF2", "KDM4A", "ALKBH4", "TET3", "DUOX1", 
#   "HR", "KDM3B")
# 
# 
# protein_of_interest <- c(
#   "ACTB", 
#   "KDM4C",  "JMJD2C",
#   "JHDM1D", "KDM7A",
#   "JMJD6",  "JMJD6",
#   "EGLN1",  "PHD2",
#   "EGLN3",  "PHD3", 
#   "P4HA2",  
#   "ALKBH5", 
#   "PLOD2",  "LH2",
#   "KDM3A",  "JMJD1A",
#   "P4HA1",  
#   "PLOD1",  "LH1",
#   "PHYH",   "PAHX",
#   "ASPH",   "BAH",
#   "P4HTM",  "PH-4", "P4H-TM", 
#   "BBOX1",  "BBOX",
#   "PLOD3",  "LH3",
#   "KDM5D",  "JARID1D",
#   "ADO",  
#   "KDM8",   "JMJD5",
#   "KDM2B",  "FBXL10",
#   "PHF2",   "KDM7C",
#   "KDM4A",  "JMJD2A",
#   "ALKBH4", 
#   "TET3", 
#   "DUOX1",  "THOX1",
#   "HR", "HAIR", "HAIRLESS", 
#   "KDM3B",  "JMJD1B"
# )
# Gene list  # KDM7A is named as JHDM1D in Ivy Gap GBM bulk RNA seq data.
gene_of_interest <- c("ACTB", "ALKBH4", "ADO", "KDM4A", "KDM3B", "BBOX1", 'PHF2', "KDM8", "HR", 'KDM5D',
                      "DUOX1", 'KDM2B', 'P4HTM', 'TET3', 'PHYH', 'PLOD3', "ASPH", "PLOD1", "EGLN1",
                      'KDM4C', 'P4HA1', 'PLOD2', 'KDM3A', "ALKBH5", 'P4HA2', 'EGLN3', "JMJD6", "JHDM1D")



# Ensure unique gene names
protein_of_interest <- unique(c("ACTB", gene_of_interest)) # 

# Subset df
df_subset <- df[df$Gene_Name %in% protein_of_interest, ]

# Optional: order Gene_Name as factor
df_subset$Gene_Name <- factor(
  df_subset$Gene_Name,
  levels = sort(unique(df_subset$Gene_Name))
)

# Preview
head(df_subset)


#################### Average data by region #############

# Define the column patterns for each tumour region
tumour_regions <- c("IT", "LE", "CT", "MVP", "PAN", "whole")

# Create a new dataframe to store the averaged values for each tumour region
# Include the Gene_Name column as the first column
df_avg <- data.frame(Gene_Name = df_subset$Gene_Name)

# Loop through each tumour region and calculate the average for each gene
for (region in tumour_regions) {
  # Get columns that match the current region pattern
  region_columns <- grep(paste0("Average ", region, "-Patient"), names(df_subset), value = TRUE)
  
  # Convert the region columns to numeric, in case they are not in the correct format
  region_data <- sapply(df_subset[, region_columns], as.numeric)
  
  # Calculate the row-wise mean for the region columns and add it to the new dataframe
  df_avg[[paste0("Average_", region)]] <- rowMeans(region_data, na.rm = TRUE)
}

# Reset the row names to numeric values
rownames(df_avg) <- NULL

# Display the first few rows of the new averaged dataframe
head(df_avg)

########### Identify outliers using the inter quartile method ################

# Define the column patterns for each tumour region
tumour_regions <- c("IT", "LE", "CT", "MVP", "PAN")

# Create a new list to store outlier information for each region
outliers_list <- list()

# Loop through each tumour region to identify outliers using the IQR method
for (region in tumour_regions) {
  # Get columns that match the current region pattern
  region_columns <- grep(paste0("Average ", region, "-Patient"), names(df_subset), value = TRUE)
  
  # Convert the region columns to numeric, in case they are not in the correct format
  region_data <- sapply(df_subset[, region_columns], as.numeric)
  
  # Calculate the 1st and 3rd quartiles for each gene
  Q1 <- apply(region_data, 1, quantile, probs = 0.25, na.rm = TRUE)
  Q3 <- apply(region_data, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # Calculate the IQR
  IQR <- Q3 - Q1
  
  # Define lower and upper bounds for outliers
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Identify outliers (TRUE if the value is less than lower bound or greater than upper bound)
  outlier_matrix <- (region_data < lower_bound) | (region_data > upper_bound)
  
  # Store the outlier matrix as a dataframe with Gene_Name as the first column
  outlier_df <- data.frame(Gene_Name = df_subset$Gene_Name, outlier_matrix)
  
  # Reset row names to numeric values
  rownames(outlier_df) <- NULL
  
  # Add the outlier dataframe to the list
  outliers_list[[region]] <- outlier_df
}

# Display the outlier data for one of the regions, e.g., CT
head(outliers_list$CT) 

######################## Average data by region (outliers removed) ##############

# Step 1: Remove the "whole" region from the outliers_list and tumour_regions
outliers_list <- outliers_list[!names(outliers_list) %in% "whole"]
tumour_regions <- tumour_regions[tumour_regions != "whole"]

# Step 2: Create a new dataframe to store the averaged values for each tumour region
df_avg <- data.frame(Gene_Name = df_subset$Gene_Name) 

# Step 3: Loop through each tumour region and calculate the average for each gene, excluding outliers
for (region in tumour_regions) {
  # Get columns that match the current region pattern
  region_columns <- grep(paste0("Average ", region, "-Patient"), names(df_subset), value = TRUE)
  
  # Convert the region columns to numeric in case they are not in the correct format
  region_data <- sapply(df_subset[, region_columns], as.numeric)
  
  # Retrieve outlier matrix for the current region from the outliers_list
  outlier_df <- outliers_list[[region]]
  
  # Ensure Gene_Name columns match before proceeding
  if (!identical(outlier_df$Gene_Name, df_subset$Gene_Name)) {
    stop("Gene names in outlier list do not match with the main dataframe.")
  }
  
  # Convert outlier matrix to logical matrix and use it to set outliers to NA in region_data
  outlier_matrix <- as.matrix(outlier_df[, -1])  # Exclude the 'Gene_Name' column
  region_data[outlier_matrix] <- NA  # Set outliers to NA
  
  # Calculate the row-wise mean for the region columns, excluding outliers (NA values)
  df_avg[[paste0("Average_", region)]] <- rowMeans(region_data, na.rm = TRUE)
  
  # ACTB-normalized average
  actb_mean <- df_avg[df_avg$Gene_Name == "ACTB", paste0("Average_", region)]
  df_avg[[paste0("Normalized_", region)]] <- df_avg[[paste0("Average_", region)]] / actb_mean
}

# Step 4: Display the final averaged dataframe
print("Final averaged dataframe:")
print(df_avg)

# write.csv(df_avg, "proteomics/IvyGapGBM/Lam2022_candidate_qc_averaged.csv", row.names=TRUE) 

# df_avg <- read.csv("proteomics/IvyGapGBM/Lam2022_candidate_qc_averaged.csv") 





################## ANOVA analysis (outliers removed) ###################

# Load necessary libraries for ANOVA and Tukey post hoc test --------------
library(multcomp)
library(dplyr)
library(tidyr)

# Define tumour regions to include in the analysis (excluding "whole")
tumour_regions <- c("IT", "LE", "CT", "MVP", "PAN")

# Initialize an empty list to store results
anova_results <- list()

# Loop through each gene in the dataframe
for (gene in unique(df_subset$Gene_Name)) {
  
  # Subset data for the current gene
  gene_data <- df_subset %>% filter(Gene_Name == gene)
  
  # Get the outliers for the current gene from the outliers_list
  # Extract outlier data for the current gene for each region
  gene_outliers <- outliers_list %>% lapply(function(outlier_df) {
    outlier_row <- outlier_df %>% filter(Gene_Name == gene)
    return(as.logical(outlier_row[1, -1]))  # Remove the Gene_Name column and keep logical values
  })
  
  # Create a long-format dataframe to organize data for ANOVA
  region_data <- data.frame()
  
  # Loop through each tumour region and collect the data
  for (region in tumour_regions) {
    # Select columns that match the current region pattern
    region_columns <- grep(paste0("Average ", region, "-Patient"), names(gene_data), value = TRUE)
    
    # Extract and reshape the data for the current region
    region_values <- gene_data[, region_columns]
    
    # Remove outliers based on the outlier matrix for the current gene and region
    if (!is.null(gene_outliers[[region]])) {
      region_values <- region_values[, !gene_outliers[[region]], drop = FALSE]
    }
    
    # Reshape the data into a long format
    region_long <- data.frame(
      Tumour_Region = rep(region, times = ncol(region_values)),  # Region label
      Value = as.numeric(unlist(region_values))  # Convert values to numeric and stack them
    )
    
    # Append to the overall region_data dataframe
    region_data <- rbind(region_data, region_long)
  }
  
  # Convert Tumour_Region to a factor for ANOVA
  region_data$Tumour_Region <- as.factor(region_data$Tumour_Region)
  
  # Perform one-way ANOVA if we have data for the gene
  if (nrow(region_data) > 0) {
    # Perform one-way ANOVA
    anova_result <- aov(Value ~ Tumour_Region, data = region_data)
    
    # Perform Tukey's post hoc test
    tukey_result <- glht(anova_result, linfct = mcp(Tumour_Region = "Tukey"))
    
    # Store results in the list
    anova_results[[gene]] <- list(
      anova = summary(anova_result),
      tukey = summary(tukey_result)
    )
  } else {
    message(paste("No data available for gene:", gene))
  }
}

# View an example result for a specific gene (replace 'Gene_Name' with a gene of interest)
anova_results[["ACTB"]]
# Define the genes you're interested in
genes_to_view <- unique(df_subset$Gene_Name)

# Loop through and print results
for (gene in genes_to_view) {
  cat("\n=====================\n")
  cat("Results for Protein:", gene, "\n")
  cat("=====================\n")
  
  if (!is.null(anova_results[[gene]])) {
    print(anova_results[[gene]]$anova)
    
    if (!is.null(anova_results[[gene]]$tukey)) {
      cat("\nTukey's HSD:\n")
      print(anova_results[[gene]]$tukey)
    } else {
      cat("\nNo Tukey result available.\n")
    }
  } else {
    cat("No ANOVA result found for this gene.\n")
  }
}





# ----------------------------
# Normalize protein abundance by ACTB
# ----------------------------
# ==========================================================

# REMOVE OUTLIERS FROM df_subset (FOR CLEAN PLOTTING)

# ==========================================================

df_clean <- df_subset

for (region in names(outliers_list)) {
  
  # Find columns for this tumour region
  
  region_columns <- grep(
    paste0("Average ", region, "-Patient"),
    names(df_clean),
    value = TRUE
  )
  
  # Ensure numeric
  
  df_clean[, region_columns] <- lapply(
    df_clean[, region_columns],
    as.numeric
  )
  
  # Get logical outlier matrix
  
  outlier_matrix <- as.matrix(outliers_list[[region]][, -1])
  
  # Replace outliers with NA
  
  df_clean[, region_columns][outlier_matrix] <- NA
}

# 1. Identify patient-region abundance columns
abundance_columns <- grep("Average .*?-Patient", names(df_subset), value = TRUE)

# 2. Extract ACTB values
actb_row <- df_subset[df_subset$Gene_Name == "ACTB", ]
if (nrow(actb_row) != 1) {
  stop("Error: ACTB gene not found or multiple entries present.")
}
actb_values <- as.numeric(actb_row[, abundance_columns])  # Ensure numeric vector

# 3. Normalize all proteins by ACTB per region/patient
df_norm <- df_clean
df_norm[, abundance_columns] <- lapply(df_norm[, abundance_columns], as.numeric)
df_norm[, abundance_columns] <- sweep(
  df_norm[, abundance_columns],
  2,
  actb_values,
  FUN = "/"
) 

# ----------------------------
# Reshape and filter for plotting
# ----------------------------

# 4. Melt the normalized data to long format
df_melt <- melt(df_norm,
                id.vars = "Gene_Name",
                measure.vars = abundance_columns,
                variable.name = "Tumour_Region_Patient",
                value.name = "Normalized_Abundance")

# 5. Extract Tumour_Region and clean up
df_melt$Tumour_Region <- sub("Average (.*?)-Patient.*", "\\1", df_melt$Tumour_Region_Patient)

# 6. Filter: exclude ACTB, exclude 'whole', and remove rows with 0 or NA abundance
df_melt <- subset(df_melt,
                  Gene_Name != "ACTB" &
                    !grepl("^whole$", Tumour_Region, ignore.case = TRUE) &
                    !is.na(Normalized_Abundance) &
                    Normalized_Abundance != 0)

# 7. Set desired Tumour_Region order
Tumour_Region_order <- c("LE", "IT", "CT", "MVP", "PAN")
df_melt$Tumour_Region <- factor(df_melt$Tumour_Region, levels = Tumour_Region_order)

protein_in_lam <- unique(as.character(df_melt$Gene_Name))


# ----------------------------
# Plot: Normalized protein abundance per Tumour_Region
# ----------------------------

fig_4b <- ggplot(df_melt, aes(x = Tumour_Region, y = Normalized_Abundance)) +
  geom_point(aes(color = Tumour_Region), alpha = 0.6, size = 1.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "black") +
  theme_minimal() +
  facet_wrap(~ Gene_Name, scales = "free_y", nrow = 1) +
  labs(
    title = "Candidate oxygen sensors in spatial proteomics of Lam et al.", 
    y = "Protein abundunce (/ACTB)",
    x = "Tumour Region"
  ) +
  theme(
    strip.text = element_text(face = "bold"), 
    axis.text.x = element_text(angle = 45),
    legend.position = "right"
  )
print(fig_4b)


# # Save
# ggsave(
#   filename = file.path(save_dir_fig_4, "fig_4b.pdf"),
#   plot = fig_4b,
#   device = grDevices::cairo_pdf,
#   width = 14,
#   height = 2.5,
#   units = "in",
#   bg = "transparent"
# )

##### Fig_4c bulk RNA expression in normal vs tumour tissue from GEPIA 3  #### 
gene_of_interest <- c("ALKBH4", "ADO", "KDM4A", "KDM3B", "BBOX1", 'PHF2', "KDM8", "HR", 'KDM5D', 
                      "DUOX1", 'KDM2B', 'P4HTM', 'TET3', 'PHYH', 'PLOD3', "ASPH", "PLOD1", "EGLN1", 
                      'KDM4C', 'P4HA1', 'PLOD2', 'KDM3A', "ALKBH5", 'P4HA2', 'EGLN3', "JMJD6", "JHDM1D")

df <- read_excel("I:/Data/Aijia/TAOF_R/data/bulk/GEPIA/GEPIA3_oxygen_sensors_20260505.xlsx")

# gene order run from 5_1_webtool_normal_vs_tumour.R 
gene_order <- rev(df$gene)

df_long <- df %>%
  pivot_longer(cols = c(tumour, normal),
               names_to = "condition",
               values_to = "expression") %>%
  mutate(gene = factor(gene, levels = gene_order))

fig_4c <- ggplot(df_long, aes(condition, gene, fill = expression)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(c("black", "white", "orangered1"))(100)) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"), 
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) + 
  labs(
    title = "GEPIA3: \nnormal vs tumor", 
    fill = "log2(TPM+1)", 
    x = NULL, 
    y = NULL) 
print(fig_4c)

# # Save
# ggsave(
#   filename = file.path(save_dir_fig_4, "fig_4c.pdf"),
#   plot = fig_4c,
#   device = grDevices::cairo_pdf,
#   width = 3,
#   height = 7.5,
#   units = "in",
#   bg = "transparent"
# )



# Keep only genes where tumour > normal
df_filtered <- df %>%
  filter(tumour > normal) 

# Optional: order genes by difference (tumour - normal)
df_filtered <- df_filtered %>%
  mutate(diff = tumour - normal) %>%
  arrange(diff)

gene_of_interest <- rev(df_filtered$gene)

# Convert to long format
df_long <- df_filtered %>%
  pivot_longer(cols = c(tumour, normal),
               names_to = "condition",
               values_to = "expression") %>%
  mutate(gene = factor(gene, levels = gene_of_interest))

# Plot
ggplot(df_long, aes(condition, gene, fill = expression)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(c("black", "white", "orangered1"))(100)) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"), 
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) + 
  labs(
    title = "GEPIA3: \nnormal < tumour", 
    fill = "log2(TPM+1)", 
    x = NULL, 
    y = NULL
  )




##### Fig_4d scRNA-seq expression in different cell types from Chloe's webtool  #### 
gene <- read.csv("scRNA_seq/webtool/expr/celltype_expr_oxygen_sensors.csv")

# Chen et al. dataset will be removed as it didn;t pass QC during data preparation. 
rownames(gene) <- gene[,1] # set gene 
gene <- gene[,-1]            # remove first column
gene <- gene[, !grepl("_Chen_", colnames(gene))] # remove Chen dataset 

gene_long <- gene %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -c(gene, group),   # exclude non-numeric column
    names_to = "variable",
    values_to = "value"
  )
gene_long <- gene_long %>%
  separate(variable, into = c("Celltype", "Dataset", "Metric"), sep = "_")
table(gene_long$Celltype)

z_df <- gene_long %>%
  dplyr::filter(Metric == "zscore") %>%
  dplyr::select(gene, Celltype, Dataset, value)

pval_df <- gene_long %>%
  dplyr::filter(Metric == "pval") %>%
  dplyr::select(gene, Celltype, Dataset, value)

pct_df <- gene_long %>%
  dplyr::filter(Metric == "pct") %>%
  dplyr::select(gene, Celltype, Dataset, value)

colnames(z_df)[4] <- "zscore"
colnames(pval_df)[4] <- "pval"
colnames(pct_df)[4] <- "pct"

df <- z_df %>%
  left_join(pval_df, by = c("gene", "Celltype", "Dataset")) %>%
  left_join(pct_df,  by = c("gene", "Celltype", "Dataset"))

df$Celltype <- factor(
  df$Celltype,
  levels = c("Malignant",
             "Myeloid",
             "Endothelial",
             "Pericyte", 
             "Oligodendrocyte", 
             "Lymphocyte")
)

df$Dataset <- factor(
  df$Dataset,
  levels = c("Abdelfattah",
             "Couturier",
             "Ebert",
             "LeBlanc",
             "Neftel")
)

gene_order <- df %>%
  filter(Celltype == "Malignant") %>%
  group_by(gene) %>%
  summarise(malignant_score = mean(zscore, na.rm = TRUE)) %>%
  arrange(malignant_score) %>%
  pull(gene) 

df$gene <- factor(df$gene, levels = gene_order)

df <- df %>%
  mutate(
    signif = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

df_subset <- df %>%
  filter(!gene %in% c("PHYH", "ASPH"))


##### Plot
fig_4d <- ggplot(df, aes(x = Dataset, y = gene, fill = zscore)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = signif), size = 5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0 
  ) +
  facet_wrap(~ Celltype, nrow=1, scales = "free_x") +
  labs(x = NULL, y = NULL, fill = "Z-score", 
       title = "ScRNA-seq: The expressions of oxygen sensor genes in different cell types") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
print(fig_4d)

# # Save
# ggsave(
#   filename = file.path(save_dir_fig_4, "fig_4d.pdf"),
#   plot = fig_4d,
#   device = grDevices::cairo_pdf,
#   width = 11,
#   height = 7.5,
#   units = "in",
#   bg = "transparent"
# )

gene_order <- df %>%
  filter(Celltype == "Malignant") %>%
  group_by(gene) %>%
  summarise(malignant_score = mean(zscore, na.rm = TRUE)) %>%
  arrange(malignant_score) %>%
  pull(gene) # This can be used for ordering genes in GEPIA3 heatmap Fig_4c. 



