
save_dir_fig_2 <- "I:/Writing/02-Meng et al_2027/Figures/Fig 2"

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
genes_OGDD <- c(unique(alias_map_OGDD$genes), "ACTB")
genes_NOX  <- c(unique(alias_map_NOX$genes), "ACTB")
gene_of_interest <- unique(c(genes_OGDD, genes_NOX, # for OGDDs + NOXs 
                      "MT-CO1", # for ETC oxygen-binding Complex IV
                      "CDO1", "ADO", # for TDOs # 69 sensor genes 
                      "ACTB")) # as house-keeping gene # 70 in total 
gene_of_interest
# > gene_of_interest
# [1] "TYW5"   "KDM4E"  "PLOD2"  "UTY"    "PHYH"   "KDM6B"  "P4HA2"  "KDM6A"  "TET3"   "HR"     "PLOD3"  "PHF2"  
# [13] "KDM4A"  "BBOX1"  "KDM4B"  "JMJD7"  "P4HA1"  "KDM5A"  "KDM5C"  "PLOD1"  "ASPH"   "ALKBH1" "P3H1"   "PHYHD1"
# [25] "KDM4D"  "TET2"   "ALKBH2" "JMJD6"  "ALKBH5" "KDM7A"  "KDM3B"  "P4HA3"  "RIOX2"  "P3H2"   "P3H3"   "KDM8"  
# [37] "OGFOD1" "TET1"   "KDM2B"  "ALKBH8" "EGLN2"  "ALKBH3" "KDM5D"  "FTO"    "EGLN1"  "KDM4C"  "RIOX1"  "EGLN3" 
# [49] "JMJD4"  "TMLHE"  "HIF1AN" "P4HTM"  "ALKBH4" "KDM5B"  "PHF8"   "KDM2A"  "KDM3A"  "RSBN1"  "ACTB"   "CYBB"  # ACTB IS INCLUDED HERE
# [61] "NOX5"   "NOX3"   "NOX4"   "DUOX2"  "DUOX1"  "NOX1"   "MT-CO1" "CDO1"   "ADO"   # 68 SENSORS + 1 ACTB = 69 IN THE GENE OF INTEREST 


##### If gene_of_interest is selected, then run from this line: #####

##### Fig 2b - c  plot heatmap ##### 
prepare_hypoxia_bins <- function(seu, cellstate_order, n_bins) {
  
  meta <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(hyp_bin = factor(hyp_bin, levels = as.character(1:n_bins)))
  
  # ---- Cell-state × bin (mean scores) ----
  state_mat <- meta %>%
    group_by(hyp_bin) %>%
    summarise(across(all_of(cellstate_order), mean, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames("hyp_bin") %>%
    as.matrix() %>%
    t()
  
  state_mat <- state_mat[, as.character(1:n_bins), drop = FALSE]
  
  # ---- Gene × bin (mean expression) -------------------------------------
  expr_df <- FetchData(seu, vars = gene_of_interest) %>% # OGDD + NOX
    # expr_df <- FetchData(seu, vars = genes_OGDD) %>% # OGDD 
    # expr_df <- FetchData(seu, vars = genes_NOX) %>% # NOX
    rownames_to_column("cell_id")
  
  gene_mat <- expr_df %>%
    left_join(meta[, c("cell_id", "hyp_bin")], by = "cell_id") %>%
    group_by(hyp_bin) %>%
    summarise(across(-cell_id, mean, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames("hyp_bin") %>%
    as.matrix() %>%
    t() 
  
  gene_mat <- gene_mat[, as.character(1:n_bins), drop = FALSE]
  
  # ---- z-score per gene ----
  gene_z <- t(scale(t(gene_mat)))
  gene_z[!is.finite(gene_z)] <- 0
  
  list(
    meta      = meta,
    state_mat = state_mat,
    gene_mat  = gene_mat, 
    gene_z    = gene_z
  )
}



plot_hypoxia_gene_program_analysis <- function(
    seurat_list,
    cellstate_order,
    housekeeping_gene = "ACTB",
    n_bins = 40
) {
  long_df_all <- list() 
  
  for (dataset_name in names(seurat_list)) {
    message("Processing dataset: ", dataset_name)
    
    # data preparation  
    seu <- seurat_list[[dataset_name]]
    prep <- prepare_hypoxia_bins(seu, cellstate_order, n_bins)
    
    # Split ACTB from gene of interest into housekeeping 
    hk <- housekeeping_gene
    
    gene_mat_hk  <- prep$gene_mat[hk, , drop = FALSE]
    gene_mat_oth <- prep$gene_mat[setdiff(rownames(prep$gene_mat), hk), , drop = FALSE]
    
    gene_z_hk  <- prep$gene_z[hk, , drop = FALSE]
    gene_z_oth <- prep$gene_z[setdiff(rownames(prep$gene_z), hk), , drop = FALSE]
    
    # average expression: sort high to low
    gene_order_avg <- order(rowMeans(gene_mat_oth, na.rm = TRUE), decreasing = TRUE)
    gene_mat_oth_avg <- gene_mat_oth[gene_order_avg, , drop = FALSE]
    
    # z-score: clustering 
    z_row_hclust <- hclust(
      dist(gene_z_oth, method = "euclidean"),
      method = "complete"
    )
    
    gene_z_oth_pheatmap_order <- gene_z_oth[z_row_hclust$order, , drop = FALSE]
    
    # Color scale
    avg_all <- rbind(gene_mat_hk, gene_mat_oth)
    avg_min <- min(avg_all, na.rm = TRUE)
    avg_max <- max(avg_all, na.rm = TRUE)
    
    avg_col_fun <- colorRamp2(
      seq(avg_min, avg_max, length.out = 100),
      viridis(100, option = "turbo")
    )
    
    z_col_fun <- colorRamp2(
      seq(-2, 2, length.out = 100),
      colorRampPalette(c("blue", "white", "red"))(100)
    )
    
    # average expression heatmap
    fig_2b <- Heatmap(
      gene_mat_oth_avg,
      name = "average expression",
      col = avg_col_fun,
      rect_gp = gpar(col = "gray", lwd = 0.05),
      
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      
      show_row_names = TRUE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 8),
      
      show_column_names = TRUE,
      column_names_rot = 45,
      column_names_gp = gpar(fontsize = 8),
      
      top_annotation = HeatmapAnnotation(
        ACTB = anno_simple(
          as.numeric(gene_mat_hk[1, ]),
          col = avg_col_fun,
          height = unit(3, "mm"),
          gp = gpar(col = "gray", lwd = 0.05)
        ),
        annotation_name_side = "right",
        annotation_name_gp = gpar(fontsize = 8)
      ),
      
      column_title = paste("69 oxygen sensors – average expression -", dataset_name),
      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
      
      heatmap_legend_param = list(
        title = "average\nexpression"
      )
    )
    
    # z-score heatmap
    fig_2c <- Heatmap(
      gene_z_oth_pheatmap_order,
      name = "z-score",
      col = z_col_fun,
      rect_gp = gpar(col = "white", lwd = 0.1),
      
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      
      show_row_names = TRUE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 8),
      
      show_column_names = TRUE,
      column_names_rot = 45,
      column_names_gp = gpar(fontsize = 8),
      
      top_annotation = HeatmapAnnotation(
        ACTB = anno_simple(
          as.numeric(gene_z_hk[1, ]),
          col = z_col_fun,
          height = unit(3, "mm"),
          gp = gpar(col = "white", lwd = 0.1)
        ),
        annotation_name_side = "right",
        annotation_name_gp = gpar(fontsize = 8)
      ),
      
      column_title = paste("69 oxygen sensors – z-score -", dataset_name),
      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
      
      heatmap_legend_param = list(
        title = "z-score",
        at = c(-2, -1, 0, 1, 2)
      )
    )
    
    # Plot and Save  
    fig_2 <- list(
      fig_2b = list(plot = fig_2b, width = 6.5, height = 8.5),  # heatmap average expression 
      fig_2c = list(plot = fig_2c, width = 6.5, height = 8.5) # heatmap z score
    )
    
    for (fig_name in names(fig_2)) {
      fig <- fig_2[[fig_name]] 
      grDevices::cairo_pdf(
        filename = file.path(save_dir_fig_2, paste0(fig_name, "_", dataset_name, ".pdf")),
        width = fig$width,
        height = fig$height,
        bg = "transparent" 
      )
      ComplexHeatmap::draw(fig$plot)
      grDevices::dev.off()
      }
  }
}

# Plot and Save 
plot_hypoxia_gene_program_analysis(
  seurat_list = seurat_list,
  cellstate_order = cellstate_order,
  housekeeping_gene = "ACTB",
  n_bins = 40
) 


##### Fig 2d - Spearman correlation analysis: gene expression and state programs across HypoScore #####
plot_gene_vs_cellstate_bin_cor <- function(
    seurat_list,
    cellstate_order,
    housekeeping_gene = "ACTB",
    n_bins = 40,
    exclude_hk = TRUE
) {
  long_df_all <- list()
  
  for (dataset_name in names(seurat_list)) {
    message("Processing dataset: ", dataset_name)
    
    # Data preparation
    seu <- seurat_list[[dataset_name]]
    prep <- prepare_hypoxia_bins(seu, cellstate_order, n_bins)
    
    gene_mat  <- prep$gene_mat
    state_mat <- prep$state_mat
    
    # Exclude housekeeping gene
    if (exclude_hk && housekeeping_gene %in% rownames(gene_mat)) {
      gene_mat <- gene_mat[
        setdiff(rownames(gene_mat), housekeeping_gene), , drop = FALSE
      ]
    }
    
    genes  <- rownames(gene_mat)
    states <- rownames(state_mat)
    
    # Spearman correlation
    rho_mat <- matrix(NA_real_, length(genes), length(states),
                      dimnames = list(genes, states))
    p_mat <- matrix(NA_real_, length(genes), length(states),
                    dimnames = list(genes, states))
    
    for (g in genes) {
      x <- as.numeric(gene_mat[g, ])
      
      for (s in states) {
        y <- as.numeric(state_mat[s, ])
        
        if (sum(complete.cases(x, y)) < 3 ||
            sd(x, na.rm = TRUE) == 0 ||
            sd(y, na.rm = TRUE) == 0) {
          next
        }
        
        ct <- suppressWarnings(
          cor.test(x, y, method = "spearman")
        )
        
        rho_mat[g, s] <- unname(ct$estimate)
        p_mat[g, s]   <- ct$p.value
      }
    }
    
    # Handle undefined correlations
    rho_mat[!is.finite(rho_mat)] <- 0
    p_mat[!is.finite(p_mat)] <- 1
    
    # Multiple testing correction (BH)
    p_adj_mat <- matrix(
      p.adjust(as.vector(p_mat), method = "BH"),
      nrow = nrow(p_mat),
      ncol = ncol(p_mat),
      dimnames = dimnames(p_mat)
    )
    
    # Significance asterisks
    star_mat <- matrix("", nrow = nrow(p_adj_mat), ncol = ncol(p_adj_mat),
                       dimnames = dimnames(p_adj_mat))
    star_mat[p_adj_mat < 0.05]  <- "*"
    star_mat[p_adj_mat < 0.01]  <- "**"
    star_mat[p_adj_mat < 0.001] <- "***"
    
    # Save long-format data (rho, p, p_adj)
    long_df <- data.frame(
      Gene = rep(rownames(rho_mat), times = ncol(rho_mat)),
      CellState = rep(colnames(rho_mat), each = nrow(rho_mat)),
      Rho = as.vector(rho_mat),
      Pval = as.vector(p_mat),
      Padj = as.vector(p_adj_mat)
    )
    long_df$Dataset <- dataset_name
    long_df_all[[dataset_name]] <- long_df
    
    # Heatmap - Fig 2d: Spearman correlation
    fig_2d <- pheatmap(
      rho_mat,
      name = "\u03C1", 
      cluster_rows = TRUE,
      treeheight_row = 60,
      cluster_cols = FALSE,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      # display_numbers = star_mat,
      number_color = "black",
      fontsize = 8, 
      angle_col = "45",
      main = paste("Spearman correlation -", dataset_name),
    )
    
    # Plot and Save
    fig_2 <- list(fig_2d = list(plot = fig_2d, width = 3.5, height = 8.5)) 
    
    for (fig_name in names(fig_2)) {
      fig <- fig_2[[fig_name]]
      grDevices::cairo_pdf(
        filename = file.path(save_dir_fig_2, paste0(fig_name, "_", dataset_name, ".pdf")),
        width = fig$width,
        height = fig$height,
        bg = "transparent"
      )
      ComplexHeatmap::draw(fig$plot)
      dev.off()
    }
  }
  return(bind_rows(long_df_all))
}

# Run
df <- plot_gene_vs_cellstate_bin_cor(
  seurat_list = seurat_list,
  cellstate_order = cellstate_order,
  housekeeping_gene = "ACTB",
  n_bins = 40,
  exclude_hk = TRUE
)


##### Fig 2e - f - Reproducible state specificity  ##### 
# Sort with Rho 
gene_order <- df %>%
  filter(!is.na(Dataset), Dataset != "NA", Dataset != "combined") %>%
  group_by(Gene) %>%
  summarise(mean_pos_rho = mean(Rho), .groups = "drop") %>%
  arrange(desc(mean_pos_rho)) %>%
  pull(Gene)

# Filter out GBO 
df <- df %>%
  mutate(
    Dataset   = factor(Dataset, levels = dataset_order),
    CellState = factor(CellState, levels = cellstate_order),
    Gene      = factor(Gene, levels = gene_order)
  ) %>%
  filter(!Dataset %in% c("GBO"))

df_sig <- df %>%
  mutate(
    sig_label = case_when(
      Padj < 0.001 ~ "***",
      Padj < 0.01  ~ "**",
      Padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )

# Check: unfiltered Heatmap of Spearman Correlation (Rho)
ggplot(df, aes(x = Dataset, y = Gene, fill = Rho)) + 
  geom_tile(color = "white") +
  facet_wrap(~ CellState, nrow = 1) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Heatmap of Spearman correlation (Rho)",
    subtitle = "Sorted by positive Rho value for each gene. ",
    x = "Cell State",
    y = "Gene",
    fill = "Rho"
  )

# Check with numerical rho values 
ggplot(df, aes(x = Dataset, y = Gene, fill = Rho)) +
  geom_tile(color = "white") +
  
  # Add Rho values on tiles
  geom_text(aes(label = sprintf("%.2f", Rho)), 
            size = 3, 
            color = "black") +
  
  facet_wrap(~ CellState, nrow = 1) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Heatmap of Spearman correlation (Rho)",
    subtitle = "Sorted by positive Rho value for each gene.",
    x = "Dataset",
    y = "Gene",
    fill = "Rho"
  )



ggplot(df_sig, aes(x = Dataset, y = Gene, fill = Rho)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sig_label),
    size = 3,
    color = "black"
  ) +
  facet_wrap(~ CellState, nrow = 1) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    subtitle = "Spearman Correlation (Rho) with Adjusted Significance - 91 substrates",
    x = "Dataset",
    y = "Gene",
    fill = "Rho"
  )


df_pos <- df %>%
  filter(!is.na(Rho), Rho > 0) 

fig_2e <- ggplot(df_pos, aes(x = Dataset, y = Gene, fill = Rho)) +
  geom_tile(
    color = "white",
    linewidth = 0.3,
    width = 0.95,
    height = 0.95
  ) +
  facet_wrap(~ CellState, nrow = 1) +
  scale_fill_gradient2(
    low = "white",
    mid = "orange",
    high = "red",
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold"), 
    plot.subtitle = element_text(size = 8), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  labs(
    title = "Positive spearman correlation between oxygen sensor \nand state program across hypoxia score",
    subtitle = "Sorted by positive Rho value for each gene. ",
    x = "ScRNA-seq datasets",
    y = "oxygen sensor genes",
    fill = "\u03C1"
  )
print(fig_2e)
ggsave(
  filename = file.path(save_dir_fig_2, "fig_2e.pdf"),
  plot = fig_2e,
  device = grDevices::cairo_pdf,
  width = 6.5,
  height = 8.5,
  units = "in",
  bg = "transparent"
)

##### Fig 2f ##### 
df_sig_pos <- df %>%
  filter(
    Padj < 0.05, # Change here for other significance ---------------- ---------------- ----------------
    Rho  > 0,
    !is.na(Padj),
    !is.na(Rho)
  )
df_B <- df_sig_pos %>%
  distinct(Gene, CellState, Dataset) %>%
  count(Gene, CellState, name = "B") %>%
  complete(Gene, CellState, fill = list(B = 0)) %>%
  mutate(
    Gene      = factor(Gene, levels = levels(df$Gene)),
    CellState = factor(CellState, levels = levels(df$CellState))
  ) 

df_B_4 <- df_B %>%
  group_by(Gene) %>%
  filter(any(B >= 4)) %>%
  ungroup() %>%
  mutate(
    B_label = ifelse(B >= 4, paste0(B, "/5"), "")
  )
cellstate_score <- c(
  NPC2 = 1,
  NPC1 = 2,
  OPC  = 3,
  AC   = 4,
  MES1 = 5,
  MES2 = 6
)
gene_order_traj <- df_B_4 %>%
  mutate(CellState_num = cellstate_score[as.character(CellState)]) %>%
  group_by(Gene) %>%
  summarize(
    traj_score = weighted.mean(CellState_num, w = B, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(traj_score) %>%
  pull(Gene)

df_B_4 <- df_B_4 %>%
  mutate(Gene = factor(Gene, levels = gene_order_traj))

fig_2f <- ggplot(df_B_4, aes(x = CellState, y = Gene, fill = B)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = B_label), size = 3) +
  scale_fill_gradientn(
    colors = c("white", "white", "red"),
    values = scales::rescale(c(0, 3, 5)),
    limits = c(0, 5),
    breaks = c(4, 5),
    labels = c("4/5", "5/5")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"), 
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title    = "ScRNA-seq: \nReproducibility of significant positive correlation",
    subtitle = paste(
      "X = number of datasets with Rho > 0 and Padj < 0.05 (out of 5). ",
      "Sorted by state enrichment trajectory. ",
      sep = "\n"
    ),
    x        = "Cellular state", 
    y        = "Gene", 
    fill    = "Significant Datasets"
  )
print(fig_2f)

ggsave(
  filename = file.path(save_dir_fig_2, "fig_2f.pdf"),
  plot = fig_2f,
  device = grDevices::cairo_pdf,
  width = 6,
  height = 7,
  units = "in",
  bg = "transparent"
)

gene_of_interest <- unique(as.character(df_B_4$Gene)) # 27 without ACTB
# > gene_of_interest
# [1] "ALKBH4" "ADO"    "KDM4A"  "KDM3B"  "BBOX1"  "PHF2"   "KDM8"  
# [8] "HR"     "KDM5D"  "DUOX1"  "KDM2B"  "P4HTM"  "TET3"   "PHYH"  
# [15] "PLOD3"  "ASPH"   "PLOD1"  "EGLN1"  "KDM4C"  "P4HA1"  "PLOD2" 
# [22] "KDM3A"  "ALKBH5" "P4HA2"  "EGLN3"  "JMJD6"  "KDM7A" 


