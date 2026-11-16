
save_dir_fig_3 <- "I:/Writing/02-Meng et al_2027/Figures/Fig 3"

# obj is prepared with neftel state scores and regional scores per sample block (orig.ident). 
# obj@meta.data$hyp_bin <- ntile(obj@meta.data$hyp_sig, 40)
region_order <- c("IT", "LE", "CT", "CTmvp", "CThbv", "CTpnz", "CTpan") # seven region type
region_order <- c("IT", "LE", "CT", "CTmvp", "CTpan") # five region signature type

obj <- readRDS("bulk/all/GSE107559_RNA.rds") # Saved IvyGAP dataset
n_bins <- 40

# > table(obj$tumour_region)  
# CT CThbv CTmvp CTpan CTpnz    IT    LE 
# 102    22    28    37    24    24    19 

region_colour <- c(
  "IT"="purple","LE"="#218fa5","CT"="green",
  "CTpan"="darkgreen","CTpnz"="#43d1f8",
  "CThbv"="orange","CTmvp"="red"
)

gene_of_interest <- c("ALKBH4", "ADO", "KDM4A", "KDM3B", "BBOX1", 'PHF2', "KDM8", "HR", 'KDM5D', 
                      "DUOX1", 'KDM2B', 'P4HTM', 'TET3', 'PHYH', 'PLOD3', "ASPH", "PLOD1", "EGLN1", 
                      'KDM4C', 'P4HA1', 'PLOD2', 'KDM3A', "ALKBH5", 'P4HA2', 'EGLN3', "JMJD6", "JHDM1D")
# JHDM1D is the another name for KDM7A/ 

## =====================================================
##                    FUNCTIONS
## =====================================================
make_bin_mat <- function(meta, cols, n_bins){
  meta %>%
    group_by(hyp_bin) %>%
    summarise(across(all_of(cols), mean, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames("hyp_bin") %>%
    as.matrix() %>%
    t()
}

make_bin_df <- function(mat, rowname, levels, n_bins){
  as.data.frame(mat) %>%
    rownames_to_column(rowname) %>%
    pivot_longer(-all_of(rowname), names_to = "hyp_bin", values_to = "Score") %>%
    mutate(
      hyp_bin = factor(hyp_bin, levels = as.character(1:n_bins)),
      !!rowname := factor(.data[[rowname]], levels = levels)
    )
}

compute_state_vs_region_bin_cor <- function(seu, cellstate_cols, region_cols, n_bins = 40){
  meta <- seu@meta.data %>%
    mutate(hyp_bin = factor(hyp_bin, levels = as.character(1:n_bins)))
  
  region_mat <- make_bin_mat(meta, region_cols, n_bins)
  state_mat  <- make_bin_mat(meta, cellstate_cols, n_bins)
  
  rho_mat <- outer(
    rownames(state_mat),
    rownames(region_mat),
    Vectorize(function(s, r)
      cor(state_mat[s, ], region_mat[r, ], method = "spearman"))
  )
  dimnames(rho_mat) <- list(rownames(state_mat), rownames(region_mat))
  
  p_mat <- outer(
    rownames(state_mat),
    rownames(region_mat),
    Vectorize(function(s, r)
      cor.test(state_mat[s, ], region_mat[r, ], method = "spearman")$p.value)
  )
  
  p_adj_mat <- matrix(
    p.adjust(as.vector(p_mat), method = "BH"),
    nrow = nrow(p_mat),
    dimnames = dimnames(p_mat)
  )
  
  expand.grid(CellState = rownames(rho_mat),
              Region    = colnames(rho_mat)) %>%
    mutate(Rho  = as.vector(rho_mat),
           Pval = as.vector(p_mat),
           Padj = as.vector(p_adj_mat))
}

## =====================================================
##              DATA PREPARATION
## =====================================================
meta <- obj@meta.data %>%
  mutate(hyp_bin = factor(hyp_bin, levels = as.character(1:n_bins)))

region_mat <- make_bin_mat(meta, region_order, n_bins)
state_mat  <- make_bin_mat(meta, cellstate_order, n_bins)

region_df <- make_bin_df(region_mat, "Region", region_order, n_bins)
state_df  <- make_bin_df(state_mat,  "CellState", cellstate_order, n_bins)

region_df$Region <- factor(region_df$Region, levels = rev(region_order))
state_df$CellState <- factor(state_df$CellState, levels = rev(cellstate_order))

hyp_df <- meta %>%
  group_by(hyp_bin) %>%
  summarise(mean_hyp = mean(hyp_sig, na.rm = TRUE), .groups = "drop")

region_strip_df <- meta %>%
  transmute(hyp_bin, tumour_region = tumour_region)

bin_n_df <- region_strip_df %>% count(hyp_bin)

rho_df <- compute_state_vs_region_bin_cor(obj, cellstate_order, region_order, n_bins) %>%
  mutate(
    asterisk = case_when(
      Padj < 0.001 ~ "***",
      Padj < 0.01  ~ "**",
      Padj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    CellState = factor(CellState, levels = cellstate_order),
    Region    = factor(Region,    levels = region_order)
  )



expr_mat <- FetchData(obj, vars = gene_of_interest) %>%
  mutate(hyp_bin = meta$hyp_bin) %>%
  group_by(hyp_bin) %>%
  summarise(across(all_of(gene_of_interest), mean), .groups = "drop") %>%
  column_to_rownames("hyp_bin") %>%
  as.matrix() %>%
  t()

expr_mat <- expr_mat[, as.character(1:n_bins)]

expr_long_bin <- FetchData(obj, vars = gene_of_interest) %>%
  mutate(hyp_bin = meta$hyp_bin) %>%
  pivot_longer(all_of(gene_of_interest),
               names_to = "Gene",
               values_to = "Expression")

label_df <- expr_long_bin %>%
  group_by(Gene, hyp_bin) %>%
  summarise(Expression = mean(Expression), .groups = "drop") %>%
  group_by(Gene) %>%
  slice_max(as.numeric(hyp_bin), n = 1)



##### Fig 3b - IvyGAP: Per-bin Region Scores ##### 
fig_3a_line <- ggplot(hyp_df, aes(hyp_bin, mean_hyp, group = 1)) +
  geom_line() +
  scale_color_manual(values = "black", name = "Metric") +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "IvyGap: Histological region annotation across HypoScore bins", 
       y = "HypoScore", x = NULL)
print(fig_3a_line) 

region_order <- c("IT", "LE", "CT", "CTmvp", "CThbv", "CTpnz", "CTpan") # seven region type
region_strip_df$tumour_region <- factor(
  region_strip_df$tumour_region,
  levels = region_order
)

fig_3a_bar <- ggplot(region_strip_df, aes(hyp_bin, 1, fill = tumour_region)) +
  geom_col(position = "fill", width = 0.8) + 
  # # Overlap HypoScore line to Region counts
  # geom_line(data = hyp_df, aes(x = hyp_bin, y = hyp_scaled, group = 1), inherit.aes = FALSE,
  #   colour = "black", linewidth = 1) +
  # geom_point(data = hyp_df, aes(x = hyp_bin, y = hyp_scaled), inherit.aes = FALSE,
  #   colour = "black", size = 2) +
  geom_text(data = bin_n_df, aes(x = hyp_bin, y = 1.05, label = n),
            inherit.aes = FALSE, vjust = 0) +
  scale_fill_manual(values = region_colour, breaks = region_order, name = "Histological region") +
  scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
  labs(x = "HypoScore bin", y = "Histological region") +
  theme_minimal() +
  theme_transparent +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), 
    panel.grid = element_blank()
  ) +
  labs(caption = "Numbers indicate sample count per bin.") 
print(fig_3a_bar)

fig_3a <- (fig_3a_line / fig_3a_bar) + 
  plot_layout(heights = c(1.5, 1), guides = "collect") & 
  theme(legend.position = "right")
print(fig_3a)

# ggsave(
#   filename = file.path(save_dir_fig_3, "fig_3a.pdf"),
#   plot = fig_3a,
#   device = grDevices::cairo_pdf,
#   width = 14,
#   height = 5,
#   units = "in",
#   bg = "transparent"
# )


##### Fig 3b - IvyGAP: Per-bin Region Scores ##### 
fig_3b <- ggplot(region_df, aes(hyp_bin, Region, fill = Score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Score)), size = 3) +
  scale_fill_viridis_c(option = "turbo", name = "Region score") +
  theme_minimal() +
  theme(panel.grid = element_blank()) + 
  labs(x= "HypoScore bin", y="Region score")
print(fig_3b) 

# ggsave(
#   filename = file.path(save_dir_fig_3, "fig_3b.pdf"),
#   plot = fig_3b,
#   device = grDevices::cairo_pdf,
#   width = 14,
#   height = 2.2,
#   units = "in",
#   bg = "transparent"
# )


##### Fig 3c - IvyGAP: Per-bin State Scores ##### 
fig_3c <- ggplot(state_df, aes(hyp_bin, CellState, fill = Score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Score)), size = 3) +
  scale_fill_viridis_c(option = "turbo", name = "Cellular state score") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(x = "HypoScore bin", y = "Cellular state score")
print(fig_3c)

# ggsave(
#   filename = file.path(save_dir_fig_3, "fig_3c.pdf"),
#   plot = fig_3c,
#   device = grDevices::cairo_pdf,
#   width = 14,
#   height = 2.5,
#   units = "in",
#   bg = "transparent"
# )

##### Fig 3d - IvyGAP: Spearman correlation state-region ##### 
fig_3d <- ggplot(rho_df, aes(Region, CellState, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = asterisk), size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45),
        panel.grid = element_blank()) +
  labs(title = "IvyGap: Spearman correlation between \ncellular state program score and tumour region signature",
       caption = "Padj < 0.05 *   Padj < 0.01 **   Padj < 0.001 ***", 
       x = "Region marker", y = "Cellular state")
fig_3d

# ggsave(
#   filename = file.path(save_dir_fig_3, "fig_3d.pdf"),
#   plot = fig_3d,
#   device = grDevices::cairo_pdf,
#   width = 5.5,
#   height = 4,
#   units = "in",
#   bg = "transparent"
# )


##### Fig_3e - sensor expression in bin-wise samples, refined in Biorender for gene label ##### 
y_gene <- range(expr_long_bin$Expression, na.rm = TRUE)
y_hyp  <- range(hyp_df$mean_hyp, na.rm = TRUE)

hyp_df2 <- hyp_df %>%
  mutate(mean_hyp_scaled =
           (mean_hyp - y_hyp[1]) / diff(y_hyp) * diff(y_gene) + y_gene[1])

p_line_bin_labeled2 <- ggplot() +
  geom_jitter(data = expr_long_bin,
              aes(as.numeric(hyp_bin), Expression, color = Gene, group = Gene),
              alpha = 0.2, size = 0.6) +
  geom_smooth(data = expr_long_bin,
              aes(as.numeric(hyp_bin), Expression, color = Gene, group = Gene),
              se = FALSE, linewidth = 1.2) +
  geom_text_repel(data = label_df,
                  aes(as.numeric(hyp_bin), Expression, label = Gene, color = Gene),
                  nudge_x = 2, direction = "y",
                  hjust = 0, size = 3, segment.color = NA) +
  # ---- ADD hypoxia line (from fig_3a_line) ----
geom_line(data = hyp_df2,
          aes(as.numeric(hyp_bin), mean_hyp_scaled, group = 1),
          color = "black", linewidth = 1.5) +
  scale_y_continuous(
    name = "Expression",
    sec.axis = sec_axis(
      ~ (. - y_gene[1]) / diff(y_gene) * diff(y_hyp) + y_hyp[1],
      name = "Mean HypoScore / bin"
    )
  ) +
  scale_fill_manual(values = region_colour, name = "Tumour region") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  guides(color = "none") +
  labs(title = "IvyGap: Candidate Gene Expression Across HypoScore bin",
       x = "HypoScore bin (low → high)")

p_line_bin_labeled2 <- p_line_bin_labeled2 +
  theme(
    axis.text.x = element_text(margin = margin(t = 2)),   # reduce gap below x ticks
    axis.text.y = element_text(margin = margin(r = 2)),   # reduce gap left of y ticks
    axis.title.x = element_text(margin = margin(t = 6)),  # bring x title closer
    axis.title.y = element_text(margin = margin(r = 6)),  # bring y title closer
    axis.title.y.right = element_text(margin = margin(l = 6)) # for secondary axis
  )
p_line_bin_labeled2
# ggsave(
#   filename = "IvyGafig_3a_lineoxia_expression.pdf",
#   plot = p_line_bin_labeled2,
#   device = cairo_pdf,   # better font rendering
#   width = 5,           # inches
#   height = 8,          # inches
#   units = "in",
#   dpi = 300
# )


# Compute gene vs state spearman correlation
compute_gene_vs_state_cor <- function(expr_mat, state_mat) {
  
  rho_mat <- outer(
    rownames(expr_mat),
    rownames(state_mat),
    Vectorize(function(g, s)
      cor(expr_mat[g, ], state_mat[s, ], method = "spearman"))
  )
  
  dimnames(rho_mat) <- list(rownames(expr_mat), rownames(state_mat))
  
  p_mat <- outer(
    rownames(expr_mat),
    rownames(state_mat),
    Vectorize(function(g, s)
      cor.test(expr_mat[g, ], state_mat[s, ], method = "spearman")$p.value)
  )
  
  p_adj_mat <- matrix(
    p.adjust(as.vector(p_mat), method = "BH"),
    nrow = nrow(p_mat),
    dimnames = dimnames(p_mat)
  )
  
  expand.grid(
    Gene = rownames(rho_mat),
    CellState = colnames(rho_mat)
  ) %>%
    mutate(
      Rho  = as.vector(rho_mat),
      Pval = as.vector(p_mat),
      Padj = as.vector(p_adj_mat)
    )
}

gene_state_df <- compute_gene_vs_state_cor(expr_mat, state_mat)

gene_state_df <- gene_state_df %>%
  mutate(
    asterisk = case_when(
      Padj < 0.001 ~ "***",
      Padj < 0.01  ~ "**",
      Padj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Gene = factor(Gene, levels = rev(gene_of_interest)),
    CellState = factor(CellState, levels = cellstate_order)
  )

##### Fig_3g - Spearman correlation gene-state in IvyGAP samples ##### 
fig_3g <- ggplot(gene_state_df, aes(CellState, Gene, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = asterisk), size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45),
        panel.grid = element_blank()) +
  labs(
    title = "IvyGap: Spearman correlation between \ngene expression and cellular state",
    # caption = "Padj < 0.05 *   Padj < 0.01 **   Padj < 0.001 ***",
    x = "Cellular state",
    y = "Gene"
  )
fig_3g 


# ggsave(
#   filename = file.path(save_dir_fig_3, "fig_3g.pdf"),
#   plot = fig_3g,
#   device = grDevices::cairo_pdf,
#   width = 4.5,
#   height = 8,
#   units = "in",
#   bg = "transparent"
# )


# # Check: rho values with significance in Fig_3g 
# gene_state_df <- gene_state_df %>%
#   mutate(
#     asterisk = case_when(
#       Padj < 0.001 ~ "***",
#       Padj < 0.01  ~ "**",
#       Padj < 0.05  ~ "*",
#       TRUE ~ ""
#     ),
#     Rho_label = sprintf("%.2f%s", Rho, asterisk),
#     Gene = factor(Gene, levels = rev(gene_of_interest)),
#     CellState = factor(CellState, levels = cellstate_order)
#   )
# 
# fig_3g <- ggplot(gene_state_df, aes(CellState, Gene, fill = Rho)) +
#   geom_tile(color = "white", linewidth = 0.4) +
#   geom_text(aes(label = Rho_label), size = 4) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red",
#                        midpoint = 0, limits = c(-1,1)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45),
#         panel.grid = element_blank()) +
#   labs(
#     title = "IvyGap vs ScRNA-seq: Spearman Correlation\nbetween Gene Expression and Cell State",
#     caption = "Padj < 0.05 *   Padj < 0.01 **   Padj < 0.001 ***",
#     x = "Cellular state",
#     y = "Gene"
#   )
# 
# fig_3g
# fig_3d / fig_3g + 
#   plot_layout(heights = c(1, 5), guides = "collect") & 
#   theme(legend.position = "right")




##### Fig 3f - IvyGAP: Spearman correlation gene-region (per sample) #####
region_scores <- c("IT", "LE", "CT", "CTmvp", "CTpan", "hyp_sig")

# Extract gene expression and region scores from the same samples
expr_mat_sample <- GetAssayData(obj, layer = "data")[gene_of_interest, , drop = FALSE]
region_score_mat <- obj@meta.data[colnames(expr_mat_sample), region_scores, drop = FALSE]

# Spearman correlation
cor_results <- data.frame()

for (score in region_scores) {
  for (gene in gene_of_interest) {
    
    x <- as.numeric(expr_mat_sample[gene, ])
    y <- as.numeric(region_score_mat[[score]])
    
    ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
    
    cor_results <- rbind(cor_results, data.frame(
      Gene = gene,
      Region = score,
      Rho = unname(ct$estimate),
      Pval = ct$p.value
    ))
  }
}

# BH correction within each region score
gene_region_df <- cor_results %>%
  group_by(Region) %>%
  mutate(Padj = p.adjust(Pval, method = "BH")) %>%
  ungroup() %>%
  mutate(
    asterisk = case_when(
      Padj < 0.001 ~ "***",
      Padj < 0.01  ~ "**",
      Padj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Gene = factor(Gene, levels = rev(gene_of_interest)),
    Region = factor(Region, levels = region_scores)
  )

# Heatmap
fig_3f <- ggplot(gene_region_df, aes(Region, Gene, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = asterisk), size = 5) +
  scale_fill_gradientn(
    colours = colorRampPalette(c("blue", "white", "red"))(101),
    limits = c(-1, 1),
    oob = scales::squish,
    name = "\u03C1"
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(
    title = "IvyGAP: Spearman correlation between\nGene Expression and Region Scores",
    caption = "Padj < 0.05 *   Padj < 0.01 **   Padj < 0.001 ***",
    x = "Region region signature",
    y = "Sensor genes"
  )

print(fig_3f) 

# # Save
# ggsave(
#   filename = file.path(save_dir_fig_3, "fig_3f.pdf"),
#   plot = fig_3f,
#   device = grDevices::cairo_pdf,
#   width = 4.5,
#   height = 8,
#   units = "in",
#   bg = "transparent"
# )



