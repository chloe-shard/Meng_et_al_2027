setwd("I:/Data/Aijia/TAOF_R/data") 

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(purrr)
library(patchwork)
library(tibble) 
library(purrr)

library(viridisLite) 
library(pheatmap)
library(readxl)
library(grid)

library(ComplexHeatmap)
library(circlize)

library(ggrepel)

library(readr)
library(reshape2)


seurat_list <- readRDS("scRNA_seq/webtool/processed/seurat_list_malignant.rds") 

cellstate_order <- c("NPC2", "NPC1", "OPC", "AC", "MES1", "MES2")
cellstate_colour <- c(
  "NPC2" = "red", 
  "NPC1" = "pink",
  "OPC"  = "lightgreen",
  "AC"   = "darkgreen",
  "MES1" = "lightblue",
  "MES2" = "darkblue"
) 

dataset_order <- c("Abdelfattah", "Couturier", "Ebert", "LeBlanc", "Neftel")

theme_transparent <- theme(
  plot.background = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.key = element_rect(fill = "transparent", colour = NA),
  strip.background = element_rect(fill = "transparent", colour = NA)
)

save_dir_fig_1 <- "I:/Writing/02-Meng et al_2027/Figures/Fig 1"

##### Figure 1 #####
# Fig 1a-1b
# Extract metadata from each Seurat object into one data frame
df <- imap_dfr(seurat_list, function(seu, dataset_name) {
  meta <- seu@meta.data
  meta$dataset <- dataset_name 
  meta$cell_id <- rownames(meta) # cells in row 
  meta
})
df$cellstate <- factor(df$cellstate, levels = cellstate_order)

# fig_1a - Scatter Plot: Cells ordered by HypoScore
# Sort by HypoScore
df <- df %>%
  group_by(dataset) %>%
  arrange(hyp_sig, .by_group = TRUE) %>%
  mutate(ranked_cell = row_number()) %>% 
  filter(Pooled_celltype == "Malignant" & cellstate != "None" 
         # & dataset != "Neftel"
  ) %>%
  ungroup()


dataset_levels <- c("Abdelfattah", "Couturier", "Ebert", "LeBlanc", "Neftel")

df$dataset <- factor(df$dataset, levels = dataset_levels) 
# show scatter plot in only one dataset 
# 
# df <- df %>%
#   filter(dataset == "Neftel") %>%
#   droplevels()

# fig 1a Plot HypoScore per cell, colored by cellstate
fig_1a <- ggplot(df, aes(x = ranked_cell, y = hyp_sig, color = cellstate)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = cellstate_colour, breaks = cellstate_order, guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    # title = "Per-Cell HypoScore",
    x = "Cells (Ranked by HypoScore)",
    y = "HypoScore",
    color = "Cellular state"
  ) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
  theme_minimal() + theme_transparent + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5)
  )
print(fig_1a)

# # Save  -----------------------------------------------------------
# ggsave(
#   filename = file.path(save_dir_fig_1, "fig_1a.pdf"),
#   plot = fig_1a,
#   device = grDevices::cairo_pdf,
#   width = 10.5,
#   height = 2.5,
#   bg = "transparent"
# )


# Mean HypoScore per dataset and cell state for vertical lines
med_hyp_by_state <- df %>%
  group_by(dataset, cellstate) %>%
  summarise(med = mean(hyp_sig), .groups = "drop")

# fig_1b - Ridgeline plot
fig_1b <- ggplot(df, aes(x = hyp_sig, y = cellstate, fill = cellstate)) +
  coord_cartesian(xlim = c(0, 0.2)) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01, alpha = 0.6, from = 0) +
  geom_segment(data = med_hyp_by_state,aes(x = med, xend = med,y = as.numeric(cellstate) - 0.25,yend = as.numeric(cellstate) + 0.25 ),
               inherit.aes = FALSE,color = "black",linewidth = 0.8) +
  scale_fill_manual(values = cellstate_colour, breaks = cellstate_order) +
  labs(
    # title = "HypoScore Distribution Across Cellular State",
    x = "HypoScore",
    y = "Cellular State",
    fill = "Cellular State"
  ) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
  theme_minimal() + theme_transparent + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5)
  ) 
print(fig_1b)

# # Save  -----------------------------------------------------------
# ggsave(
#   filename = file.path(save_dir_fig_1, "fig_1b.pdf"),
#   plot = fig_1b,
#   device = grDevices::cairo_pdf,
#   width = 10.5,
#   height = 2.5,
#   bg = "transparent"
# )


# Fig 1c-1f
# Main function
plot_hyp_bin_summary <- function(seurat_obj, name) {

  # Shared metadata and bin factor levels
  meta <- seurat_obj@meta.data %>%
    mutate(
      dataset = name, 
      cellstate = factor(cellstate, levels = rev(cellstate_order)),
      hyp_bin = factor(hyp_bin, levels = as.character(1:40))
    )
  

  # Summary statistics used across multiple plots
  bin_state_df <- meta %>%
    group_by(hyp_bin, cellstate) %>%
    summarise(
      cell_count     = n(),
      mean_hyp_state = mean(hyp_sig),
      .groups = "drop"
    ) %>%
    complete(
      hyp_bin = factor(as.character(1:40), levels = as.character(1:40)),
      cellstate     = factor(cellstate_order, levels = cellstate_order),
      fill = list(cell_count = 0)
    )
  
  # Bin-level summaries derived from bin_state_df
  bin_counts <- bin_state_df %>%
    group_by(hyp_bin) %>%
    summarise(bin_size = sum(cell_count), .groups = "drop") # cell counts 
  
  # Mean HypoScore per bin (global)
  mean_hyp_bin <- bin_state_df %>%
    group_by(hyp_bin) %>%
    summarise(
      mean_hyp_sig = weighted.mean(mean_hyp_state, cell_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Mean HypoScore per bin and cell state
  mean_hyp_bin_state <- bin_state_df %>%
    select(hyp_bin, cellstate, mean_hyp_state)
  
  # Cell count per bin and cellstate
  cell_counts_per_bin <- bin_state_df %>%
    select(hyp_bin, cellstate, cell_count)
  plot_df <- cell_counts_per_bin %>%
    left_join(mean_hyp_bin_state, by = c("hyp_bin", "cellstate")) %>%
    left_join(mean_hyp_bin, by = "hyp_bin") %>%
    mutate(delta_hyp = mean_hyp_state - mean_hyp_sig)
  
  # Program scores reshaped long 
  score_summary <- meta %>%
    select(hyp_bin, any_of(cellstate_order)) %>%
    pivot_longer(-hyp_bin, names_to = "cell_state", values_to = "score") %>%
    group_by(hyp_bin, cell_state) %>%
    summarise(mean_score = mean(score), .groups = "drop") %>%
    complete(
      hyp_bin = factor(as.character(1:40), levels = as.character(1:40)),
      cell_state    = factor(cellstate_order, levels = cellstate_order)
    )
  score_summary$cell_state <- factor(score_summary$cell_state,
                                     levels = cellstate_order)
  
  # Join score_summary with bin-level HypoScore and count
  score_summary_ext <- score_summary %>%
    left_join(bin_counts,  by = "hyp_bin") %>%
    left_join(mean_hyp_bin, by = "hyp_bin")
  # --- compute rho per state using score_summary (mean_score across bins) 
  state_mat <- score_summary %>%
    select(hyp_bin, cell_state, mean_score) %>%
    mutate(hyp_bin_num = as.numeric(as.character(hyp_bin))) %>%
    arrange(cell_state, hyp_bin_num) %>%
    select(-hyp_bin) %>%
    pivot_wider(names_from = hyp_bin_num, values_from = mean_score) %>%
    as.data.frame()
  
  rownames(state_mat) <- state_mat$cell_state
  state_mat$cell_state <- NULL
  state_mat <- as.matrix(state_mat)
  
  hyp_index <- as.numeric(colnames(state_mat))
  
  rho_vec <- numeric(nrow(state_mat))
  p_vec   <- numeric(nrow(state_mat))
  names(rho_vec) <- rownames(state_mat)
  names(p_vec)   <- rownames(state_mat)
  
  for (s in rownames(state_mat)) {
    x <- as.numeric(state_mat[s, ])
    ct <- suppressWarnings(cor.test(x, hyp_index, method = "spearman"))
    rho_vec[s] <- unname(ct$estimate)
    p_vec[s]   <- ct$p.value
  }
  
  # handle constants / NA
  rho_vec[!is.finite(rho_vec)] <- 0
  p_vec[!is.finite(p_vec)]     <- 1
  
  p_adj <- p.adjust(p_vec, method = "BH")
  
  star_vec <- ifelse(p_adj < 0.001, "***",
                     ifelse(p_adj < 0.01,  "**",
                            ifelse(p_adj < 0.05,  "*", "")))
  
  rho_df <- tibble(
    cell_state = factor(names(rho_vec), levels = cellstate_order),
    metric = factor("rho", levels = "rho"),
    rho = as.numeric(rho_vec),
    star = star_vec
  ) 

  # Scaling factors (used by plots as-is)
  scale_factor <- max(score_summary$mean_score) /
    max(mean_hyp_bin$mean_hyp_sig)
  
  global_scale <- max(cell_counts_per_bin$cell_count) /
    max(mean_hyp_bin$mean_hyp_sig)  
  
  delta_scale <- max(plot_df$cell_count, na.rm = TRUE) /
    max(abs(plot_df$delta_hyp), na.rm = TRUE)
  
  # Fig 1 Plots -----------------------------------------------------------------------------------------
  # fig_1c — Total cells per bin with HypoScore overlay
  fig_1c <- bin_counts %>%
    left_join(mean_hyp_bin, by = "hyp_bin") %>%
    ggplot(aes(x = hyp_bin)) +
    geom_col(aes(y = bin_size), fill = "grey60") +
    geom_text(aes(y = bin_size, label = bin_size), vjust = 0.5 , size = 3, angle = 45) +
    geom_line(aes(y = mean_hyp_sig * max(bin_size) / max(mean_hyp_sig),
                  color = "Mean HypoScore",
                  group = 1),
              linewidth = 1) +
    geom_point(aes(y = mean_hyp_sig * max(bin_size) / max(mean_hyp_sig),
                   color = "Mean HypoScore"),
               size = 1.5) +
    
    scale_color_manual(
      values = c("Mean HypoScore" = "black"),
      name = NULL
    ) +
    
    labs(
      title = paste("Cell count and mean HypoScore per bin -", name),
      x = "HypoScore bin (ranked by HypoScore from low to high)",
      y = "Number of cells"
    ) +
    theme_minimal() + theme_transparent + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  
  
  # fig_1d — Mean program score trend + mean HypoScore overlay
  fig_1d <- ggplot(score_summary_ext, aes(x = hyp_bin)) +
    geom_line(aes(y = mean_score, color = cell_state, group = cell_state), linewidth = 1) +
    geom_point(aes(y = mean_score, color = cell_state), size = 1.5) +
    
    geom_line(
      aes(
        y = mean_hyp_sig * scale_factor,
        color = "Mean HypoScore",
        group = 1
      ),
      linewidth = 1
    ) +
    geom_point(
      aes(
        y = mean_hyp_sig * scale_factor,
        color = "Mean HypoScore"
      ),
      size = 1.5
    ) +
    
    scale_color_manual(
      values = c(cellstate_colour, "Mean HypoScore" = "black"),
      breaks = c(cellstate_order, "Mean HypoScore")
    ) +
    scale_y_continuous(
      name = "Mean state score",
      sec.axis = sec_axis(~ . / scale_factor, name = "Mean HypoScore")
    ) +
    labs(
      title = paste("Cellular state program transcriptomic score across HypoScore bin -", name),
      x = "HypoScore bin (ranked by HypoScore from low to high)",
      color = "Cellular state"
    ) +
    theme_minimal() + theme_transparent + 
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      axis.title.y.right = element_text(color = "black"),
      axis.text.y.right  = element_text(color = "black")
    )


  # Fig 1e — Heatmap of mean program scores per bin/state
  fig_1e_heatmap <- ggplot(score_summary, aes(hyp_bin, cell_state, fill = mean_score)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "turbo", name = "Score") +
    labs(
      title = paste("Cellular state transcriptomic program score per HypoScore bin -", name),
      x = "HypoScore bin (ranked by HypoScore from low to high)",
      y = "Cellular state"
    ) +
    theme_minimal() + theme_transparent + 
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      panel.grid = element_blank()
    )
  # Fig 1e_rho — Spearman's correlation between cellstate vs HypoScore per bin. 
  fig_1e_rho <- ggplot(rho_df, aes(x = metric, y = cell_state, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = star), size = 4) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Rho"
    ) +
    labs(x = NULL, y = NULL, title = "State vs Bin") +
    theme_minimal() + theme_transparent + 
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10)
    )
  fig_1e <- (fig_1e_heatmap + fig_1e_rho) + 
               plot_layout(widths = c(1, 0.08), guides = "collect") & 
    theme(legend.position = "right", legend.box = "horizontal")
  
  
  # fig_1f — Faceted cell counts + state-specific and global HypoScore
  fig_1f <- ggplot(plot_df, aes(x = hyp_bin)) +
    geom_col(aes(y = cell_count, fill = cellstate, width = 0.8), show.legend = FALSE) +
    # geom_line(aes(y = mean_hyp_state * global_scale, group = 1), color = "grey", linewidth = 1) +
    # geom_point(aes(y = mean_hyp_state * global_scale), color = "grey", size = 1) +
    geom_line(aes(y = mean_hyp_sig * global_scale, group = 1), color = "black", linewidth = 1) +
    facet_wrap(~ cellstate, ncol = 1, scales = "fixed") +
    scale_fill_manual(values = cellstate_colour) +
    scale_y_continuous(
      name = "Number of cells",
      sec.axis = sec_axis(~ . / global_scale, name = "Mean HypoScore")
    ) +
    labs(
      title = paste("Cell counts per cellular state across HypoScore bins -", name),
      x = "HypoScore bin (ranked by HypoScore from low to high)",
    ) +
    theme_minimal() + theme_transparent + 
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y.right = element_text(color = "black"),
      axis.text.y.right  = element_text(color = "black")
    )
  
  # Plot 
  # print(fig_1c)
  # print(fig_1d) 
  # print(fig_1e)
  # print(fig_1f)
  
  
  # # Combined plot panel display
  #   layout <- "
  # AA
  # BB
  # CD
  # "
  # 
  #   fig_1c_to_1e <- (fig_1c / fig_1d / fig_1e) +
  #       plot_layout(guides = "collect") &
  #       theme(legend.position = "right")
  #   fig_1c_to_1f <- (
  #     wrap_elements(full = fig_1c_to_1e) | fig_1f) +
  #     plot_layout(widths = c(2.5, 1))
  # 
  # print(fig_1c_to_1f)
  # 
  
  ##### Not in paper but checked    ##### 
  # fig_1c_stacked — Stacked bar plot of cell state composition per bin 
  # fig_1c_stacked <- ggplot(cell_counts_per_bin, aes(hyp_bin, cell_count, fill = cellstate)) +
  #   geom_col() +
  #   scale_fill_manual(values = cellstate_colour) +
  #   labs(
  #     title = paste("Cellular State Composition -", name),
  #     x = "HypoScore bin (low → high)", y = "Number of cells", fill = "Cellular state"
  #   ) +
  #   theme_minimal() + theme_transparent + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  # fig_1f_delta — Cell counts with Δ-hypoxia (state − global)
  # plot_df <- cell_counts_per_bin %>%
  #   left_join(mean_hyp_bin_state, by = c("hyp_bin", "cellstate")) %>%
  #   left_join(mean_hyp_bin, by = "hyp_bin") %>%
  #   mutate(delta_hyp = mean_hyp_state - mean_hyp_sig)
  # fig_1f_delta <- ggplot(plot_df, aes(x = hyp_bin)) +
  #   geom_col(aes(y = cell_count, fill = cellstate), show.legend = FALSE) +
  #   geom_hline(yintercept = 0, linetype = "dashed") +
  #   geom_line(aes(y = delta_hyp * delta_scale, group = 1),
  #             color = "orange", linewidth = 1, na.rm = TRUE) +
  #   geom_point(aes(y = delta_hyp * delta_scale),
  #              color = "orange", size = 1, na.rm = TRUE) +
  #   facet_wrap(~ cellstate, ncol = 1) +
  #   scale_fill_manual(values = cellstate_colour) +
  #   scale_y_continuous(
  #     name = "Number of cells",
  #     sec.axis = sec_axis(~ . / delta_scale,
  #                         name = expression(Delta * " Mean HypoScore (state − global)"))
  #   ) +
  #   labs(
  #     title = paste("State-specific Deviation from Global Hypoxia Across Bins -", name),
  #     subtitle = "Orange line = deviation from global mean hypoxia within each bin",
  #     x = "HypoScore bin (low → high)"
  #   ) +
  #   theme_minimal() + theme_transparent + 
  #   theme(
  #     strip.text = element_text(face = "bold"),
  #     axis.text.x = element_text(angle = 45, hjust = 1)
  #   )
  # print(fig_1f_delta) # cell counts per state across HypoScore with deviation from global mean hypxoai scores
  
  
  # # Save -----------------------------------------------------------
  # fig_1 <- list(
  #   fig_1c = list(plot = fig_1c, width = 10, height = 3),
  #   fig_1d = list(plot = fig_1d, width = 10, height = 3),
  #   fig_1e = list(plot = fig_1e, width = 10, height = 3),
  #   fig_1f = list(plot = fig_1f, width = 6, height = 9)
  # )
  # 
  # for(fig_name in names(fig_1)) {
  #   fig <- fig_1[[fig_name]]
  #   ggsave(
  #     filename = file.path(save_dir_fig_1, paste0(fig_name, "_", name, ".pdf")),
  #     device = grDevices::cairo_pdf,
  #     plot = fig$plot,
  #     width = fig$width,
  #     height = fig$height,
  #     bg = "transparent"
  # )}
} # -----------------

# Plot and Save
imap(seurat_list, plot_hyp_bin_summary) 



cellstate_order <- c("NPC2", "NPC1", "OPC", "AC", "MES1", "MES2")

##### Check: dataset, state, rho and FDR for fig_1e_rho ##### 
rho_fdr_summary <- imap_dfr(seurat_list, function(seurat_obj, dataset_name) {
  
  bin_scores <- seurat_obj@meta.data %>%
    mutate(hyp_bin = as.numeric(as.character(hyp_bin))) %>%
    filter(hyp_bin %in% 1:40) %>%
    select(hyp_bin, all_of(cellstate_order)) %>%
    pivot_longer(
      cols = all_of(cellstate_order),
      names_to = "state",
      values_to = "score"
    ) %>%
    group_by(hyp_bin, state) %>%
    summarise(mean_score = mean(score), .groups = "drop")
  
  map_dfr(cellstate_order, function(s) {
    
    dat <- bin_scores %>%
      filter(state == s, is.finite(mean_score))
    
    # Retain undefined correlations as NA
    if (nrow(dat) < 3L || n_distinct(dat$mean_score) < 2L) {
      return(tibble(
        dataset = dataset_name,
        state = s,
        rho = NA_real_,
        p_value = NA_real_
      ))
    }
    
    ct <- cor.test(
      dat$mean_score,
      dat$hyp_bin,
      method = "spearman",
      exact = FALSE
    )
    
    tibble(
      dataset = dataset_name,
      state = s,
      rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  }) %>%
    # BH adjustment across the six states within this dataset
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    select(dataset, state, rho, FDR)
})

# Show rho to 3 decimal places and FDR in scientific notation
rho_fdr_display <- rho_fdr_summary %>%
  mutate(
    rho = sprintf("%.3f", rho),
    FDR = sprintf("%.3e", FDR)
  )

print(rho_fdr_display, n = Inf) 


##### Check fig_1b_patient_wise: Patient-wise HypoScore distribution #####
# Mean HypoScore per patient and cellular state
med_hyp_by_patient <- df %>%
  group_by(dataset, orig.ident, cellstate) %>%
  summarise(med = mean(hyp_sig, na.rm = TRUE), .groups = "drop")


# Patient-wise ridgeline plot
plot_patient_ridge <- function(data, name) {
  
  patient_df <- data %>%
    filter(dataset == name)
  
  patient_med <- med_hyp_by_patient %>%
    filter(dataset == name)
  
  ggplot(patient_df, aes(x = hyp_sig, y = cellstate, fill = cellstate)) +
    coord_cartesian(xlim = c(0, 0.2)) +
    geom_density_ridges(
      scale = 1.2, rel_min_height = 0.01,
      alpha = 0.6, from = 0
    ) +
    geom_segment(
      data = patient_med,
      aes(x = med, xend = med,
          y = as.numeric(cellstate) - 0.25,
          yend = as.numeric(cellstate) + 0.25),
      inherit.aes = FALSE,
      color = "black", linewidth = 0.8
    ) +
    scale_fill_manual(values = cellstate_colour, breaks = cellstate_order) +
    labs(
      title = name,
      x = "HypoScore",
      y = "Cellular State",
      fill = "Cellular State"
    ) +
    facet_wrap(~ orig.ident, ncol = 5) +
    theme_minimal() + theme_transparent +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5)
    ) 
  }

fig_1b_patient_wise <-
  wrap_plots(
    plot_patient_ridge(df, "Neftel"),
    plot_patient_ridge(df, "Abdelfattah"),
    plot_patient_ridge(df, "Couturier"),
    plot_patient_ridge(df, "Ebert"),
    plot_patient_ridge(df, "LeBlanc"),
    ncol = 1, 
    heights = c( 4, 3, 3, 1, 2)
  )
print(fig_1b_patient_wise)

# # Save - Supp Figure 1 patient-wise check for distribution of per-state HypoScore
# ggsave(
#   path = "I:/Writing/02-Meng et al_2027/Figures/Figure v2 for manuscript v4 + comments reply/Supplementary Figure 1. Patient-wise HypoScore distribution per cellular state",
#   filename = "fig_1b_patient_wise.pdf", 
#   plot = fig_1b_patient_wise,
#   width = 10,
#   height = 16,
#   units = "in",
#   bg = "transparent", 
#   device = cairo_pdf
# )



