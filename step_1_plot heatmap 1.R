
# heatmap on log2 values and z-scored log2 values -------------------------

heatmap_input_log2_values <- as.matrix(post_peeling_cutoff)

pheatmap(heatmap_input_log2_values, scale = "none", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,
         show_rownames = FALSE,       
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         clustering_method = "complete",
         main = "input = log2value")

pheatmap(heatmap_input_log2_values, scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,         
         show_rownames = FALSE,       
         clustering_method = "complete",
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored log2value")

post_peeling_cutoff_subset <- post_peeling_cutoff %>% 
  dplyr::select(1:32)

heatmap_input_log2_values_subset <- as.matrix(post_peeling_cutoff_subset)

pheatmap(heatmap_input_log2_values_subset, scale = "none", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,         
         show_rownames = FALSE,       
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         clustering_method = "complete",
         main = "input = log2value, only Labeled samples")

pheatmap(heatmap_input_log2_values_subset, scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,          
         show_rownames = FALSE,       
         clustering_method = "complete",
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored log2value, only Labeled samples")


# heatmap on labeled / non-labeled ratios (8 ratios) ----------------------

# Define each L condition and its matching nL condition
lnl_pairs <- list(
  SL     = list(L = paste0("SL_",     1:4), nL = paste0("SNL_",     1:2)),
  FLAL   = list(L = paste0("FLAL_",   1:4), nL = paste0("FLANL_",   1:2)),
  LPSL   = list(L = paste0("LPSL_",   1:4), nL = paste0("LPSNL_",   1:2)),
  PAML   = list(L = paste0("PAML_",   1:4), nL = paste0("PAMNL_",   1:2)),
  R848L  = list(L = paste0("R848L_",  1:4), nL = paste0("R848NL_",  1:2)),
  TNFL   = list(L = paste0("TNFL_",   1:4), nL = paste0("TNFNL_",   1:2)),
  POLYL  = list(L = paste0("POLYL_",  1:4), nL = paste0("POLYNL_",  1:2)),
  GAMPL  = list(L = paste0("GAMPL_",  1:4), nL = paste0("GAMPNL_",  1:2))
)
# Initialize empty list to store results
log2_fc_list <- list()

# For each condition (e.g., SL, FLAL, etc.)
for (condition in names(lnl_pairs)) {
  L_samples  <- lnl_pairs[[condition]]$L
  nL_samples <- lnl_pairs[[condition]]$nL
  
  # For each L and each nL combination
  for (L in L_samples) {
    for (nL in nL_samples) {
      col_name <- paste(L, "-", nL, sep = "")
      log2_fc_list[[col_name]] <- post_peeling_cutoff[, L] - post_peeling_cutoff[, nL]
    }
  }
}

# Combine all into a data frame
log2_fc_df <- as.data.frame(log2_fc_list)
rownames(log2_fc_df) <- rownames(post_peeling_cutoff)

pheatmap(log2_fc_df, scale = "none", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,         
         show_rownames = FALSE,       
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         clustering_method = "complete",
         main = "input = L/nL")

pheatmap(log2_fc_df, scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,             
         show_rownames = FALSE,       
         clustering_method = "complete",
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored L/nL")

# heatmap on labeled / non-labeled ratios (4 ratios) ----------------------

log2_fc_list_nL_mean <- list()

for (cond in names(lnl_pairs)) {
  L_cols <- lnl_pairs[[cond]]$L
  nL_cols <- lnl_pairs[[cond]]$nL
  
  # Compute mean of nL columns (vector of length = # of genes)
  nL_mean <- rowMeans(post_peeling_cutoff[, nL_cols, drop = FALSE])
  
  for (L in L_cols) {
    col_name <- paste0(L, "_vs_", cond, "_nLmean")
    
    # Calculate log2FC and store as named vector
    log2_fc_list_nL_mean[[col_name]] <- post_peeling_cutoff[, L] - nL_mean
  }
}

log2_fc_df_nL_mean <- as.data.frame(do.call(cbind, log2_fc_list_nL_mean))
rownames(log2_fc_df_nL_mean) <- rownames(post_peeling_cutoff)

pheatmap(log2_fc_df_nL_mean, scale = "none", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,         
         show_rownames = FALSE,       
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         clustering_method = "complete",
         main = "input = L/nL_mean")

pheatmap(log2_fc_df_nL_mean, scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,        
         show_rownames = FALSE,       
         clustering_method = "complete",
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored L/nL_mean")


# same heatmaps with locked order -----------------------------------------

heatmap_input_log2_values <- as.matrix(post_peeling_cutoff)
group_labels_1 <- gsub("_.*", "", colnames(heatmap_input_log2_values))
annotation_col_1 <- data.frame(Group = group_labels_1)
rownames(annotation_col_1) <- colnames(heatmap_input_log2_values)

pheatmap(heatmap_input_log2_values, scale = "row", 
         cluster_cols = FALSE, 
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,        
         show_rownames = FALSE,       
         clustering_method = "complete",
         annotation_col = annotation_col_1,
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored log2value_locked_order")

group_labels_2 <- heatmap_input_log2_values_subset
group_labels_2 <- gsub("_.*", "", colnames(heatmap_input_log2_values))
annotation_col_2 <- data.frame(Group = group_labels_2)
rownames(annotation_col_2) <- colnames(heatmap_input_log2_values)

pheatmap(heatmap_input_log2_values_subset, scale = "row", 
         cluster_cols = FALSE,
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,          
         show_rownames = FALSE,       
         clustering_method = "complete",
         annotation_col = annotation_col_2,
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored log2value_locked_order")

group_labels_3 <- log2_fc_df
group_labels_3 <- gsub("_.*", "", colnames(log2_fc_df))
annotation_col_3 <- data.frame(Group = group_labels_3)
rownames(annotation_col_3) <- colnames(log2_fc_df)

pheatmap(log2_fc_df, scale = "row", 
         cluster_cols = FALSE,
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,             
         show_rownames = FALSE,       
         clustering_method = "complete",
         annotation_col = annotation_col_3,
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored L/nL_locked_order")

group_labels_4 <- log2_fc_df_nL_mean
group_labels_4 <- gsub("_.*", "", colnames(log2_fc_df_nL_mean))
annotation_col_4 <- data.frame(Group = group_labels_4)
rownames(annotation_col_4) <- colnames(log2_fc_df_nL_mean)

pheatmap(log2_fc_df_nL_mean, scale = "row", 
         cluster_cols = FALSE,
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0,        
         show_rownames = FALSE,       
         clustering_method = "complete",
         annotation_col = annotation_col_4,
         color = colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101),
         breaks = seq(-3, 3, length.out = 101),
         main = "input = z_scored L/nL_mean_locked_order")
