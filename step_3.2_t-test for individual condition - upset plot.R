head(post_peeling_cutoff_subset)
get_col_indices(post_peeling_cutoff_subset)

FLA_stats <- perform_stats_analysis_by_t_test(post_peeling_cutoff_subset, "FLA", 5:8, 1:4)
LPS_stats <- perform_stats_analysis_by_t_test(post_peeling_cutoff_subset, "LPS", 9:12, 1:4)
PAM_stats <- perform_stats_analysis_by_t_test(post_peeling_cutoff_subset, "PAM", 13:16, 1:4)
R848_stats <- perform_stats_analysis_by_t_test(post_peeling_cutoff_subset, "R848", 17:20, 1:4)
TNF_stats <- perform_stats_analysis_by_t_test(post_peeling_cutoff_subset, "TNF", 21:24, 1:4)
POLY_stats <- perform_stats_analysis_by_t_test(post_peeling_cutoff_subset, "POLY", 25:28, 1:4)
GAMP_stats <- perform_stats_analysis_by_t_test(post_peeling_cutoff_subset, "GAMP", 29:32, 1:4)

# put all stats files into one large file --------------------------------

post_peeling_log2FC_mean_1SD <- data.frame(
  Genes = FLA_stats$Genes,
  
  # FLA
  FLA_p_values = FLA_stats$FLA_p_values,
  FLA_Log2FC   = FLA_stats$FLA_Log2FC,
  FLA_vs_saline = ifelse(
    FLA_stats$FLA_significance == "significant" & grepl("beyond", FLA_stats$FLA_SD_category),
    FLA_stats$FLA_Direction, ""),
  
  # GAMP
  GAMP_p_values = GAMP_stats$GAMP_p_values,
  GAMP_Log2FC   = GAMP_stats$GAMP_Log2FC,
  GAMP_vs_saline = ifelse(
    GAMP_stats$GAMP_significance == "significant" & grepl("beyond", GAMP_stats$GAMP_SD_category),
    GAMP_stats$GAMP_Direction, ""),
  
  # LPS
  LPS_p_values = LPS_stats$LPS_p_values,
  LPS_Log2FC   = LPS_stats$LPS_Log2FC,
  LPS_vs_saline = ifelse(
    LPS_stats$LPS_significance == "significant" & grepl("beyond", LPS_stats$LPS_SD_category),
    LPS_stats$LPS_Direction, ""),
  
  # PAM
  PAM_p_values = PAM_stats$PAM_p_values,
  PAM_Log2FC   = PAM_stats$PAM_Log2FC,
  PAM_vs_saline = ifelse(
    PAM_stats$PAM_significance == "significant" & grepl("beyond", PAM_stats$PAM_SD_category),
    PAM_stats$PAM_Direction, ""),
  
  # POLY
  POLY_p_values = POLY_stats$POLY_p_values,
  POLY_Log2FC   = POLY_stats$POLY_Log2FC,
  POLY_vs_saline = ifelse(
    POLY_stats$POLY_significance == "significant" & grepl("beyond", POLY_stats$POLY_SD_category),
    POLY_stats$POLY_Direction, ""),
  
  # R848
  R848_p_values = R848_stats$R848_p_values,
  R848_Log2FC   = R848_stats$R848_Log2FC,
  R848_vs_saline = ifelse(
    R848_stats$R848_significance == "significant" & grepl("beyond", R848_stats$R848_SD_category),
    R848_stats$R848_Direction, ""),
  
  # TNF
  TNF_p_values = TNF_stats$TNF_p_values,
  TNF_Log2FC   = TNF_stats$TNF_Log2FC,
  TNF_vs_saline = ifelse(
    TNF_stats$TNF_significance == "significant" & grepl("beyond", TNF_stats$TNF_SD_category),
    TNF_stats$TNF_Direction, "")
)

# make UpSet input for UP and DOWN
UP_df_t_test <- data.frame(
  Gene = post_peeling_log2FC_mean_1SD$Genes,
  sapply(conditions, function(drug) {
    ifelse(post_peeling_log2FC_mean_1SD[[paste0(drug, "_vs_saline")]] == "UP", 1, 0)
  })
)
UP_df_t_test <- UP_df_t_test[rowSums(UP_df_t_test[,-1]) > 0, ]
colnames(UP_df_t_test)[-1] <- conditions
rownames(UP_df_t_test) <- UP_df_t_test$Gene
UP_df_t_test$Gene <- NULL

DOWN_df_t_test <- data.frame(
  Gene = post_peeling_log2FC_mean_1SD$Genes,
  sapply(conditions, function(drug) {
    ifelse(post_peeling_log2FC_mean_1SD[[paste0(drug, "_vs_saline")]] == "DOWN", 1, 0)
  })
)

DOWN_df_t_test <- DOWN_df_t_test[rowSums(DOWN_df_t_test[,-1]) > 0, ]
colnames(DOWN_df_t_test)[-1] <- conditions
rownames(DOWN_df_t_test) <- DOWN_df_t_test$Gene
DOWN_df_t_test$Gene <- NULL

upset(UP_df_t_test, 
      sets = colnames(UP_df_t_test),
      keep.order = TRUE,
      sets.bar.color = "steelblue",
      order.by = "degree",
      text.scale = 1.8,
      nintersects = NA)

upset(DOWN_df_t_test, 
      sets = colnames(DOWN_df_t_test),
      keep.order = TRUE,
      sets.bar.color = "steelblue",
      order.by = "degree",
      text.scale = 1.8,
      nintersects = NA)

# -------- Step 1: Convert to long format --------
UP_long_ttest <- UP_df_t_test %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Condition", values_to = "Present") %>%
  filter(Present == 1)

# -------- Step 2: Aggregate conditions per gene --------
UP_combined_ttest <- UP_long_ttest %>%
  group_by(Gene) %>%
  summarise(Combination = paste(sort(Condition), collapse = " / "), .groups = "drop")

# -------- Step 3: Group genes by unique intersection --------
UP_intersections_ttest <- UP_combined_ttest %>%
  group_by(Combination) %>%
  summarise(GeneList = list(Gene), .groups = "drop")

# -------- Step 4: Format to wide for export --------
max_genes_UP <- max(lengths(UP_intersections_ttest$GeneList))

UP_wide_ttest <- t(sapply(1:nrow(UP_intersections_ttest), function(i) {
  c(UP_intersections_ttest$Combination[i],
    UP_intersections_ttest$GeneList[[i]],
    rep(NA, max_genes_UP - length(UP_intersections_ttest$GeneList[[i]])))
}))

UP_wide_df_ttest <- as.data.frame(UP_wide_ttest, stringsAsFactors = FALSE)
colnames(UP_wide_df_ttest)[1] <- "Condition"

# -------- Step 5: Export to Excel --------
write.xlsx(UP_wide_df_ttest, file = "UP_gene_intersect_from_t_test.xlsx", rowNames = FALSE)

# do the same DOWN df
DOWN_long_ttest <- DOWN_df_t_test %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Condition", values_to = "Present") %>%
  filter(Present == 1)

DOWN_combined_ttest <- DOWN_long_ttest %>%
  group_by(Gene) %>%
  summarise(Combination = paste(sort(Condition), collapse = " / "), .groups = "drop")

DOWN_intersections_ttest <- DOWN_combined_ttest %>%
  group_by(Combination) %>%
  summarise(GeneList = list(Gene), .groups = "drop")

max_genes_DOWN <- max(lengths(DOWN_intersections_ttest$GeneList))

DOWN_wide_ttest <- t(sapply(1:nrow(DOWN_intersections_ttest), function(i) {
  c(DOWN_intersections_ttest$Combination[i],
    DOWN_intersections_ttest$GeneList[[i]],
    rep(NA, max_genes_DOWN - length(DOWN_intersections_ttest$GeneList[[i]])))
}))

DOWN_wide_df_ttest <- as.data.frame(DOWN_wide_ttest, stringsAsFactors = FALSE)
colnames(DOWN_wide_df_ttest)[1] <- "Condition"

write.xlsx(DOWN_wide_df_ttest, file = "DOWN_gene_intersect_from_t_test.xlsx", rowNames = FALSE)
