sample_cols <- colnames(post_peeling_cutoff_subset)

meta <- bind_rows(lapply(sample_cols, generate_meta_data))

expr <- as.matrix(post_peeling_cutoff_subset)

# relevel groups (baseline = Saline)
group <- fct_inorder(meta$Condition)            
group <- relevel(group, ref = "Saline")

# design without intercept
design <- model.matrix(~ 0 + group)                  
colnames(design) <- levels(group)

fit  <- lmFit(expr, design)

contrasts_matrix <- makeContrasts(
  FLA_vs_saline    = FLA - Saline,
  LPS_vs_saline    = LPS - Saline,
  PAM_vs_saline    = PAM - Saline,
  R848_vs_saline   = R848 - Saline,
  TNF_vs_saline    = TNF - Saline,
  POLY_vs_saline   = POLY - Saline,
  GAMP_vs_saline   = GAMP - Saline,
  levels = design
)

fit_individual_comparison <- eBayes(contrasts.fit(fit, contrasts_matrix))

get_coef_and_sig <- function(fit, condition, sig_cutoff, log2fc_cutoff){
  
  coef_name <- paste0(condition, "_vs_saline")
  df <- topTable(fit, coef = coef_name,    number = Inf) %>%
    rownames_to_column("Genes") %>%
    mutate(
      sig = ifelse(adj.P.Val < sig_cutoff, "Yes", "No"),
      direction = case_when(
        logFC >  (mean(logFC) + sd(logFC))  & sig == "Yes" ~ "UP",       # (mean(logFC) + sd(logFC))    OR      log2fc_cutoff
        logFC <  (mean(logFC) - sd(logFC)) & sig == "Yes" ~ "DOWN",   # check careful, for down, it is mean-1sd
        TRUE ~ "NS"   # Non-significant
      )
    )
  return(df)
}

conditions <- c("FLA", "LPS", "PAM", "R848", "TNF", "POLY", "GAMP")

res_list <- lapply(conditions, function(cond) {
  get_coef_and_sig(fit_individual_comparison, cond, 0.05, 1)
})
names(res_list) <- conditions
res_list$FLA
up_list <- lapply(res_list, function(df) df$Genes[df$direction == "UP"])
down_list <- lapply(res_list, function(df) df$Genes[df$direction == "DOWN"])

all_UP_genes <- unique(unlist(up_list))
UP_df <- data.frame(
  row.names = all_UP_genes,
  FLA  = as.integer(all_UP_genes %in% up_list$FLA),
  LPS  = as.integer(all_UP_genes %in% up_list$LPS),
  PAM  = as.integer(all_UP_genes %in% up_list$PAM),
  R848 = as.integer(all_UP_genes %in% up_list$R848),
  TNF  = as.integer(all_UP_genes %in% up_list$TNF),
  POLY = as.integer(all_UP_genes %in% up_list$POLY),
  GAMP = as.integer(all_UP_genes %in% up_list$GAMP)
)

all_down_genes <- unique(unlist(down_list))
DOWN_df <- data.frame(
  row.names = all_down_genes,
  FLA  = as.integer(all_down_genes %in% down_list$FLA),
  LPS  = as.integer(all_down_genes %in% down_list$LPS),
  PAM  = as.integer(all_down_genes %in% down_list$PAM),
  R848 = as.integer(all_down_genes %in% down_list$R848),
  TNF  = as.integer(all_down_genes %in% down_list$TNF),
  POLY = as.integer(all_down_genes %in% down_list$POLY),
  GAMP = as.integer(all_down_genes %in% down_list$GAMP)
)

upset(UP_df, 
      sets = colnames(UP_df),
      keep.order = TRUE,
      sets.bar.color = "steelblue",
      order.by = "degree",
      text.scale = 1.8,
      nintersects = NA)

upset(DOWN_df, 
      sets = colnames(DOWN_df),
      keep.order = TRUE,
      sets.bar.color = "steelblue",
      order.by = "degree",
      text.scale = 1.8,
      nintersects = NA)

# Step 1: Convert binary matrix to long format
UP_long <- UP_df %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Condition", values_to = "Present") %>%
  filter(Present == 1)

# Step 2: Re-aggregate into a list of conditions per gene
UP_combined <- UP_long %>%
  group_by(Gene) %>%
  summarise(Combination = paste(sort(Condition), collapse = " / "), .groups = "drop")

# Step 3: Group genes by each unique intersection
UP_intersections <- UP_combined %>%
  group_by(Combination) %>%
  summarise(GeneList = list(Gene), .groups = "drop")

max_genes_UP <- max(lengths(UP_intersections$GeneList))

UP_wide <- t(sapply(1:nrow(UP_intersections), function(i) {
  c(UP_intersections$Combination[i], UP_intersections$GeneList[[i]], rep(NA, max_genes_UP - length(UP_intersections$GeneList[[i]])))
}))

UP_wide_df <- as.data.frame(UP_wide, stringsAsFactors = FALSE)
colnames(UP_wide_df)[1] <- "Condition"
write.xlsx(UP_wide_df, file = "UP_gene_intersect_from_limma.xlsx", rowNames = FALSE)

# Now make the same output for DOWN
DOWN_long <- DOWN_df %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Condition", values_to = "Present") %>%
  filter(Present == 1)

DOWN_combined <- DOWN_long %>%
  group_by(Gene) %>%
  summarise(Combination = paste(sort(Condition), collapse = " / "), .groups = "drop")

DOWN_intersections <- DOWN_combined %>%
  group_by(Combination) %>%
  summarise(GeneList = list(Gene), .groups = "drop")

max_genes_DOWN <- max(lengths(DOWN_intersections$GeneList))

DOWN_wide <- t(sapply(1:nrow(DOWN_intersections), function(i) {
  c(DOWN_intersections$Combination[i], DOWN_intersections$GeneList[[i]], rep(NA, max_genes_DOWN - length(DOWN_intersections$GeneList[[i]])))
}))

DOWN_wide_df <- as.data.frame(DOWN_wide, stringsAsFactors = FALSE)
colnames(DOWN_wide_df)[1] <- "Condition"
write.xlsx(DOWN_wide_df, file = "DOWN_gene_intersect_from_limma.xlsx", rowNames = FALSE)
