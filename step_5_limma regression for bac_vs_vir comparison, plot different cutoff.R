sample_cols <- colnames(post_peeling_cutoff_subset)

meta <- bind_rows(lapply(sample_cols, generate_meta_data))

bacterial <- c("FLA","LPS","PAM")
viral     <- c("TNF","POLY","R848")

meta_bac_vir <- meta %>%
  mutate(
    Category = case_when(
      Condition %in% bacterial ~ "Bacterial",
      Condition %in% viral     ~ "Viral",
      Condition %in% "GAMP"     ~ "GAMP",
      Condition == "Saline"    ~ "Saline"
    )
  )

Category <- factor(meta_bac_vir$Category, levels = c("Saline","Bacterial","Viral", "GAMP"))
design_bac_vir   <- model.matrix(~ 0 + Category)
colnames(design_bac_vir) <- levels(Category)

fit_bac_vir  <- lmFit(as.matrix(post_peeling_cutoff_subset), design_bac_vir)

contrasts_bac_vir <- makeContrasts(
  Bac_vs_Vir    = Bacterial - Viral,
  Bac_vs_Saline = Bacterial - Saline,
  Vir_vs_Saline = Viral     - Saline,
  levels = design_bac_vir
)

fit_bac_vir_2 <- eBayes(contrasts.fit(fit_bac_vir, contrasts_bac_vir))

res_bac_vir <- topTable(fit_bac_vir_2, coef = "Bac_vs_Vir",    number = Inf) %>%
  rownames_to_column("Genes")

res_bac_saline <- topTable(fit_bac_vir_2, coef = "Bac_vs_Saline", number = Inf) %>%
  rownames_to_column("Genes")

res_vir_saline <- topTable(fit_bac_vir_2, coef = "Vir_vs_Saline", number = Inf) %>%
  rownames_to_column("Genes")

# re-organize output
bac_vir_vs_saline <- res_bac_saline %>%
  dplyr::select(Genes, logFC_bac = logFC, p_value_bac = P.Value, fdr_bac = adj.P.Val) %>%
  dplyr::inner_join(res_vir_saline %>%
                      dplyr::select(Genes, logFC_vir = logFC, p_value_viral = P.Value, fdr_vir = adj.P.Val),
                    by = "Genes")

bac_vs_vir <- res_bac_vir %>%
  dplyr::mutate(dir = ifelse(logFC > 0, "Significant Higher in Bacterial", "Significant Higher in Viral"),
                p_value_bac_vir = P.Value,
                fdr_bac_vir  = adj.P.Val, sig = adj.P.Val < 0.05, logFC_bac_vir = logFC) %>%
  dplyr::select(Genes, logFC_bac_vir, p_value_bac_vir, fdr_bac_vir,  sig, dir)

df_fc <- bac_vir_vs_saline %>%  dplyr::left_join(bac_vs_vir, by = "Genes")

lab_df <- df_fc %>% filter(sig)

df_fc_filtered_by_vs_saline <- df_fc %>%
  filter(p_value_bac <= 0.05 | p_value_viral <= 0.05)

# plot fixed log2fc cutoff
tau_fc   <- 1.0   # strength vs Saline (on log2 scale) 
delta_fc <- 0.5   # difference between categories (on log2 scale)

df_fc_hl_fixed_cutoff <- df_fc_filtered_by_vs_saline %>%
  mutate(
    strong_vs_saline   = pmax(abs(logFC_bac), abs(logFC_vir)) >= tau_fc,          # Only highlight genes has vs_saline log2fc > tau_fc
    strong_across_condition = abs(logFC_bac - logFC_vir) >= delta_fc,             # Only highlight genes has |bac-vir| log2fc > delta_fc
    keep         = sig & strong_across_condition & strong_vs_saline,              # combine (tweak as desired)
    size_var = abs(logFC_bac_vir),
    color_var = -log10(pmax(p_value_bac_vir, 1e-300))
)

ggplot(df_fc_hl_fixed_cutoff, aes(x = logFC_bac, y = logFC_vir)) +
  # all genes (grey)
  geom_point(color = "grey80", alpha = 0.6, size = 1.4) +
  # diagonal guide and the ±delta band
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  
  geom_abline(slope = 1, intercept =  delta_fc, linetype = "dashed", color = "royalblue") +
  geom_abline(slope = 1, intercept = -delta_fc, linetype = "dashed", color = "firebrick") +
  # highlighted genes after filters
  geom_point(data = subset(df_fc_hl_fixed_cutoff, keep),
             aes(color = color_var, size = size_var), show.legend = TRUE) +
  # labels for highlighted genes (top N by Bac vs Viral FDR)
  geom_text_repel(
    data = df_fc_hl_fixed_cutoff %>% filter(keep),
    aes(label = Genes, color = color_var),
    size = 3, box.padding = 0.25, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_viridis_c(
    option = "C",
    name   = expression(-log[10]~"p (Bac vs Viral)"),
    limits = c(0, 10), oob = scales::squish
  ) +
  scale_size_continuous(
    range  = c(2, 5.0), 
    breaks = c(0.5, 1, 2, 3, 4),
    name   = "|Δ log2FC| (Bac − Vir)"
  ) +
  scale_y_continuous(limits = c(-5,10)) + 
  scale_x_continuous(limits = c(-5,10)) +
  theme_classic(base_size = 14) +
  labs(
    x = "log2FC (Bacterial vs Saline)",
    y = "log2FC (Viral vs Saline)",
    title = "Effect–effect scatter with strength (vs_Saline log2FC>|1|) and difference filters (bac-vir log2FC>|0.5|)"
  )

# plot floating log2fc cutoff for both vs_saline_log2fc and bac_vir_log2fc
df_fc_hl_floating_cutoff <- df_fc_filtered_by_vs_saline %>%
  mutate(
    strong_vs_saline   = abs(logFC_bac - mean(logFC_bac)) >= sd(logFC_bac) |  abs(logFC_vir - mean(logFC_vir)) >= sd(logFC_vir),    # Only highlight genes has vs_saline log2fc larger than mean +- 1SD
    strong_across_condition = abs(logFC_bac_vir - mean(logFC_bac_vir)) >= sd(logFC_bac_vir),                                        # Only highlight genes has bac-vir log2fc larger than mean +- 1SD 
    keep         = sig & strong_across_condition & strong_vs_saline,          
    size_var = abs(logFC_bac_vir),
    color_var = -log10(pmax(p_value_bac_vir, 1e-300))
)

ggplot(df_fc_hl_floating_cutoff, aes(x = logFC_bac, y = logFC_vir)) +
  # all genes (grey)
  geom_point(color = "grey80", alpha = 0.6, size = 1.4) +
  # diagonal guide and the ±delta band
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  
  geom_abline(slope = 1, intercept =  delta_fc, linetype = "dashed", color = "royalblue") +
  geom_abline(slope = 1, intercept = -delta_fc, linetype = "dashed", color = "firebrick") +
  # highlighted genes after filters
  geom_point(data = subset(df_fc_hl_floating_cutoff, keep),
             aes(color = color_var, size = size_var), show.legend = TRUE) +
  # labels for highlighted genes (top N by Bac vs Viral FDR)
  geom_text_repel(
    data = df_fc_hl_floating_cutoff %>% filter(keep),
    aes(label = Genes, color = color_var),
    size = 3, box.padding = 0.25, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_viridis_c(
    option = "C",
    name   = expression(-log[10]~"p (Bac vs Viral)"),
    limits = c(0, 10), oob = scales::squish
  ) +
  scale_size_continuous(
    range  = c(2, 5.0), 
    breaks = c(0.5, 1, 2, 3, 4),
    name   = "|Δ log2FC| (Bac − Vir)"
  ) +
  scale_y_continuous(limits = c(-5,10)) + 
  scale_x_continuous(limits = c(-5,10)) +
  theme_classic(base_size = 14) +
  labs(
    x = "log2FC (Bacterial vs Saline)",
    y = "log2FC (Viral vs Saline)",
    title = "Effect–effect scatter with strength (vs_Saline |log2FC - mean| >= 1SD) and difference filters (|log2FC - mean| >= 1SD)"
  )

# floating vs_saline_log2fc, ignore bac_vir_log2fc
df_fc_hl_floating_cutoff_only_for_saline <- df_fc_filtered_by_vs_saline %>%
  mutate(
    strong_vs_saline   = abs(logFC_bac - mean(logFC_bac)) >= sd(logFC_bac) |  abs(logFC_vir - mean(logFC_vir)) >= sd(logFC_vir),    # Only highlight genes has vs_saline log2fc larger than mean +- 1SD
    keep         = sig & strong_vs_saline,          
    size_var = abs(logFC_bac_vir),
    color_var = -log10(pmax(p_value_bac_vir, 1e-300))
  )

ggplot(df_fc_hl_floating_cutoff_only_for_saline, aes(x = logFC_bac, y = logFC_vir)) +
  # all genes (grey)
  geom_point(color = "grey80", alpha = 0.6, size = 1.4) +
  # diagonal guide and the ±delta band
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  
  geom_abline(slope = 1, intercept =  delta_fc, linetype = "dashed", color = "royalblue") +
  geom_abline(slope = 1, intercept = -delta_fc, linetype = "dashed", color = "firebrick") +
  # highlighted genes after filters
  geom_point(data = subset(df_fc_hl_floating_cutoff_only_for_saline, keep),
             aes(color = color_var, size = size_var), show.legend = TRUE) +
  # labels for highlighted genes (top N by Bac vs Viral FDR)
  geom_text_repel(
    data = df_fc_hl_floating_cutoff_only_for_saline %>% filter(keep),
    aes(label = Genes, color = color_var),
    size = 3, box.padding = 0.25, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_viridis_c(
    option = "C",
    name   = expression(-log[10]~"p (Bac vs Viral)"),
    limits = c(0, 10), oob = scales::squish
  ) +
  scale_size_continuous(
    range  = c(2, 5.0), 
    breaks = c(0.5, 1, 2, 3, 4),
    name   = "|Δ log2FC| (Bac − Vir)"
  ) +
  scale_y_continuous(limits = c(-5,10)) + 
  scale_x_continuous(limits = c(-5,10)) +
  theme_classic(base_size = 14) +
  labs(
    x = "log2FC (Bacterial vs Saline)",
    y = "log2FC (Viral vs Saline)",
    title = "Effect–effect scatter with strength (vs_Saline |log2FC - mean| >= 1SD)"
  )

df_fc_hl_no_cutoff <- df_fc_filtered_by_vs_saline %>%
  mutate(
    keep         = sig,          
    size_var = abs(logFC_bac_vir),
    color_var = -log10(pmax(p_value_bac_vir, 1e-300))
  )

ggplot(df_fc_hl_no_cutoff, aes(x = logFC_bac, y = logFC_vir)) +
  # all genes (grey)
  geom_point(color = "grey80", alpha = 0.6, size = 1.4) +
  # diagonal guide and the ±delta band
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_hline(yintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = tau_fc, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = -tau_fc, linetype = "dotted", color = "grey80") +
  
  geom_abline(slope = 1, intercept =  delta_fc, linetype = "dashed", color = "royalblue") +
  geom_abline(slope = 1, intercept = -delta_fc, linetype = "dashed", color = "firebrick") +
  # highlighted genes after filters
  geom_point(data = subset(df_fc_hl_no_cutoff, keep),
             aes(color = color_var, size = size_var), show.legend = TRUE) +
  # labels for highlighted genes (top N by Bac vs Viral FDR)
  # geom_text_repel(
  #   data = df_fc_hl_no_cutoff %>% filter(keep),
  #   aes(label = Genes, color = color_var),
  #   size = 3, box.padding = 0.25, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_viridis_c(
    option = "C",
    name   = expression(-log[10]~"p (Bac vs Viral)"),
    limits = c(0, 10), oob = scales::squish
  ) +
  scale_size_continuous(
    range  = c(2, 5.0), 
    breaks = c(0.5, 1, 2, 3, 4),
    name   = "|Δ log2FC| (Bac − Vir)"
  ) +
  scale_y_continuous(limits = c(-5,10)) + 
  scale_x_continuous(limits = c(-5,10)) +
  theme_classic(base_size = 14) +
  labs(
    x = "log2FC (Bacterial vs Saline)",
    y = "log2FC (Viral vs Saline)",
    title = "Effect–effect scatter plotting all genes significantly different from saline, highlight genes significant different between bac and vir"
  )
