# For step 3 - already have list from limma output
# Perform GO enrichment for UP-regulated genes ----------------------------

ego_list_UP <- lapply(names(up_list), function(cond) {
  enrichGO(
    gene          = up_list[[cond]],
    OrgDb         = org.Mm.eg.db,     # Change to org.Hs.eg.db if human
    keyType       = "SYMBOL",         # Change if your IDs are Entrez, etc.
    ont           = "BP",             # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
})
names(ego_list_UP) <- names(up_list)


# Combine enrichment results into a single data frame
go_df <- bind_rows(lapply(names(ego_list_UP), function(cond) {
  ego <- ego_list_UP[[cond]]
  print(ego)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  ego@result %>%
    mutate(Condition = cond) %>%
    select(Condition, ID, Description, pvalue, p.adjust, Count, GeneRatio)
}))

# Select top 20 GO terms per condition based on adjusted p-value
top_go_df <- go_df %>%
  group_by(Condition) %>%
  slice_min(order_by = p.adjust, n = 20) %>%
  ungroup()

# Plot the enrichment results
ggplot(top_go_df, aes(x = Condition, y = fct_reorder(Description, p.adjust))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  # Color gradient: from blue (low p) to red (high p)
  scale_color_gradientn(
    colors = c("#d73027", "#4575b4"),
    name = "p.adjust",
    limits = c(0, max(top_go_df$p.adjust, na.rm = TRUE)),
    oob = scales::squish
  ) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "GO Enrichment (BP) for Upregulated Genes"
  )


# Perform GO enrichment for DOWN-regulated genes --------------------------

ego_list_DOWN <- lapply(names(down_list), function(cond) {
  enrichGO(
    gene          = down_list[[cond]],
    OrgDb         = org.Mm.eg.db,     # Use org.Hs.eg.db if human
    keyType       = "SYMBOL",         # Change if your IDs are Entrez, etc.
    ont           = "BP",             # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
})
names(ego_list_DOWN) <- names(down_list)

# Combine enrichment results into a single data frame
go_df_down <- bind_rows(lapply(names(ego_list_DOWN), function(cond) {
  ego <- ego_list_DOWN[[cond]]
  print(ego)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  ego@result %>%
    mutate(Condition = cond) %>%
    select(Condition, ID, Description, p.adjust, Count, GeneRatio)
}))

# Select top 20 GO terms per condition based on adjusted p-value
top_go_df_down <- go_df_down %>%
  group_by(Condition) %>%
  slice_min(order_by = p.adjust, n = 20) %>%
  ungroup()

# Plot the enrichment results for DOWN-regulated genes
ggplot(top_go_df_down, aes(x = Condition, y = fct_reorder(Description, p.adjust))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradientn(
    colors = c("#4575b4",  "#d73027"),
    name = "p.adjust",
    limits = c(0, max(top_go_df_down$p.adjust, na.rm = TRUE)),
    oob = scales::squish
  ) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "GO Enrichment (BP) for Downregulated Genes"
  )




# after step 3.2 - need to get genes in list from t-test output -------------

# --- From UP_long_ttest ---
up_list <- UP_long_ttest %>%
  group_by(Condition) %>%
  summarise(GeneList = list(Gene), .groups = "drop") %>%
  deframe()  # Convert to named list

# --- From DOWN_long_ttest ---
down_list <- DOWN_long_ttest %>%
  group_by(Condition) %>%
  summarise(GeneList = list(Gene), .groups = "drop") %>%
  deframe()  # Convert to named list

# Perform GO enrichment for UP-regulated genes ----------------------------

ego_list_UP <- lapply(names(up_list), function(cond) {
  enrichGO(
    gene          = up_list[[cond]],
    OrgDb         = org.Mm.eg.db,     # Change to org.Hs.eg.db if human
    keyType       = "SYMBOL",         # Change if your IDs are Entrez, etc.
    ont           = "BP",             # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
})
names(ego_list_UP) <- names(up_list)


# Combine enrichment results into a single data frame
go_df <- bind_rows(lapply(names(ego_list_UP), function(cond) {
  ego <- ego_list_UP[[cond]]
  print(ego)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  ego@result %>%
    mutate(Condition = cond) %>%
    select(Condition, ID, Description, pvalue, p.adjust, Count, GeneRatio)
}))

# Select top 20 GO terms per condition based on adjusted p-value
top_go_df <- go_df %>%
  group_by(Condition) %>%
  slice_min(order_by = p.adjust, n = 20) %>%
  ungroup()

# Plot the enrichment results
ggplot(top_go_df, aes(x = Condition, y = fct_reorder(Description, p.adjust))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  # Color gradient: from blue (low p) to red (high p)
  scale_color_gradientn(
    colors = c("#d73027", "#4575b4"),
    name = "p.adjust",
    limits = c(0, max(top_go_df$p.adjust, na.rm = TRUE)),
    oob = scales::squish
  ) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "GO Enrichment (BP) for Upregulated Genes"
  )
# Perform GO enrichment for DOWN-regulated genes --------------------------

ego_list_DOWN <- lapply(names(down_list), function(cond) {
  enrichGO(
    gene          = down_list[[cond]],
    OrgDb         = org.Mm.eg.db,     # Use org.Hs.eg.db if human
    keyType       = "SYMBOL",         # Change if your IDs are Entrez, etc.
    ont           = "BP",             # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
})
names(ego_list_DOWN) <- names(down_list)

# Combine enrichment results into a single data frame
go_df_down <- bind_rows(lapply(names(ego_list_DOWN), function(cond) {
  ego <- ego_list_DOWN[[cond]]
  print(ego)
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  ego@result %>%
    mutate(Condition = cond) %>%
    select(Condition, ID, Description, p.adjust, Count, GeneRatio)
}))

# Select top 20 GO terms per condition based on adjusted p-value
top_go_df_down <- go_df_down %>%
  group_by(Condition) %>%
  slice_min(order_by = p.adjust, n = 20) %>%
  ungroup()

# Plot the enrichment results for DOWN-regulated genes
ggplot(top_go_df_down, aes(x = Condition, y = fct_reorder(Description, p.adjust))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradientn(
    colors = c("#4575b4",  "#d73027"),
    name = "p.adjust",
    limits = c(0, max(top_go_df_down$p.adjust, na.rm = TRUE)),
    oob = scales::squish
  ) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "GO Enrichment (BP) for Downregulated Genes"
  )
