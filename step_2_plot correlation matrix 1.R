
# prepare the color and condition/category --------------------------------
cols <- colnames(dat_numeric_reordered)
base <- sub("_[0-9]+$", "", cols)
cond_code <- sub("(NL|L)$", "", base)
cond_lookup <- c(
  "S"   = "Saline",
  "FLA" = "FLA",
  "LPS" = "LPS",
  "PAM" = "PAM",
  "R848"= "R848",
  "TNF" = "TNF",
  "POLY"= "POLY",
  "GAMP"= "GAMP"
)
cat_lab <- ifelse(grepl("NL$", base), "non-labeled", "labeled")

ann_col <- data.frame(
  Condition = factor(cond_lookup[cond_code], levels = c("Saline", "FLA", "LPS", "PAM", "R848", "TNF", "POLY", "GAMP")),
  Category  = cat_lab,
  row.names = cols,
  check.names = FALSE
)

my_cols <- c(Saline = "grey80", FLA = "#A51C36", LPS = "#7ABBDB", PAM = "#84BA42", R848 = "#682487", TNF = "#DBB428", POLY = "#D4562E", GAMP = "#4485C7")

cat_cols <- c("labeled" = "black", "non-labeled" = "grey60")

ann_colors <- list(
  Condition = my_cols,
  Category  = cat_cols
)

pal <- colorRampPalette(c("#2b8cbe", "grey90", "#d7301f"))(101)  
gaps_col <- c(4, 16, 28, 32)


# plot P1-4 before peeling as input --------------------------------------------------------------

cors1 <- cor(
  dat_numeric_reordered,     
  method = "pearson"
)

pheatmap(
  cors1,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P1"
)

cors2 <- cor(
  dat_numeric_reordered,    
  method = "spearman"
)

pheatmap(
  cors2,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P2"
)

cors3 <- cor(
  dat_numeric_reordered_scaled,  
  method = "pearson"
)

pheatmap(
  cors3,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P3"
)

cors4 <- cor(
  dat_numeric_reordered_scaled,  
  method = "spearman"
)

pheatmap(
  cors4,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P4"
)

# plot P5-8 after peeling as input -----------------------------------------

cors5 <- cor(
  post_peeling_cutoff,     
  method = "pearson"
)

pheatmap(
  cors5,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P5"
)

cors6 <- cor(
  post_peeling_cutoff,     
  method = "spearman"
)

pheatmap(
  cors6,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P6"
)

cors7 <- cor(
  post_peeling_cutoff_scaled,     
  method = "pearson"
)

pheatmap(
  cors7,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P7"
)

cors8 <- cor(
  post_peeling_cutoff_scaled,     
  method = "spearman"
)

pheatmap(
  cors8,
  color             = pal,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = gaps_col,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_row      = 8,
  legend            = TRUE,
  main              = "P8"
)
