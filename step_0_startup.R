library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(GSEABase)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggsankey)
library(rio)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(pheatmap)
library(limma)
library(UpSetR)
library(forcats)
library(openxlsx)
library(patchwork)


# set diUpSetR# set dir
output_dir <- "C:/Users/wangy9/Desktop/lung_data_analysis/output/"
setwd("C:/Users/wangy9/Desktop/lung_data_analysis/")

# load functions
reorder_conditions <- function(df, order_prefix) {
  cols <- colnames(df)
  
  # Handle columns without "_"
  prefixes <- ifelse(grepl("_", cols), sub("_.*", "", cols), cols)
  
  # Get ordering index
  ordering <- match(prefixes, order_prefix)
  
  # Check for unmatched prefixes
  if (any(is.na(ordering))) {
    warning("Some column prefixes not found in 'order_prefix'. They will be placed at the end.")
  }
  
  # Build a sortable data frame
  sorting_df <- data.frame(
    colname = cols,
    prefix = prefixes,
    order_index = ordering,
    stringsAsFactors = FALSE
  )
  
  # Sort by custom prefix order, then by column name
  sorting_df <- sorting_df[order(sorting_df$order_index, sorting_df$colname), ]
  
  # Reorder and return the data frame
  df[, sorting_df$colname]
}

get_col_indices <- function(df) {
  # Ensure input is a data frame
  if (!is.data.frame(df)) {
    stop("Input must be a data frame.")
  }
  
  # Create matrix with column numbers and names
  result <- cbind(
    Column_Number = seq_along(df),
    Column_Name = names(df)
  )
  
  return(result)
}

generate_pass_cutoff_pool <- function(df, X) {
  # Initialize an empty data frame to store results
  Enrichment_pool <- data.frame(UniprotID = character(), stringsAsFactors = FALSE)
  
  for (i in X) {
    # Extract and split the UniprotID column for the current row
    Condition_i <- data.frame(
      UniprotID = unlist(strsplit(df[i, 2], split = ",", fixed = TRUE)), 
      stringsAsFactors = FALSE
    )
    
    # Store each Condition in the global environment with the condition name
    new_df_name <- as.character(df[i, 1])
    assign(new_df_name, Condition_i)
    
    # Merge to maintain a unique list of UniprotIDs
    Enrichment_pool <- merge(Enrichment_pool, Condition_i, by = "UniprotID", all = TRUE)
  }
  
  return(Enrichment_pool)
}

generate_meta_data <- function(s) {
  tok <- sub("_.*", "", s)                 # strip trailing _1/_2...
  # detect label flag
  if (grepl("NL$", tok)) {
    Label <- "NL";  CondTok <- sub("NL$", "", tok)
  } else if (grepl("L$", tok) || tok == "SL") { # SL is saline-labeled alias
    Label <- "L";   CondTok <- sub("L$", "", tok)
  } else {
    stop(paste("Cannot parse L/NL from sample name:", s))
  }
  # normalize saline aliases to "Saline"
  if (CondTok %in% c("S","Sa","Sal","SL","Saline")) CondTok <- "Saline"
  data.frame(SampleID = s, Condition = CondTok, LabelStatus = Label, stringsAsFactors = FALSE)
}

perform_stats_analysis_by_t_test <- function(df, drug_name, drug_cols, control_cols) {
  
  p_values <- apply(df, 1, function(row) {
    t.test(as.numeric(row[drug_cols]), as.numeric(row[control_cols]), var.equal = FALSE)$p.value
  })
  
  log2FC <- apply(df, 1, function(row) {
    mean(as.numeric(row[drug_cols])) - mean(as.numeric(row[control_cols]))
  })
  
  q_values <- p.adjust(p_values, method = "BH")
  
  results_df <- data.frame(
    Genes = rownames(df))
  
  results_df[[paste0(drug_name, "_p_values")]] <- p_values
  results_df[[paste0(drug_name, "_Log2FC")]] <- log2FC
  results_df[[paste0(drug_name, "_q_values")]] <- q_values
  
  results_df[[paste0(drug_name, "_significance")]] <- ifelse(p_values < 0.05, "significant", "")        #if want significance by p_value, change it to p_values
  
  results_df[[paste0(drug_name, "_vs_Saline_threshold=1")]] <- ifelse(q_values < 0.05 & log2FC > 1, "UP",
                                                                      ifelse(q_values < 0.05 & log2FC < -1, "DOWN", ""))
  
  mean_fc <- mean(log2FC)
  sd_fc <- sd(log2FC)
  
  # Assign SD category symmetrically (both sides of the mean)
  results_df[[paste0(drug_name, "_SD_category")]] <- ifelse(abs(log2FC - mean_fc) > 3*sd_fc, "beyond 3SD",
                                                            ifelse(abs(log2FC - mean_fc) > 2*sd_fc, "beyond 2SD",
                                                                   ifelse(abs(log2FC - mean_fc) > 1*sd_fc, "beyond 1SD", "")))
  
  # Assign Direction only if it is beyond threshold
  results_df[[paste0(drug_name, "_Direction")]] <- ifelse(log2FC > mean_fc & results_df[[paste0(drug_name, "_SD_category")]] != "" 
                                                          & results_df[[paste0(drug_name, "_significance")]] == "significant", "UP",
                                                          ifelse(log2FC < mean_fc & results_df[[paste0(drug_name, "_SD_category")]] != ""
                                                                 & results_df[[paste0(drug_name, "_significance")]] == "significant", "DOWN", ""))
  
  return(results_df)
}
