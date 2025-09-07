# load data
# all imputated values
dat <- import("C:/Users/wangy9/Desktop/MS_data/20250820_oneSaline/perseus_output.xlsx", sheet = "imputation")
dat_numeric <- dat%>%
  dplyr::select(1:48)
rownames(dat_numeric) <- dat$Genes

# re-order the data
my_order <- c("SL", "FLAL", "LPSL", "PAML", "R848L", "TNFL", "POLYL", "GAMPL", "SNL","FLANL", "LPSNL", "PAMNL", "R848NL", "TNFNL", "POLYNL", "GAMPNL")
dat_numeric_reordered <- as.matrix(reorder_conditions(dat_numeric, my_order))

dat_numeric_reordered_scaled <- t(scale(t(dat_numeric_reordered)))



# post peeling cutoff list
peeling_output <- import("C:/Users/wangy9/Desktop/MS_data/20250820_oneSaline/peeling input_NL mean/peeling_output_tol=0.xlsx", header = FALSE)
Enrichment_pool <- generate_pass_cutoff_pool(peeling_output, 1:8)
post_peeling_cutoff <- dat %>%
  filter(UniprotID %in% unlist(Enrichment_pool))%>%
  dplyr::select(c(50,51), c(1:48))
rownames(post_peeling_cutoff) <- post_peeling_cutoff$Genes
post_peeling_cutoff <- post_peeling_cutoff[,c(-1:-2)]

# re-order the data
my_order <- c("SL", "FLAL", "LPSL", "PAML", "R848L", "TNFL", "POLYL", "GAMPL", "SNL","FLANL", "LPSNL", "PAMNL", "R848NL", "TNFNL", "POLYNL", "GAMPNL")
post_peeling_cutoff <- reorder_conditions(post_peeling_cutoff, my_order)

post_peeling_cutoff_scaled <- t(scale(t(dat_numeric_reordered)))

post_peeling_cutoff_subset <- post_peeling_cutoff %>% 
  dplyr::select(1:32)
