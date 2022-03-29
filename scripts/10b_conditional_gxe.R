library(tidyverse)

ca_df <- readRDS("../data/processed/ewis/conditional_analysis/conditional_analysis_dataset.rds")

sig_combos <- read_csv("../data/processed/ewis/significant_variant_exposure_combos.csv")
sig_threshold <- scan("../data/processed/ewis/significance_threshold.txt", what=double())

test_gxe <- function(exposure, bm, SNP, df, adj_exposures=c()) {
  bm_adj <- paste0(bm, "_adj")
  g <- df[[SNP]]
  e <- df[[paste0("e", exposure)]]  # Need proper variable name starting with letter
  lm_str <- paste0(bm_adj, " ~ g * e")
  if (length(adj_exposures) > 0) {
    adj_exposure_fields <- paste0("e", adj_exposures)
    lm_str <- paste(lm_str, "+", paste0("g * ", adj_exposure_fields, collapse=" + "))
  }
  pval <- tryCatch({
    lm_res <- broom::tidy(lm(as.formula(lm_str), data=df))
    filter(lm_res, term == "g:e")$p.value
  }, error = function(e) as.numeric(NA))
  if (length(pval) == 1) pval else as.numeric(NA) 
}

run_conditional_tests <- function(all_exposures, bm, SNP, p_threshold) {
  remaining_exposures <- all_exposures
  adj_exposures <- c()
  while(length(remaining_exposures) > 0) {
    pvals <- map_dbl(remaining_exposures, test_gxe, bm, SNP, ca_df, adj_exposures)
    if (any(pvals < p_threshold, na.rm=TRUE)) {
      adj_exposures <- c(adj_exposures, remaining_exposures[which.min(pvals)])
    }
    nonsig_exposures <- remaining_exposures[which(
      (pvals > p_threshold) | is.na(pvals)
    )]
    remaining_exposures <- setdiff(remaining_exposures, 
                                   c(adj_exposures, nonsig_exposures))
  }
  if (length(adj_exposures) > 0) adj_exposures else "none"
}

indep_res <- sig_combos %>%
  rowwise() %>%
  mutate(independent = list(run_conditional_tests(
    str_split(sig_exposures, ",")[[1]], bm, SNP, sig_threshold)
  )) %>%
  mutate(independent_nom = list(run_conditional_tests(
    str_split(sig_exposures, ",")[[1]], bm, SNP, 0.05)
  )) %>%
  ungroup() %>%
  mutate(n_sig = lengths(independent),
         n_sig_nom = lengths(independent_nom))

saveRDS(indep_res, "../data/processed/ewis/conditional_analysis/independent_hits.rds")
indep_res %>%
  mutate(indep_exposures = map_chr(independent, paste, collapse=","),
         indep_exposures_nom = map_chr(independent_nom, paste, collapse=",")) %>%
  select(bm, SNP, sig_exposures, 
         indep_exposures, n_sig, 
         indep_exposures_nom, n_sig_nom) %>%
  write_csv("../data/processed/ewis/conditional_analysis/independent_hits.csv")