library(data.table)
library(tidyverse)


biomarkers <- scan("../data/processed/metabolic_biomarkers.txt", what=character())
exposures <- c("PC1", "PC2")


# META-ANALYSIS

read_gxe_anthroPC <- function(bm, e) {
  fn <- paste0("../data/processed/ewis/", bm, "/anthro", e, "_EUR")
  df <- tryCatch({
    fread(fn, data.table=F, stringsAsFactors=F, na.strings=c("NA", "nan", "nane-nan", "")) %>%
      select(SNP=SNPID, Non_Effect_Allele, Effect_Allele, 
	     Beta=paste0("Beta_G-", e), Variance=paste0("Var_Beta_G-", e), P="P_Value_Interaction") %>%
      mutate(SE=sqrt(Variance), P=as.numeric(P))
  }, error=function (e) tibble())
  if (e == exposures[1]) print(bm)
  if (nrow(df) > 1) df else NULL
}


ewis_df <- expand_grid(
  bm = biomarkers,
  exposure = exposures
) %>%
  mutate(res = pmap(list(.$bm, .$exposure), read_gxe_anthroPC)) %>%
  unnest(res)

ewis_df %>%
  write_csv("../data/processed/ewis/ewis_ma_results_anthroPCs.csv")

ewis_df %>%
  filter(P < 0.05) %>%
  write_csv("../data/processed/ewis/ewis_ma_results_anthroPCs_nom.csv")
