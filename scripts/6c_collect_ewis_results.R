library(data.table)
library(tidyverse)


biomarkers <- scan("../data/processed/metabolic_biomarkers.txt", what=character())
exposures <- scan("../data/processed/ewis/ewis_phenotype_list.txt", what=character())


# META-ANALYSIS

read_gxe_ma <- function(bm, e) {
  fn <- paste0("../data/processed/ewis/", bm, "/", bm, "_", e, "_MA_1.tbl")
  df <- tryCatch({
    fread(fn, data.table=F, stringsAsFactors=F, na.strings=c("NA", "nan", "nane-nan", "")) %>%
      select(SNP=MarkerName, Allele1, Allele2, 
	     Beta=Effect, SE=StdErr, P="P-value", 
             Dir=Direction, ISq=HetISq, HetP=HetPVal) %>%
      mutate(P=as.numeric(P))
  }, error=function (e) tibble())
  if (e == exposures[1]) print(bm)
  if (nrow(df) > 1) df else NULL
}


ewis_df <- expand_grid(
  bm = biomarkers,
  exposure = exposures
) %>%
  mutate(res = pmap(list(.$bm, .$exposure), read_gxe_ma)) %>%
  unnest(res)

ewis_df %>%
  filter(P < 0.05) %>%
  write_csv("../data/processed/ewis/ewis_ma_results_nom.csv")

ewis_df %>%
  write_csv("../data/processed/ewis/ewis_ma_results.csv")


## ANCESTRY-SPECIFIC

read_gxe_anc <- function(bm, e, anc) {
  fn <- paste0("../data/processed/ewis/", bm, "/", e, "_", anc)
  df <- tryCatch({
    fread(fn, data.table=F, stringsAsFactors=F, na.strings=c("NA", "nan", "nane-nan", "")) %>%
      select(SNP=ID, Allele1=Non_Effect_Allele, Allele2=Effect_Allele, 
	     Beta=Beta_Interaction_1, P=P_Value_Interaction) %>%
      mutate(P=as.numeric(P)) %>%
      filter(P < 0.05)
  }, error=function (e) tibble())
  if (e == exposures[1]) print(paste(bm, anc))
  if (nrow(df) > 1) df else NULL
}

ancestries <- c("EUR", "AFR", "EAS", "SAS")

anc_ewis_df <- expand_grid(
  bm = biomarkers,
  exposure = exposures,
  anc = ancestries
) %>%
  mutate(res = pmap(list(.$bm, .$exposure, .$anc), read_gxe_anc)) %>%
  unnest(res)

anc_ewis_df %>%
  filter(P < 0.05) %>%  # Redundant
  write_csv("../data/processed/ewis/ewis_anc_results_nom.csv")
