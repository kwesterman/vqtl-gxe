library(data.table)
library(tidyverse)


biomarkers <- scan("../data/processed/metabolic_biomarkers.txt", what=character())
exposures <- scan("../data/processed/diet_ewis/24hr_exposures.txt", what=character())


## EUROPEAN-SPECIFIC

read_gxe_anc <- function(bm, e, anc) {
  fn <- paste0("../data/processed/diet_ewis/", bm, "/", e, "_", anc)
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

ancestries <- c("EUR")

anc_ewis_df <- expand_grid(
  bm = biomarkers,
  exposure = exposures,
  anc = ancestries
) %>%
  mutate(res = pmap(list(.$bm, .$exposure, .$anc), read_gxe_anc)) %>%
  unnest(res)

anc_ewis_df %>%
  filter(P < 0.05) %>%  # Redundant
  write_csv("../data/processed/diet_ewis/diet_ewis_EUR_results_nom.csv")
