#!/bin/sh


biomarker=$1


source /broad/software/scripts/useuse
use .r-3.6.0

cd kw/ukbb-vqtl/scripts

R --vanilla <<EOF
library(tidyverse)
ewis_genotypes_pvar <- read_tsv("../data/processed/ewis/ewis_genotypes.pvar") %>%
  select(1:3) %>%
  setNames(c("CHR", "POS", "ID"))
read_csv("../data/processed/ewis/ewis_variant_data.csv") %>%
  filter(bm == "${biomarker}", ancestry == "MA", P < 5e-8 / 16.2) %>%
  left_join(ewis_genotypes_pvar, by=c("CHR", "POS")) %>%
  select(ID, ALT, BETA) %>%
  write_tsv("${biomarker}_vPRS_weights.tmp")
EOF

plink2=../../opt/plink2
${plink2} \
	--pfile ../data/processed/ewis/ewis_genotypes \
	--score ${biomarker}_vPRS_weights.tmp \
	--out ../data/processed/vprs/${biomarker}_vprs

rm ${biomarker}_vPRS_weights.tmp

R --vanilla <<EOF
library(tidyverse)
phenos <- data.table::fread("../data/processed/ewis/ewis_phenos_EUR.csv", stringsAsFactors=F)
score_df <- read_tsv(paste0("../data/processed/vprs/${biomarker}_vprs.sscore")) %>%
  select(id=IID, score=SCORE1_AVG) %>%
  inner_join(phenos, by="id")
all_exposures <- names(phenos)[!(names(phenos) == "id" | grepl("_adj", names(phenos)))]
vprs_pvals <- map_dbl(all_exposures, function(e) tryCatch({
  score_df[["e"]] <- score_df[[e]]
  vprs_int_lm <- lm(${biomarker}_adj ~ score * e, data=score_df)
  tidy_lm <- broom::tidy(vprs_int_lm)
  tidy_lm[["p.value"]][nrow(tidy_lm)]
}, error = function (e) as.numeric(NA)))
tibble(exposure = all_exposures, P_vprs_int = vprs_pvals) %>%
  write_csv("../data/processed/vprs/${biomarker}_vprs_interactions.csv")
EOF
