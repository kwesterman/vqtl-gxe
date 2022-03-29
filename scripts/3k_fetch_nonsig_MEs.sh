#!/bin/bash


source /broad/software/scripts/useuse
use .r-3.6.0


R --vanilla <<EOF
library(tidyverse)
primary_vqtl_df <- read_csv("../data/processed/ewis/ewis_variant_data.csv")
ME_res <- map_dfr(1:nrow(primary_vqtl_df), function(i) {
  fn <- paste0("../data/processed/main_effect_ss/metal/", primary_vqtl_df[["bm"]][i], "_MA_1.tbl")
  me_res_str <- system(paste0("grep -w ", primary_vqtl_df[["index_var"]][i], " ", fn, " | cut -f 4,6"), intern=TRUE)
  me_res <- as.numeric(strsplit(me_res_str, split="\t")[[1]])
  tibble(beta_ME = me_res[1], P_ME = me_res[2])	
})
primary_vqtl_df <- bind_cols(primary_vqtl_df, ME_res)
write_csv(primary_vqtl_df, "../data/processed/sensitivity/primary_vqtl_df_withME.csv")
EOF
