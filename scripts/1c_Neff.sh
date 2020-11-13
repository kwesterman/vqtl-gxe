#!/bin/bash


cd ~/kw/ukbb-vqtl/scripts
source /broad/software/scripts/useuse
use .r-3.6.0

R --vanilla <<EOF

library(tidyverse)
library(data.table)


eur_phenos <- fread("../data/processed/vqtl_phenos_EUR.csv",
                    data.table=F, stringsAsFactors=F)

bm_pca <- eur_phenos %>%
  select(contains("_adj")) %>%
  mutate_all(~ifelse(is.na(.), mean(., na.rm=T), .)) %>%
  prcomp(scale.=T)
n_eff_bm <- sum(bm_pca[["sdev"]]) ** 2 / sum(bm_pca[["sdev"]] ** 2)
print(paste0("There are ", round(n_eff_bm, 1), " effective biomarkers."))

eur_exposures <- fread("../data/processed/ewis/ewis_phenos_EUR.csv", data.table=F, stringsAsFactors=F)

exp_pca <- eur_exposures %>%
  select(-id, -contains("_adj")) %>%
  mutate_all(~ifelse(is.na(.), mean(., na.rm=T), .)) %>% 
  as.matrix() %>%
  prcomp(scale.=T)
n_eff_exp <- sum(exp_pca[["sdev"]]) ** 2 / sum(exp_pca[["sdev"]] ** 2)
print(paste0("There are ", round(n_eff_exp, 1), " effective exposures."))

EOF
