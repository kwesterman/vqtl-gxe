#!/bin/bash


cd ~/kw/ukbb-vqtl/scripts
source /broad/software/scripts/useuse
use .r-3.6.0

R --vanilla <<EOF

library(tidyverse)
library(data.table)


afr_phenos <- fread("../data/processed/vqtl_phenos_AFR.csv",
                    data.table=F, stringsAsFactors=F)
eas_phenos <- fread("../data/processed/vqtl_phenos_EAS.csv",
                    data.table=F, stringsAsFactors=F)
eur_phenos <- fread("../data/processed/vqtl_phenos_EUR.csv",
                    data.table=F, stringsAsFactors=F)
sas_phenos <- fread("../data/processed/vqtl_phenos_SAS.csv",
                    data.table=F, stringsAsFactors=F)
all_phenos <- bind_rows(afr_phenos, eas_phenos, eur_phenos, sas_phenos)

metabolic_biomarkers <- scan("../data/processed/metabolic_biomarkers.txt", what=character())

bm_pca <- all_phenos %>%
  select(all_of(paste0(metabolic_biomarkers, "_adj"))) %>%
  mutate_all(~ifelse(is.na(.), mean(., na.rm=T), .)) %>%
  prcomp(scale.=T)
bm_eigenvals <- bm_pca[["sdev"]] ** 2
n_eff_bm <- sum(bm_eigenvals) ** 2 / sum(bm_eigenvals ** 2)
print(paste0("There are ", round(n_eff_bm, 1), " effective biomarkers."))

afr_exposures <- fread("../data/processed/ewis/ewis_phenos_AFR.csv", data.table=F, stringsAsFactors=F)
eas_exposures <- fread("../data/processed/ewis/ewis_phenos_EAS.csv", data.table=F, stringsAsFactors=F)
eur_exposures <- fread("../data/processed/ewis/ewis_phenos_EUR.csv", data.table=F, stringsAsFactors=F)
sas_exposures <- fread("../data/processed/ewis/ewis_phenos_SAS.csv", data.table=F, stringsAsFactors=F)

all_exposures <- bind_rows(afr_exposures, eas_exposures, eur_exposures, sas_exposures)

exp_pca <- all_exposures %>%
  select(-id, -contains("_adj")) %>%
  mutate_all(~ifelse(is.na(.), mean(., na.rm=T), .)) %>% 
  as.matrix() %>%
  prcomp(scale.=T)
exp_eigenvals <- exp_pca[["sdev"]] ** 2
n_eff_exp <- sum(exp_eigenvals) ** 2 / sum(exp_eigenvals ** 2)
print(paste0("There are ", round(n_eff_exp, 1), " effective exposures."))

EOF
