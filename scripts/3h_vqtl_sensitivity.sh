#!/bin/bash


pheno=$1


cd ~/kw/ukbb-vqtl/scripts


ewis_dir=../data/processed/ewis


## Run European vQTL test using log-transformed biomarker values
source /broad/software/scripts/useuse
use .r-3.6.0
ancestry=EUR

R --slave <<EOF
library(tidyverse)
INT <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))
read_csv("../data/processed/vqtl_phenos_${ancestry}.csv") %>%
  mutate(FID=id, IID=id) %>% 
  #mutate(${pheno}_adj_log = log(${pheno}_adj + 20)) %>%
  mutate(${pheno}_adj_int = INT(${pheno}_adj)) %>%
  select(FID, IID, ${pheno}_adj_log) %>% 
  write_tsv("${pheno}_${ancestry}_vqtl_sensitivity.pheno")
EOF

osca=../../opt/osca_Linux
bfile=${ewis_dir}/ewis_genotypes
${osca} \
	--vqtl \
	--bfile ${bfile} \
	--pheno ${pheno}_${ancestry}_vqtl_sensitivity.pheno \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--out ../data/processed/sensitivity/${pheno}_${ancestry}_INT_vqtl

rm ${pheno}_${ancestry}_vqtl_sensitivity.pheno
