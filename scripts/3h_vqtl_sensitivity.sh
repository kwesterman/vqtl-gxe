#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=01:00:00
#$ -cwd
#$ -j y


bm=$1


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
  mutate(${bm}_adj_INT = INT(${bm}_adj)) %>%
  select(FID, IID, ${bm}_adj_INT) %>% 
  write_tsv("${bm}_${ancestry}_vqtl_sensitivity.pheno")
EOF

osca=../../opt/osca_Linux
bfile=${ewis_dir}/ewis_genotypes
${osca} \
	--vqtl \
	--bfile ${bfile} \
	--pheno ${bm}_${ancestry}_vqtl_sensitivity.pheno \
	--vqtl-mtd 2 \
	--out ../data/processed/sensitivity/${bm}_${ancestry}_INT_vqtl

rm ${bm}_${ancestry}_vqtl_sensitivity.pheno
