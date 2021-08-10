#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00

#$ -cwd
#$ -j y


bm=$1


source /broad/software/scripts/useuse
use .r-3.6.0

sensdir=../data/processed/sensitivity

head -1 ${sensdir}/${bm}_chr1_allExpAdj.vqtl > ${sensdir}/${bm}_allExpAdj_vqtl_merged
for chr in {1..22}; do
	echo "Chromosome " ${chr}
	tail -n +2 ${sensdir}/${bm}_chr${chr}_allExpAdj.vqtl >> ${sensdir}/${bm}_allExpAdj_vqtl_merged
done

# Join primary and exposure-adjusted vQTL tests for comparison
R --vanilla <<EOF
library(tidyverse)
primary_vqtl <- read_tsv("../data/processed/vqtl_ss/${bm}_EUR_vqtl_merged", col_types=cols_only(SNP="c", P="d"))
adj_vqtl <- read_tsv("${sensdir}/${bm}_allExpAdj_vqtl_merged", col_types=cols_only(SNP="c", P="d"))
merged <- inner_join(primary_vqtl, adj_vqtl, by="SNP", suffix=c(".primary", ".allExpAdj"))
write_tsv(merged, "${sensdir}/${bm}_allExpAdj_comparison.txt")
EOF
