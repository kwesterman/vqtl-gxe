#!/bin/bash


biomarker=$1
exposure=$2


#$ -l h_vmem=15G
#$ -l h_rt=01:00:00
#$ -cwd
#$ -j y


ewis_dir=../data/processed/ewis


## Run European vQTL test using log-transformed biomarker values
source /broad/software/scripts/useuse
use .r-3.6.0
ancestry=EUR

R --slave <<EOF
library(tidyverse)
p <- data.table::fread("../data/processed/ewis/ewis_phenos_${ancestry}.csv", 
		       select=c("id", "${biomarker}_adj", "${exposure}"), 
		       data.table=F, stringsAsFactors=F)
p[["biomarker"]] <- p[["${biomarker}_adj"]]
p[["e"]] <- p[["${exposure}"]]
p[["biomarker_resid"]] <- resid(lm(biomarker ~ e, data=p, na.action=na.exclude))
select(p, FID=id, IID=id, biomarker_resid) %>%
  write_tsv("${biomarker}_${ancestry}_vqtl_sensitivity_${exposure}adj.pheno")
EOF

osca=../../opt/osca_Linux
bfile=${ewis_dir}/ewis_genotypes
${osca} \
	--vqtl \
	--bfile ${bfile} \
	--pheno ${biomarker}_${ancestry}_vqtl_sensitivity_${exposure}adj.pheno \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--out ../data/processed/sensitivity/${biomarker}_${ancestry}_${exposure}adj_vqtl

rm ${biomarker}_${ancestry}_vqtl_sensitivity_${exposure}adj.pheno
