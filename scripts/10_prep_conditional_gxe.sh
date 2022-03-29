#!/bin/bash


#$ -l h_vmem=40G
#$ -l h_rt=1:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use .r-3.6.0


plink2=../../opt/plink2

ewis_dir=../data/processed/ewis
ca_dir=${ewis_dir}/conditional_analysis

mkdir -p ${ca_dir}


# Create genotype matrix with only significant variants
${plink2} \
	--pfile ${ewis_dir}/ewis_genotypes \
	--extract ${ewis_dir}/significant_variants.txt \
	--export A ref-first \
	--out ${ca_dir}/sig_ewis_genotypes

# Create phenotype file with all necessary variants and phenotypes
R --vanilla <<EOF
library(tidyverse)

geno_df <- read_tsv("${ca_dir}/sig_ewis_genotypes.raw") %>%
  select(id=IID, 7:ncol(.))
names(geno_df) <- sub("(_[ACGT]*$)", "", names(geno_df))

metabolic_biomarkers <- scan("../data/processed/metabolic_biomarkers.txt", what=character())
pheno_df <- read_csv("../data/processed/vqtl_phenos_EUR.csv") %>%
  select(id, paste0(metabolic_biomarkers, "_adj"))

sig_exposures <- scan("${ewis_dir}/significant_exposures.txt", what=character())
exposure_df <- data.table::fread("${ewis_dir}/ewis_phenos_EUR.csv", 
  select=c("id", sig_exposures),
  data.table=F, stringsAsFactors=F) %>%
  rename_at(vars(-id), ~paste0("e", .))

ca_df <- geno_df %>%
  inner_join(pheno_df, by="id") %>%
  inner_join(exposure_df, by="id")
saveRDS(ca_df, "${ca_dir}/conditional_analysis_dataset.rds")
EOF
