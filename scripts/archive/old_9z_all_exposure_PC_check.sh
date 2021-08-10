#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=12:00:00

#$ -cwd
#$ -j y

#$ -pe smp 8
#$ -binding linear:8
#$ -R y


bm=$1

source /broad/software/scripts/useuse
use .r-3.6.0

ewisdir=../data/processed/ewis
sensdir=../data/processed/sensitivity
osca=../../opt/osca_Linux


# Create outcome adjusted for all exposures
R --vanilla <<EOF
library(tidyverse)
p <- data.table::fread("${ewisdir}/ewis_phenos_EUR.csv", data.table=F, stringsAsFactors=F)
all_exposures <- scan("${ewisdir}/ewis_phenotype_list.txt", what=character())
all_exposures <- setdiff(all_exposures, c("age", "sex"))
#p <- dplyr::mutate_at(p, all_exposures, ~ifelse(is.na(.), median(., na.rm=TRUE), .))
exp_pca <- p %>%
  select(seq(which(names(p) == "age"), ncol(p))) %>%
  mutate_all(~ifelse(is.na(.), mean(., na.rm=T), .)) %>% 
  as.matrix() %>%
  prcomp(scale.=T)
saveRDS(exp_pca, "${ewisdir}/exposure_PCA_EUR.rds")
exp_pca <- readRDS("${ewisdir}/exposure_PCA_EUR.rds")
names(exp_pca[["x"]]) <- paste0("e", names(exp_pca[["x"]]))
p <- bind_cols(p, as.data.frame(exp_pca[["x"]][, 1:200]))
lm_formula <- as.formula(paste0("${bm}_adj ~", paste(names(exp_pca[["x"]])[1:200], collapse=" + ")))
p[["${bm}_allExpPcAdj"]] <- resid(lm(lm_formula, data=p, na.action=na.exclude))
readr::write_tsv(p[, c("id", "id", "${bm}_allExpPcAdj")], "${sensdir}/${bm}_allExpPcAdj_vqtl_phenos.csv")
EOF


# Test vQTL with all-exposure-adjusted outcome
${osca} \
	--vqtl \
	--bfile /broad/hptmp/kwesterm/ukb_bfiles/chr22 \
	--pheno ${sensdir}/${bm}_allExpPcAdj_vqtl_phenos.csv \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--thread-num 4 \
	--out ${sensdir}/${bm}_allExpPcAdj


# Join primary and exposure-adjusted vQTL tests for comparison
R --vanilla <<EOF
library(tidyverse)
primary_vqtl <- read_tsv("../data/processed/vqtl_ss/${bm}_EUR_chr22.vqtl", col_types=cols_only(SNP="c", P="d"))
adj_vqtl <- read_tsv("${sensdir}/${bm}_allExpPcAdj.vqtl", col_types=cols_only(SNP="c", P="d"))
merged <- inner_join(primary_vqtl, adj_vqtl, by="SNP", suffix=c(".primary", ".allExpPcAdj"))
write_tsv(merged, "${sensdir}/${bm}_allExpPcAdj_comparison.txt")
EOF

rm ../data/processed/sensitivity/${bm}_allExpPcAdj_vqtl_phenos.csv
