#!/bin/bash


#$ -l h_vmem=30G
#$ -l h_rt=2:00:00

#$ -cwd
#$ -j y


bm=$1

source /broad/software/scripts/useuse
use .r-3.6.0

ewisdir=../data/processed/ewis
sensdir=../data/processed/sensitivity
osca=../../opt/osca_Linux


# Create outcome adjusted for all exposures
R --vanilla <<EOF
p <- data.table::fread("${ewisdir}/ewis_phenos_EUR.csv", data.table=F, stringsAsFactors=F)
all_exposures <- scan("${ewisdir}/ewis_phenotype_list.txt", what=character())
all_exposures <- setdiff(all_exposures, c("age", "sex"))
p <- dplyr::mutate_at(p, all_exposures, ~ifelse(is.na(.), median(., na.rm=TRUE), .))
names(p)[names(p) %in% all_exposures] <- paste0("e", names(p)[names(p) %in% all_exposures])
lm_formula <- as.formula(paste0("${bm}_adj ~", paste0("e", all_exposures, collapse=" + ")))
p[["${bm}_allExpAdj"]] <- resid(lm(lm_formula, data=p, na.action=na.exclude))
readr::write_tsv(p[, c("id", "id", "${bm}_allExpAdj")], "${sensdir}/${bm}_allExpAdj_vqtl_phenos.csv")
EOF
