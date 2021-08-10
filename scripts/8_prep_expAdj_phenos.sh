#!/bin/bash


#$ -l h_vmem=30G
#$ -l h_rt=2:00:00

#$ -cwd
#$ -j y


bm=$1
anc=$2

source /broad/software/scripts/useuse
use .r-3.6.0

ewisdir=../data/processed/ewis
sensdir=../data/processed/sensitivity
osca=../../opt/osca_Linux


# Create outcome adjusted for all exposures
R --vanilla <<EOF
p <- data.table::fread("${ewisdir}/ewis_phenos_${anc}.csv", data.table=F, stringsAsFactors=F)
sig_exposures <- scan("${ewisdir}/significant_exposures.txt", what=character())
sig_exposures <- setdiff(sig_exposures, c("age", "sex"))  # Because these have already been adjusted
p <- dplyr::mutate_at(p, sig_exposures, ~ifelse(is.na(.), median(., na.rm=TRUE), .))
names(p)[names(p) %in% sig_exposures] <- paste0("e", names(p)[names(p) %in% sig_exposures])
lm_formula <- as.formula(paste0("${bm}_adj ~", paste0("e", sig_exposures, collapse=" + ")))
p[["${bm}_adj_expAdj"]] <- resid(lm(lm_formula, data=p, na.action=na.exclude))
readr::write_tsv(p[, c("id", "id", "${bm}_adj_expAdj")], "${sensdir}/${bm}_expAdj_vqtl_phenos_${anc}.csv")
EOF
