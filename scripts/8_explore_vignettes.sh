#!/bin/bash


bm=$1  # Ex. tg
snp=$2  # Ex. rs12345
exposure=$3  # Ex. 21001 (corresponds to BMI)
ancestry=$4  # Ex. EUR
raw_exposure=$5  # Ex. f.21001.0.0
raw_phenofile=$6  # Ex. ukb10528.tab


#$ -l os=RedHat7
#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -j y


source /broad/software/scripts/useuse
use .r-3.6.0

cd kw/ukbb-vqtl/scripts

plink2=../../opt/plink2
osca=../../opt/osca_Linux

ewisdir=../data/processed/ewis
tmpfile=${bm}_${snp}_${exposure}_${ancestry}


# Create single-SNP files
echo ${snp} > ${tmpfile}.snp
  # bed/bim/fam (for vQTL test)
${plink2} \
	--pfile ${ewisdir}/ewis_genotypes \
	--extract ${tmpfile}.snp \
	--make-bed \
	--out ${tmpfile}
  # flat file (for reading into R)
${plink2} \
	--pfile ${ewisdir}/ewis_genotypes \
	--extract ${tmpfile}.snp \
	--export A \
	--out ${tmpfile}

# Create exposure-adjusted outcome
R --vanilla <<EOF
p <- data.table::fread("${ewisdir}/ewis_phenos_${ancestry}.csv", 
		       select=c("id", "${bm}_adj", "${exposure}"),
		       data.table=F, stringsAsFactors=F)
p[["exposure"]] <- p[["${exposure}"]]
p[["${bm}_resid"]] <- resid(lm(${bm}_adj ~ exposure, data=p, na.action=na.exclude))
readr::write_tsv(p[, c("id", "id", "${bm}_resid")], "${tmpfile}_vqtl_phenos.csv")
readr::write_csv(p[, c("id", "${bm}_adj", "${bm}_resid", "exposure")], "${tmpfile}_phenos.csv")
EOF

# Test vQTL with exposure-adjusted outcome
${osca} \
	--vqtl \
	--bfile ${tmpfile} \
	--pheno ${tmpfile}_vqtl_phenos.csv \
	--vqtl-mtd 2 \
	--out ../data/processed/vignettes/${tmpfile}

# Various follow-up tests using R
R --vanilla <<EOF
library(tidyverse)
p <- read_csv("${tmpfile}_phenos.csv")
snp_df <- read_table2("${tmpfile}.raw") %>%
  select(id=IID, snp=7)
raw_bm_df <- data.table::fread("../data/processed/vqtl_phenos_${ancestry}.csv", 
		       select=c("id", "${bm}"),
		       data.table=F, stringsAsFactors=F)
raw_exp_df <- data.table::fread("../data/raw/ukb_phenofiles/${raw_phenofile}", 
		       select=c("f.eid", "${raw_exposure}"),
		       data.table=F, stringsAsFactors=F)
raw_exp_df[["raw_exposure"]] <- raw_exp_df[["${raw_exposure}"]]
raw_exp_df <- select(raw_exp_df, id=f.eid, raw_exposure)
df <- inner_join(p, snp_df, by="id") %>%
  inner_join(raw_bm_df, by="id") %>%
  inner_join(raw_exp_df, by="id")
saveRDS(df, "../data/processed/vignettes/${tmpfile}_df.rds")
EOF

rm ${tmpfile}.snp
rm ${tmpfile}_vqtl_phenos.csv
rm ${tmpfile}_phenos.csv
rm ${tmpfile}.log
rm ${tmpfile}.bed
rm ${tmpfile}.bim
rm ${tmpfile}.fam
rm ${tmpfile}.raw
