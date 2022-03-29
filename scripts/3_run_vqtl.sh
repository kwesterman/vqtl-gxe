#!/bin/sh


#$ -l h_vmem=22G
#$ -l h_rt=28:00:00

#$ -pe smp 6
#$ -binding linear:6
#$ -R y

#$ -cwd
#$ -j y


pheno=$1
ancestry=$2
chr=$3
threads=$4

genodir=/broad/hptmp/kwesterm/ukb_bfiles


cd ~/kw/ukbb-vqtl/scripts


# Create phenotype file with OSCA-specific format
source /broad/software/scripts/useuse
use .r-3.6.0
R --slave <<EOF
library(tidyverse)
read_csv("../data/processed/vqtl_phenos_${ancestry}.csv") %>%
  mutate(FID=id, IID=id) %>% 
  select(FID, IID, ${pheno}_adj) %>% 
  write_tsv("${pheno}_${ancestry}_chr${chr}_tmp.pheno")
EOF

## Generate list of unique SNPs from .bim files
#awk '{print $2}' ${genodir}/chr${chr}.bim | sort | uniq > ${genodir}/chr${chr}_unique_snps.txt

bfile=${genodir}/chr${chr}_${ancestry}
if [ ${ancestry} == "EUR" ] || [ ${ancestry} == "all" ]
then
        bfile=${genodir}/chr${chr}
fi

# Run vQTL
osca=../../opt/osca_Linux
${osca} \
	--vqtl \
	--bfile ${bfile} \
	--pheno ${pheno}_${ancestry}_chr${chr}_tmp.pheno \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--thread-num $threads \
	--out ../data/processed/vqtl_ss/${pheno}_${ancestry}_chr${chr}
	#--extract-snp ${genodir}/chr${chr}_unique_snps.txt \

# Run analogous main-effect GWAS
plink2=../../opt/plink2
${plink2} \
	--bfile ${bfile} \
	--pheno ${pheno}_${ancestry}_chr${chr}_tmp.pheno \
	--pheno-name ${pheno}_adj \
	--maf 0.01 \
	--glm omit-ref cols=+a1freq \
	--out ../data/processed/main_effect_ss/${pheno}_${ancestry}_chr${chr}

rm ${pheno}_${ancestry}_chr${chr}_tmp.pheno
