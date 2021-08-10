#!/bin/bash


bm=$1
snp=$2
ancestry=$4


#$ -l os=RedHat7
#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -j y


source /broad/software/scripts/useuse
use .r-3.6.0

cd kw/ukbb-vqtl/scripts

plink2=../../opt/plink2

ewisdir=../data/processed/ewis
tmpfile=${bm}_${snp}_${ancestry}


# Create single-SNP file
echo ${snp} > ${tmpfile}.snp
${plink2} \
	--pfile ${ewisdir}/ewis_genotypes \
	--extract ${tmpfile}.snp \
	--export A \
	--out ${tmpfile}

rm ${tmpfile}.snp
