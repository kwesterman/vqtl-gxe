#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use .r-3.6.0


biomarker=$1

cd ~/florez_ukb_projects/ukbb-vqtl/scripts


# Concatenate chromosomes

for anc in EUR AFR EAS SAS; do
sumstats_root=../data/processed/sensitivity/${biomarker}
head -1 ${sumstats_root}_chr1_${anc}_expAdj.vqtl > ${sumstats_root}_${anc}_expAdj_merged
for chr in {1..22}; do
	echo "Chromosome " ${chr}
	tail -n +2 ${sumstats_root}_chr${chr}_${anc}_expAdj.vqtl >> ${sumstats_root}_${anc}_expAdj_merged
done
done


# vQTL meta-analysis

METAL=../opt/generic-metal/metal

SENS_DIR=../data/processed/sensitivity

echo "Performing vQTL meta-analysis for ${biomarker}..."

echo """
SCHEME SAMPLESIZE
MARKER SNP 
ALLELE A2 A1 
WEIGHT NMISS
EFFECT beta 
PVALUE P 

ADDFILTER P != 0
ADDFILTER freq > 0.01
ADDFILTER freq < 0.99

PROCESS ${SENS_DIR}/${biomarker}_EUR_expAdj_merged 
PROCESS ${SENS_DIR}/${biomarker}_AFR_expAdj_merged 
PROCESS ${SENS_DIR}/${biomarker}_EAS_expAdj_merged 
PROCESS ${SENS_DIR}/${biomarker}_SAS_expAdj_merged 

OUTFILE ${SENS_DIR}/metal/${biomarker}_expAdj_MA_ .tbl
ANALYZE HETEROGENEITY

QUIT
""" > tmp_metal_script_${biomarker}.txt

$METAL tmp_metal_script_${biomarker}.txt
rm tmp_metal_script_${biomarker}.txt

RES_FILE=${SENS_DIR}/metal/${biomarker}_MA_1.tbl
Rscript 3b_postprocess_vqtl.R ${RES_FILE}
