#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y

biomarker=$1

METAL=../opt/generic-metal/metal

source /broad/software/scripts/useuse
use .r-3.6.0


# vQTL meta-analysis
# Note: in practice, the P-value != 0 METAL command only removes variants where beta/SE are missing

VQTL_DIR=../data/processed/vqtl_ss
mkdir -p ${VQTL_DIR}

echo "Performing vQTL meta-analysis for ${biomarker}..."

echo """
SCHEME SAMPLESIZE
MARKER SNP 
ALLELE A1 A2
WEIGHT NMISS
EFFECT beta 
STDERR se
PVALUE P 

ADDFILTER P != 0
ADDFILTER freq > 0.01
ADDFILTER freq < 0.99

PROCESS ${VQTL_DIR}/${biomarker}_EUR_vqtl_merged 
PROCESS ${VQTL_DIR}/${biomarker}_AFR_vqtl_merged 
PROCESS ${VQTL_DIR}/${biomarker}_EAS_vqtl_merged 
PROCESS ${VQTL_DIR}/${biomarker}_SAS_vqtl_merged 

OUTFILE ${VQTL_DIR}/metal/${biomarker}_MA_ .tbl
ANALYZE HETEROGENEITY

QUIT
""" > tmp_metal_script_${biomarker}.txt

$METAL tmp_metal_script_${biomarker}.txt
rm tmp_metal_script_${biomarker}.txt

RES_FILE=${VQTL_DIR}/metal/${biomarker}_MA_1.tbl
Rscript 3b_postprocess_vqtl.R ${RES_FILE}


# Main-effect meta-analysis

ME_DIR=../data/processed/main_effect_ss
mkdir -p ${ME_DIR}

echo "Performing main-effect meta-analysis for ${biomarker}..."

echo """
SCHEME STDERR
MARKER ID 
ALLELE ALT REF
EFFECT BETA 
STDERR SE
PVALUE P 

PROCESS ${ME_DIR}/${biomarker}_EUR_ME_merged 
PROCESS ${ME_DIR}/${biomarker}_AFR_ME_merged 
PROCESS ${ME_DIR}/${biomarker}_EAS_ME_merged 
PROCESS ${ME_DIR}/${biomarker}_SAS_ME_merged 

OUTFILE ${ME_DIR}/metal/${biomarker}_MA_ .tbl
ANALYZE HETEROGENEITY

QUIT
""" > tmp_metal_script_${biomarker}.txt

$METAL tmp_metal_script_${biomarker}.txt
rm tmp_metal_script_${biomarker}.txt

RES_FILE=${ME_DIR}/metal/${biomarker}_MA_1.tbl
Rscript 3b_postprocess_vqtl.R ${RES_FILE}
