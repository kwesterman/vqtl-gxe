#!/bin/bash


biomarker=$1

cd ~/kw/ukbb-vqtl/scripts

METAL=../opt/generic-metal/metal

source /broad/software/scripts/useuse
use .r-3.6.0


# vQTL meta-analysis

VQTL_DIR=../data/processed/vqtl_ss

echo "Performing vQTL meta-analysis for ${biomarker}..."

echo """
SCHEME STDERR
MARKER SNP 
ALLELE A2 A1 
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

echo "Performing main-effect meta-analysis for ${biomarker}..."

echo """
SCHEME STDERR
MARKER ID 
ALLELE REF ALT
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
