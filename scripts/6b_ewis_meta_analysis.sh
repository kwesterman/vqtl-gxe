#!/bin/bash


biomarker=$1

cd ~/kw/ukbb-vqtl/scripts

METAL=../opt/generic-metal/metal


# EWIS meta-analysis for each phenotype

EWIS_DIR=../data/processed/ewis

cat ${EWIS_DIR}/ewis_phenotype_list.txt | while read pheno; do 

# Create new temporary file with STDERR column
for ancestry in EUR AFR EAS SAS; do
	echo "SNP REF ALT BETA STDERR" > ${EWIS_DIR}/${biomarker}/${ancestry}.tmp
	tail -n +2 ${EWIS_DIR}/${biomarker}/${pheno}_${ancestry} \
	| grep -v nan \
	| awk '{print $1,$4,$5,$10,sqrt($11)}'  \
	>> ${EWIS_DIR}/${biomarker}/${ancestry}.tmp
done

echo """
SCHEME STDERR

MARKER SNP
ALLELE REF ALT
EFFECT BETA 
STDERR STDERR

PROCESS ${EWIS_DIR}/${biomarker}/EUR.tmp
PROCESS ${EWIS_DIR}/${biomarker}/AFR.tmp
PROCESS ${EWIS_DIR}/${biomarker}/EAS.tmp
PROCESS ${EWIS_DIR}/${biomarker}/SAS.tmp

OUTFILE ${EWIS_DIR}/${biomarker}/${biomarker}_${pheno}_MA_ .tbl
ANALYZE HETEROGENEITY

QUIT
""" > tmp_metal_script_${biomarker}.txt

$METAL tmp_metal_script_${biomarker}.txt
rm tmp_metal_script_${biomarker}.txt ${EWIS_DIR}/${biomarker}/*.tmp
done
