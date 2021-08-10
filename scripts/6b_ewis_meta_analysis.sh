#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=24:00:00

#$ -cwd
#$ -j y


biomarker=$1

METAL=../opt/generic-metal/metal


# EWIS meta-analysis for each phenotype

EWIS_DIR=../data/processed/ewis

cat ${EWIS_DIR}/ewis_phenotype_list.txt | while read exp; do 

# Create new temporary file with STDERR column
for ancestry in EUR AFR EAS SAS; do
	echo "SNP REF ALT BETA STDERR" > ${EWIS_DIR}/${biomarker}/${ancestry}.tmp
	tail -n +2 ${EWIS_DIR}/${biomarker}/${exp}_${ancestry} \
	| grep -v nan \
	| awk '{print $1,$4,$5,$11,sqrt($13)}'  \
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

OUTFILE ${EWIS_DIR}/${biomarker}/${biomarker}_${exp}_MA_ .tbl
ANALYZE HETEROGENEITY

QUIT
""" > tmp_metal_script_${biomarker}.txt

$METAL tmp_metal_script_${biomarker}.txt
rm tmp_metal_script_${biomarker}.txt ${EWIS_DIR}/${biomarker}/*.tmp
done
