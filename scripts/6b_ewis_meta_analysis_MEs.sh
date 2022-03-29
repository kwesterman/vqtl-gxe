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
	echo "SNP NONEFF EFF BETA STDERR" > ${EWIS_DIR}/${biomarker}/${ancestry}_ME.tmp
	tail -n +2 ${EWIS_DIR}/${biomarker}/${exp}_${ancestry}_ME \
	| grep -v nan \
	| awk '{print $1,$4,$5,$11,sqrt($13)}'  \
	>> ${EWIS_DIR}/${biomarker}/${ancestry}_ME.tmp
done

echo """
SCHEME STDERR

MARKER SNP
ALLELE EFF NONEFF
EFFECT BETA 
STDERR STDERR

PROCESS ${EWIS_DIR}/${biomarker}/EUR_ME.tmp
PROCESS ${EWIS_DIR}/${biomarker}/AFR_ME.tmp
PROCESS ${EWIS_DIR}/${biomarker}/EAS_ME.tmp
PROCESS ${EWIS_DIR}/${biomarker}/SAS_ME.tmp

OUTFILE ${EWIS_DIR}/${biomarker}/${biomarker}_${exp}_ME_MA_ .tbl
ANALYZE HETEROGENEITY

QUIT
""" > tmp_ewis_ME_metal_script_${biomarker}.txt

$METAL tmp_ewis_ME_metal_script_${biomarker}.txt
rm tmp_ewis_ME_metal_script_${biomarker}.txt ${EWIS_DIR}/${biomarker}/*.tmp
done
