#!/bin/bash


#$ -l h_vmem=23G
#$ -l h_rt=24:00:00

#$ -cwd
#$ -j y

#$ -pe smp 6
#$ -binding linear:6
#$ -R y


bm=$1
chr=$2

sensdir=../data/processed/sensitivity
osca=../../opt/osca_Linux


# Test vQTL with all-exposure-adjusted outcome
${osca} \
	--vqtl \
	--bfile /broad/hptmp/kwesterm/ukb_bfiles/chr${chr} \
	--pheno ${sensdir}/${bm}_allExpAdj_vqtl_phenos.csv \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--thread-num 6 \
	--out ${sensdir}/${bm}_chr${chr}_allExpAdj
