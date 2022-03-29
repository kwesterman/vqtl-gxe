#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00

#$ -cwd
#$ -j y


bm=$1
anc=$2

sensdir=../data/processed/sensitivity
osca=../../opt/osca_Linux


# Test vQTL with all-exposure-adjusted outcome
${osca} \
	--vqtl \
	--bfile ../data/processed/ewis/ewis_genotypes \
	--pheno ${sensdir}/${bm}_expAdj_vqtl_phenos_${anc}.csv \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--thread-num 6 \
	--out ${sensdir}/${bm}_${anc}_expAdj
