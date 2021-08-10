#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=00:20:00

#$ -cwd
#$ -j y


bm=$1

sensdir=../data/processed/sensitivity
osca=../../opt/osca_Linux


# Test vQTL with all-exposure-adjusted outcome
${osca} \
	--vqtl \
	--bfile ../data/processed/ewis/ewis_genotypes \
	--pheno ${sensdir}/${bm}_allExpAdj_vqtl_phenos.csv \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--out ${sensdir}/${bm}_allExpAdj_gw_snps
