#!/bin/bash


#$ -l h_vmem=22G
#$ -l h_rt=28:00:00

#$ -pe smp 6
#$ -binding linear:6
#$ -R y

#$ -cwd
#$ -j y


bm=$1
anc=$2
chr=$3

sensdir=../data/processed/sensitivity
osca=../../opt/osca_Linux


# Test vQTL with all-exposure-adjusted outcome
${osca} \
	--vqtl \
	--bfile /broad/hptmp/kwesterm/ukb_bfiles/chr${chr} \
	--pheno ${sensdir}/${bm}_expAdj_vqtl_phenos_${anc}.csv \
	--vqtl-mtd 2 \
	--maf 0.01 \
	--thread-num 6 \
	--out ${sensdir}/${bm}_chr${chr}_${anc}_expAdj
