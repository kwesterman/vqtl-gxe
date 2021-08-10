#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=00:30:00

#$ -cwd
#$ -j y


pheno=$1
anc=$2
p_col_vqtl=12
p_col_ME=12

source /broad/software/scripts/useuse
use .r-3.6.0


echo "Creating summary files for ${anc} ${pheno}..."
# NOTE: the below "_merged" file is not INFO filtered (that is applied in R script)

sumstats_root=../data/processed/vqtl_ss/${pheno}_${anc}
head -1 ${sumstats_root}_chr1.vqtl > ${sumstats_root}_vqtl_merged
for chr in {1..22}; do
	echo "Chromosome " ${chr}
	tail -n +2 ${sumstats_root}_chr${chr}.vqtl >> ${sumstats_root}_vqtl_merged
done


ME_root=../data/processed/main_effect_ss/${pheno}_${anc}
head -1 ${ME_root}_chr1.${pheno}_adj.glm.linear > ${ME_root}_ME_merged
for chr in {1..22}; do
	echo "Chromosome " ${chr}
	tail -n +2 ${ME_root}_chr${chr}.${pheno}_adj.glm.linear >> ${ME_root}_ME_merged
done

Rscript 3b_postprocess_vqtl.R ${sumstats_root}_vqtl_merged
Rscript 3b_postprocess_vqtl.R ${ME_root}_ME_merged
