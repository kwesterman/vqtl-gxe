#!/bin/bash


pheno=$1
anc=$2
p_col_vqtl=12
p_col_ME=12

cd kw/ukbb-vqtl/scripts

source /broad/software/scripts/useuse
use .r-3.6.0


echo "Creating summary files for ${anc} ${pheno}..."

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

#echo "Creating nominal subset..."
#head -1 ${sumstats_root}_vqtl_merged > ${sumstats_root}_vqtl_merged_nom
#awk -v col=$p_col_vqtl '$col < 1e-2' ${sumstats_root}_vqtl_merged >> ${sumstats_root}_vqtl_merged_nom
#
#head -1 ${ME_root}_ME_merged > ${ME_root}_ME_merged_nom
#awk -v col=$p_col_ME '$col < 1e-2' ${ME_root}_ME_merged >> ${ME_root}_ME_merged_nom
#
#echo "Creating suggestive subset..."
#head -1 ${sumstats_root}_vqtl_merged > ${sumstats_root}_vqtl_merged_sugg
#awk -v col=$p_col_vqtl '$col < 1e-5' ${sumstats_root}_vqtl_merged >> ${sumstats_root}_vqtl_merged_sugg
#
#head -1 ${ME_root}_ME_merged > ${ME_root}_ME_merged_sugg
#awk -v col=$p_col_ME '$col < 1e-5' ${ME_root}_ME_merged >> ${ME_root}_ME_merged_sugg
#
#echo "Writing p-value file..."
#tail -n +2 ${sumstats_root}_merged | awk '{print $1,$12}' > ${sumstats_root}_merged_pvals.txt
#tail -n +2 ${ME_root}_merged | awk '{print $1,$12}' > ${ME_root}_merged_pvals.txt
