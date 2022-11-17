#!/bin/bash


datadir=../data/processed

rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vqtl_phenos*csv ../data/processed/

vqtldir=../data/processed/vqtl_ss
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vqtl_ss/*gw ${vqtldir}/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vqtl_ss/*nom ${vqtldir}/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vqtl_ss/metal/*gw ${vqtldir}/metal/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vqtl_ss/metal/*nom ${vqtldir}/metal/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vqtl_ss/qq_plots ${vqtldir}/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vqtl_ss/ldsc/*log ${vqtldir}/ldsc/

medir=../data/processed/main_effect_ss
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/main_effect_ss/*gw ${medir}/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/main_effect_ss/*nom ${medir}/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/main_effect_ss/metal/*gw ${medir}/metal/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/main_effect_ss/metal/*nom ${medir}/metal/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/main_effect_ss/qq_plots ${medir}/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/main_effect_ss/ldsc/*log ${medir}/ldsc/

ewisdir=../data/processed/ewis
#rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/*_MA_1.tbl ${ewisdir}/
#rsync -avP --include='*tbl' --include='*/' --exclude='*' uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/ ${ewisdir}/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/ewis_genotypes.pvar ${ewisdir}/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/ewis_ma_results_nom.csv ${ewisdir}/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/ewis_anc_results_nom.csv ${ewisdir}/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/ewis_ma_results_anthroPCs.csv ${ewisdir}/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/*ewis_ma_results_qq.pdf ${ewisdir}/
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/ewis_ME_ma_results_nom.csv ${ewisdir}/

sensitivitydir=../data/processed/sensitivity
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/sensitivity/*INT_vqtl.vqtl ${sensitivitydir}/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/sensitivity/*adj_vqtl.vqtl ${sensitivitydir}/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/sensitivity/metal/*tbl_* ${sensitivitydir}/metal/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/sensitivity/metal/*tbl ${sensitivitydir}/metal/
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/sensitivity/primary_vqtl_df_withME.csv ${sensitivitydir}/

ancestryspecdir=../data/processed/ancestry_specific
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ancestry_specific/*raw ${ancestryspecdir}/

vignettesdir=../data/processed/vignettes
rsync -ravP uger:florez_ukb_projects/ukbb-vqtl/data/processed/vignettes ../data/processed/

conditionaldir=../data/processed/ewis/conditional_analysis
rsync -avP uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/conditional_analysis/independent_hits.csv ${conditionaldir}/
