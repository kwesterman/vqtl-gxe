#!/bin/bash


vqtldir=../data/processed/vqtl_ss

rsync -avP uger:kw/ukbb-vqtl/data/processed/vqtl_ss/*gw ${vqtldir}/
rsync -avP uger:kw/ukbb-vqtl/data/processed/vqtl_ss/*nom ${vqtldir}/
rsync -avP uger:kw/ukbb-vqtl/data/processed/vqtl_ss/metal/*gw ${vqtldir}/metal/
rsync -avP uger:kw/ukbb-vqtl/data/processed/vqtl_ss/metal/*nom ${vqtldir}/metal/
rsync -ravP uger:kw/ukbb-vqtl/data/processed/vqtl_ss/qq_plots ${vqtldir}/


medir=../data/processed/main_effect_ss

rsync -avP uger:kw/ukbb-vqtl/data/processed/main_effect_ss/*gw ${medir}/
rsync -avP uger:kw/ukbb-vqtl/data/processed/main_effect_ss/*nom ${medir}/
rsync -avP uger:kw/ukbb-vqtl/data/processed/main_effect_ss/metal/*gw ${medir}/metal/
rsync -avP uger:kw/ukbb-vqtl/data/processed/main_effect_ss/metal/*nom ${medir}/metal/
rsync -ravP uger:kw/ukbb-vqtl/data/processed/main_effect_ss/qq_plots ${medir}/


ewisdir=../data/processed/ewis

#rsync -ravP uger:kw/ukbb-vqtl/data/processed/ewis/*_MA_1.tbl ${ewisdir}/
#rsync -avP --include='*tbl' --include='*/' --exclude='*' uger:kw/ukbb-vqtl/data/processed/ewis/ ${ewisdir}/
rsync -ravP uger:kw/ukbb-vqtl/data/processed/ewis/ewis_ma_results_nom.csv ${ewisdir}/
rsync -ravP uger:kw/ukbb-vqtl/data/processed/ewis/ewis_anc_results_nom.csv ${ewisdir}/
