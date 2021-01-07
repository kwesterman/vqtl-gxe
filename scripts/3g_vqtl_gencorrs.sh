#!/bin/sh


#$ -j y


bm1=$1
bm2=$2
ancestry=EUR


source /broad/software/scripts/useuse
use Anaconda3

cd ~/kw/ukbb-vqtl/scripts
vqtl_dir=../data/processed/vqtl_ss
me_dir=../data/processed/main_effect_ss
source activate ldsc


# Run LDSC to calculate genetic correlations

# For vQTLs
../opt/ldsc/ldsc.py \
	--rg ${vqtl_dir}/ldsc/${bm1}_${ancestry}_ldsc.sumstats.gz,${vqtl_dir}/ldsc/${bm2}_${ancestry}_ldsc.sumstats.gz \
	--ref-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--w-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--out ${vqtl_dir}/ldsc/${bm1}_${bm2}_${ancestry}

# ...and for main effects
../opt/ldsc/ldsc.py \
	--rg ${me_dir}/ldsc/${bm1}_${ancestry}_ldsc.sumstats.gz,${me_dir}/ldsc/${bm2}_${ancestry}_ldsc.sumstats.gz \
	--ref-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--w-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--out ${me_dir}/ldsc/${bm1}_${bm2}_${ancestry}

conda deactivate
