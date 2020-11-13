#!/bin/sh


##$ -l h_vmem=25G
##$ -l h_rt=28:00:00
##$ -j y
#
##$ -pe smp 4
##$ -binding linear:4
##$ -R y

bm1=$1
bm2=$2
ancestry=EUR


source /broad/software/scripts/useuse
use Anaconda3

cd ~/kw/ukbb-vqtl/scripts
vqtl_dir=../data/processed/vqtl_ss
source activate ldsc


# Run LDSC to calculate the genetic correlation
../opt/ldsc/ldsc.py \
	--rg ${vqtl_dir}/ldsc/${bm1}_${ancestry}_ldsc.sumstats.gz,${vqtl_dir}/ldsc/${bm2}_${ancestry}_ldsc.sumstats.gz \
	--ref-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--w-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--out ${vqtl_dir}/ldsc/${bm1}_${bm2}_${ancestry}

conda deactivate
