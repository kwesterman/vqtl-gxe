#!/bin/sh


#$ -l h_vmem=15G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y


bm=$1
ancestry=EUR


source /broad/software/scripts/useuse
use Anaconda3
use .r-3.6.0


vqtl_dir=../data/processed/vqtl_ss
mkdir -p ${vqtl_dir}/ldsc	
me_dir=../data/processed/main_effect_ss
mkdir -p ${me_dir}/ldsc	

source activate ldsc


# Download necessary reference LD files
if [ ! -e ../data/raw/ldsc/eur_w_ld_chr ]; then
	wget -P ../data/raw/ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
	bzip2 -d ../data/raw/ldsc/eur_w_ld_chr.tar.bz2
	tar xvf ../data/raw/ldsc/eur_w_ld_chr.tar -C ../data/raw/ldsc/ 
fi
if [ ! -f ../data/raw/ldsc/w_hm3.snplist ]; then
	wget -P ../data/raw/ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
	bzip2 -d ../data/raw/ldsc/w_hm3.snplist.bz2
fi

# Associate vQTL summary stats files with rsIDs
R --vanilla <<EOF
  library(tidyverse)
  hm3_vars <- scan("../data/raw/ldsc/w_hm3.snplist", what=character())
  var_annot <- data.table::fread("../data/processed/ukb_variants_maf0.005_info0.5.txt", data.table=F, stringsAsFactors=F,
				 col.names=c("CHR", "SNPID", "rsID", "POS", "REF", "ALT", "MAF", "Q", "INFO")) %>%
    filter(rsID %in% hm3_vars) %>%
    mutate(CHR = gsub(":.*", "", SNPID))
  sumstats <- data.table::fread("${vqtl_dir}/${bm}_${ancestry}_vqtl_merged", data.table=F, stringsAsFactors=F) %>%
    dplyr::rename(CHR=Chr, POS=bp) %>%
    mutate(CHR = as.character(CHR)) %>%
    inner_join(var_annot, by=c("CHR", "POS")) %>%
    select(rsID, A1, A2, freq, NMISS, BETA=beta, SE=se, P) %>%
    write_tsv("${vqtl_dir}/ldsc/${bm}_${ancestry}_hm3")
  sumstats <- data.table::fread("${me_dir}/${bm}_${ancestry}_ME_merged", data.table=F, stringsAsFactors=F) %>%
    setNames(c("CHR", "POS", "ID", "A2", "ALT", "A1", "A1_FREQ", "TEST", "N", "BETA", "SE", "T_STAT", "P")) %>%
    mutate(CHR = as.character(CHR),
	   P = as.numeric(P)) %>%
    inner_join(var_annot, by=c("CHR", "POS")) %>%
    select(rsID, A1, A2, N, BETA, SE, P) %>%
    write_tsv("${me_dir}/ldsc/${bm}_${ancestry}_hm3")
EOF

# Munge sumstats for vQTL
../opt/ldsc/munge_sumstats.py \
	--sumstats ${vqtl_dir}/ldsc/${bm}_${ancestry}_hm3 \
	--merge-alleles ../data/raw/ldsc/w_hm3.snplist \
	--snp rsID \
	--N-col NMISS \
	--a1 A1 \
	--a2 A2 \
	--p P \
	--frq freq \
	--signed-sumstats BETA,0 \
	--out ${vqtl_dir}/ldsc/${bm}_${ancestry}_ldsc

# Munge sumstats for main effect
../opt/ldsc/munge_sumstats.py \
	--sumstats ${me_dir}/ldsc/${bm}_${ancestry}_hm3 \
	--merge-alleles ../data/raw/ldsc/w_hm3.snplist \
	--snp rsID \
	--N-col N \
	--a1 A1 \
	--a2 A2 \
	--p P \
	--signed-sumstats BETA,0 \
	--out ${me_dir}/ldsc/${bm}_${ancestry}_ldsc
