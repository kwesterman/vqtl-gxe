 #!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=10:00:00

#$ -cwd
#$ -j y


vqtl_dir=../data/processed/vqtl_ss
cmdkp_dir=../data/processed/cmdkp


source /broad/software/scripts/useuse
use .r-3.6.0


# Prepare vQTL summary stats
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do
	# Ancestry-specific
	echo "Preparing ${bm} summary stats for CMDKP."
	echo "${anc}..."
	ss=${vqtl_dir}/${bm}_${anc}_vqtl_merged
	output_fn=${cmdkp_dir}/${bm}_${anc}_vqtls_full.csv
	R --vanilla <<EOF
	library(tidyverse)
	ss_df <- read_tsv("${ss}") %>%
		mutate(P = ifelse(P == 0, 1e-300, P)) %>%
		select(CHR=Chr, POS=bp, SNP, REF=A2, ALT=A1, P) %>%
		write_csv("${output_fn}")
EOF
	gzip ${output_fn}
done

	# Meta-analysis
	echo "Meta-analysis..."
	ss=${vqtl_dir}/metal/${bm}_MA_1.tbl
	output_fn=${cmdkp_dir}/${bm}_metaanalysis_vqtls_full.csv
	cp ${ss} ${output_fn}
	gzip ${output_fn}
done

# Prepare GWIS meta-analysis summary stats
ss=../data/processed/ewis/ewis_ma_results.csv
output_fn=${cmdkp_dir}/gwis_metaanalysis_full.csv
R --vanilla <<EOF
library(tidyverse)
ss_df <- read_csv("${ss}") %>%
	select(Biomarker=bm, SNP, Exposure=exposure, REF=Allele2, ALT=Allele1, Beta, SE, P, Dir, ISq, HetP) %>%
	write_csv("${output_fn}")
EOF
gzip ${output_fn}
