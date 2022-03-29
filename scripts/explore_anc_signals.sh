#!/bin/bash


plink2=../../opt/plink2

ewis_dir=../data/processed/ewis
anc_spec_dir=../data/processed/ancestry_specific

mkdir -p ${anc_spec_dir}

# Create genotype matrix with only significant variants
${plink2} \
	--pfile ${ewis_dir}/ewis_genotypes \
	--extract ${anc_spec_dir}/ancestry_specific_variants.txt \
	--export A ref-first \
	--out ${anc_spec_dir}/ancestry_specific_genotypes
