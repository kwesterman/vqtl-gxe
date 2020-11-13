#!/bin/bash


cd ~/kw/ukbb-vqtl/scripts/

ewis_dir=../data/processed/ewis

# Use qctool v2 to subset 
source /broad/software/scripts/useuse
use GCC-5.2
qctool="../../opt/qctool_v2.0.6-CentOS\ Linux7.3.1611-x86_64/qctool"

for chr in {1..22}; do
	eval "${qctool}" \
		-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
		-s /humgen/diabetes/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
		-incl-positions ${ewis_dir}/ewis_variants.txt \
		-og ${ewis_dir}/ewis_genotypes_chr${chr}.bgen \
		-os ${ewis_dir}/ewis_genotypes_chr${chr}.sample
done

## Merge subsets
cp ${ewis_dir}/ewis_genotypes_chr22.bgen ${ewis_dir}/ewis_genotypes_base.bgen  # Initialize "base" .bgen file as chr22 .bgen
cp ${ewis_dir}/ewis_genotypes_chr22.sample ${ewis_dir}/ewis_genotypes_base.sample
for chr in {21..1}; do
	echo "Merging chromosome ${chr}..."
	eval "${qctool}" \
		-g ${ewis_dir}/ewis_genotypes_base.bgen \
		-s ${ewis_dir}/ewis_genotypes_base.sample \
		-merge-in ${ewis_dir}/ewis_genotypes_chr${chr}.bgen ${ewis_dir}/ewis_genotypes_chr${chr}.sample \
		-og ${ewis_dir}/ewis_genotypes.bgen \
		-os ${ewis_dir}/ewis_genotypes.sample
	cp ${ewis_dir}/ewis_genotypes.bgen ${ewis_dir}/ewis_genotypes_base.bgen  # Merged file becomes the new base for next round
	cp ${ewis_dir}/ewis_genotypes.sample ${ewis_dir}/ewis_genotypes_base.sample
done
rm ${ewis_dir}/ewis_genotypes_base.bgen ${ewis_dir}/ewis_genotypes_base.sample

## Convert to PLINK2 format
plink2=../../opt/plink2
${plink2} \
	--bgen ${ewis_dir}/ewis_genotypes.bgen ref-first \
	--sample ${ewis_dir}/ewis_genotypes.sample \
	--make-pgen \
	--out ${ewis_dir}/ewis_genotypes

## Convert PLINK2 format to PLINK1 format
${plink2} \
	--pfile ${ewis_dir}/ewis_genotypes \
	--make-bed \
	--hard-call-threshold 0.1 \
	--rm-dup force-first \
	--out ${ewis_dir}/ewis_genotypes
