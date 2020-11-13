#!/bin/bash


cd ~/kw/ukbb-vqtl/scripts

ukb_dir=/broad/ukbb/imputed_v3


echo "" > ../data/processed/ukb_variants_maf0.005_info0.5.txt
for chr in {1..22} X; do
	echo "Filtering chromosome ${chr}..."
	awk -v chr=${chr} '$6 > 0.005 && $8 > 0.5 {print chr"\t"$0}' ${ukb_dir}/ukb_mfi_chr${chr}_v3.txt >> ../data/processed/ukb_variants_maf0.005_info0.5.txt
done
