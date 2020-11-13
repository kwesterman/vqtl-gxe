#!/bin/sh


chr=22
dir=/broad/hptmp/kwesterm/ukb_bfiles
plink2=../../opt/plink2

source /broad/software/scripts/useuse

cd kw/ukbb-vqtl/scripts

samplefile=/humgen/diabetes/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
if [ chr = "X" ]
then
	samplefile=/humgen/diabetes/UKBB_app27892/ukb27892_imp_chrX_v3_s486743.sample
fi

# Read in UKB .bgen file and convert to bed/bim/fam fileset
${plink2} \
	--bgen /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen ref-first \
	--sample $samplefile \
	--maf 0.495 \
	--hard-call-threshold 0.1 \
	--make-bed \
	--out test_genos

#cp ${dir}/chr${chr}.bim ${dir}/chr${chr}.bim_raw

# Use PLINK to format .bim file for the OSCA program
${plink2} \
	--bfile ${dir}/chr${chr} \
	--set-all-var-ids @:#_\$r_\$a \
	--new-id-max-allele-len 650 \
	--make-just-bim \
	--out ${dir}/chr${chr}

if [ chr = "12" ]  # Keep the first instance of any variants with duplicate IDs
then
	${plink2} \
		--bfile ${dir}/chr${chr} \
		--rm-dup force-first \
		--make-bed \
		--out ${dir}/chr${chr}
fi

if [ chr = "X" ]
then
	cp ${dir}/chr${chr}.bim ${dir}/chr${chr}.bim_tmp
	sed 's/X/23/g' ${dir}/chr${chr}.bim_tmp > ${dir}/chr${chr}.bim
