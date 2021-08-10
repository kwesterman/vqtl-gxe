#!/bin/sh


#$ -l h_vmem=70G
#$ -l h_rt=20:00:00
#$ -j y

chr=$1
ukb_dir=/broad/ukbb/imputed_v3
dir=/broad/hptmp/kwesterm/ukb_bfiles
plink2=../../opt/plink2

source /broad/software/scripts/useuse

cd kw/ukbb-vqtl/scripts

samplefile=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
if [ chr = "X" ]
then
	samplefile=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrX_v3_s486743.sample
fi

# Generate list of variants w/ MAF > 0.05 in full sample
awk '$6 > 0.005 {print $2}' ${ukb_dir}/ukb_mfi_chr${chr}_v3.txt > ${dir}/maf005_variants_chr${chr}.txt

# Read in UKB .bgen file and convert to pgen/pvar/psam fileset
${plink2} \
	--bgen ${ukb_dir}/ukb_imp_chr${chr}_v3.bgen ref-first \
	--sample $samplefile \
	--extract ${dir}/maf005_variants_chr${chr}.txt \
	--memory 50000 \
	--make-pgen \
	--out ${dir}/chr${chr}

# Read in pgen/pvar/psam fileset and convert to bed/bim/fam fileset
${plink2} \
	--pfile ${dir}/chr${chr} \
	--rm-dup force-first \
	--hard-call-threshold 0.1 \
	--make-bed \
	--out ${dir}/chr${chr}

#cp ${dir}/chr${chr}.bim ${dir}/chr${chr}.bim_raw

# Use PLINK to format .bim file for the OSCA program
${plink2} \
	--bfile ${dir}/chr${chr} \
	--new-id-max-allele-len 650 \
	--make-just-bim \
	--out ${dir}/chr${chr}
	#--set-all-var-ids @:#_\$r_\$a \

#if [ chr = "12" ]  # Keep the first instance of any variants with duplicate IDs
#then
#	${plink2} \
#		--bfile ${dir}/chr${chr} \
#		--rm-dup force-first \
#		--make-bed \
#		--out ${dir}/chr${chr}
#fi

if [ chr = "X" ]
then
	cp ${dir}/chr${chr}.bim ${dir}/chr${chr}.bim_tmp
	sed 's/X/23/g' ${dir}/chr${chr}.bim_tmp > ${dir}/chr${chr}.bim
fi

# Create separate (much smaller) ancestry-specific plinksets
for anc in AFR EAS SAS; do
	awk -F, '{print $1,$1}' ../data/processed/vqtl_phenos_${anc}.csv > keep_IDs_chr${chr}_${anc}.tmp
	${plink2} \
		--bfile ${dir}/chr${chr} \
		--keep keep_IDs_chr${chr}_${anc}.tmp \
		--make-bed \
		--out ${dir}/chr${chr}_${anc}
	rm keep_IDs_chr${chr}_${anc}.tmp
done
