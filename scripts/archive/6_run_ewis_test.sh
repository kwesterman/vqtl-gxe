#!/bin/sh


biomarker=$1
ancestry=$2
start_pheno_idx=$3
end_pheno_idx=$4


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00
#$ -j y


source /broad/software/scripts/useuse
use .r-3.6.0

cd kw/ukbb-vqtl/scripts

mkdir -p ../data/processed/ewis/${biomarker}

cat ../data/processed/ewis/ewis_phenotype_list.txt \
	| head -${end_pheno_idx} \
	| tail -n +${start_pheno_idx} \
	| while read exposure; do

R --vanilla <<EOF
p <- data.table::fread("../data/processed/ewis/ewis_phenos_${ancestry}.csv", 
		       select=c("id", "${biomarker}_adj", "${exposure}"), 
		       data.table=F, stringsAsFactors=F)
p[["pheno"]] <- p[["${biomarker}_adj"]]
p[["e"]] <- p[["${exposure}"]]
p[["pheno"]] <- resid(lm(pheno ~ e, data=p, na.action=na.exclude))
if (
  sum(!is.na(p[["e"]])) < 100 |  # Less than 100 non-missing exposure values
  (all(p[["e"]] %in% c(0, 1)) & (min(table(p[["e"]])) < 50))  # Binary w/ either exposure having <50 instances 
) p$pheno <- NA
readr::write_csv(dplyr::select(p, id, e, pheno), "ewis_phenos_${ancestry}_${exposure}.csv")
EOF

singularity exec \
	-B ../data/processed/ewis:/data \
	../../singularity/gem-v1.2-workflow.simg \
	/bin/bash <<EOF
/GEM/GEM \
	--pfile /data/ewis_genotypes \
	--pheno-file ewis_phenos_${ancestry}_${exposure}.csv \
	--sampleid-name id \
	--pheno-name pheno \
	--pheno-type 0 \
	--exposure-names e \
	--delim , \
	--maf 0.01 \
	--missing-value NA \
	--robust 1 \
	--threads 1 \
	--out /data/${biomarker}/TEST_${biomarker}_${exposure}_${ancestry}
EOF

rm ewis_phenos_${ancestry}_${exposure}.csv

done
