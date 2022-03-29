#!/bin/sh


biomarker=$1
ancestry=$2
start_pheno_idx=$3
end_pheno_idx=$4


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00

#$ -pe smp 4
#$ -binding linear:4
#$ -R y

#$ -j y
#$ -cwd


source /broad/software/scripts/useuse
use .r-3.6.0

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
if (
  sum(!is.na(p[["e"]])) < 100 |  # Less than 100 non-missing exposure values
  (all(p[["e"]] %in% c(0, 1, NA)) & (min(table(p[["e"]])) < 50))  # Binary w/ either exposure having <50 instances 
) p[["pheno"]] <- NA
readr::write_csv(dplyr::select(p, id, e, pheno), "../data/processed/ewis/${biomarker}/ewis_phenos_${ancestry}_${exposure}.csv")
EOF

singularity exec \
	-B ../data/processed/ewis:/data \
	../../singularity/gem-v1.3-workflow.simg \
	/bin/bash <<EOF
/GEM/GEM \
	--pfile /data/ME_genotypes \
	--pheno-file /data/${biomarker}/ewis_phenos_${ancestry}_${exposure}.csv \
	--sampleid-name id \
	--pheno-name pheno \
	--pheno-type 0 \
	--exposure-names e \
	--delim , \
	--maf 0.01 \
	--missing-value NA \
	--robust 1 \
	--threads 4 \
	--output-style meta \
	--out /data/${biomarker}/${exposure}_${ancestry}_ME
EOF

rm ../data/processed/ewis/${biomarker}/ewis_phenos_${ancestry}_${exposure}.csv

done
