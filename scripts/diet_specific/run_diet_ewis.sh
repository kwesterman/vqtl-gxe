#!/bin/sh


biomarker=$1
ancestry=$2


#$ -l os=RedHat7
#$ -l h_vmem=5G
#$ -l h_rt=5:00:00
#$ -j y

#$ -pe smp 4
#$ -binding linear:4

source /broad/software/scripts/useuse
use .r-3.6.0

cd kw/ukbb-vqtl/scripts

mkdir -p ../data/processed/diet_ewis/${biomarker}

diet_jcole=/humgen/diabetes2/users/jcole/UKBB/diet/UKB_24HRdiet_genoQCEUR_adjusted_INT_2019Nov
#head -1 $jcole | tr '\t' '\n' | grep '_average\b' > ../data/processed/diet_ewis/24hr_exposures.txt

R --vanilla <<EOF
library(tidyverse)
diet <- data.table::fread("${diet_jcole}", stringsAsFactors=F, data.table=F) %>%
  select(id=Florez_IID, everything())
p <- data.table::fread("../data/processed/ewis/ewis_phenos_${ancestry}.csv", 
		       data.table=F, stringsAsFactors=F) %>%
  select(id, contains("_adj")) %>%
  inner_join(diet, by="id")
p[["tcals"]] <- p[["energy_yesterday_total.100002_average"]]
readr::write_csv(p, "../data/processed/diet_ewis/diet_ewis_phenos_${ancestry}.csv")
EOF
#p[["pheno"]] <- p[["${biomarker}_adj"]]
#p[["e"]] <- p[["${exposure}"]]

cat ../data/processed/diet_ewis/24hr_exposures.txt \
	| while read exposure; do

echo ${exposure}

singularity exec \
	-B ../data/processed:/data \
	../../singularity/gem-v1.2-workflow.simg \
	/bin/bash <<EOF
/GEM/GEM \
	--pfile /data/ewis/ewis_genotypes \
	--pheno-file /data/diet_ewis/diet_ewis_phenos_${ancestry}.csv \
	--sampleid-name id \
	--pheno-name ${biomarker}_adj \
	--pheno-type 0 \
	--exposure-names ${exposure} \
	--int-covar-names tcals \
	--delim , \
	--maf 0.01 \
	--missing-value NA \
	--robust 1 \
	--threads 2 \
	--out /data/diet_ewis/${biomarker}/${exposure}_${ancestry}
EOF

done
