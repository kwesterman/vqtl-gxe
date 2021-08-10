#!/bin/sh


biomarker=$1
ancestry=EUR


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00
#$ -j y
#$ -cwd


source /broad/software/scripts/useuse
use .r-3.6.0

R --vanilla <<EOF
p <- data.table::fread("../data/processed/ewis/ewis_phenos_EUR_PCA_wbiomarkers", 
                       select=c("id", "${biomarker}_adj", "PC1", "PC2"), 
                       data.table=F, stringsAsFactors=F)
p[["pheno"]] <- p[["${biomarker}_adj"]]
readr::write_csv(dplyr::select(p, id, PC1, PC2, pheno), "../data/processed/ewis/${biomarker}/ewis_phenos_${ancestry}_anthroPCs.csv")
EOF


singularity exec \
    -B ../data/processed/ewis:/data \
    ../../singularity/gem-v1.3-workflow.simg \
    /bin/bash <<EOF

/GEM/GEM \
    --pfile /data/ewis_genotypes \
    --pheno-file /data/${biomarker}/ewis_phenos_${ancestry}_anthroPCs.csv \
    --sampleid-name id \
    --pheno-name pheno \
    --pheno-type 0 \
    --exposure-names PC1 \
    --delim , \
    --maf 0.01 \
    --missing-value NA \
    --robust 1 \
    --threads 1 \
    --output-style meta \
    --out /data/${biomarker}/anthroPC1_${ancestry}

/GEM/GEM \
    --pfile /data/ewis_genotypes \
    --pheno-file /data/${biomarker}/ewis_phenos_${ancestry}_anthroPCs.csv \
    --sampleid-name id \
    --pheno-name pheno \
    --pheno-type 0 \
    --exposure-names PC2 \
    --delim , \
    --maf 0.01 \
    --missing-value NA \
    --robust 1 \
    --threads 1 \
    --output-style meta \
    --out /data/${biomarker}/anthroPC2_${ancestry}
EOF

rm ../data/processed/ewis/${biomarker}/ewis_phenos_${ancestry}_anthroPCs.csv
