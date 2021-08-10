#!/bin/bash


#$ -l h_vmem=90G
#$ -l h_rt=30:00:00

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use .r-3.6.0

Rscript 1_prep_vqtl_phenos.R
