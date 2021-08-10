#!/bin/bash


#$ -l h_vmem=40G
#$ -l h_rt=6:00:00

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use .r-3.6.0

Rscript 6c_collect_ewis_results.R
