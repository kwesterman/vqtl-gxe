#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use .r-3.6.0

Rscript 9c_collect_ewis_anthroPC_results.R
