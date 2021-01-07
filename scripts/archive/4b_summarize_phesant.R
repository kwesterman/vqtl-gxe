#!/bin/bash


source /broad/software/scripts/useuse
use .r-3.6.0

cd ~/kw/ukbb-vqtl/scripts

library(data.table)

source("../opt/PHESANT/summarise_phenotypes.r")

phenofile <- "../data/processed/phesant_phenos/output..tsv"
logfile <- "../data/processed/phesant_phenos/output..log"

outcome_infofile <- "../data/processed/phesant_phenos/exposure_info.tsv"
data_codingfile <- "../variable-info/data-coding-ordinal-info-nov2019-update.txt"

hist_filename <- gsub("\\.tsv", "_hist", phenofile)
phenosummary_filename <- gsub("\\.tsv", "_phenosummary.tsv", phenofile)

phenos <- fread(phenofile, data.table=F, stringsAsFactors=F)
names(phenos)[2:ncol(phenos)] <- paste0("x", names(phenos)[2:ncol(phenos)])
outcome_info <- read.table(outcome_infofile, sep='\t', quote="", comment.char="", header=TRUE)
summ <-  get_hists_and_notes(hist_filename, phenos, logfile, outcome_info, data_codingfile, 
			     samples_for_inclusion=TRUE, check=FALSE, start_column=4)
write.table(summ, file=phenosummary_filename, col.names=TRUE, row.names=TRUE, 
	    sep='\t', quote=FALSE)
