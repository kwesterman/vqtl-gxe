#$ -l h_vmem=100G
#$ -l h_rt=24:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use .r-3.6.0


R --vanilla <<EOF
library(data.table)
library(tidyverse)

main_phenos <- fread("../data/raw/ukb_phenofiles/ukb10528.tab", stringsAsFactors=F, data.table=F)
diet_24HR <- fread("../data/raw/ukb_phenofiles/ukb22861.tab", stringsAsFactors=F, data.table=F)
biochemistry <- fread("../data/raw/ukb_phenofiles/ukb28679.tab", stringsAsFactors=F, data.table=F)
meds <- fread("../data/raw/ukb_phenofiles/ukb29937.tab", stringsAsFactors=F, data.table=F)
clinical_various <- fread("../data/raw/ukb_phenofiles/ukb38040.tab", stringsAsFactors=F, data.table=F)
body_comp <- fread("../data/raw/ukb_phenofiles/ukb38381.tab", stringsAsFactors=F, data.table=F)
alc_plus_various <- fread("../data/raw/ukb_phenofiles/ukb40167.tab", stringsAsFactors=F, data.table=F)
pa_edu <- fread("../data/raw/ukb_phenofiles/ukb42628.tab", stringsAsFactors=F, data.table=F)
phenos <- main_phenos %>%
  full_join(diet_24HR, by="f.eid") %>%
  full_join(biochemistry, by="f.eid") %>%
  full_join(meds, by="f.eid") %>%
  full_join(clinical_various, by="f.eid") %>%
  full_join(body_comp, by="f.eid") %>%
  full_join(alc_plus_various, by="f.eid") %>%
  full_join(pa_edu, by="f.eid")
phenos[grepl("\\\.y", names(phenos))] <- NULL  # Remove duplicated fields from merge
names(phenos) <- gsub("\\\.x", "", names(phenos))  # Rename other instance of duplicated fields
stopifnot(!any(duplicated(names(phenos))))  # Confirm no remaining duplicated fields
names(phenos) <- gsub("f\\\.", "x", names(phenos))
names(phenos) <- gsub("\\\.", "_", names(phenos))
names(phenos)[1] <- "id"
fwrite(phenos, "../data/processed/phesant_phenos/input_phenos.csv")

read_tsv("../opt/PHESANT/variable-info/outcome_info_final_multi_ancestry_jan2020.tsv", na="NA") %>%
  #mutate(EXCLUDED = case_when(
  #  EXCLUDED == "REMOVE" ~ "YES-ORIGINAL",  # Accept all removals from original PHESANT approach
  #  EXCLUDED_NEALE %in% c("YES-SEX", "YES-AGE") ~ "",  # Keep sex and age as exposures
  #  EXCLUDED_NEALE != "" ~ EXCLUDED_NEALE,  # Accept all other removals from Neale Mega-GWAS
  #  grepl("4[0-9]{4}", FieldID) ~ "YES-OUTCOMES",  # Remove disease outcomes
  #  FieldID == 20003 ~ "YES-MEDICATION",  # Remove medications (for now)
  #  TRUE ~ ""
  #)) %>%
  #select(-EXCLUDED_NEALE, -EXCLUDED_PHARMA) %>%
  #write_tsv("../data/processed/phesant_phenos/exposure_info.tsv")
  mutate(EXCLUDED = case_when(
    EXCLUDED_NEALE %in% c("YES-SEX", "Yes-AGE", "YES-ASSESSMENT-CENTRE") ~ "",  # Neale exclusions to be kept
    grepl("Diet by 24-hour recall", Path) ~ "",  # Also keep all 24HR diet data
    EXCLUDED_NEALE != "" ~ "REMOVE_NEALE",  # Otherwise, remove all Neale exclusions
    Category %in% c(  # Additionally remove the following...
      2000:2025,  # Hospital records
      100092:100093,  # Cancer and death registers
      9081,
      100081, 17518, 51428  # Blood-based assays
    ) ~ "REMOVE",
    TRUE ~ ""
  )) %>%
  select(FieldID, CAT_MULT_INDICATOR_FIELDS, CAT_SINGLE_TO_CAT_MULT,
         DATA_CODING, Path, Category, Field, ValueType, Units, Notes,
         Link, EXCLUDED) %>%
  write_tsv("../data/processed/phesant_phenos/ewis_variable_info.tsv")
EOF

cd ../opt/PHESANT/WAS

Rscript phenomeScan.r \
	--phenofile="../../../data/processed/phesant_phenos/input_phenos.csv" \
	--variablelistfile="../../../data/processed/phesant_phenos/ewis_variable_info.tsv" \
	--datacodingfile="../variable-info/data-coding-ordinal-info-nov2019-update.txt" \
	--resDir="../../../data/processed/phesant_phenos/" \
	--userId="id" \
	--partIdx=1 \
	--numParts=1 \
	--out="ewis_phenotypes"

cd ../../../scripts

#head -1 ewis_phenotypes.1.tsv | tr '\t' '\n' | tail -n +2 > ewis_phenotype_list.txt  # Write all phenotype names to a file
#cp ewis_phenotypes.1.tsv ewis_phenotypes_bkp.1.tsv  # Replace TRUE/FALSE with 1/0
#sed -e 's/FALSE/0/g; s/TRUE/1/g' ewis_phenotypes_bkp.1.tsv > ewis_phenotypes.1.tsv

R --vanilla <<EOF
library(data.table)
library(tidyverse)

ancestries <- c("EUR", "AFR", "EAS", "SAS") 
biomarkers <- scan("../data/processed/metabolic_biomarkers.txt", what=character())

phesant_phenos <- fread("../data/processed/phesant_phenos/ewis_phenotypes.1.tsv", data.table=F, stringsAsFactors=F)
phesant_phenos <- phesant_phenos[, -which(duplicated(names(phesant_phenos)))]
phesant_phenos <- mutate_if(phesant_phenos, is.logical, as.numeric)
names(phesant_phenos)[1] <- "id"

write(names(phesant_phenos)[-1], file="../data/processed/ewis/ewis_phenotype_list.txt")

INT <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))

for (anc in ancestries) {
	fread(paste0("../data/processed/vqtl_phenos_", anc, ".csv")) %>%
		select(id, all_of(biomarkers),
			all_of(paste0(biomarkers, "_adj"))) %>%
    		#mutate_at(biomarkers, list(INT = INT)) %>%
		inner_join(phesant_phenos, by="id") %>%
		write_csv(paste0("../data/processed/ewis/ewis_phenos_", anc, ".csv"))
}
EOF
