library(tidyverse)


### Define each of the ancestry groups --------------------

eur_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/UKBiobank_genoQC_reportedANDgeneticEUR_N455146_FLOREZ_EUR_PCA_covariates_40dim.txt") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY=genotyping.array) %>%
  mutate(id = as.character(id))

afr_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_AFR_N7447_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY) %>%
  mutate(id = as.character(id))

eas_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_EAS_N2264_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY) %>%
  mutate(id = as.character(id))

sas_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_SAS_N8669_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY) %>%
  mutate(id = as.character(id))


### Basic data fields ---------------------------------

batch <- data.table::fread("/humgen/florezlab/UKBB_app27892/ukb27892_cal_chrAUT_v2_s488363.fam", 
                           header=FALSE, data.table=FALSE, stringsAsFactors=FALSE,
                           col.names=c("FID", "IID", "P1", "P2", "S", "Batch")) %>% 
  mutate(id = as.character(format(IID, scientific=FALSE))) %>%
  select(id, Batch)

# 4 black to 4003 "any other black background" 
# 2 mixed to 2004 any other mixed background"
# 1 white to 1003 "any other white background"
# 3 asian or asian british to 3004 "any other asian background"
# -1 do not know to 6 "other ethnic group"

fix.ethnicity = c("4"=4003, "2"=2004, "1"=1003, "3"=3004, "-1" =6)

eth <- read_tsv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/ethnicity.phe") %>%
  select(id=f.eid, Ethnicity=f.21000.0.0) %>%
  mutate(id = as.character(format(id, scientific=FALSE)),
         Ethnicity = ifelse(as.character(Ethnicity) %in% names(fix.ethnicity), fix.ethnicity[as.character(Ethnicity)], Ethnicity),
         Ethnicity = as.factor(Ethnicity))  # now a factor with 17 levels

bmi <- read_tsv("/humgen/florezlab/UKBB_app27892/ukb10528.tab.gz") %>%
  select(id=f.eid,
         bmi=f.21001.0.0) %>%
  mutate(id = as.character(format(id, scientific=FALSE)))

print("Reading covariates")
base_phenos <- read_tsv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/merged.phe") %>%
  dplyr::rename(id=IID) %>%
  mutate(id = as.character(format(id, scientific=FALSE)),
         ageIndicator = as.factor(ifelse(age < 50, 50, ifelse(age > 78, 78, age))),
         ageBin = as.factor(ntile(age, 5)),
         FastingTime = as.factor(ifelse(FastingTime > 18, 18, ifelse(FastingTime == 0, 1, FastingTime))),
         DrawTime=as.factor(ntile(DrawTime, 20)),
         DilutionFactor=as.factor(ntile(DilutionFactor, 20))) %>%
  left_join(batch, by="id") %>% 
  left_join(eth, by="id") %>%
  left_join(bmi, by="id")

print("Unfiltered")
print(dim(base_phenos))

print("Not NA")
print(dim(na.omit(base_phenos)))

print("Generating covariate formula")
kept.terms <- setdiff(names(base_phenos), c("id", "age", "ac_date", "bmi"))
covariates <- c(kept.terms, "sex*ageIndicator", "sex*DrawTime", "sex*FastingTime", "sex*Ethnicity", "ageBin*FastingTime")
print(covariates)

### Biomarker data ----------------------------------

ukb_biomarker_fields <- c(
  alt = 30620, alb = 30600, alp = 30610, apoA = 30630, apoB = 30640,
  ast = 30650, hscrp = 30710, Ca = 30680, chol = 30690, creatinine = 30700,
  cysC = 30720, bilirubin_dir = 30660, ggt = 30730, glu = 30740, hba1c = 30750,
  hdl = 30760, igf1 = 30770, ldl = 30780, lipA = 30790, oestradiol = 30800,
  phos = 30810, rheum_factor = 30820, shbg = 30830, tes = 30850,
  bilirubin_tot = 30840, protein_tot = 30860, tg = 30870, urate = 30880,
  urea = 30670, vitD = 30890
)

# This table includes the codes for extra fields for adjustment
# i.e. 30621 = alanine aminotransferase assay date
# i.e. 30622 = alanine aminotransferase aliquot
biomarker2day <- read_tsv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariate_correction/biomarker2day_JBC.txt") %>%
  select(field=biomarker, name, extra_fields) %>%
  mutate(biomarker_KWname = names(ukb_biomarker_fields)[match(field, ukb_biomarker_fields)])

biomarker_phenos <- read_tsv("/humgen/florezlab/UKBB_app27892/ukb28679.tab.gz") %>%
  select(id=f.eid,
         all_of(setNames(paste0("f.", ukb_biomarker_fields, ".0.0"), 
                         names(ukb_biomarker_fields)))) %>%  # "Everything" in order to keep assay data and aliquot covariates as well 
  mutate(id = as.character(format(id, scientific=FALSE)))

extra = read_tsv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/biochem_fields_florez.phe", 
                 col_types=cols(.default="f", IID="c")) %>%
  dplyr::rename(id=IID)
str(extra)  # They are all factors in the regression model


### Outcomes for sample exclusions --------------------------------

outcomes <- read_tsv("/humgen/florezlab/UKBB_app27892/ukb38040.tab.gz") %>%
  select(id=f.eid, everything()) %>%
  mutate(id = as.character(format(id, scientific=FALSE)))
outcomes$diabetes <- rowSums(outcomes[, grepl("f\\.2443\\.0\\.", names(outcomes))] == 1, na.rm=T) > 0
outcomes$MI <- rowSums(outcomes[, grepl("f\\.6150\\.0\\.", names(outcomes))] == 1, na.rm=T) > 0
outcomes$angina <- rowSums(outcomes[, grepl("f\\.6150\\.0\\.", names(outcomes))] == 2, na.rm=T) > 0
cirrhosis_codes <- c(paste0("K70", 2:4), "K717", paste0("K74", 0:6))
cirrhosis_primary_ids <- c()
for (f in grep("f\\.41202\\.0\\.", names(outcomes), value=T)) {
  cirrhosis_primary_ids <- c(cirrhosis_primary_ids, outcomes$id[outcomes[[f]] %in% cirrhosis_codes])
}
cirrhosis_secondary_ids <- c()
for (f in grep("f\\.41204\\.0\\.", names(outcomes), value=T)) {
  cirrhosis_secondary_ids <- c(cirrhosis_secondary_ids, outcomes$id[outcomes[[f]] %in% cirrhosis_codes])
}
outcomes$pregnant <- rowSums(outcomes[, grepl("f\\.3140\\.0\\.", names(outcomes))] == 1, na.rm=T) > 0
cancer_tmp <- select(outcomes, 1, contains("f.40005."))
cancer <- inner_join(cancer_tmp[, 1:7], select(base_phenos, id, ac_date), by="id")  # Add assessment center dates
cancer$ac_date = as.Date(cancer$ac_date)
for (i in 2:7) {
  x <- ifelse(abs(difftime(cancer[, i, drop=TRUE], cancer$ac_date, units="days")) <= 365, TRUE, FALSE)  # TRUE if cancer diagnosis within a year of assessment center visit
  cancer <- cbind(cancer, x)
}
cancer$cancer_within_1yearac = apply(cancer[, 9:14], 1, function(x) {
  ifelse(any(x == TRUE, na.rm=TRUE), TRUE, FALSE)
})
cancer[names(cancer) == "x"] <- NULL

outcomes <- outcomes %>%
  inner_join(base_phenos, by="id") %>%
  left_join(select(biomarker_phenos, id, creatinine), by="id") %>%
  left_join(select(cancer, id, cancer_within_1yearac), by="id") %>%
  mutate(CHD = MI | angina,
         cirrhosis = id %in% c(cirrhosis_primary_ids, cirrhosis_secondary_ids),
         esrd_algorithm = !is.na(f.42026.0.0) & f.42026.0.0 <= ac_date,
         creatinine_conv = creatinine / 88.42,  # Convert from umol/L to mg/dL
         black_indicator = ifelse(id %in% afr_samples$id, 1, 0),
         eGFR = nephro::CKDEpi.creat(creatinine_conv, sex, age, black_indicator),
         esrd = esrd_algorithm | (!is.na(eGFR) & eGFR < 15)) %>%
  select(id, diabetes, CHD, cirrhosis, esrd_algorithm, eGFR, esrd, pregnant, cancer_within_1yearac)

### Statin adjustment ---------------------------------------

# field 20003 - drugs taken at baseline
drugs.taken <- read_tsv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/drugs_taken.phe") %>%
  mutate(id = as.character(format(f.eid, scientific=FALSE)))
# The list of statin drug codes we want to correct for.
drugs.adjusted <- scan("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariate_correction/statins_ids.txt", what=character())
# LDL, cholesterol, and Apolipoprotein B - adjustment factors 
statin.adj = read_tsv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/statin_adjustment.phe")

# For each statin drug, for each person, if they are on any of these drugs, make the Drugs=1, then add to the previous Drug value
# This means individuals who are on more than one statin will have values > 1?
num.drugs <- select(drugs.taken, id) # data frame of ids
num.drugs$Drugs <- 0 # set Drugs to 0 for everyone
for (i in drugs.adjusted) {
    print(i)
    num.drugs$Drugs <- num.drugs$Drugs + rowSums(drugs.taken[, 2:ncol(drugs.taken)] == i, na.rm=T)  # Add drug if 
}
on.drugs <- filter(num.drugs, Drugs > 0)$id  # A list of f.eid of individuals on drugs

for (bm in statin.adj$biomarker_kw){  # For each biomarker needing statin adjustment...
  # If the person is on a drug, then divide their biomarker value by the adjustment factor into a new biomarker adjusted value, if not, give the original value
  adj_factor = statin.adj$statin_adj_factor[statin.adj$biomarker_kw == bm]
  print(paste0(bm, " adj factor: ",adj_factor))
  biomarker_phenos[[paste0(bm, "_statinadj")]] <- ifelse(
    biomarker_phenos$id %in% on.drugs,
    biomarker_phenos[[bm]] / adj_factor,
    biomarker_phenos[[bm]]
  )
}

### Merge all phenotype datasets and write raw phenotype file ----------------------

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/w27892_20210201.csv", what=character())

phenos <- base_phenos %>%
  inner_join(biomarker_phenos, by="id") %>%
  left_join(outcomes, by="id") %>%
  left_join(extra, by="id") %>%
  filter(!(id %in% withdrawn_consent))

write_csv(phenos, "../data/processed/vqtl_phenos_raw.csv")
saveRDS(phenos, "../data/processed/vqtl_phenos_raw.rds")

### Preprocess biomarkers -------------------------------------

preprocess_vqtl_biomarker <- function(bm, df, ancestry) {
  
  # Given a biomarker vector and a set of relevant IDs, preprocess as follows:
  # 1) Set sample exclusions to missing
  # 2) Adjust the biomarker for age and 10 genetic principal components ... and a bunch of other covariates
  # 3) Remove outliers (> 5 SDs from the mean)
  # 4) Standardize residuals to z-scores (mean = 0, variance = 1)
  
  print(paste0("Preprocessing for ", bm, " in ", ancestry))
  
  # Sample exclusions
  df[[bm]][df$diabetes | df$CHD | df$cirrhosis | df$esrd | df$pregnant | df$cancer_within_1yearac] <- NA  # set to missing - those with chronic cardiometabolic disease, pregnancy, cancer
  
  # Run linear regression and print summary metrics
  raw_bm_name <- gsub("_log||_statinadj","", bm)
  extra_fields.regex <- biomarker2day$extra_fields[(biomarker2day$biomarker_KWname == raw_bm_name & 
                                                      !is.na(biomarker2day$biomarker_KWname))]
  extra_covar.columns = paste0("f.", unlist(strsplit(extra_fields.regex,"\\|")), ".0.0")
  covariates_complete <- paste(c(covariates, paste0("PC", 1:10), extra_covar.columns),
                               collapse=" + ")
  adj_lm_str <- paste0(bm, " ~ ", covariates_complete)
  adj_lm <- lm(as.formula(adj_lm_str), data=df, na.action=na.exclude)
  adj_pheno <- resid(adj_lm)
  all_estimates[[bm]][[ancestry]] <- broom::tidy(adj_lm)
  print(broom::tidy(adj_lm))
  all_r2[[bm]][[ancestry]] <- summary(adj_lm)$adj.r.squared
  print(paste("Adjusted R-squared:", summary(adj_lm)$adj.r.squared))
  
  # Remove outliers
  adj_pheno <- ifelse(  # Remove outliers (outside of 5 SDs from the mean)
    findInterval(adj_pheno, mean(adj_pheno, na.rm = T) +
                   5 * c(-1, 1) * sd(adj_pheno, na.rm = T)) == 1,
    adj_pheno, NA
  )
  
  # Scale phenotype
  scale(adj_pheno)
}


raw_biomarkers <- names(ukb_biomarker_fields)
biomarkers_withStatinAdj <- c(raw_biomarkers, "chol_statinadj", "ldl_statinadj", "apoB_statinadj")
log_biomarkers <- paste0(biomarkers_withStatinAdj, "_log")
all_biomarkers <- c(biomarkers_withStatinAdj, log_biomarkers)

phenos <- phenos %>%
  mutate_at(biomarkers_withStatinAdj, list(log = log))

eur_phenos <- eur_samples %>%
  inner_join(phenos, by="id")

afr_phenos <- afr_samples %>%
  inner_join(phenos, by="id")

eas_phenos <- eas_samples %>%
  inner_join(phenos, by="id")

sas_phenos <- sas_samples %>%
  inner_join(phenos, by="id")

joint_PCs <- read_csv("/humgen/florezlab/UKBB_app27892/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_app27892.csv") %>%
  mutate(id = as.character(f.eid)) %>%
  select(id, contains("PC")) %>%
	 # contains("related"))  # Ensure unrelatedness by membership in ancestry-specific groups instead of using this value
  filter(!duplicated(id))
names(joint_PCs) <- gsub("_return2442", "", names(joint_PCs))
all_phenos <- phenos %>%
  inner_join(joint_PCs, by="id") %>%
  filter(id %in% c(eur_phenos$id, afr_phenos$id, eas_phenos$id, sas_phenos$id))

all_estimates <- list()
all_r2 <- list()

for (bm in all_biomarkers) {
  eur_phenos[[paste0(bm, "_adj")]] <- preprocess_vqtl_biomarker(bm, eur_phenos, "eur")
  afr_phenos[[paste0(bm, "_adj")]] <- preprocess_vqtl_biomarker(bm, afr_phenos, "afr")
  eas_phenos[[paste0(bm, "_adj")]] <- preprocess_vqtl_biomarker(bm, eas_phenos, "eas")
  sas_phenos[[paste0(bm, "_adj")]] <- preprocess_vqtl_biomarker(bm, sas_phenos, "sas")
  all_phenos[[paste0(bm, "_adj")]] <- preprocess_vqtl_biomarker(bm, all_phenos, "all")
}

write_csv(eur_phenos, "../data/processed/vqtl_phenos_EUR.csv")
write_csv(afr_phenos, "../data/processed/vqtl_phenos_AFR.csv")
write_csv(eas_phenos, "../data/processed/vqtl_phenos_EAS.csv")
write_csv(sas_phenos, "../data/processed/vqtl_phenos_SAS.csv")
write_csv(all_phenos, "../data/processed/vqtl_phenos_all.csv")

lapply(all_estimates, bind_rows, .id="ancestry") %>%
  bind_rows(.id="biomarker") %>%
  write_tsv("../output/phenotyping_model_fit_estimates.txt")

lapply(all_r2, bind_rows, .id="ancestry") %>%
  bind_rows(.id="biomarker") %>%
  write_tsv("../output/phenotyping_model_fit_adjustedR2.txt")
