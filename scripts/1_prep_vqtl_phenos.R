library(tidyverse)


admin <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/UKBiobank_assessment_center_f.54_birthplace_f.1647.txt") %>%
  select(id = Florez_FID, assessment_centre=f.54.0.0_categorical, birthplace=f.1647.0.0_categorical)

eur_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/UKBiobank_genoQC_reportedANDgeneticEUR_N455146_FLOREZ_EUR_PCA_covariates_40dim.txt") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY=genotyping.array)

afr_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_AFR_N7447_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY)

eas_samples<- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_EAS_N2264_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY)

sas_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_SAS_N8669_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(id=Florez_IID, contains("PC"), cov_GENO_ARRAY)

base_phenos <- read_tsv("../data/raw/ukb_phenofiles/ukb10528.tab") %>%
  select(id = f.eid,
         sex = f.31.0.0,
         age = f.21022.0.0,
         fasting_hrs = f.74.0.0,
         bmi = f.21001.0.0,
         smoking = f.20116.0.0,
         ac_date = f.53.0.0) %>%
  mutate(age_squared = age ** 2)

ukb_biomarker_fields <- c(
  alt = 30620, alb = 30600, alp = 30610, apoA = 30630, apoB = 30640, 
  ast = 30650, hscrp = 30710, Ca = 30680, chol = 30690, creatinine = 30700, 
  cysC = 30720, bilirubin_dir = 30660, ggt = 30730, glu = 30740, hba1c = 30750, 
  hdl = 30760, igf1 = 30770, ldl = 30780, lipA = 30790, oestradiol = 30800, 
  phos = 30810, rheum_factor = 30820, shbg = 30830, tes = 30850, 
  bilirubin_tot = 30840, protein_tot = 30860, tg = 30870, urate = 30880, 
  urea = 30670, vitD = 30890 
)
ukb_biomarker_fields <- setNames(paste0("f.", ukb_biomarker_fields, ".0.0"), 
                                 names(ukb_biomarker_fields))
  
biomarker_phenos <- read_tsv("/humgen/diabetes/UKBB_app27892/ukb28679.tab.gz") %>%
  select(id = f.eid, 
         all_of(ukb_biomarker_fields))

# DM_phenos <- read_tsv("/humgen/diabetes/users/jcole/UKBB/pheno/round2_T2D_updateMay2019/UKBiobank_ALLethnicities_diabetes_complete_2020Feb.txt") %>%
#   select(id = f.eid,
#          t1d = prob_poss_t1dm_all_plus,
#          t2d = prob_poss_t2dm_all_plus_t2d_controls_strict_hba1c) %>%
#   replace_na(list(t1d = F, t2d = F))

outcomes <- read_tsv("/humgen/diabetes/UKBB_app27892/ukb38040.tab.gz") %>%
  select(id = f.eid, everything())
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
outcomes <- outcomes %>%
  inner_join(base_phenos, by="id") %>%
  left_join(select(biomarker_phenos, id, creatinine), by="id") %>%
  mutate(CHD = MI | angina,
         cirrhosis = id %in% c(cirrhosis_primary_ids, cirrhosis_secondary_ids),
         esrd_algorithm = !is.na(f.42026.0.0) & f.42026.0.0 <= ac_date,
         creatinine_conv = creatinine / 88.42,  # Convert from umol/L to mg/dL
         black_indicator = ifelse(id %in% afr_samples$id, 1, 0),
         eGFR = nephro::CKDEpi.creat(creatinine_conv, sex, age, rep(0, nrow(.))),
         esrd = esrd_algorithm | (!is.na(eGFR) & eGFR < 15)) %>%
  select(id, diabetes, CHD, cirrhosis, esrd_algorithm, eGFR, esrd)

withdrawn_consent <- scan("/humgen/diabetes/UKBB_app27892/w27892_20200204.csv", what=character())

phenos <- base_phenos %>%
  inner_join(biomarker_phenos, by="id") %>%
  left_join(outcomes, by="id") %>%
  filter(!(id %in% withdrawn_consent))

write_csv(phenos, "../data/processed/vqtl_phenos_raw.csv")


preprocess_vqtl_biomarker <- function(bm, df) {
  
  # Given a biomarker vector and a set of relevant IDs, preprocess as follows:
  # 1) Set individuals with diabetes to missing
  # 2) Adjust the biomarker for age and 10 genetic principal components
  # 3) Remove outliers (> 5 SDs from the mean)
  # 4) Standardize residuals to z-scores within each sex
  
  bm[df$diabetes | df$CHD | df$cirrhosis | df$esrd] <- NA  # Remove those with chronic cardiometabolic diseases
  adj_lm_str <- paste0("bm ~ age + ", paste0("PC", 1:10, collapse = " + "))  # Adjust for age and 10 PCs
  adj_lm <- lm(as.formula(adj_lm_str), data = df, na.action = na.exclude)
  adj_pheno <- resid(adj_lm)
  adj_pheno <- ifelse(  # Remove outliers (outside of 5 SDs from the mean)
    findInterval(adj_pheno, mean(adj_pheno, na.rm = T) + 
                   5 * c(-1, 1) * sd(adj_pheno, na.rm = T)) == 1,
    adj_pheno, NA
  )
  scale(adj_pheno)
}


biomarkers <- names(ukb_biomarker_fields)
write(biomarkers, "../data/processed/30biomarkers.txt")


eur_phenos <- inner_join(eur_samples, phenos, by="id") 
eur_phenos %>%
  nest(data = -sex) %>%  # Split into male and female datasets
  mutate(data = map(data, function(d) {  # Stratified preprocessing
    mutate_at(d, biomarkers, list(adj = preprocess_vqtl_biomarker), d)
  })) %>%
  unnest(data) %>% 
  select(id, everything()) %>%
  write_csv("../data/processed/vqtl_phenos_EUR.csv")

afr_phenos <- inner_join(afr_samples, phenos, by="id")
afr_phenos %>%
  nest(data = -sex) %>%
  mutate(data = map(data, function(d) {
    mutate_at(d, biomarkers, list(adj = preprocess_vqtl_biomarker), d)
  })) %>%
  unnest(data) %>% 
  select(id, everything()) %>%
  write_csv("../data/processed/vqtl_phenos_AFR.csv")

eas_phenos <- inner_join(eas_samples, phenos, by="id")
eas_phenos %>%
  nest(data = -sex) %>%
  mutate(data = map(data, function(d) {
    mutate_at(d, biomarkers, list(adj = preprocess_vqtl_biomarker), d)
  })) %>%
  unnest(data) %>% 
  select(id, everything()) %>%
  write_csv("../data/processed/vqtl_phenos_EAS.csv")

sas_phenos <- inner_join(sas_samples, phenos, by="id")
sas_phenos %>%
  nest(data = -sex) %>%
  mutate(data = map(data, function(d) {
    mutate_at(d, biomarkers, list(adj = preprocess_vqtl_biomarker), d)
  })) %>%
  unnest(data) %>% 
  select(id, everything()) %>%
  write_csv("../data/processed/vqtl_phenos_SAS.csv")
