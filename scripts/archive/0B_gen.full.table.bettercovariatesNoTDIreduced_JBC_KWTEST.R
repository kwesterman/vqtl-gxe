library(dplyr)
library(broom)
library(RNOmni)

batch <- data.table::fread("/humgen/florezlab/UKBB_app27892/ukb27892_cal_chrAUT_v2_s488363.fam", header=F, col.names=c("FID", "IID", "P1", "P2", "S", "Batch")) %>% select(IID, Batch)

# 4 black to 4003 "any other black background" 
# 2 mixed to 2004 any other mixed background"
# 1 white to 1003 "any other white background"
# 3 asian or asian british to 3004 "any other asian background"
# -1 do not know to 6 "other ethnic group"

fix.ethnicity = c("4"=4003, "2"=2004, "1"=1003, "3"=3004, "-1" =6)

eth <- read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/ethnicity.phe", header=T, col.names=c("IID", "Ethnicity")) %>%
    mutate(Ethnicity = ifelse(as.character(Ethnicity) %in% names(fix.ethnicity), fix.ethnicity[as.character(Ethnicity)], Ethnicity)) %>%
    mutate(Ethnicity = as.factor(Ethnicity))
## now a factor with 17 levels
## table(eth$Ethnicity)

print("Reading covariates")
base_phenos <- read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/merged.phe", header=T, stringsAsFactors=F) %>%
    mutate(ageIndicator = as.factor(ifelse(age < 50, 50, ifelse(age > 78, 78, age)))) %>%
    mutate(ageBin = as.factor(ntile(age, 5))) %>%
    mutate(FastingTime = as.factor(ifelse(FastingTime > 18, 18, ifelse(FastingTime == 0, 1, FastingTime)))) %>%
    #mutate(UrineSampleMinutes=as.factor(ifelse(UrineSampleMinutes == 0, 0, ntile(UrineSampleMinutes, 20)))) %>%
    mutate(DrawTime=as.factor(ntile(DrawTime, 20))) %>%
    mutate(DilutionFactor=as.factor(ntile(DilutionFactor, 20))) %>%
    left_join(batch) %>% left_join(eth)

print("Unfiltered")
print(dim(base_phenos))

print("not na")
print(dim(na.omit(base_phenos)))

print("Generating covariate formula")
#kept.terms <- colnames(base_phenos %>% select(-IID, -FID, -age, -Array, -TDI))
kept.terms <- colnames(base_phenos %>% select(-IID, -age, -ac_date))
#covariates <- paste0(paste(kept.terms, collapse=" + "),
#     " + ageIndicator * sex + sex * DrawTime + sex * UrineSampleMinutes + sex * FastingTime + ageBin * FastingTime + Ethnicity * sex")
covariates <- paste0(paste(kept.terms, collapse=" + "),
     " + ageIndicator * sex + sex * DrawTime + sex * FastingTime + ageBin * FastingTime + Ethnicity * sex")

print(covariates)

#regression.data <- base_phenos %>% mutate(IID=as.numeric(IID)) %>% na.omit
#print("regression data")
#print(dim(regression.data))

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

library(tidyverse)

eur_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/UKBiobank_genoQC_reportedANDgeneticEUR_N455146_FLOREZ_EUR_PCA_covariates_40dim.txt") %>%
  filter(unrelateds == 1) %>%
  select(IID=Florez_IID, contains("PC"), cov_GENO_ARRAY=genotyping.array)

afr_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_AFR_N7447_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(IID=Florez_IID, contains("PC"), cov_GENO_ARRAY)

eas_samples<- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_EAS_N2264_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(IID=Florez_IID, contains("PC"), cov_GENO_ARRAY)

sas_samples <- read_tsv("/humgen/diabetes2/users/jcole/UKBB/pheno/multiethnic/UKB_genoQC_SAS_N8669_geneticPCs_genoarray") %>%
  filter(unrelateds == 1) %>%
  select(IID=Florez_IID, contains("PC"), cov_GENO_ARRAY)

ukb_biomarker_fields <- c(
  alt = 30620, alb = 30600, alp = 30610, apoA = 30630, apoB = 30640,
  ast = 30650, hscrp = 30710, Ca = 30680, chol = 30690, creatinine = 30700,
  cysC = 30720, bilirubin_dir = 30660, ggt = 30730, glu = 30740, hba1c = 30750,
  hdl = 30760, igf1 = 30770, ldl = 30780, lipA = 30790, oestradiol = 30800,
  phos = 30810, rheum_factor = 30820, shbg = 30830, tes = 30850,
  bilirubin_tot = 30840, protein_tot = 30860, tg = 30870, urate = 30880,
  urea = 30670, vitD = 30890
)

######### This table includes the codes for extra fields for adjustment
# i.e. 30621 = alanine aminotransferase assay date
# i.e. 30622 = alanine aminotransferase alliquot
biomarker2day <- read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariate_correction/biomarker2day_JBC.txt", header=T, stringsAsFactors=F)
biomarker2day$biomarker_KWname = names(ukb_biomarker_fields)[match(biomarker2day$biomarker, ukb_biomarker_fields)]

### Now reformat these fields so we can pull them from UKB data
ukb_biomarker_fields <- setNames(paste0("f.", ukb_biomarker_fields, ".0.0"),
                                 names(ukb_biomarker_fields))

biomarker_phenos <- read_tsv("/humgen/florezlab/UKBB_app27892/ukb28679.tab.gz") %>%
  select(IID = f.eid,
         all_of(ukb_biomarker_fields))

######### This table includes the actual extra fields for adjustment
########## create a fake adjustment file
##
##cut -d$'\t' -f3 /humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariate_correction/biomarker2day.txt | tail -n+2 | tr '|' '\n'  > \
##/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariate_correction/biomarker_extra.txt
##
#########
##df = read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariate_correction/biomarker_extra.txt", header=FALSE, col.names="extra")
##length(df$extra)
##x = replicate(n=length(df$extra), sample(1:4, nrow(base_phenos), replace=TRUE))
##colnames(x) = paste0("f.",df$extra,".0.0")
##y = cbind(base_phenos$IID,x)
##colnames(y)[1] ="IID"
##
##write.table(y,"/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/biomarker_extra_random.phe", row.names=FALSE, col.names=TRUE, quote=FALSE)
##### oops dup columns
##cut -d$'\t' -f1-65 /humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/biomarker_extra_random.phe > ./tmp
##mv ./tmp /humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/biomarker_extra_random.phe
##columns /humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/biomarker_extra_random.phe

####### please load the real one...
#extra = read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/biomarker_extra_random.phe", header=T, stringsAsFactors=F)
extra = read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/biochem_fields_florez.phe", header=T, colClasses = c("numeric",rep("factor",179)))
# they are all factors in the regression model
str(extra)

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

#############################
################## Sample exclusions

outcomes <- read_tsv("/humgen/florezlab/UKBB_app27892/ukb38040.tab.gz") %>%
  select(IID = f.eid, everything())
outcomes$diabetes <- rowSums(outcomes[, grepl("f\\.2443\\.0\\.", names(outcomes))] == 1, na.rm=T) > 0
outcomes$MI <- rowSums(outcomes[, grepl("f\\.6150\\.0\\.", names(outcomes))] == 1, na.rm=T) > 0
outcomes$angina <- rowSums(outcomes[, grepl("f\\.6150\\.0\\.", names(outcomes))] == 2, na.rm=T) > 0
cirrhosis_codes <- c(paste0("K70", 2:4), "K717", paste0("K74", 0:6))
cirrhosis_primary_ids <- c()
for (f in grep("f\\.41202\\.0\\.", names(outcomes), value=T)) {
  cirrhosis_primary_ids <- c(cirrhosis_primary_ids, outcomes$IID[outcomes[[f]] %in% cirrhosis_codes])
}
cirrhosis_secondary_ids <- c()
for (f in grep("f\\.41204\\.0\\.", names(outcomes), value=T)) {
  cirrhosis_secondary_ids <- c(cirrhosis_secondary_ids, outcomes$IID[outcomes[[f]] %in% cirrhosis_codes])
}
######### pregnant
outcomes$pregnant <- rowSums(outcomes[, grepl("f\\.3140\\.0\\.", names(outcomes))] == 1, na.rm=T) > 0
table(outcomes$pregnant, useNA="always")

######## cancer diagnosis in last year
## if any of their cancer diagnosis dates were in the last year -> cancer=1
## date of cancer diagnosis
cancer_tmp = as.data.frame(outcomes[,c(1,grep("f.40005.",colnames(outcomes)))])
# extra cancer fields that are empty
cancer = merge(cancer_tmp[,1:7],base_phenos[,c("IID","ac_date")], by="IID",all.y=T)
cancer$ac_date = as.Date(cancer$ac_date)
str(cancer)

for (i in 2:7) {
x = ifelse(abs(difftime(cancer[,i],cancer$ac_date, units="days"))<=365, TRUE, FALSE)
cancer=cbind(cancer,x)
}

cancer$cancer_within_1yearac = apply(cancer[,9:14], 1, function(x) ifelse(any(x==TRUE,na.rm=T),TRUE ,FALSE))
#head(cancer,30)
table(cancer$cancer_within_1yearac)
#############################

outcomes <- outcomes %>%
  inner_join(base_phenos, by="IID") %>%
  left_join(select(biomarker_phenos, IID, creatinine), by="IID") %>%
  left_join(cancer[,c("IID","cancer_within_1yearac")], by="IID") %>%
  mutate(CHD = MI | angina,
         cirrhosis = IID %in% c(cirrhosis_primary_ids, cirrhosis_secondary_ids),
         esrd_algorithm = !is.na(f.42026.0.0) & f.42026.0.0 <= ac_date,
         creatinine_conv = creatinine / 88.42,  # Convert from umol/L to mg/dL
         black_indicator = ifelse(IID %in% afr_samples$IID, 1, 0),
         eGFR = nephro::CKDEpi.creat(creatinine_conv, sex, age, black_indicator),
         esrd = esrd_algorithm | (!is.na(eGFR) & eGFR < 15)) %>%
  select(IID, diabetes, CHD, cirrhosis, esrd_algorithm, eGFR, esrd, pregnant, cancer_within_1yearac)

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/w27892_20210201.csv", what=character())

############################ Add statin adjustment here:

## LDL, cholesterol, and Apolipoprotein B - adjustment factors 
statin.adj = read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/statin_adjustment.phe", header=T)
# field 20003 - drugs taken at baseline
drugs.taken <- read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariates/drugs_taken.phe", header=T)
## The list of statin drug codes we want to correct for.
drugs.adjusted <- read.table("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/src/covariate_correction/statins_ids.txt", header=F)$V1

## for each statin drug, for each person, if they are on any of these drugs, make the Drugs=1, then add to the previous Drug value
## This means individuals who are on more than one statin will have values > 1?
num.drugs <- select(drugs.taken, f.eid) # data frame of ids
num.drugs$Drugs <- 0 # set Drugs to 0 for everyone
for (i in drugs.adjusted) {
    print(i)
    num.drugs$Drugs <- num.drugs$Drugs + rowSums(drugs.taken[,2:ncol(drugs.taken)] == i, na.rm=T) #
}
## a list of f.eid of individuals on drugs
on.drugs <- num.drugs %>% filter(Drugs > 0) %>% .$f.eid
table(on.drugs>0)

## for each drug in bm_adj
## if the person is on a drug, then divide their biomarker value by the adjustment factor into a new biomarker adjusted value, if not, give the original value
bm_adj = c("ldl","chol","apoB")

biomarker_phenos = as.data.frame(biomarker_phenos)
for (bm in bm_adj){
  adj_factor = statin.adj[statin.adj$biomarker_kw==bm,"statin_adj_factor"]
  print(paste0(bm, " adj factor: ",adj_factor))
  biomarker_phenos$new_column = ifelse((biomarker_phenos$IID %in% on.drugs), biomarker_phenos[,bm]/adj_factor, biomarker_phenos[,bm])
  colnames(biomarker_phenos)[ncol(biomarker_phenos)] = paste0(bm,"_statinadj")
}

## checkme
#head(biomarker_phenos)
#head(num.drugs)
# this id is on statins
#biomarker_phenos[biomarker_phenos$IID=="1000069",c("IID","ldl","chol","apoB","ldl_statinadj","chol_statinadj","apoB_statinadj")]

############################

phenos <- base_phenos %>%
  inner_join(biomarker_phenos, by="IID") %>%
  left_join(outcomes, by="IID") %>%
  left_join(extra, by="IID") %>%
  filter(!(IID %in% withdrawn_consent))

# write_csv(phenos, "/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/data/vqtl_phenos_raw_rivas2.csv")

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

## KW Function
preprocess_vqtl_biomarker <- function(bm, df) {

  # Given a biomarker vector and a set of relevant IDs, preprocess as follows:
  # 1) Set sample exclusions to missing
  # 2) Adjust the biomarker for age and 10 genetic principal components
  # 3) Remove outliers (> 5 SDs from the mean)
  # 4) Standardize residuals to z-scores within each sex

# Sample exlcusions
  bm[df$diabetes | df$CHD | df$cirrhosis | df$esrd | df$pregnant | df$cancer_within_1yearac ] <- NA  # set to missing -  those with chronic cardiometabolic , pregnancy, cancer
# Run Linear Regression
  adj_lm <- lm(as.formula(adj_lm_str), data = df, na.action = na.exclude)
  adj_pheno <- resid(adj_lm)
# Remove outliers
  adj_pheno <- ifelse(  # Remove outliers (outside of 5 SDs from the mean)
    findInterval(adj_pheno, mean(adj_pheno, na.rm = T) +
                   5 * c(-1, 1) * sd(adj_pheno, na.rm = T)) == 1,
    adj_pheno, NA
  )
 
# Scale phenotype  
  scale(adj_pheno)

}


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

biomarkers <- names(ukb_biomarker_fields)
statinadj_biomarkers <- paste0(bm_adj, "_statinadj")
biomarkers2 = c(biomarkers, statinadj_biomarkers)
log_biomarkers <- paste0(biomarkers2, "_log")
all_biomarkers <- c(biomarkers2, log_biomarkers)

######## log transform biomarkers 

eur_phenos <- inner_join(eur_samples, phenos, by="IID")
eur_phenos = eur_phenos %>%
  mutate_at(biomarkers2, list(log = log)) 

afr_phenos <- inner_join(afr_samples, phenos, by="IID")
afr_phenos = afr_phenos %>%
  mutate_at(biomarkers2, list(log = log)) 

sas_phenos <- inner_join(sas_samples, phenos, by="IID")
sas_phenos = sas_phenos %>%
  mutate_at(biomarkers2, list(log = log)) 

eas_phenos <- inner_join(eas_samples, phenos, by="IID")
eas_phenos = eas_phenos %>%
  mutate_at(biomarkers2, list(log = log)) 
  
######## Remove sample exclusions, adjust, and scale
estimates_all=NULL
r2_all=NULL

for (i in 1:length(all_biomarkers)){

######################################################################################
## the biomarker
######################################################################################

bm_name=all_biomarkers[i]
print(paste0("biomarker: ",bm_name))

######################################################################################
## covariate formula - ancestry and biomarker-specific 
######################################################################################

bm_name_for_covariates=gsub("_log||_statinadj","",bm_name)
extra_fields.regex <- biomarker2day[(biomarker2day$biomarker_KWname == bm_name_for_covariates & !is.na(biomarker2day$biomarker_KWname)), "extra_fields"]
 
extra.covar.columns = unlist(strsplit(extra_fields.regex,"\\|")) %>%
  gsub("^","f.",.) %>%
  gsub("$",".0.0",.)
print(paste0("extra covariates: ", paste0(extra.covar.columns, collapse=" ") ))

covariates_complete <- paste0(paste0(extra.covar.columns, collapse=" + "), " + ", covariates, " + ", paste0("PC", 1:10, collapse = " + "), collapse=" + ")
adj_lm_str <- paste0("bm ~ ", covariates_complete)  

# ######################################################################################
# ## run function and print our regression summaries
# ######################################################################################
# 
# eur_phenos = eur_phenos %>% 
#   mutate_at(., .vars=bm_name, list(adj = preprocess_vqtl_biomarker), .) %>%
#   select(IID, everything())
# colnames(eur_phenos)[ncol(eur_phenos)] = paste0(bm_name,"_adj")
# 
# ###########################################
# ## print out tables of regression summaries
# test = eur_phenos
# bm = test[,bm_name]
# # Sample exlcusions
# test[test$diabetes | test$CHD | test$cirrhosis | test$esrd | test$pregnant | test$cancer_within_1yearac ,] <- NA  # set to missing -  those with chronic cardiometabolic , pregnancy, cancer
# # Run Linear Regression
# #formula2 = gsub("bm ~ ","",adj_lm_str)
# adj_lm <- lm(as.formula(paste0(bm_name, " ~ ", covariates_complete) ), data = test, na.action = na.exclude)
# # write out estimates
# eur_estimates=tidy(adj_lm)
# colnames(eur_estimates)[2:5] = paste0("eur_",colnames(eur_estimates)[2:5])
# eur_estimates
# # adjusted R2
# eur_r2=summary(adj_lm)$adj.r.squared
# eur_r2
# 
# ######################################################################################
# 
# afr_phenos = afr_phenos %>% 
#   mutate_at(., .vars=bm_name, list(adj = preprocess_vqtl_biomarker), .) %>%
#   select(IID, everything())
# colnames(afr_phenos)[ncol(afr_phenos)] = paste0(bm_name,"_adj")
# 
# ###########################################
# ## print out tables of regression summaries
# test = afr_phenos
# bm = test[,bm_name]
# # Sample exlcusions
# test[test$diabetes | test$CHD | test$cirrhosis | test$esrd | test$pregnant | test$cancer_within_1yearac ,] <- NA  # set to missing -  those with chronic cardiometabolic , pregnancy, cancer
# # Run Linear Regression
# #formula2 = gsub("bm ~ ","",adj_lm_str)
# adj_lm <- lm(as.formula(paste0(bm_name, " ~ ", covariates_complete) ), data = test, na.action = na.exclude)
# # write out estimates
# afr_estimates=tidy(adj_lm)
# colnames(afr_estimates)[2:5] = paste0("afr_",colnames(afr_estimates)[2:5])
# afr_estimates
# # adjusted R2
# afr_r2=summary(adj_lm)$adj.r.squared
# afr_r2
# 
# ######################################################################################

sas_phenos = sas_phenos %>% 
  mutate_at(., .vars=bm_name, list(adj = preprocess_vqtl_biomarker), .) %>%
  select(IID, everything())
colnames(sas_phenos)[ncol(sas_phenos)] = paste0(bm_name,"_adj")

###########################################
## print out tables of regression summaries
test = sas_phenos
bm = test[,bm_name]
# Sample exlcusions
test[test$diabetes | test$CHD | test$cirrhosis | test$esrd | test$pregnant | test$cancer_within_1yearac ,] <- NA  # set to missing -  those with chronic cardiometabolic , pregnancy, cancer
# Run Linear Regression
#formula2 = gsub("bm ~ ","",adj_lm_str)
adj_lm <- lm(as.formula(paste0(bm_name, " ~ ", covariates_complete) ), data = test, na.action = na.exclude)
# write out estimates
sas_estimates=tidy(adj_lm)
colnames(sas_estimates)[2:5] = paste0("sas_",colnames(sas_estimates)[2:5])
sas_estimates
# adjusted R2
sas_r2=summary(adj_lm)$adj.r.squared
sas_r2

# ######################################################################################
# eas_phenos = eas_phenos %>% 
#   mutate_at(., .vars=bm_name, list(adj = preprocess_vqtl_biomarker), .) %>%
#   select(IID, everything())
# colnames(eas_phenos)[ncol(eas_phenos)] = paste0(bm_name,"_adj")
# 
# ###########################################
# ## print out tables of regression summaries
# test = eas_phenos
# bm = test[,bm_name]
# # Sample exlcusions
# test[test$diabetes | test$CHD | test$cirrhosis | test$esrd | test$pregnant | test$cancer_within_1yearac ,] <- NA  # set to missing -  those with chronic cardiometabolic , pregnancy, cancer
# # Run Linear Regression
# #formula2 = gsub("bm ~ ","",adj_lm_str)
# adj_lm <- lm(as.formula(paste0(bm_name, " ~ ", covariates_complete) ), data = test, na.action = na.exclude)
# # write out estimates
# eas_estimates=tidy(adj_lm)
# colnames(eas_estimates)[2:5] = paste0("eas_",colnames(eas_estimates)[2:5])
# eas_estimates
# # adjusted R2
# eas_r2=summary(adj_lm)$adj.r.squared
# eas_r2
# ######################################################################################

############# combine ancestries for model fit output
estimates = eur_estimates %>%
  full_join(afr_estimates, by="term") %>%
  full_join(sas_estimates, by="term") %>%
  full_join(eas_estimates, by="term") %>%
  add_column(., .before=1, biomarker = bm_name)
estimates_all = rbind(estimates_all, estimates)
# print.data.frame(estimates[grepl("Ethnicity",estimates$term),]) # ancestry-specific ethnicity indicators 

#estimates$biomarker = bm_name
#estimates[,c(ncol(estimates),1,2:(ncol(estimates)-1))]
r2 = c(bm_name, eur_r2, afr_r2, sas_r2, eas_r2)
r2_all = rbind(r2_all, r2)

}
########
# I need to test this mutate: mutate(FID=id, IID=id). I previously just sed the ID correction 
# eur_phenos %>% mutate(id=IID) %>% write_csv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/data/vqtl_phenos_EUR_modelrivas2.csv")
# afr_phenos %>% mutate(id=IID) %>% write_csv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/data/vqtl_phenos_AFR_modelrivas2.csv")
sas_phenos %>% mutate(id=IID) %>% write_csv("../data/processed/vqtl_phenos_SAS_modelrivas2.csv")
# eas_phenos %>% mutate(id=IID) %>% write_csv("/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/data/vqtl_phenos_EAS_modelrivas2.csv")

############################# Print model fit estimates and adjusted r2 for all
# estimates_all$term = gsub("\\(Intercept\\)","Intercept",estimates_all$term)
# write.table(estimates_all, "/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/data/model_fit/ALL_rivas2_modelfit_estimates.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# colnames(r2_all) = c("biomarker","eur_r2","afr_r2","sas_r2","eas_eur")
# write.table(r2_all, "/humgen/florezlab/users/jcole/UKBB/diet/vqtl/pheno_processing/data/model_fit/ALL_rivas2_modelfit_adjustedRsquared.txt", col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

