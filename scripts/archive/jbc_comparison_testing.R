test = sas_phenos
bm = test[,bm_name]
# Sample exlcusions
test[test$diabetes | test$CHD | test$cirrhosis | test$esrd | test$pregnant | test$cancer_within_1yearac ,] <- NA  # set to missing -  those with chronic cardiometabolic , pregnancy, cancer
# Run Linear Regression
#formula2 = gsub("bm ~ ","",adj_lm_str)
adj_lm <- lm(as.formula(paste0(bm_name, " ~ ", covariates_complete) ), data = test, na.action = na.exclude)
# write out estimates
eur_estimates=tidy(adj_lm)
colnames(eur_estimates)[2:5] = paste0("eur_",colnames(eur_estimates)[2:5])
eur_estimates
# adjusted R2
eur_r2=summary(adj_lm)$adj.r.squared
eur_r2



#############

## KW portion ##

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

eur_phenos$alt_log_adj <- preprocess_vqtl_biomarker("alt_log", eur_phenos, "eur")

## JBC portion ##

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

jep <- fread("/humgen/diabetes2/users/jcole/UKBB/diet/vqtl/pheno_processing/data/vqtl_phenos_EUR_modelrivas2.csv", data.table=F, stringsAsFactors=F)

bm_name <- "alt_log"
bm_name_for_covariates <- "alt"
biomarker2day <- as.data.frame(biomarker2day)  # Makes next line extract string rather than data frame
extra_fields.regex <- biomarker2day[(biomarker2day$biomarker_KWname == bm_name_for_covariates & !is.na(biomarker2day$biomarker_KWname)), "extra_fields"]
extra.covar.columns = unlist(strsplit(extra_fields.regex,"\\|")) %>%
  gsub("^","f.",.) %>%
  gsub("$",".0.0",.)
covariates_complete <- paste0(paste0(extra.covar.columns, collapse=" + "), " + ", covariates, " + ", paste0("PC", 1:10, collapse = " + "), collapse=" + ")
# covariates_complete <- paste0(setdiff(covariates, c(grep("^f", covariates, value=T), 
#                                                     grep("sex|FastingTime|Ethnicity|Batch", covariates, value=T))), 
#                               collapse=" + ")
covariates_complete
adj_lm_str <- paste0("bm ~ ", covariates_complete)  

eur_phenos$alt_log_adj_jbcLike <- preprocess_vqtl_biomarker(eur_phenos$alt_log, eur_phenos)
jep$alt_log_adj_jbcLike <- preprocess_vqtl_biomarker(jep$alt_log, jep)

cor.test(jep$alt_log_adj, eur_phenos$alt_log_adj)
cor.test(jep$alt_log_adj, eur_phenos$alt_log_adj_jbcLike)


#########

sas_phenos$mytest <- preprocess_vqtl_biomarker(sas_phenos$alt_log, sas_phenos)
a$mytest <- preprocess_vqtl_biomarker(a$alt_log, a)
cor.test(sas_phenos$mytest, a$mytest)

