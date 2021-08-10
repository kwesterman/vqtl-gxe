#STEP0A - JBC

setwd("/humgen/diabetes2/users/jcole/UKBB/diet/vqtl/pheno_processing/src")

library(data.table)
library(dplyr)

## read in covariate files
age = fread("covariates/age_sex.phe", data.table=FALSE,verbose=T)
drawtime <- fread("covariates/draw_time_minutes.phe", data.table=FALSE,verbose=T)
fast <- fread("covariates/fastingtime.phe", data.table=FALSE,verbose=T)
dilution <- read.table("covariates/DilutionFactor_30897.phe",header=T)
colnames(dilution)[1] = "f.eid"

## merge 
age %>% inner_join(drawtime) %>% inner_join(fast) %>% inner_join(dilution) -> full.cov
colnames(full.cov) = c("IID","sex","age","DrawTime","FastingTime","DilutionFactor")
dim(full.cov)
# an inner join removes missing info, but not too bad

## format assessment month and day
assess <- fread("covariates/assessment.phe", data.table=FALSE,verbose=T)
assess$ac_date = assess$f.53.0.0
assess$f.53.0.0 <- sub("-[^-]*$", "", assess$f.53.0.0)
library(mltools)
library(data.table)
assess$f.53.0.0 <- as.factor(assess$f.53.0.0)
assess$f.54.0.0 <- as.factor(assess$f.54.0.0)
assess$f.53.0.0 <- as.factor(gsub("^2006-.*", "2006", assess$f.53.0.0))
assess$f.53.0.0 <- as.factor(gsub("^2010-08", "2010-0810", assess$f.53.0.0))
assess$f.53.0.0 <- as.factor(gsub("^2010-09", "2010-0810", assess$f.53.0.0))
assess$f.53.0.0 <- as.factor(gsub("^2010-10", "2010-0810", assess$f.53.0.0))
oneassess <- one_hot(as.data.table(assess)) # default columns are all unordered factored columns
colnames(oneassess)[1] <- "IID"
colnames(oneassess) = gsub("-",".", colnames(oneassess))
colnames(oneassess)
full.cov %>% inner_join(oneassess) -> full.cov

colnames(full.cov)

write.table(full.cov, "covariates/merged.phe", quote=F, sep="\t", row.names=F, col.names=T)
