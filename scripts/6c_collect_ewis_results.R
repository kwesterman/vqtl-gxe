library(data.table)
library(tidyverse)


biomarkers <- scan("../data/processed/metabolic_biomarkers.txt", what=character())
exposures <- scan("../data/processed/ewis/ewis_phenotype_list.txt", what=character())


# META-ANALYSIS

read_gxe_ma <- function(bm, e) {
  fn <- paste0("../data/processed/ewis/", bm, "/", bm, "_", e, "_MA_1.tbl")
  df <- tryCatch({
    fread(fn, data.table=F, stringsAsFactors=F, na.strings=c("NA", "nan", "nane-nan", "")) %>%
      select(SNP=MarkerName, Allele1, Allele2, 
	     Beta=Effect, SE=StdErr, P="P-value", 
             Dir=Direction, ISq=HetISq, HetP=HetPVal) %>%
      mutate(P=as.numeric(P))
  }, error=function (e) tibble())
  if (e == exposures[1]) print(bm)
  if (nrow(df) > 1) df else NULL
}


ewis_df <- expand_grid(
  bm = biomarkers,
  exposure = exposures
) %>%
  mutate(res = pmap(list(.$bm, .$exposure), read_gxe_ma)) %>%
  unnest(res)

ewis_df %>%
  filter(P < 0.05) %>%
  write_csv("../data/processed/ewis/ewis_ma_results_nom.csv")

ewis_df %>%
  write_csv("../data/processed/ewis/ewis_ma_results.csv")


## ANCESTRY-SPECIFIC

read_gxe_anc <- function(bm, e, anc) {
  fn <- paste0("../data/processed/ewis/", bm, "/", e, "_", anc)
  df <- tryCatch({
    fread(fn, data.table=F, stringsAsFactors=F, na.strings=c("NA", "nan", "nane-nan", "")) %>%
      select(SNP=SNPID, Allele1=Non_Effect_Allele, Allele2=Effect_Allele,
	     Beta=`Beta_G-e`, P=P_Value_Interaction) %>%
      mutate(P=as.numeric(P)) %>%
      filter(P < 0.05)
  }, error=function (e) tibble())
  if (e == exposures[1]) print(paste(bm, anc))
  if (nrow(df) > 1) df else NULL
}

ancestries <- c("EUR", "AFR", "EAS", "SAS")

anc_ewis_df <- expand_grid(
  bm = biomarkers,
  exposure = exposures,
  anc = ancestries
) %>%
  mutate(res = pmap(list(.$bm, .$exposure, .$anc), read_gxe_anc)) %>%
  unnest(res)

anc_ewis_df %>%
  write_csv("../data/processed/ewis/ewis_anc_results_nom.csv")

## QQ PLOT

calc_lambda <- function(x, p=0.5){
  # Calculate genomic inflation lambda value
  x = x[!is.na(x)]
  x.quantile <- quantile(x, p)
  round(qchisq(1 - x.quantile, 1) / qchisq(p, 1), 2)
}

make_qq <- function(data, pval_col, main=""){
  # Make a quantile-quantile plot
  data <- filter(data, data[[pval_col]] > 0)  # In case extremely low p-values are stored as zero
  
  # Process p-values
  y_vals <- sort(-log10(data[[pval_col]]))
  x_vals <- -log10(rev(ppoints(length(y_vals))))  # ppoints generates a uniform probability distribution
  
  # Trim points at higher p-values (credit to RaMWAS package for code snippet)
  levels = as.integer((x_vals - x_vals[1]) / (tail(x_vals, 1) - x_vals[1]) * 2000)
  keep = c(TRUE, diff(levels) != 0)
  levels = as.integer((y_vals - y_vals[1])/(tail(y_vals, 1) - y_vals[1]) * 2000)
  keep = keep | c(TRUE, diff(levels) != 0)
  keep = which(keep)
  
  par(ps=18)
  plot(x=x_vals[keep], y=y_vals[keep], 
       xlab=expression(-log[10](italic(p)) * " (Expected)"), 
       ylab=expression(-log[10](italic(p)) * " (Observed)"),
       main=main, cex=0.8, 
       cex.lab=0.8, cex.main=0.9, 
       pch=16, ylim=c(0, ceiling(max(y_vals))))
  abline(0, 1, lty=2)
  legend(x='topleft', y='topleft',
         bquote(lambda == .(calc_lambda(data[[pval_col]]))), 
         cex=0.9, bty="n")
}

pdf(file="../data/processed/ewis/all_ewis_ma_results_qq.pdf")
ma_qq <- make_qq(ewis_df, "P")
dev.off()

primary_vqtl_df <- read_csv("../data/processed/ewis/ewis_variant_data.csv") %>%
  filter(ancestry == "MA", P < 2.38e-7) %>%
  select(bm, SNP=index_var)
pdf(file="../data/processed/ewis/stage2_ewis_ma_results_qq.pdf")
ma_qq <- ewis_df %>%
  inner_join(primary_vqtl_df, by=c("bm", "SNP")) %>%
  make_qq("P")
dev.off()

