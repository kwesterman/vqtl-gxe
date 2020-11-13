# Perform initial processing of vQTL and main-effect summary stats
# (formatting, INFO filter, subsetting) and create a basic QQ plot


library(tidyverse)
library(data.table)


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

args <- commandArgs(trailingOnly=T)
filepath <- args[1]


vqtl_cols <- c(
  CHR="Chr", SNP="SNP", POS="bp",
  ALT="A1", REF="A2", AF="freq",
  N="NMISS", BETA="beta", SE="se", P="P"
)
me_cols <- c(
  CHR="#CHROM", SNP="ID", POS="POS",
  ALT="ALT", REF="REF",
  N="OBS_CT", BETA="BETA", SE="SE", P="P"
)
use_cols <- if (grepl("vqtl_ss", filepath)) vqtl_cols else me_cols

high_qual_variants <- scan("../data/processed/info0.5_variants_alleleSwitch.txt", 
                           what=character())

ss_df <- fread(filepath, stringsAsFactors = F, data.table = F) %>%
  select(use_cols) %>%
  filter(SNP %in% high_qual_variants)

ss_df %>%
  filter(P < 0.01) %>%
  write_tsv(paste0(filepath, "_nom"))

ss_df %>%
  filter(P < 5e-8) %>%
  write_tsv(paste0(filepath, "_gw"))
  
file_base <- basename(filepath)

plot_filepath <- gsub(basename)
png(gsub(""))
make_qq(ss_df, "P")
dev.off()
