# Run phenotyping script
qsub run_vqtl_phenos.sh

# Run PHESANT and attach phenotypes
qsub 4_prep_ewis_phenos.sh

# Determine number of effective biomarkers and exposures
qsub Neff.sh

# Run vQTLs
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for chr in {1..22}; do qsub run_vqtl.sh $bm EUR $chr 6; done; done

# Postprocess
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub postprocess_vqtl.sh $bm $anc; done; done

# Meta-analyze vQTLs
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub vqtl_meta_analysis.sh $bm; done

# Copy new ewis variant file to cluster 
rsync -avP ../data/processed/ewis/ewis_variants.txt uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/
rsync -avP ../data/processed/ewis/ewis_variant_data.csv uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/
rsync -avP ../data/processed/ewis/ME_variants.txt uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/
rsync -avP ../data/processed/ewis/ME_variant_data.csv uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/

# Subset to vQTL SNPs
qsub subset_vqtls.sh
qsub subset_MEs.sh

# Run EWIS
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis.sh $bm $anc 1 500; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis.sh $bm $anc 501 1000; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis.sh $bm $anc 1001 1500; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis.sh $bm $anc 1501 2000; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis.sh $bm $anc 2001 3000; done; done

cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis_MEs.sh $bm $anc 1 500; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis_MEs.sh $bm $anc 501 1000; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis_MEs.sh $bm $anc 1001 1500; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis_MEs.sh $bm $anc 1501 2000; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do qsub run_ewis_MEs.sh $bm $anc 2001 3000; done; done

# Meta-analyze EWIS
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub ewis_meta_analysis.sh $bm; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub ewis_meta_analysis_MEs.sh $bm; done

# Collect EWIS results
qsub collect_ewis_results.sh
qsub collect_ewis_results_MEs.sh

# vQTL follow-up stuff
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub prep_vqtl_gencorrs.sh $bm; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm1; do cat ../data/processed/metabolic_biomarkers.txt | while read bm2; do qsub vqtl_gencorrs.sh $bm1 $bm2; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub vqtl_sensitivity.sh $bm; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub vqtl_zhang_check.sh $bm 21001; done
./3j_run_annovar.sh  # Locally
./3k_fetch_nonsig_MEs.sh  # Interactively on HPC
./7b_explore_anc_signals.sh  # Interactively on HPC
# Run FORGE2 tool for epigenomic enrichment on vQTL and ME variants: e.g. pbcopy < ../data/processed/ewis/vqtl_rsids.txt --> mv ~/Downloads/Unnamed.GWAS.erc2-DHS.chart.tsv.gz ../data/processed/forge2/vqtl_variants_DHS.chart.tsv.gz

# All-exposure adjustment sensitivity analyses genome-wide
rsync -avP ../data/processed/ewis/significant_exposures.txt uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do ./8_prep_expAdj_phenos.sh $bm $anc; done; done
# Genome-wide version....
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for chr in {1..22}; do qsub run_vqtl_expAdj.sh $bm EUR $chr; done; done  # Then repeat with lower resources for other ancestries
for chr in {1..22}; do echo "chr${chr}"; ls -lh ../data/processed/sensitivity/*_chr${chr}_EUR_expAdj.vqtl | wc -l; echo "\n"; done  # Quality check for EUR vQTL runs
for chr in {1..22}; do ls -lh ../data/processed/sensitivity/*_chr${chr}_EUR_expAdj.vqtl; echo "\n"; done  # Quality check for EUR vQTL runs
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub expAdj_meta_analysis.sh $bm; done
# Sensitivity version...
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do for anc in EUR AFR EAS SAS; do ./8b_run_vqtl_expAdj_sigOnly.sh $bm $anc; done; done
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do ./8c_expAdj_meta_analysis_sigOnly.sh $bm; done

# Prepare and run conditional GEI tests
rsync -avP ../data/processed/ewis/signif* uger:florez_ukb_projects/ukbb-vqtl/data/processed/ewis/
qsub prep_conditional_gxe.sh
qsub conditional_gxe.sh

# Follow-up for vignettes
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub run_ewis_anthroPC.sh $bm; done
qsub collect_ewis_anthroPC_results.sh

# Prepare all summary stats for upload to CMDKP
cat ../data/processed/metabolic_biomarkers.txt | while read bm; do qsub cmdkp_exports.sh $bm; done
