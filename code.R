###code for Study 1 and 3

library(TwoSampleMR)
library(data.table)
exposure <- read_exposure_data(filename = "xxx.csv",sep = ",",snp_col = "SNP",
                               beta_col = "Beta",se_col = "Se",effect_allele_col = "Effect_allele",
                               other_allele_col = "Other_allele",eaf_col = "Effect_allele_frequency",
                               pval_col = "P",phenotype_col ="Gene_transcript",samplesize_col = "N",
                               clump=TRUE)
outcome <- read_outcome_data(snps = exposure$SNP,filename = "xxx.txt",sep = "\t",
                             snp_col = "SNP",beta_col = "Beta",se_col = "SE",
                             effect_allele_col = "alt",other_allele_col = "ref",
                             eaf_col = "af_alt", pval_col = "pval")
outcome <- extract_outcome_data(snps = exposure$SNP,outcomes = 'GWAS ID')

dat <- harmonise_data(exposure_dat = exposure,outcome_dat = outcome)
result <- mr(hm,method_list = c('mr_wald_ratio',"mr_ivw","mr_egger_regression","mr_weighted_median"))
or <- generate_odds_ratios(mr_res = result)



###code for Study 2
#sensitivity analysis
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(hm))
h <- mr_heterogeneity(dat)
p <- mr_pleiotropy_test(dat)

#colocalization analysis
library(coloc)
library(dplyr)
library(tidyr)

coloc<-coloc.abf(
  dataset1 = list(snp=input$SNP,
                  pvalues=input$pval.exposure,
                  type="quant",
                  N=input$samplesize.exposure),
  dataset2 = list(snp=input$SNP,
                  pvalues=input$pval.outcome,
                  type="quant",
                  N=input$samplesize.outcome),
  MAF = input$eaf.exposure)%>%.$results%>%filter(SNP.PP.H4>0.7)