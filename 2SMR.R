# Running Two Sample Mendelian Randomization using SNP ~ CpG associations (EWAS results) and SNP ~ T1D associations (GWAS summary statistics Data 1)
#
#
library(TwoSampleMR)
#
######################################################################
# generate the exposure data file from the previous EWAS analysis
######################################################################
#
#
mQTLIV<-read.csv("mQTLproxy.csv",header=T,sep=",", stringsAsFactors = F)
#
# remove trans-associations between mQTL and CpGs, Trans.association was coded as 'Y' 
#
mQTLIV<-subset(mQTLIV, Trans.association=='N')
#
# re-format the table
names(mQTLIV)[names(mQTLIV) =='beta_se']<-"se_exposure"
names(mQTLIV)[names(mQTLIV) =='beta']<-'beta_exposure'
names(mQTLIV)[names(mQTLIV) =='R2']<-'R2_exposure'
names(mQTLIV)[names(mQTLIV) =='pvalue']<-'pvalue_exposure'
names(mQTLIV)[names(mQTLIV) =='FDR']<-'FDR_exposure'
#
#
#
write.table(mQTLIV, "T1D_mQTL_expo.txt",sep=",",row.names=F)
#
#
#################################################################
# load exposure data, obtained from the previous EWAS analysis
#################################################################
#
T1D_mQTL_expo<-read_exposure_data(filename="T1D_mQTL_expo.txt", sep=",", 
                                  snp_col = "SNP",beta_col = "beta_exposure",
                                  se_col="se_exposure",effect_allele_col = "effect_allele_exposure",
                                  other_allele_col='other_allele_exposure', eaf_col = "eat_exposure",
                                  pval_col = "pvalue_exposure", phenotype_col = "gene")
#
# perform LD clumping 
# T1D_mQTL_expo_clump<-clump_data(T1D_mQTL_expo)
# 
# ##########################
# load outcome data (Data1)
############################
#
# outcome_data<-read_outcome_data(snps=T1D_mQTL_expo_clump$SNP, filename="T1D_mQTL_expo.txt", 
                                sep=",",snp_col = "SNP",beta_col="beta.meta",se_col="se.meta",
                                effect_allele_col="a1",other_allele_col = "a0",pval_col="p.meta",
                                eaf_col = "eaf_outcome") # Note that the the odds ratio in Data 1 is the addition of A1 thus A1 is the effect allele
# 
# the outcome data file cannot contain rows with empty values  in the 'effect_allele.outcome' column, to remove empty rows
outcome_data<-subset(outcome_data,mr_keep.outcome=='TRUE')
#
#
#
# harmoziation
data<-harmonise_data(exposure_dat = T1D_mQTL_expo_clump,outcome_dat = outcome_data)
#
#
###########################################
# running 2SMR using the Wald Ratio method
###########################################
#
result<-mr(data)






