# Running Two Sample Mendelian Randomization using SNP ~ CpG associations (EWAS results) and SNP ~ T1D associations (GWAS summary statistics Data 1)
#
#
library(TwoSampleMR)
#
######################################################################
# load mQTLs as instruments, obtained from the previous EWAS analysis
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
