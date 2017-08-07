# SNP  ~ CpG association analysis (EWAS) for T1D GWAS variants and ALSPAC DNA methylation 450k array data
#
# BlueCystal cluster: Linux
# 
# 
module add languages/R-3.2.2-ATLAS
R
# 
#
# library(MatrixEQTL)
#
# create a txt file containing 67 T1D GWAS variants to test, these SNPs have been previously clumped to make sure they are not in LD (r2 <0.1) 
#
cat > t1d_snp_proxy.txt
rs2476601
rs6691977
rs3024493
rs478222
rs9653442
rs4849135
rs1990760
rs3087243
rs113010081
rs10517086
rs17388568
rs2611215
rs11954020
rs924043
rs6920220
rs3104636
rs597325
rs1538171
rs1738074
rs7804356
rs62447205
rs4948088
rs10758593
rs61839660
rs10795791
rs41295121
rs11258747
rs722988
rs12416116
rs689
rs694739
rs11065979
rs10492166
rs11170466
rs11171710
rs9585056
rs1456988
rs7140939
rs7149271
rs56994090
rs72727394
rs12908309
rs3825932
rs12927355
rs193778
rs9924471
rs8056814
rs2290400
rs7221109
rs1893217
rs763361
rs34536443
rs12720356
rs402072
rs602662
rs6043409
rs11203202
rs5753037
rs229533
rs2664170
rs2269242
rs78037997
rs12068671
rs10801121
rs118000057
#
# use plink to extract individual level genotype data from ASPAC population
plink --bfile /projects/ALSPAC/studies/originals/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/filtered/bestguess/maf0.01_info0.8/combined/data \
--extract t1d_snpslist_proxy.txt \
--recode A \
--out t1d_snps 
#
#
# plink returned 65 SNPs found, rs9268645 and rs1052553 were not found in the ALSPAC data
# load the returned genotype file into R
#
geno<-read.table("t1d_snps.raw", sep="", header=T)
#
str(geno) # returns 17825 individuals’ genotype on the 65 SNPs
#
# re-format genotype file
rownames(geno)<-geno$FID
#
geno<-geno[,grep("^rs", colnames(geno))]
#
# load samplesheet for ALSPAC 450k methylation data
load("samplesheet.Robj",verbose=TRUE)
#
#
# load ALSPAC 450k DNA methylation array data 
#
load(“meth.Robj”,verbose=T)
#
# load estimated cell counts for 450k array data
# 
counts<-read.table("cell-counts.txt",row.names=1,header=T,sep="\t")
#
# 
# subset samplesheet, methylation data and cell counts to participants aged 15 years and older
#
idx<-which(samplesheet$time_point=="15up")
#
samplesheet<-samplesheet[idx,]
#
meth<-norm.beta.random[,idx]
#
counts<-counts[idx,]
#
geno<-geno[common.ids,]
#
#after matching, there are n=907 individuals with their genotype available 
save(geno,file="genotype.Robj")
#
#
#
#
#
#
################################################
# prepare genotype data for MatrixEQTL
################################################
# convert 'geno' from data frame into matrix in order to transform the matrix, then convert it back to data frame
#
geno1<-as.matrix(geno)
geno1<-t(geno1)
geno1<-as.data.frame(geno1)
#
geno1<-cbind(row.names(geno1),geno1)
row.names(geno1)<-NULL
#
#
geno1<-rbind(colnames(geno1),geno1)
colnames(geno1)<-NULL
#
#
save(geno1,file="genotype.Robj")
# write it into a text file for matrixEQTL
write.table(geno1,file="genotype.txt",row.names=F,col.names=F,sep=”\t”,quote=FALSE )
#
#
#
#
#################################################
# prepare 450k methylation data for MatrixEQTL
#################################################
# 
# transform beta values to M values
#
# first check the distribution of beta values
hist(meth[,10])
#
# define a function to transform beta values to M values
M.value<-function(x){log2(x/(1-x))}
#plot a histogram of a transformed M values 
hist(M.value(meth[,10]))
#
#
# to remove outliers, rank transform M values from bimodal distribution to normal distribution
#
rntransformedM=apply(M.value(meth),2,FUN=function(y) qnorm(p=(rank (y,na.last='keep')-.5)/length(y)))
#
#
# reformat M values
#                    
#
rntransformedM<-rbind(colnames(rntransformedM),rntransformedM) 
colnames(rntransformedM)<-NULL
rntransformedM<-cbind(row.names(rntransformedM),rntransformedM)										 
row.names(rntransformedM)<-NULL                     
# save the rank transformed M values                    
save(rntransformedM,file="rntransformedM.Robj")
#
# write a text file for rank transformed M values
write.table(rntransformedM,file="rntransformedM.txt",row.names=F,col.names=F,sep="\t",quote=FALSE)										 



