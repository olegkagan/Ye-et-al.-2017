# mQTL enrichment analysis for Data1

#########################################################
# load GWAS summary statistics (Data1) and merge LD score
#########################################################

library(data.table)
gwasdata<-fread("meta_all.txt",fill=TRUE)

# reformat Data1

setnames(gwasdata,old=c("chromosome","position","a0", "a1", "EUR_MAF_1kG",
                        "rsid","majmin","aff.pvalue","aff.OR","aff.se","ill.pvalue",
                        "ill.OR","ill.se","OR.meta","beta.meta","se.meta","z.meta",
                        "p.meta","max_posterior.A","info_score.A","cases_maf.A","controls_maf.A",
                        "max_posterior.I","info_score.I","cases_maf.I","controls_maf.I","qc.check",
                        "warnings","V29"),new=c("chromosome","chr", "positition","a0", "a1", 
                                                "EUR_MAF_1kG","rsid","majmin","aff.pvalue","aff.OR",
                                                "aff.se","ill.pvalue","ill.OR","ill.se","OR.meta","beta.meta",
                                                "se.meta","z.meta","p.meta","max_posterior.A","info_score.A",
                                                "cases_maf.A","controls_maf.A","max_posterior.I","info_score.I",
                                                "cases_maf.I","controls_maf.I","qc.check","warnings"))


# extraction LD score for CEU panel and save in local drive
# LD score can be downloaded from https://github.com/bulik/ldsc


chr1<-read.table(gzfile("1.l2.ldscore.gz"),header=T)
chr2<-read.table(gzfile("2.l2.ldscore.gz"),header=T)
chr3<-read.table(gzfile("3.l2.ldscore.gz"),header=T)
chr4<-read.table(gzfile("4.l2.ldscore.gz"),header=T)
chr5<-read.table(gzfile("5.l2.ldscore.gz"),header=T)
chr6<-read.table(gzfile("6.l2.ldscore.gz"),header=T)
chr7<-read.table(gzfile("7.l2.ldscore.gz"),header=T)
chr8<-read.table(gzfile("8.l2.ldscore.gz"),header=T)
chr9<-read.table(gzfile("9.l2.ldscore.gz"),header=T)
chr10<-read.table(gzfile("10.l2.ldscore.gz"),header=T)
chr11<-read.table(gzfile("11.l2.ldscore.gz"),header=T)
chr12<-read.table(gzfile("12.l2.ldscore.gz"),header=T)
chr13<-read.table(gzfile("13.l2.ldscore.gz"),header=T)
chr14<-read.table(gzfile("14.l2.ldscore.gz"),header=T)
chr15<-read.table(gzfile("15.l2.ldscore.gz"),header=T)
chr16<-read.table(gzfile("16.l2.ldscore.gz"),header=T)
chr17<-read.table(gzfile("17.l2.ldscore.gz"),header=T)
chr18<-read.table(gzfile("18.l2.ldscore.gz"),header=T)
chr19<-read.table(gzfile("19.l2.ldscore.gz"),header=T)
chr20<-read.table(gzfile("20.l2.ldscore.gz"),header=T)
chr21<-read.table(gzfile("21.l2.ldscore.gz"),header=T)
chr22<-read.table(gzfile("22.l2.ldscore.gz"),header=T)

merge<-rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)


# merge the ldscore with GWAS summary statistics Data1

Data1_hg19<-merge(x=merge,y=gwasdata,by.x='SNP',by.y='rsid')


# after merging, 1,167,927 SNPs are avaialble

Data1_hg19_ldscore<-subset(Data1_hg19,select=c("SNP","CHR","BP","CM","MAF","L2","p.meta"))

write.table(Data1_hg19_ldscore,"Data1_hg19_ldscore.csv",sep='\t')



##############################################################################
# retrieve mQTL from adolescent and perform Fisher's combined probability test
##############################################################################



library(TwoSampleMR)
library(MRInstruments)
data("aries_mqtl")
cpglist<-subset(aries_mqtl,age=='Adolescence')

# remove duplicated mQTLs from cpglist

cpglist<-subset(cpglist,!duplicated(SNP))

# merge Aries mQTLs with Data1

available_mQTL<-merge(x=Data1_hg19_ldscore,y=cpglist,by.x='SNP',by.y='SNP')

# obtained 6062 mQTLs

# select mQTLs that have strong exposure effect 
available_mQTL_test<-subset(available_mQTL,pval<1e-14)
str(available_mQTL_test)

# this gives n=4562 mQTLs to test
# perform fisher's combined test
combined.p.for.mQTL<-fishers_combined_test(available_mQTL_test$p.meta)
# this gives pval = 2.312334e-42


#######################################
#matching null SNPs via SNP properties
#######################################

# null SNPs must be further than 1000kb away from the mQTL
# max MAF devidation from reference SNP is 0.02
# LDscore L2 need to be in the same bin

# add a new column 'SNPcat'

Data1_hg19_ldscore$SNPcat<-'null'

#define mQTLs with weak effect and mQTLs with strong effect on CpGs



#weak effect 
Data1_hg19_ldscore$SNPcat[Data1_hg19_ldscore$SNP %in% available_mQTL$SNP]<-'mQTL'

#strong effect
Data1_hg19_ldscore$SNPcat[Data1_hg19_ldscore$SNP %in% available_mQTL_test$SNP]<-'mQTL_test'




#find maximum LDscore
max(Data1_hg19_ldscore$L2)
#305.1991


min(Data1_hg19_ldscore$L2)
#-0.2644381

#find percentiles of LDscore
duration=Data1_hg19_ldscore$L2
quantile(duration,c(.10,.20,.30,.40,.50,.60,.70,.80,.90))



#10%       20%       30%       40%       50%       60%       70%       80%       90% 
#7.012412  9.967493 12.704047 15.535856 18.644254 22.269589 26.798721 33.192342 44.460718



#create 10 LDscore bins according to percentiles
Data1_hg19_ldscore$L2bin<-cut(Data1_hg19_ldscore$L2,
                              breaks=c(0,7.012412,9.967493,12.704047,15.535856,18.644254,22.269589,26.798721,33.192342,44.460718,305.1991),
                              labels=c("group1","group2","group3","group4","group5","group6","group7","group8","group9","group10"))



#define mQTL and null SNPs
mQTL<-subset(Data1_hg19_ldscore,Data1_hg19_ldscore$SNPcat=='mQTL_test')

null<-subset(Data1_hg19_ldscore,Data1_hg19_ldscore$SNPcat=='null')



#define location of null SNPs


location<-(null$CHR==mQTL$CHR) & (null$BP-mQTL$BP>1000000) | (mQTL$BP-null$BP>1000000)

#define MAF range

mafrange<-(null$MAF-mQTL$MAF<0.02) | (mQTL$MAF-null$MAF<0.02)

#define LDscore
ldscore<-null$L2bin == mQTL$L2bin


#select null SNPs

p.values<-Data1_hg19_ldscore$p.meta[location&mafrange&ldscore]


sample.size<-length(mQTL$SNP)


#define permuation function

set.seed(200)
mQTL_enrichment<-function(x){
  p.val.sample<-sample(p.values,sample.size)
  fishers_combined_test(p.val.sample)$pval
}

results<-sapply(1:10000,mQTL_enrichment)

results2<--log10(results)

hist(results2,xlim=c(0,60),xlab='-log10(Combined probability)')
empirical.pval<-sum(results2>=41.63595)/10000
empirical.pval



###########################################
#matching null SNPs via genomic annotation
###########################################

#required packages

source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
biocLite("GenomicFeatures")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
biocLite("org.Hs.eg.db")
library(devtools)
install_github("MRCIEU/TwoSampleMR")
install.packages("qqman")
library(org.Hs.eg.db)

library(VariantAnnotation)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TwoSampleMR)
library(MRInstruments)

#open GWAS Data1 file with merged LD score

Data1_hg19_ldscore<-read.csv("Data1_hg19_ldscore.csv",he=T,sep="")

data("aries_mqtl")
cpglist<-subset(aries_mqtl,age=="Adolescence")
cpglist<-subset(cpglist,!duplicated(SNP))
#
#
# cpglist contained 30087 mQTLs available
#
#
#
#
#



