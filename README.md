# mQTL-enrichment analysis

This page describes how to test whether mQTLs are enriched in GWAS SNPs with small P-values against type 1 diabetes

Requires R.3.3.3

Required packages for this analysis:

	source("https://bioconductor.org/biocLite.R")
	biocLite("VariantAnnotation")
	biocLite("GenomicFeatures")
	biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
	biocLite("org.Hs.eg.db")
	
Also install package "TwoSampleMR" which is avaiable here https://github.com/MRCIEU/TwoSampleMR
```
library(devtools)
install_github("MRCIEU/TwoSampleMR")
```
Load packages
```
library(org.Hs.eg.db)
library(VariantAnnotation)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TwoSampleMR)
library(MRInstruments)
```

Load GWAS summary statistics containing rsid, GWAS P-values, SNP position, MAF, and pre-calculated LDscore for the CEU population. LDscore can be downloaded from https://github.com/bulik/ldsc
```
JohnToddgwas_hg19_ldscore<-read.csv("JohnToddgwas_hg19_ldscore.csv",he=T,sep="")
```
Find existing ALSPAC adolescent mQTLs from the 'TwoSampleMR' package
```
data("aries_mqtl")
cpglist<-subset(aries_mqtl,age=="Adolescence")
cpglist<-subset(cpglist,!duplicated(SNP))
```
Create another column and split SNPs into ALSPAC mQTLs and null SNPs
```
JohnToddgwas_hg19_ldscore$group<-"null"
JohnToddgwas_hg19_ldscore$group[JohnToddgwas_hg19_ldscore$SNP %in% cpglist$SNP]<-"mQTL"
```
We run the enrichment analysis only using mQTLs with strong genetic effects on DNA methylation, we thus set the Bonferroni threshold to P<1e-14 and we name this group of mQTLs "mQTL_test"
```
JohnToddgwas_hg19_ldscore$group[JohnToddgwas_hg19_ldscore$SNP %in% cpglist$SNP[which(cpglist$pval<1e-14)]]<-"mQTL_test"
```
this gives n=4562 target mQTLs for enrichment analysis

We match null SNPs to mQTLs in two ways. First, we match them via similar SNP features.  



#Match null SNPs to mQTLs via same gennomic annotations







