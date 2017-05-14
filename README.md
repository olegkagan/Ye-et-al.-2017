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

#Need to download and merge LDscore

Load GWAS summary statistics containing GWAS P-values, rsid, SNP position, MAF, and pre-calculated LDscore for the CEU population
```
JohnToddgwas_hg19_ldscore<-read.csv("JohnToddgwas_hg19_ldscore.csv",he=T,sep="")
```
Find existing ALSPAC adolescent mQTLs from the 'TwoSampleMR' package, these mQTLs have strong genetic effect on DNA methylation (filtered with Bonferroni threshold p<1e-14)
```
data("aries_mqtl")
cpglist<-subset(aries_mqtl,age=="Adolescence")
cpglist<-subset(cpglist,!duplicated(SNP))
```
