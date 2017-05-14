# mQTL-enrichment analysis

This page describes how to test whether mQTLs are enriched in GWAS SNPs with small P-values against type 1 diabetes

Requires R.3.3.3

Required packages to load:

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

