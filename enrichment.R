# mQTL enrichment analysis for Data1

#load GWAS summary statistics (Data1)


library(data.table)
gwasdata<-fread("meta_all.txt",fill=TRUE)
