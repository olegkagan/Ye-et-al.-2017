# SNP  ~ CpG association analysis (EWAS)
#
# BlueCystal cluster: linux
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
