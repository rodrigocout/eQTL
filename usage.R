setwd("/Users/Rodrigo/Dropbox/R_scripts/functions4eqtls/")
########################################
###  run eQTL
########################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(annotate)
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(IRanges)
source("./lmEQTL.R")
source("http://bioconductor.org/biocLite.R")
biocLite("S4Vectors")
library(S4Vectors)


#---get markers
lMarkers<-read.table("./mirsnps_proxies_breastcancer_affy_snp6_coord.txt", 
                     header = TRUE, stringsAsFactors = FALSE)

#---genomic annotatio for snps
annotsnps <- as(lMarkers, "GRanges")

#---genomic annotation for genes
entrez<-c("5567","148418","58511","260425","64858")
annotgene<-getTxDb(entrez, getTSS=TRUE, removeDuplicates=TRUE)

#---get cis gene-to-snp candicate map
geneSnpList<-gene2snpCandidates(annotgene,annotsnps)

#---load  gxdata, gtypedata and gxids
#apenas um exemplo para demonstrar dados de entrada, nao gera resultado significativo.
load(file="sampledata.RData")

#---run eQTL
eqtls<-lmEQTL(geneSnpList=geneSnpList, gxdata=gxdata, gtypedata=gtypedata, gxids=gxids, cutoff=1)
#write.table(eqtls,file="eqtlsCohort.txt", quote=TRUE)

########################################
###  plot eQTL
########################################
ploteqtl(geneid="5567", snpid="rs2134646", gxdata=gxdata, gtypedata=gtypedata,
         eqtls=eqtls, plotpdf=FALSE)


