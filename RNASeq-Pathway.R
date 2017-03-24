## Install Bioconductor

source("http://bioconductor.org/biocLite.R")

## Install Required Packages
biocLite(c("pathview", "gage", "gageData", "GenomicAlignments", "TxDb.Hsapiens.UCSC.hg19.knownGene"))

## Load all necessary libraries
library(GenomicAlignments)
library(pathview)

## Load hg19 gene and add read counts
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exonByGene <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")

## Load up Gene Alignments
tophat_fls <- list.files("tophat_all/", pattern="bam$", full.names =T)
bamfls <- BamFileList(tophat_fls)
flag <- scanBamFlag(isSecondaryAlignment=FALSE, isProperPair=TRUE)
param <- ScanBamParam(flag=flag)

## Give a summary of Overlaps 
gnCnt <- summarizeOverlaps(exonByGene, bamfls, mode="Union",ignore.strand=TRUE, singleEnd=FALSE, param=param)
hnrnp.cnts=assay(gnCnt)

## Start Preprocessing 

require(gageData)
data(hnrnp.cnts) 
cnts=hnrnp.cnts
sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
cnts.norm=log2(cnts.norm+8)

## Pathway Visualization

#differential expression: log2 ratio or fold change, uppaired samples
cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
#up-regulated pathways (top 3) visualized by pathview
sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 & !is.na(cnts.kegg.p$greater[,"q.val"])
path.ids <- rownames(cnts.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = cnts.d, pathway.id = pid, species = "hsa"))

