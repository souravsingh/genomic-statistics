
install.packages(c("devtools","MatrixEQTL"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","goseq","DESeq2"))

library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)

## Check for supported genomes and Gene ID.
head(supportedGenomes())
head(supportedGeneIDs())

##  Read data.
temp_data =read.table(system.file("extdata","Li_sum.txt",
                                     package="goseq"),sep="\t",
                                     header=TRUE,
                                     stringsAsFactors=FALSE)
## Filter dataset based on Row Means and Clean them.

expression= temp_data[,-1]
rownames(expression) = temp_data[,1]
expression = expression[rowMeans(expression) > 5,]
grp=factor(rep(c("Control","Treated"),times=c(4,3)))
pdata  = data.frame(grp)

## USe DESeq package for Differential Expression Analysis.

de = DESeqDataSetFromMatrix(expr, pdata, ~grp)
de_fit = DESeq(de)
de_results = results(de_fit)

## remove NA values
genes = as.integer(de_results$padj < 0.05)
not_na = !is.na(genes)
names(genes) = rownames(expr)
genes = genes[not_na]

## Show the suppported genomes in dataset

head(supportedGenomes(),n=12)[,1:4]

## ------------------------------------------------------------------------
pwf=nullp(genes,"hg19","ensGene")
head(pwf)

## ------------------------------------------------------------------------
GO.wall=goseq(pwf,"hg19","ensGene")
head(GO.wall)

## ------------------------------------------------------------------------
GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
head(GO.MF)
