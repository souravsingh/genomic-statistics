install.packages(c("devtools","gplots"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","org.Hs.eg.db","AnnotationDbi"))
biocLite("alyssafrazee/RSkittleBrewer")

library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)



# Load Bodymap dataset from Bowtie archives

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
p_data=pData(bm)
e_data=exprs(bm)
f_data = fData(bm)
ls()

# Load Tables on Gender and Race
table(p_data$gender)
table(p_data$gender,p_data$race)

# Show a summary of the dataset
summary(e_data)

# Display Missing Age Values

table(p_data$age,useNA="ifany")
table(is.na(p_data$age))
sum(p_data$age==" ")

# Check genomic data for NAs
sum(is.na(e_data))

# Make the distribution of NA's by genes

gene_na = rowSums(is.na(e_data))
table(gene_na)

# Make the distribution of NA's by samples

sample_na = rowSums(is.na(e_data))
table(sample_na)


# Display Dimensions of 
dim(f_data)
dim(p_data)
dim(e_data)

# Make a Boxplot of Expression data
boxplot(log2(e_data+1),col=2,range=0)

# Create Histograms
par(mfrow=c(1,2))
hist(log2(e_data[,1]+1),col=2)
hist(log2(e_data[,2]+1),col=2)

# Display log densities
plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=3)

# Make qqplot
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3)

# Convert Expression Data to Dataframe, filter by rowmeans and create boxplot
e_data = as.data.frame(e_data)
filt_edata = filter(edata,rowMeans(edata) > 1)
boxplot(as.matrix(log2(filt_edata+1)),col=2)

aeid = as.character(fdata[,1])
chr = AnnotationDbi::select(org.Hs.eg.db,keys=aeid,keytype="ENSEMBL",columns="CHR")
head(chr)

dim(chr)
dim(edata)
chr = chr[!duplicated(chr[,1]),]

all(chr[,1] == rownames(edata))

edatay = dplyr::filter(edata,chr$CHR=="Y")

boxplot(colSums(edatay) ~ pdata$gender)
points(colSums(edatay) ~ jitter(as.numeric(pdata$gender)),
        col=as.numeric(pdata$gender),
        pch=19)


ematrix = as.matrix(edata)[rowMeans(edata) > 10000,]
heatmap(ematrix)

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(ematrix,col=colramp)

heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA)

heatmap.2(ematrix,col=colramp,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none")
