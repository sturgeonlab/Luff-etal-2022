library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("limma")

#Reads in a count table/matrix of genes in rows with samples in their own columns
#Sometimes this has problems, see below for table manipulation to fix duplicate gene names and other nonsense
#bckCountTable <- read.table("all.gene_counts-meso.txt", header=TRUE, row =1)
#bckCountTable <- read.table("all.gene_counts-WNTd.txt", header=TRUE, row =1)
#head(bckCountTable)

#Sometimes this manipulation is necessary if R or DEseq2 complains about non unquie genes
#forces uniqueness of first column and sets as row names
bckCountTable <- read.table("all.gene_counts-WNTd.txt", header=TRUE) #A (Fig. S3D)
bckCountTable <- read.table("all.gene_counts_HE.txt", header=TRUE) #B (Fig. S8Aii)
x=bckCountTable[,1]
loss=bckCountTable[,-1]
rownames(loss)=make.names(x,unique=TRUE)
head(loss)
bckCountTable <- loss


#making data frame to hold the sample names and the condition types/number of different conditions
colnames_bck <- colnames(bckCountTable)

samples <- data.frame(row.names = colnames_bck,
condition=as.factor(c(rep("CXCR4neg",3),rep("CXCR4pos",3)))) #A

samples <- data.frame(row.names = colnames_bck,
condition=as.factor(c(rep("IWP HE",3),rep("RAi HE",3),rep("RAd HE",3)))) #B


#Creates a Summarized Experiment (SE) object using the matrix and design objects
bckCDS <- DESeqDataSetFromMatrix(countData = bckCountTable, colData= samples, design= ~ condition )

#performs the differential expression analysis
bckCDS_1 <- DESeq(bckCDS)

# calculate the varaiance stabilizing transformation, for sample to sample distance heatmap/PCA
vsd <- vst(bckCDS_1, blind=FALSE)

#generate PCA plot
pdf()
plotPCA(vsd, intgroup=c("condition"))
dev.off()
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
write.csv(pcaData,file="pcaData.csv")
