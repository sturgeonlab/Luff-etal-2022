library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("limma")

#Reads in a count table/matrix of genes in rows with samples in their own columns

bckCountTable <- read.table("RAi-Pr1.txt", header=TRUE) #A
bckCountTable <- read.table("RAi-SP.txt", header=TRUE) #B
bckCountTable <- read.table("RAd-Pr1.txt", header=TRUE) #C
bckCountTable <- read.table("RAd-SP.txt", header=TRUE) #D
bckCountTable <- read.table("HPC-Pr1-SP.txt", header=TRUE) #E
x=bckCountTable[,1]
loss=bckCountTable[,-1]
rownames(loss)=make.names(x,unique=TRUE)
head(loss)
bckCountTable <- loss


#making data frame to hold the sample names and the condition types/number of different conditions
colnames_bck <- colnames(bckCountTable)

samples <- data.frame(row.names = colnames_bck,
condition=as.factor(c(rep("RAi-HPC",3),rep("AGM-Pr1-90neg",1)))) #A

samples <- data.frame(row.names = colnames_bck,
condition=as.factor(c(rep("RAi-HPC",3),rep("AGM-Pr1-90pos",1)))) #B

samples <- data.frame(row.names = colnames_bck,
condition=as.factor(c(rep("RAd-HPC",3),rep("AGM-Pr1-90neg",1)))) #C

samples <- data.frame(row.names = colnames_bck,
condition=as.factor(c(rep("RAd-HPC",3),rep("AGM-Pr1-90pos",1)))) #D

samples <- data.frame(row.names = colnames_bck,
condition=as.factor(c(rep("RAi-HPC",3),rep("RAd-HPC",3),rep("AGM-Pr1-90neg",1),rep("AGM-SP-90pos",1)))) #E


#Creates a Summarized Experiment (SE) object using the matrix and design objects
bckCDS <- DESeqDataSetFromMatrix(countData = bckCountTable, colData= samples, design= ~ condition )

#performs the differential expression analysis
bckCDS_1 <- DESeq(bckCDS)

# generates differential gene expression list with log2 fold change, standard error, p value, and p adjusted
bck_res <- results(bckCDS_1)

#examine the combinations of conditions/genotype/batch that can be examined
resultsNames(bckCDS_1)

#makes a differential gene expression list using the condition names (for multiple samples sets)
bck_res <- results(bckCDS_1)

#order results by p value
res_ordered <- bck_res[order(bck_res$padj),]

#displays summary of genes with padj from the results function (default .1)
summary(res_ordered)

#Save the differential gene expression list to table
write.csv(res_ordered, file="RAiHPC-Pr1.csv")
write.csv(res_ordered, file="RAiHPC-SP.csv")
write.csv(res_ordered, file="RAdHPC-Pr1.csv")
write.csv(res_ordered, file="RAdHPC-SP.csv")



# calculate the rlog dispersion for data set, alternative to VSD but slower
rld <- rlog(bckCDS_1, blind=FALSE)

# calculate the varaiance stabilizing transformation, for sample to sample distance heatmap/PCA
vsd <- vst(bckCDS_1, blind=FALSE)

#generate PCA plot
pdf()
plotPCA(vsd, intgroup=c("condition"))
dev.off()
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
write.csv(pcaData, file="HPC-PCA.csv")
