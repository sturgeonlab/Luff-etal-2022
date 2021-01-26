library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("limma")

#Reads in a count table/matrix of genes in rows with samples in their own columns

bckCountTable <- read.table("RAiHE-AGM34.txt", header=TRUE) #A
bckCountTable <- read.table("RAdHE-AGM34.txt", header=TRUE) #B

x=bckCountTable[,1]
loss=bckCountTable[,-1]
rownames(loss)=make.names(x,unique=TRUE)
head(loss)
bckCountTable <- loss


#making data frame to hold the sample names and the condition types/number of different conditions
colnames_bck <- colnames(bckCountTable)

samples <- data.frame(row.names= colnames_bck,
condition=as.factor(c(rep("RAi",3),rep("AGM",1)))) #A

samples <- data.frame(row.names= colnames_bck,
condition=as.factor(c(rep("RAd",3),rep("AGM",1)))) #B


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
write.csv(res_ordered, file="RAi-AGM34.csv") #A
write.csv(res_ordered, file="RAd-AGM34.csv") #B
