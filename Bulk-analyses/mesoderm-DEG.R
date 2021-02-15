library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("limma")


#forces uniqueness of first column and sets as row names
bckCountTable <- read.table("mesodermcounts.txt", header=TRUE)
x=bckCountTable[,1]
loss=bckCountTable[,-1]
rownames(loss)=make.names(x,unique=TRUE)
head(loss)
bckCountTable <- loss

table1 <- bckCountTable[1:7]
table2a <- bckCountTable[1:4]
table2b <- bckCountTable[8:10]
table2 <- cbind(table2a,table2b)
table3 <- bckCountTable[5:10]

#making data frame to hold the sample names and the condition types/number of different conditions
colnames_bck1 <- colnames(table1)
colnames_bck2 <- colnames(table2)
colnames_bck3 <- colnames(table3)
samples1 <- data.frame(row.names = colnames_bck1, condition=as.factor(c(rep("CD235a.pos",4),rep("CXCR4.neg",3))))
samples2 <- data.frame(row.names = colnames_bck2, condition=as.factor(c(rep("CD235a.pos",4),rep("CXCR4.pos",3))))
samples3 <- data.frame(row.names = colnames_bck3, condition=as.factor(c(rep("CXCR4.neg",3),rep("CXCR4.pos",3))))

#Creates a Summarized Experiment (SE) object using the matrix and design objects
bckCDS1 <- DESeqDataSetFromMatrix(countData = table1, colData = samples1, design = ~condition)
bckCDS2 <- DESeqDataSetFromMatrix(countData = table2, colData = samples2, design = ~condition)
bckCDS3 <- DESeqDataSetFromMatrix(countData = table3, colData = samples3, design = ~condition)

#performs the differential expression analysis
bckCDS_1 <- DESeq(bckCDS1)
bckCDS_2 <- DESeq(bckCDS2)
bckCDS_3 <- DESeq(bckCDS3)

# generates differential gene expression list with log2 fold change, standard error, p value, and p adjusted
bck_res1 <- results(bckCDS_1)
bck_res2 <- results(bckCDS_2)
bck_res3 <- results(bckCDS_3)

#order results by p value
res_ordered1 <- bck_res1[order(bck_res1$padj),]
res_ordered2 <- bck_res2[order(bck_res2$padj),]
res_ordered3 <- bck_res3[order(bck_res3$padj),]

#Save the differential gene expression list to table
write.csv(res_ordered1, file="CD235-CXCR4neg.csv")
write.csv(res_ordered2, file="CD235-CXCR4pos.csv")
write.csv(res_ordered3, file="CXCR4neg-CXCR4pos.csv")


### generating PCA plot of WNTd populations and export coordinates
rld <- rlog(bckCDS_3, blind=FALSE)
vsd <- vst(bckCDS_3, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
write.csv(pcaData, file="WNTd-mesoderm.csv")
