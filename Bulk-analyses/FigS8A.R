library("DESeq2")
library("ggplot2")
library("pheatmap")
library("limma")
library("sva")

# Supplementary Fig. 8Ai
bckCountTable <- read.table("FigS8A.txt", header = TRUE)
x = bckCountTable[,1]
loss = bckCountTable[,-1]
rownames(loss) = make.names(x,unique = TRUE)
head(loss)
bckCountTable <- loss
colnames_bck <- colnames(bckCountTable)
samples <- data.frame(row.names = colnames_bck,
condition = as.factor(c(rep("WNTi HE",3),rep("RAi HE",3),rep("RAd HE",3),rep("RAi HPC",3),rep("RAd HPC",3))),
batch = as.factor(c(rep("4",1),rep("4",1),rep("4",1),rep("1",1),rep("2",1),rep("3",1),rep("1",1),rep("2",1),rep("3",1),
rep("1",1),rep("2",1),rep("3",1),rep("1",1),rep("2",1),rep("3",1))))
samples
bckCountTable_corr <- ComBat_seq(as.matrix(bckCountTable),samples$batch)
bckCDS <- DESeqDataSetFromMatrix(countData = bckCountTable_corr, colData = samples, design = ~ condition)
DESeq.ds <- DESeq(bckCDS)
DESeq.rlog <- rlogTransformation(DESeq.ds, blind = TRUE)
pdf()
plotPCA(DESeq.rlog,ntop=300)
dev.off()
pcaData <- plotPCA(DESeq.rlog,ntop=300, returnData = TRUE)
write.csv(pcaData,file = "pcaData.csv")


# Supplementary Fig. 8Aii
wntd <- bckCountTable[,4:9]
colnames_bck <- colnames(wntd)
samples <- data.frame(row.names = colnames_bck,
condition = as.factor(c(rep("RAi HE",3),rep("RAd HE",3))),
batch = as.factor(c(rep("1",1),rep("2",1),rep("3",1),rep("1",1),rep("2",1),rep("3",1))))
wntd_corr <- ComBat_seq(as.matrix(wntd),samples$batch)
bckCDS <- DESeqDataSetFromMatrix(countData = wntd_corr, colData = samples, design = ~ condition)
DESeq.ds <- DESeq(bckCDS)
DESeq.rlog <- rlogTransformation(DESeq.ds, blind = TRUE)
mat <- assay(DESeq.rlog)
sampleDists_corr <- dist(t(mat))
sampleDistMatrix_corr <- as.matrix(sampleDists_corr)
rownames(sampleDistMatrix_corr) <- paste(rownames(sampleDistMatrix_corr))
colnames(sampleDistMatrix_corr) <- paste(colnames(sampleDistMatrix_corr))
pdf(width = 3.5, height = 2.75)
pheatmap(sampleDistMatrix_corr, clustering_distance_rows=sampleDists_corr,clustering_distance_cols=sampleDists_corr, display_numbers = TRUE, treeheight_row = 10, treeheight_col = 10)
dev.off()


# Supplementary Table 5B
bck_res <- results(DESeq.ds)
res_ordered <- bck_res[order(bck_res$padj),]
write.csv(res_ordered, file="Table5B.csv")
