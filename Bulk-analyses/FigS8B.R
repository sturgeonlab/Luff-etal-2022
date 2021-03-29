library("DESeq2")
library("ggplot2")
library("pheatmap")
library("limma")

bckCountTable <- read.table("FigS8B.txt", header = TRUE)
x = bckCountTable[,1]
loss = bckCountTable[,-1]
rownames(loss) = make.names(x,unique = TRUE)
head(loss)
bckCountTable <- loss

# RAi HE vs AGM CD34
RAi.34a <- bckCountTable[,1:3]
RAi.34b <- bckCountTable[,7]
RAi.34 <- cbind(RAi.34a,RAi.34b)
colnames_bck <- colnames(RAi.34)
samples <- data.frame(row.names = colnames_bck, condition = as.factor(c(rep("RAi",3),rep("CD34",1))))
bckCDS <- DESeqDataSetFromMatrix(countData = RAi.34, colData = samples, design = ~ condition)
DESeq.ds <- DESeq(bckCDS)
bck_res <- results(DESeq.ds)
res_ordered <- bck_res[order(bck_res$padj),]
write.csv(res_ordered, file="RAi-CD34.csv")

# RAd HE vs AGM CD34
RAd.34 <- bckCountTable[,4:7]
colnames_bck <- colnames(RAd.34)
samples <- data.frame(row.names = colnames_bck, condition = as.factor(c(rep("RAd",3),rep("CD34",1))))
bckCDS <- DESeqDataSetFromMatrix(countData = RAd.34, colData = samples, design = ~ condition)
DESeq.ds <- DESeq(bckCDS)
bck_res <- results(DESeq.ds)
res_ordered <- bck_res[order(bck_res$padj),]
write.csv(res_ordered, file="RAd-CD34.csv")

len <- read.table("len.txt", header = TRUE)
mcols(DESeq.ds) <- cbind(mcols(DESeq.ds), len)
fpkm <- fpkm(DESeq.ds)
write.csv(fpkm, file = "fpkm.csv")
