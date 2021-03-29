library("DESeq2")
library("ggplot2")
library("limma")
library("pheatmap")

install.packages("sva_3.39.0.tar.gz", repos = NULL, type="source")
library("sva")

bckCountTable <- read.table("FigS3D.txt", header = TRUE)
x = bckCountTable[,1]
loss = bckCountTable[,-1]
rownames(loss) = make.names(x, unique = TRUE)
head(loss)
bckCountTable <- loss
colnames_bck1 <- colnames(bckCountTable)
samples <- data.frame(row.names = colnames_bck1, condition = as.factor(c(rep("CXCR4.neg",3),rep("CXCR4.pos",3))), 
                      batch = as.factor(c(rep("1",1),rep("2",1),rep("3",1),rep("1",1),rep("2",1),rep("3",1))))

bckCountTable_corr <- ComBat_seq(as.matrix(bckCountTable),samples$batch)
bckCDS <- DESeqDataSetFromMatrix(countData = bckCountTable_corr, colData = samples, design = ~ batch + condition)
DESeq.ds <- DESeq(bckCDS)
DESeq.rlog <- rlogTransformation(DESeq.ds, blind = TRUE)
pdf()
plotPCA(DESeq.rlog,ntop=300)
dev.off()
pcaData <- plotPCA(DESeq.rlog,ntop=300, returnData = TRUE)
write.csv(pcaData,file = "pcaData.csv")

# Supplementary Table 3B
bck_res <- results(DESeq.ds)
res_ordered <- bck_res[order(bck_res$padj),]
write.csv(res_ordered, file="Table3B.csv")
