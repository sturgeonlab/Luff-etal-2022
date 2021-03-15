library("DESeq2")

bckCountTable <- read.table("fCount-master-coding.txt", header = TRUE)
x = bckCountTable[,1]
loss = bckCountTable[,-1]
rownames(loss) = make.names(x, unique = TRUE)
head(loss)
bckCountTable <- loss
colnames_bck <- colnames(bckCountTable)
samples <- data.frame(row.names = colnames_bck, condition = as.factor(c(rep("WNTi meso",4),rep("WNTd meso",4),
rep("NEG",3),rep("POS",3),rep("WNTi HE",3),rep("WNTd HE",3),rep("RAi HE",3),rep("RAd HE",3),rep("AGM PR",1),
rep("AGM CD34",1),rep("AGM HSPC",1),rep("RAi HPC",3),rep("RAd HPC",3),rep("AEC",3))))
bckCDS <- DESeqDataSetFromMatrix(countData = bckCountTable, colData = samples, design = ~ condition)
DESeq.ds <- DESeq(bckCDS)
len <- read.table("lengths.txt", header = TRUE) 
mcols(DESeq.ds) <- cbind(mcols(DESeq.ds), len) #adding gene lengths
fpkm <- fpkm(DESeq.ds)
write.csv(fpkm, file = "fpkm.csv")
