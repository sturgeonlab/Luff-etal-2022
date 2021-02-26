#load in the counts file
readcountsRaw <- read.delim("fCount.txt", header = TRUE, row= 1, comment.char="#")
fcounts <- readcountsRaw[,-(1:5)]
#produce trimmed counts file
write.csv(fcounts[-1,],"fCount_genes_only.csv")

#load in the counts file
readcountsRaw <- read.delim("fCount.txt", header = TRUE, comment.char="#")
GeneNamesRaw <- readcountsRaw[,1]
RPKM <- readcountsRaw
RPKM <- RPKM[,-c(1:5)]
lengths <- RPKM[,1]
RPKM <- RPKM[,-1]
RPKM <- as.matrix(RPKM)
#calculate RPKMS
RPKM <- apply(RPKM, 2, function(x) x/sum(x))
RPKM <- RPKM/lengths
RPKM <- RPKM*1000000000
row.names(RPKM) <- GeneNamesRaw
write.csv(RPKM,file=paste("RPKM.csv",sep=""))
