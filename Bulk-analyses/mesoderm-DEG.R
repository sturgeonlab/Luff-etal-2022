library("DESeq2")
library("ggplot2")
library("limma")

#forces uniqueness of first column and sets as row names
bckCountTable <- read.table("mesodermcounts.txt", header=TRUE)
x=bckCountTable[,1]
loss=bckCountTable[,-1]
rownames(loss)=make.names(x,unique=TRUE)
head(loss)
bckCountTable <- loss

table1 <- bckCountTable[,5:10]
table2 <- bckCountTable[,1:7]
table3a <- bckCountTable[,1:4]
table3b <- bckCountTable[,8:10]
table3 <- cbind(table3a,table3b)

#making data frame to hold the sample names and the condition types/number of different conditions
colnames_bck1 <- colnames(table1)
colnames_bck2 <- colnames(table2)
colnames_bck3 <- colnames(table3)
samples1 <- data.frame(row.names = colnames_bck1, condition=as.factor(c(rep("CXCR4.pos",3),rep("CXCR4.neg",3))))
samples2 <- data.frame(row.names = colnames_bck2, condition=as.factor(c(rep("CD235a.pos",4),rep("CXCR4.pos",3))))
samples3 <- data.frame(row.names = colnames_bck3, condition=as.factor(c(rep("CD235a.pos",4),rep("CXCR4.neg",3))))


#Creates a Summarized Experiment (SE) object using the matrix and design objects
bckCDS1 <- DESeqDataSetFromMatrix(countData = table1, colData = samples1, design = ~condition)
bckCDS2 <- DESeqDataSetFromMatrix(countData = table2, colData = samples2, design = ~condition)
bckCDS3 <- DESeqDataSetFromMatrix(countData = table3, colData = samples3, design = ~condition)

#performs the differential expression analysis
bckCDS1_1 <- DESeq(bckCDS1)
bckCDS2_2 <- DESeq(bckCDS2)
bckCDS3_3 <- DESeq(bckCDS3)

# generates differential gene expression list with log2 fold change, standard error, p value, and p adjusted
bck_res1 <- results(bckCDS1_1)
bck_res2 <- results(bckCDS2_2)
bck_res3 <- results(bckCDS3_3)

#order results by p value
res_ordered1 <- bck_res1[order(bck_res1$padj),]
res_ordered2 <- bck_res2[order(bck_res2$padj),]
res_ordered3 <- bck_res3[order(bck_res3$padj),]

#Save the differential gene expression list to table
write.csv(res_ordered1, file="CXCR4pos-CXCR4neg.csv")
write.csv(res_ordered2, file="CXCR4pos-CD235a.csv")
write.csv(res_ordered3, file="CXCR4neg-CD235a.csv")


### generating PCA plot of WNTd populations and export coordinates
rld <- rlog(bckCDS_1, blind=FALSE)
vsd <- vst(bckCDS_1, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
write.csv(pcaData, file="WNTd-PCA.csv")




> sessionInfo()
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] limma_3.40.6                ggplot2_3.3.2              
 [3] DESeq2_1.24.0               SummarizedExperiment_1.14.1
 [5] DelayedArray_0.10.0         BiocParallel_1.18.1        
 [7] matrixStats_0.55.0          Biobase_2.44.0             
 [9] GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
[11] IRanges_2.18.3              S4Vectors_0.22.1           
[13] BiocGenerics_0.30.0        

loaded via a namespace (and not attached):
 [1] bit64_0.9-7            splines_3.6.2          Formula_1.2-3         
 [4] latticeExtra_0.6-29    blob_1.2.1             GenomeInfoDbData_1.2.1
 [7] pillar_1.4.3           RSQLite_2.2.0          backports_1.1.5       
[10] lattice_0.20-40        glue_1.4.2             digest_0.6.24         
[13] RColorBrewer_1.1-2     XVector_0.24.0         checkmate_2.0.0       
[16] colorspace_1.4-1       htmltools_0.4.0        Matrix_1.2-18         
[19] XML_3.99-0.3           pkgconfig_2.0.3        genefilter_1.66.0     
[22] zlibbioc_1.30.0        purrr_0.3.3            xtable_1.8-4          
[25] scales_1.1.0           jpeg_0.1-8.1           htmlTable_1.13.3      
[28] tibble_3.0.4           annotate_1.62.0        generics_0.0.2        
[31] ellipsis_0.3.1         withr_2.4.0            nnet_7.3-12           
[34] survival_3.1-8         magrittr_1.5           crayon_1.3.4          
[37] memoise_1.1.0          foreign_0.8-75         tools_3.6.2           
[40] data.table_1.12.8      lifecycle_0.2.0        stringr_1.4.0         
[43] locfit_1.5-9.1         munsell_0.5.0          cluster_2.1.0         
[46] AnnotationDbi_1.46.1   compiler_3.6.2         rlang_0.4.10          
[49] grid_3.6.2             RCurl_1.98-1.1         rstudioapi_0.11       
[52] htmlwidgets_1.5.1      bitops_1.0-6           base64enc_0.1-3       
[55] gtable_0.3.0           DBI_1.1.0              R6_2.4.1              
[58] gridExtra_2.3          knitr_1.28             dplyr_1.0.2           
[61] bit_1.1-15.2           Hmisc_4.3-1            stringi_1.4.6         
[64] Rcpp_1.0.3             geneplotter_1.62.0     vctrs_0.3.4           
[67] rpart_4.1-15           acepack_1.4.1          png_0.1-7             
[70] tidyselect_1.1.0       xfun_0.12             
