library("DESeq2")
library("ggplot2")
library("limma")

bckCountTable <- read.table("FigS3.txt", header = TRUE)
x = bckCountTable[,1]
loss = bckCountTable[,-1]
rownames(loss) = make.names(x, unique = TRUE)
head(loss)
bckCountTable <- loss
colnames_bck1 <- colnames(bckCountTable)
samples <- data.frame(row.names = colnames_bck1, condition = as.factor(c(rep("CXCR4.neg",3),rep("CXCR4.pos",3))))
bckCDS <- DESeqDataSetFromMatrix(countData = bckCountTable, colData = samples, design = ~ condition)
DESeq.ds <- DESeq(bckCDS)
bck_res <- results(DESeq.ds)
res_ordered <- bck_res[order(bck_res$padj),]
write.csv(res_ordered, file="CXCR4neg-pos.csv")

vsd <- vst(DESeq.ds, blind = FALSE)
plotPCA(vsd, intgroup = c("condition"))
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
write.csv(pcaData,file = "pcaData.csv")




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
[28] tibble_3.0.4           annotate_1.62.0        farver_2.0.3          
[31] generics_0.0.2         ellipsis_0.3.1         withr_2.4.0           
[34] nnet_7.3-12            survival_3.1-8         magrittr_1.5          
[37] crayon_1.3.4           memoise_1.1.0          foreign_0.8-75        
[40] tools_3.6.2            data.table_1.12.8      lifecycle_0.2.0       
[43] stringr_1.4.0          locfit_1.5-9.1         munsell_0.5.0         
[46] cluster_2.1.0          AnnotationDbi_1.46.1   compiler_3.6.2        
[49] rlang_0.4.10           grid_3.6.2             RCurl_1.98-1.1        
[52] rstudioapi_0.11        htmlwidgets_1.5.1      labeling_0.3          
[55] bitops_1.0-6           base64enc_0.1-3        gtable_0.3.0          
[58] DBI_1.1.0              R6_2.4.1               gridExtra_2.3         
[61] knitr_1.28             dplyr_1.0.2            bit_1.1-15.2          
[64] Hmisc_4.3-1            stringi_1.4.6          Rcpp_1.0.3            
[67] geneplotter_1.62.0     vctrs_0.3.4            rpart_4.1-15          
[70] acepack_1.4.1          png_0.1-7              tidyselect_1.1.0      
[73] xfun_0.12  
