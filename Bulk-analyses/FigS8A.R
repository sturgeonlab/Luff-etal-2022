library("DESeq2")
library("ggplot2")
library("pheatmap")
library("limma")

install.packages("sva_3.39.0.tar.gz", repos = NULL, type="source")
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


# GSEA parameters for Supplementary Fig. 8Aiii
gsea-cli.sh GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/c5.go.bp.v7.2.symbols.gmt -collapse No_Collapse -mode Max_probe -norm meandiv 
-nperm 1000 -rnk /Table5C.rnk.txt -scoring_scheme weighted -rpt_label rad.rai -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 
-rnd_seed timestamp -set_max 500 -set_min 5 -zip_report false -out /output



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
 [1] sva_3.39.0                  genefilter_1.66.0          
 [3] mgcv_1.8-31                 nlme_3.1-144               
 [5] limma_3.40.6                pheatmap_1.0.12            
 [7] ggplot2_3.3.2               DESeq2_1.24.0              
 [9] SummarizedExperiment_1.14.1 DelayedArray_0.10.0        
[11] BiocParallel_1.18.1         matrixStats_0.55.0         
[13] Biobase_2.44.0              GenomicRanges_1.36.1       
[15] GenomeInfoDb_1.20.0         IRanges_2.18.3             
[17] S4Vectors_0.22.1            BiocGenerics_0.30.0        

loaded via a namespace (and not attached):
 [1] edgeR_3.28.1           bit64_0.9-7            splines_3.6.2         
 [4] Formula_1.2-3          latticeExtra_0.6-29    blob_1.2.1            
 [7] GenomeInfoDbData_1.2.1 pillar_1.4.3           RSQLite_2.2.0         
[10] backports_1.1.5        lattice_0.20-40        glue_1.4.2            
[13] digest_0.6.24          RColorBrewer_1.1-2     XVector_0.24.0        
[16] checkmate_2.0.0        colorspace_1.4-1       htmltools_0.4.0       
[19] Matrix_1.2-18          XML_3.99-0.3           pkgconfig_2.0.3       
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
