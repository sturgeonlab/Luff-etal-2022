library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(viridis)
library(matrixStats)
library(outliers)
library(R.utils)
library(pheatmap)

library(devtools)
install_github("sturgeonlab/Luff-etal-2021/SingleR/versioncontrolscripts")
library(versioncontrolscripts)

# set up Zeng dataset
#obtain matrices from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135202
cs10 <-read.table("GSM3993420_CS10_Body_rawdata.txt", header = T, row.names=1, as.is=T) #10X
cs11 <-read.table("GSM3993421_CS11_CH_rawdata.txt", header = T, row.names=1, as.is=T) #10X
cs12 <-read.table("GSM3993423_CS12_CH_UMI_raw.txt", header = T, row.names=1, as.is=T) #dropseq
cs13X <-read.table("GSM3993422_CS13_DA_rawdata.txt", header = T, row.names=1, as.is=T) #10X
cs13D <-read.table("GSM3993424_CS13_DA_UMI_raw.txt", header = T, row.names=1, as.is=T) #dropseq
cs14 <-read.table("GSM3993425_CS14_DA_UMI_raw.txt", header = T, row.names=1, as.is=T) #dropseq
cs15 <-read.table("GSM3993426_CS15_DA_UMI_raw.txt", header = T, row.names=1, as.is=T) #dropseq
data <- CreateSeuratObject(cs10, project = "zeng", min.cells = 1)
data2 <- CreateSeuratObject(cs11, project = "zeng", min.cells = 1)
data3 <- CreateSeuratObject(cs12, project = "zeng", min.cells = 1)
data4 <- CreateSeuratObject(cs13X, project = "zeng", min.cells = 1)
data5 <- CreateSeuratObject(cs13D, project = "zeng", min.cells = 1)
data6 <- CreateSeuratObject(cs14, project = "zeng", min.cells = 1)
data7 <- CreateSeuratObject(cs15, project = "zeng", min.cells = 1)
data$sample <- "CS10"
data2$sample <- "CS11"
data3$sample <- "CS12"
data4$sample <- "CS13X"
data5$sample <- "CS13D"
data6$sample <- "CS14"
data7$sample <- "CS15"
zeng <- merge(data, y = c(data2,data3,data4,data5,data6,data7), add.cell.ids = c("CS10","CS11","CS12","CS13X","CS13D","CS14","CS15"), project = "zeng")
zeng <- NormalizeData(zeng, normalization.method = "LogNormalize", scale.factor = 10000)
zeng <- FindVariableFeatures(zeng, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(zeng)
zeng <- ScaleData(zeng, features = all.genes)

# identify HE and arterial endothelial cells for SingleR
pdf()
  plot(density(zeng@assays$RNA@data['CDH5',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['CXCR4',]))
  abline(v=0.3)
  plot(density(zeng@assays$RNA@data['GJA5',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['DLL4',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['PTPRC',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['SPN',]))
  abline(v=0.05)
  plot(density(zeng@assays$RNA@data['HEY2',]))
  abline(v=0.15)
  plot(density(zeng@assays$RNA@data['RUNX1',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['HOXA1',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['HOXA2',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['HOXA3',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['HOXA4',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['HOXA5',]))
  abline(v=0.15)
  plot(density(zeng@assays$RNA@data['HOXA6',]))
  abline(v=0.05)
  plot(density(zeng@assays$RNA@data['HOXA7',]))
  abline(v=0.15)
  plot(density(zeng@assays$RNA@data['HOXA9',]))
  abline(v=0.05)
  plot(density(zeng@assays$RNA@data['HOXA10',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['HOXA11',]))
  abline(v=0.02)
  plot(density(zeng@assays$RNA@data['HOXA13',]))
  abline(v=0.01)
  plot(density(zeng@assays$RNA@data['ITGA2B',]))
  abline(v=0.1)
  dev.off()
aec <- WhichCells(zeng, expression = CDH5 > 0.2 & CXCR4 > 0.3 & GJA5 > 0.2 & DLL4 > 0.1 & PTPRC < 0.1 & SPN < 0.05 & HEY2 > 0.15)
he1 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA1 > 0.2 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he2 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA2 > 0.1 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he3 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA3 > 0.2 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he4 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA4 > 0.1 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he5 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA5 > 0.15 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he6 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA6 > 0.05 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he7 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA7 > 0.15 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he9 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA9 > 0.05 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he10 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA10 > 0.2 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he11 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA11 > 0.02 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
#he13 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA13 > 0.01 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05) #no cells
allhe <- union(he1,c(he2,he3,he4,he5,he6,he7,he9,he10,he11))
Idents(zeng) <- "exp"
zeng <- SetIdent(zeng, cells = aec, value = "AEC")
zeng <- SetIdent(zeng, cells = allhe, value = "HE")
zeng[["newlabels"]] <- Idents(zeng)
zeng <- RunPCA(zeng, npcs = 50)
zeng <- JackStraw(zeng, num.replicate = 100, dims = 50)
zeng <- ScoreJackStraw(zeng, dims = 1:50)
pdf()
JackStrawPlot(zeng, dims = 1:50)
dev.off()
zeng <- FindNeighbors(zeng, reduction = "pca", dims = 1:50)
zeng <- FindClusters(zeng, resolution = 1)
zeng <- RunUMAP(zeng, reduction = "pca", dims = 1:50)

#add Zeng metadata
pheno <- read.table("phenotype.txt", header = TRUE, sep ="\t")
phenos <- pheno$Phenotype
zeng <- AddMetaData(zeng, metadata = phenos, col.name = "pheno")

# Supplementary Fig. 9A,B
pdf()
DimPlot(zeng, reduction = "umap", group.by = "sample", pt.size = 2.25) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20))
DimPlot(zeng, reduction = "umap", group.by = "sample", pt.size = 2.25) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20)) + NoLegend()
DimPlot(zeng, reduction = "umap", group.by = "pheno", pt.size = 2.25) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20))
DimPlot(zeng, reduction = "umap", group.by = "pheno", pt.size = 2.25) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20)) + NoLegend()
DimPlot(zeng, reduction = "umap", group.by = "newlabels", order = c("AEC","HE"),
cols = c("#c7c7c7","#7cae00","#00b4f0"), pt.size = 2.25) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20))
DimPlot(zeng, reduction = "umap", group.by = "newlabels", order = c("AEC","HE"),
cols = c("#c7c7c7","#7cae00","#00b4f0"), pt.size = 2.25) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20)) + NoLegend()
dev.off()
save(zeng, file="zeng.new.RData")


# set up reference dataset
name = "HE_Endo"
bulk <- read.table("bulkseq-input.txt", header=T, sep="\t", row.names = 1)
bulk = as.matrix(bulk)
types = list("IWP","RAi","RAd","CD184")
types = as.character(types)
main_types = list("HE","HE","HE","Endo")
main_types = as.character(main_types)
myref = list(data = bulk, main_types=main_types, types=types)
myref[["de.genes"]] = CreateVariableGeneSet(bulk,types,200)
myref[["de.genes.main"]] = CreateVariableGeneSet(bulk,main_types,300)

# assign categories to cells in cs dataset for "clusters" when plotted later
sub <- SubsetData(zeng, subset.name = "newlabels", accept.value = c("HE", "AEC"))
labels <- sub@meta.data$newlabels
stage <- sub@meta.data$sample
stage2 <- data.frame(stage)
stage2$stage <- as.character(stage2$stage)
stage2$stage[stage2$stage == "CS13D"] <- "CS13"
stage2$stage[stage2$stage == "CS13X"] <- "CS13"
stage2 <- stage2$stage

# create singlecell dataset for comparison to bulk RNAseq
data <- GetAssayData(sub, assay = "RNA", slot = "data")
data <- as.matrix(data)

# run main SingleR function & add names for cluster
obj <- SingleR.CreateObject(data, myref, clusters = NULL, variable.genes = "de", fine.tune = F)
obj[["names"]] <- stage2
obj[["names2"]] <- labels

# plot results by hierarchial clustering and "names"
pdf(width = 8, height = 2)
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names, order.by.clusters = TRUE)
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names2, order.by.clusters = TRUE)
dev.off()
#fine tune bypassed because it does not change scores with current type setup


# export relative correlation scores and names for each cell in cs dataset
scores = obj$SingleR.single.main$r
datascores <- as.matrix(scores)
mmax = rowMaxs(datascores)
mmin = rowMins(datascores)
datascoresrel = (datascores-mmin)/(mmax-mmin)
datascoresrel = datascoresrel^3
stage <- data.frame(obj$names)
scores <- data.frame(datascoresrel)
scores$stage <- stage$obj.names
type <- data.frame(obj$names2)
scores$type <- type$obj.names2


#successive sorting orders cells in aesthetically-pleasing manner
earlyCS <- c("CS10","CS11")
early <- filter(scores, stage %in% earlyCS)
early.HEC <- filter(early, type == "HE")
early.HEC <- early.HEC[order(early.HEC$IWP),]
early.HEC <- early.HEC[order(early.HEC$HER),]
early.HEC <- early.HEC[sort.list(early.HEC$CD184,decreasing = TRUE),]
early.HEC <- early.HEC[order(early.HEC$HED),]
early.AEC <- filter(early, type == "AEC")
early.AEC <- early.AEC[order(early.AEC$CD184),]
early.AEC <- early.AEC[order(early.AEC$IWP),]
early.AEC <- early.AEC[order(early.AEC$HER),]
early.AEC <- early.AEC[order(early.AEC$HED),]

lateCS <- c("CS12","CS13","CS14")
late <- filter(scores, stage %in% lateCS)
late.HEC <- filter(late, type == "HE")
late.HEC <- late.HEC[order(late.HEC$IWP),]
late.HEC <- late.HEC[order(late.HEC$HED),]
late.HEC <- late.HEC[order(late.HEC$CD184),]
late.HEC <- late.HEC[order(late.HEC$HER),]
late.AEC <- filter(late, type == "AEC")
late.AEC <- late.AEC[order(late.AEC$CD184),]
late.AEC <- late.AEC[order(late.AEC$IWP),]
late.AEC <- late.AEC[order(late.AEC$HED),]
late.AEC <- late.AEC[order(late.AEC$HER),]

early.HEC.cells <- rownames(early.HEC)
early.AEC.cells <- rownames(early.AEC)
late.HEC.cells <- rownames(late.HEC)
late.AEC.cells <- rownames(late.AEC)
order <- c(early.AEC.cells,early.HEC.cells,late.AEC.cells,late.HEC.cells)

# Fig. 3Bi
pdf(width = 8, height = 2)
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names, order.by.clusters = FALSE, cells_order = order)
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names2, order.by.clusters = FALSE, cells_order = order)
dev.off()

# Fig. 3Bii
write.table(datascoresrel, file="datascores.txt", sep="\t")
write.table(order, file="cellorder.txt", sep="\t")
write.table(type, file="type.txt", sep="\t")
write.table(stage, file="stage.txt", sep="\t")
# import into excel and use vlookup to match ordered barcodes with scores and type label
# see cellorder-matched.txt for final result


# Fig. 3Ci
genes <- c("CDH5","CD34","PECAM1","TEK","ENG","ESAM","KDR","THY1","SOX7","SOX17","CXCR4","GJA4",
"GJA5","HEY2","DLL4","NOTCH1","EFNB2","RUNX1")
zeng <- SetIdent(zeng, cells = early.AEC.cells, value = "early.AEC.cells")
zeng <- SetIdent(zeng, cells = early.HEC.cells, value = "early.HEC.cells")
zeng <- SetIdent(zeng, cells = late.AEC.cells, value = "late.AEC.cells")
zeng <- SetIdent(zeng, cells = late.HEC.cells, value = "late.HEC.cells")
zeng[["newestlabels"]] <- Idents(zeng)

#group cells based on scores
Idents(sub) <- "newlabels"
sub <- SetIdent(sub, cells = early.AEC.cells, value = "group1")
sub <- SetIdent(sub, cells = late.AEC.cells, value = "group4")
group2 <- early.HEC[1:46,]
group2.cells <- rownames(group2)
sub <- SetIdent(sub, cells = group2.cells, value = "group2")
group3 <- early.HEC[47:56,]
group3.cells <- rownames(group3)
sub <- SetIdent(sub, cells = group3.cells, value = "group3")
group5 <- late.HEC[1:41,]
group5.cells <- rownames(group5)
sub <- SetIdent(sub, cells = group5.cells, value = "group5")
group6 <- late.HEC[42:99,]
group6.cells <- rownames(group6)
sub <- SetIdent(sub, cells = group6.cells, value = "group6")
sub[["groups"]] <- Idents(sub)
Idents(sub) <- "groups"
my_levels <- c("group1","group2","group3","group4","group5","group6")
levels(sub) <- my_levels

# Fig. S9C
zeng <- SetIdent(zeng, cells = early.AEC.cells, value = "group1")
zeng <- SetIdent(zeng, cells = late.AEC.cells, value = "group4")
zeng <- SetIdent(zeng, cells = group2.cells, value = "group2")
zeng <- SetIdent(zeng, cells = group3.cells, value = "group3")
zeng <- SetIdent(zeng, cells = group5.cells, value = "group5")
zeng <- SetIdent(zeng, cells = group6.cells, value = "group6")
zeng[["groups"]] <- Idents(zeng)
pdf()
DimPlot(zeng, reduction = "umap", group.by = "groups", pt.size = 2.25, cols = c("#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7",
"#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7",
"#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#00b6eb","#f8766d","#c49a00","#a58aff","#53b400","#fb61d7"),
order = c("group5","group6","group4","group3","group2","group1")) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20))
DimPlot(zeng, reduction = "umap", group.by = "groups", pt.size = 2.25, cols = c("#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7",
"#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7",
"#c7c7c7","#c7c7c7","#c7c7c7","#c7c7c7","#00b6eb","#f8766d","#c49a00","#a58aff","#53b400","#fb61d7"),
order = c("group5","group6","group4","group3","group2","group1")) + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20)) + NoLegend()
dev.off()

# Supplementary Fig. 9Di
sub <- ScaleData(sub, features = rownames(sub))
tiff(width = 6500, height = 2000)
DoHeatmap(sub, draw.lines = FALSE, disp.min = 0, raster = FALSE, size = 3, features = genes, cells = order) + theme(axis.text.y = element_text(size = 6)) + scale_fill_viridis()
dev.off()

# Supplementary Fig. 9Dii
pdf(width = 2.5, height = 10)
VlnPlot(sub, features = c("CXCR4","GJA4","GJA5","HEY2","DLL4"), log = TRUE, pt.size = 0.25, ncol = 1)
dev.off()

# Supplementary Table 6A
Idents(sub) <- "groups"
markers <- FindAllMarkers(sub, logfc = 0.176)
write.table(markers, file = "group.markers.txt", sep = "\t")

# Supplementary Table 6B
avg <- AverageExpression(sub, assay = "RNA")
write.table(avg, file = "avg.txt", sep = "\t")

save(zeng,sub, file="zeng.sub.new.RData")
