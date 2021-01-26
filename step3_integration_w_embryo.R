library(Seurat)
library(viridis)
library(ggplot2)

## create embryo (Tyser) dataset based on dataset downloaded from http://www.human-gastrula.net/ on 22 July 2020
exp <- readRDS(file = "expression_values.rds") #obtained from "download gene expression matrix" button
expt <- as.data.frame(t(exp))
rownames(expt)[rownames(expt) == "TBXT"] <- "T" #hPSC data uses T instead of TBXT
annot <- readRDS(file = "annot_umap.rds") #obtained from "download annotation and UMAP" button
ty <- CreateSeuratObject(expt, project = "tyser", assay = "RNA")
id <- read.table("id.txt", header = TRUE, sep ="\t")
ids <- id$cluster_id
ty <- AddMetaData(ty, metadata = ids, col.name = "cluster_id")


# process data
ty <- NormalizeData(ty)
ty <- FindVariableFeatures(ty, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ty)
ty <- ScaleData(ty, features = all.genes)
ty <- RunPCA(ty, features = VariableFeatures(object = ty))
ty <- JackStraw(ty, num.replicate = 100, dims = 50)
ty <- ScoreJackStraw(ty, dims = 1:50)
pdf()
JackStrawPlot(ty, dims = 1:50)
dev.off()
ty <- FindNeighbors(ty, dims = 1:27)
ty <- FindClusters(ty, resolution = 0.5)
ty <- RunUMAP(ty, dims = 1:27)
save(ty, file = "ty.only.RData")


# Integrate hPSC WNTd with Tyser dataset
load("WNTd.RData")

ty$exp <- "ty"
CH$exp <- "hPSC"

int.list = list(ty,CH)
all_genes <- Reduce(intersect, lapply(int.list, rownames))
aligned.anchors <- FindIntegrationAnchors(object.list = int.list)
tymeso <- IntegrateData(anchorset = aligned.anchors, dims = 1:50, features.to.integrate = all_genes)
DefaultAssay(tymeso) <- "integrated"
tymeso <- ScaleData(tymeso, verbose = FALSE)
tymeso <- RunPCA(tymeso, npcs = 50)
tymeso <- JackStraw(tymeso, num.replicate = 100, dims = 50)
tymeso <- ScoreJackStraw(tymeso, dims = 1:50)
pdf()
JackStrawPlot(tymeso, dims = 1:50)
dev.off()
tymeso <- FindNeighbors(tymeso, reduction = "pca", dims = 1:50)
tymeso <- FindClusters(tymeso, resolution = 1)
tymeso <- RunUMAP(tymeso, reduction = "pca", dims = 1:50)


# replacing NAs so plotting order can be manipulated
library(tidyr)
ID <- tymeso@meta.data$cluster_id
ID.replaced <- replace_na(ID, "WNTd")
tymeso <- AddMetaData(tymeso, metadata = ID.replaced, col.name = "cluster_id")

# refine clustering resolution
library('clustree')
tymeso <- FindClusters(tymeso, resolution = 0.0)
  tymeso <- FindClusters(tymeso, resolution = 0.1)
  tymeso <- FindClusters(tymeso, resolution = 0.2)
  tymeso <- FindClusters(tymeso, resolution = 0.3)
  tymeso <- FindClusters(tymeso, resolution = 0.4)
  tymeso <- FindClusters(tymeso, resolution = 0.5)
  tymeso <- FindClusters(tymeso, resolution = 0.6)
  tymeso <- FindClusters(tymeso, resolution = 0.7)
  tymeso <- FindClusters(tymeso, resolution = 0.8)
  tymeso <- FindClusters(tymeso, resolution = 0.9)
  tymeso <- FindClusters(tymeso, resolution = 1)
  tymeso <- FindClusters(tymeso, resolution = 1.1)
  tymeso <- FindClusters(tymeso, resolution = 1.2)
  tymeso <- FindClusters(tymeso, resolution = 1.3)
  tymeso <- FindClusters(tymeso, resolution = 1.4)
pdf(width = 8, height = 14)
clustree(tymeso, prefix = "integrated_snn_res.", layout = "sugiyama", use_core_edges = FALSE)
clustree(tymeso, prefix = "integrated_snn_res.", layout = "sugiyama", use_core_edges = FALSE, node_colour = "sc3_stability") + scale_color_viridis(option="plasma")
dev.off()

save(tymeso, file = "ty-WNTd-integrated-v16.RData")



# Fig. S5A
pdf()
DimPlot(tymeso, reduction = "umap", group.by = "exp", pt.size = 2.25, order = "ty", cols = c("#3c3c3c","#bbbbbb")) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20))
DimPlot(tymeso, reduction = "umap", group.by = "exp", pt.size = 2.25, order = "ty", cols = c("#3c3c3c","#bbbbbb")) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20)) + NoLegend()
dev.off()

# Fig. S5Bi
pdf()
DimPlot(tymeso, reduction = "umap", pt.size = 2.25, group.by = "cluster_id", cols = c("#c7c7c7","#02c08b","#c77cff","#ff64b0","#7cae00","#01bfc4","#619cff","#ff6464",
"#b79f00","#de8c00","#f564e3","#00b4f0"), order = c("Nascent Mesoderm","Axial Mesoderm","Emergent Mesoderm","Endoderm","Erythrocytes","Hemogenic Endothelial Progenitors",
"Epiblast","Primitive Streak","YS Mesoderm","Ectoderm","Advanced Mesoderm")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
DimPlot(tymeso, reduction = "umap", pt.size = 2.25, group.by = "cluster_id", cols = c("#c7c7c7","#02c08b","#c77cff","#ff64b0","#7cae00","#01bfc4","#619cff","#ff6464",
"#b79f00","#de8c00","#f564e3","#00b4f0"), order = c("Nascent Mesoderm","Axial Mesoderm","Emergent Mesoderm","Endoderm","Erythrocytes","Hemogenic Endothelial Progenitors",
"Epiblast","Primitive Streak","YS Mesoderm","Ectoderm","Advanced Mesoderm")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()

# Fig. S5Bii
DefaultAssay(tymeso) <- "integrated"
ty.list <- SplitObject(tymeso, split.by = "exp")
ty.query <- ty.list[["hPSC"]]
ty.ref <- ty.list[["ty"]]
ty.anchors <- FindTransferAnchors(reference = ty.ref, query = ty.query, dims = 1:30)
predictions <- TransferData(anchorset = ty.anchors, refdata = ty.ref$cluster_id, dims = 1:30)
ty.query <- AddMetaData(ty.query, metadata = predictions)
ty.query$prediction.match <- ty.query$predicted.id == ty.query$cluster_id
table(ty.query$prediction.match)
table(ty.query$predicted.id)
pt1 <- ty@meta.data$orig.ident
pt2 <- ty.query$predicted.id
transfer <- c(pt1,pt2)
transfer <- data.frame(transfer)
transfers <- transfer$transfer
tymeso <- AddMetaData(tymeso, metadata = transfers, col.name = "transfer")
pdf()
DimPlot(tymeso, reduction = "umap", pt.size = 2.25, group.by = "transfer", cols = c("#c7c7c7","#02c08b","#c77cff","#ff64b0","#7cae00","#01bfc4","#619cff","#ff6464",
"#b79f00","#de8c00","#f564e3","#00b4f0"), order = c("Nascent Mesoderm","Axial Mesoderm","Emergent Mesoderm","Endoderm","Erythrocytes","Hemogenic Endothelial Progenitors",
"Epiblast","Primitive Streak","YS Mesoderm","Ectoderm","Advanced Mesoderm")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
DimPlot(tymeso, reduction = "umap", pt.size = 2.25, group.by = "transfer", cols = c("#c7c7c7","#02c08b","#c77cff","#ff64b0","#7cae00","#01bfc4","#619cff","#ff6464",
"#b79f00","#de8c00","#f564e3","#00b4f0"), order = c("Nascent Mesoderm","Axial Mesoderm","Emergent Mesoderm","Endoderm","Erythrocytes","Hemogenic Endothelial Progenitors",
"Epiblast","Primitive Streak","YS Mesoderm","Ectoderm","Advanced Mesoderm")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()


# Fig. S5C
DefaultAssay(tymeso) <- "RNA"
Idents(tymeso) <- "exp"
embryo <- WhichCells(tymeso, idents = "ty")
vitro <- WhichCells(tymeso, idents = "hPSC")
pdf()
FeaturePlot(tymeso, cells = embryo, features = "ALDH1A2", reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = vitro, features = c("ALDH1A2"), reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = embryo, features = "CDX4", reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = vitro, features = c("CDX4"), reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = embryo, features = "CYP26A1", reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = vitro, features = c("CYP26A1"), reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = embryo, features = "GYPA", reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = vitro, features = c("GYPA"), reduction = "umap", order = TRUE, pt.size = 2.25) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()


### Embryo dataset
embryo <- SubsetData(tymeso, subset.name = "orig.ident", accept.value = "tyser")

# Table 4A
Idents(embryo) <- "cluster_id"
avg <- AverageExpression(embryo, assay = "RNA")
write.table(avg, file = "embryo-avg.txt", sep = "\t")

# Table 4B
markers <- FindAllMarkers(embryo, logfc = 0.176)
write.table(markers, file = "markers.txt", sep = "\t")

# Table 4C
pdf()
plot(density(embryo@assays$RNA@data['CDX1',]))
  abline(v=0.25)
  plot(density(embryo@assays$RNA@data['CDX2',]))
  abline(v=0.25)
  plot(density(embryo@assays$RNA@data['CDX4',]))
  abline(v=0.2)
  plot(density(embryo@assays$RNA@data['T',]))
  abline(v=0.4)
  plot(density(embryo@assays$RNA@data['MIXL1',]))
  abline(v=0.25)
  plot(density(embryo@assays$RNA@data['MESP1',]))
  abline(v=0.3)
  plot(density(embryo@assays$RNA@data['PDGFRA',]))
  abline(v=0.4)
  plot(density(embryo@assays$RNA@data['HAND1',]))
  abline(v=0.4)
  plot(density(embryo@assays$RNA@data['KDR',]))
  abline(v=0.3)
  plot(density(embryo@assays$RNA@data['TEK',]))
  abline(v=0.3)
  plot(density(embryo@assays$RNA@data['FLT1',]))
  abline(v=0.25)
  plot(density(embryo@assays$RNA@data['CD34',]))
  abline(v=0.2)
  plot(density(embryo@assays$RNA@data['CDH5',]))
  abline(v=0.15)
  plot(density(embryo@assays$RNA@data['GYPA',]))
  abline(v=0.25)

  plot(density(embryo@assays$RNA@data['GYPB',]))
  abline(v=0.3)
  plot(density(embryo@assays$RNA@data['GYPC',]))
  abline(v=0.45)
  plot(density(embryo@assays$RNA@data['GYPE',]))
  abline(v=0.25)
  plot(density(embryo@assays$RNA@data['NOTO',]))
  abline(v=0.1)
  plot(density(embryo@assays$RNA@data['CHRD',]))
  abline(v=0.15)
  plot(density(embryo@assays$RNA@data['FOXA2',]))
  abline(v=0.25)
  plot(density(embryo@assays$RNA@data['SOX17',]))
  abline(v=0.2)
  dev.off()


### Embryo subset
embryo.sub <- SubsetData(tymeso, subset.name = "cluster_id", accept.value = c("Nascent Mesoderm","Primitive Streak"))
DefaultAssay(embryo.sub) <- "RNA"
embryo.sub <- FindVariableFeatures(embryo.sub, selection.method = "vst", nfeatures = 2000)
embryo.sub <- ScaleData(embryo.sub, features = rownames(embryo.sub))
embryo.sub <- RunPCA(embryo.sub, features = VariableFeatures(embryo.sub), dims = 50)
embryo.sub <- JackStraw(embryo.sub, num.replicate = 100, dims = 50)
embryo.sub <- ScoreJackStraw(embryo.sub, dims = 1:50)
pdf()
JackStrawPlot(embryo.sub, dims = 1:50)
ElbowPlot(embryo.sub, ndims = 50)
dev.off()
embryo.sub <- FindNeighbors(embryo.sub, dims = 1:10)
embryo.sub <- RunUMAP(embryo.sub, dims = 1:10)
embryo.sub <- FindClusters(embryo.sub, resolution = 1)
library('clustree')
embryo.sub <- FindClusters(embryo.sub, resolution = 0.0)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.1)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.2)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.3)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.4)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.5)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.6)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.7)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.8)
  embryo.sub <- FindClusters(embryo.sub, resolution = 0.9)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.1)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.2)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.3)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.4)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.5)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.6)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.7)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.8)
  embryo.sub <- FindClusters(embryo.sub, resolution = 1.9)
  embryo.sub <- FindClusters(embryo.sub, resolution = 2)
  pdf(width = 8, height = 14)
  clustree(embryo.sub, prefix = "RNA_snn_res.", layout = "sugiyama", use_core_edges = FALSE)
  clustree(embryo.sub, prefix = "RNA_snn_res.", layout = "sugiyama", use_core_edges = FALSE, node_colour = "sc3_stability") + scale_color_viridis(option="plasma")
  dev.off()


# Fig. S5Di
pdf()
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "cluster_id", cols = c("#6f9414","#dfc51d")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "cluster_id", cols = c("#6f9414","#dfc51d")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
dev.off()

# Fig. S5Dii
Idents(embryo.sub) <- "RNA_snn_res.1"
embryo.sub <- RenameIdents(embryo.sub,'0'='5','1'='3','2'='2','3'='4','4'='6','5'='1')
embryo.sub[["newclusters"]] <- Idents(embryo.sub)
my_levels <- c("1","2","3","4","5","6")
levels(embryo.sub) <- my_levels
pdf()
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "newclusters", order = c("6","5","4","3","2","1")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "newclusters", order = c("6","5","4","3","2","1")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
dev.off()

# Fig. S5Diii
pdf(width = 6, height = 2.3)
VlnPlot(embryo.sub, features = c("ALDH1A2","CXCR4"), log = TRUE, pt.size = 0.25, ncol = 3)
dev.off()


# Fig. S5E
Idents(embryo.sub) <- "newclusters"
markers.all <- FindAllMarkers(embryo.sub, logfc = 0, return.thresh = 1) # Table 4D

library(clusterProfiler)
library(dplyr)

pos.up <- read.table("newpos05.txt", header = TRUE, sep ="\t") # Table 4E
neg.up <- read.table("newneg05.txt", header = TRUE, sep ="\t") # Table 4E

cluster1 <- filter(markers.all, cluster == "1")
  genes.1 <- as.matrix(cluster1[,7])
  FC.1 <- as.matrix(cluster1[,2])
  genes.FC.1 <- matrix(c(genes.1,FC.1), ncol = 2)
  c1.geneList = as.numeric(genes.FC.1[,2])
  names(c1.geneList) = as.character(genes.FC.1[,1])
  c1.geneList.sorted = sort(c1.geneList, decreasing = TRUE)
cluster2 <- filter(markers.all, cluster == "2")
  genes.2 <- as.matrix(cluster2[,7])
  FC.2 <- as.matrix(cluster2[,2])
  genes.FC.2 <- matrix(c(genes.2,FC.2), ncol = 2)
  c2.geneList = as.numeric(genes.FC.2[,2])
  names(c2.geneList) = as.character(genes.FC.2[,1])
  c2.geneList.sorted = sort(c2.geneList, decreasing = TRUE)
cluster3 <- filter(markers.all, cluster == "3")
  genes.3 <- as.matrix(cluster3[,7])
  FC.3 <- as.matrix(cluster3[,2])
  genes.FC.3 <- matrix(c(genes.3,FC.3), ncol = 2)
  c3.geneList = as.numeric(genes.FC.3[,2])
  names(c3.geneList) = as.character(genes.FC.3[,1])
  c3.geneList.sorted = sort(c3.geneList, decreasing = TRUE)
cluster4 <- filter(markers.all, cluster == "4")
  genes.4 <- as.matrix(cluster4[,7])
  FC.4 <- as.matrix(cluster4[,2])
  genes.FC.4 <- matrix(c(genes.4,FC.4), ncol = 2)
  c4.geneList = as.numeric(genes.FC.4[,2])
  names(c4.geneList) = as.character(genes.FC.4[,1])
  c4.geneList.sorted = sort(c4.geneList, decreasing = TRUE)
cluster5 <- filter(markers.all, cluster == "5")
  genes.5 <- as.matrix(cluster5[,7])
  FC.5 <- as.matrix(cluster5[,2])
  genes.FC.5 <- matrix(c(genes.5,FC.5), ncol = 2)
  c5.geneList = as.numeric(genes.FC.5[,2])
  names(c5.geneList) = as.character(genes.FC.5[,1])
  c5.geneList.sorted = sort(c5.geneList, decreasing = TRUE)
cluster6 <- filter(markers.all, cluster == "6")
  genes.6 <- as.matrix(cluster6[,7])
  FC.6 <- as.matrix(cluster6[,2])
  genes.FC.6 <- matrix(c(genes.6,FC.6), ncol = 2)
  c6.geneList = as.numeric(genes.FC.6[,2])
  names(c6.geneList) = as.character(genes.FC.6[,1])
  c6.geneList.sorted = sort(c6.geneList, decreasing = TRUE)
c1.posup <- GSEA(c1.geneList.sorted, TERM2GENE = pos.up, verbose = FALSE, pvalueCutoff = 1)
  c2.posup <- GSEA(c2.geneList.sorted, TERM2GENE = pos.up, verbose = FALSE, pvalueCutoff = 1)
  c3.posup <- GSEA(c3.geneList.sorted, TERM2GENE = pos.up, verbose = FALSE, pvalueCutoff = 1)
  c4.posup <- GSEA(c4.geneList.sorted, TERM2GENE = pos.up, verbose = FALSE, pvalueCutoff = 1)
  c5.posup <- GSEA(c5.geneList.sorted, TERM2GENE = pos.up, verbose = FALSE, pvalueCutoff = 1)
  c6.posup <- GSEA(c6.geneList.sorted, TERM2GENE = pos.up, verbose = FALSE, pvalueCutoff = 1)
  c1.negup <- GSEA(c1.geneList.sorted, TERM2GENE = neg.up, verbose = FALSE, pvalueCutoff = 1)
  c2.negup <- GSEA(c2.geneList.sorted, TERM2GENE = neg.up, verbose = FALSE, pvalueCutoff = 1)
  c3.negup <- GSEA(c3.geneList.sorted, TERM2GENE = neg.up, verbose = FALSE, pvalueCutoff = 1)
  c4.negup <- GSEA(c4.geneList.sorted, TERM2GENE = neg.up, verbose = FALSE, pvalueCutoff = 1)
  c5.negup <- GSEA(c5.geneList.sorted, TERM2GENE = neg.up, verbose = FALSE, pvalueCutoff = 1)
  c6.negup <- GSEA(c6.geneList.sorted, TERM2GENE = neg.up, verbose = FALSE, pvalueCutoff = 1)
write.table(c1.posup, file = "c1.posup.txt", sep = "\t")
  write.table(c2.posup, file = "c2.posup.txt", sep = "\t")
  write.table(c3.posup, file = "c3.posup.txt", sep = "\t")
  write.table(c4.posup, file = "c4.posup.txt", sep = "\t")
  write.table(c5.posup, file = "c5.posup.txt", sep = "\t")
  write.table(c6.posup, file = "c6.posup.txt", sep = "\t")
  write.table(c1.negup, file = "c1.negup.txt", sep = "\t")
  write.table(c2.negup, file = "c2.negup.txt", sep = "\t")
  write.table(c3.negup, file = "c3.negup.txt", sep = "\t")
  write.table(c4.negup, file = "c4.negup.txt", sep = "\t")
  write.table(c5.negup, file = "c5.negup.txt", sep = "\t")
  write.table(c6.negup, file = "c6.negup.txt", sep = "\t")
 
