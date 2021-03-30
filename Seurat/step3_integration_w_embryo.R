library(Seurat)
library(viridis)
library(ggplot2)
library(tidyr)
library(clustree)
library(dplyr)

## create embryo (Tyser) dataset based on dataset downloaded from http://www.human-gastrula.net/ on 22 July 2020
exp <- readRDS(file = "expression_values.rds") #obtained from "download gene expression matrix" button
expt <- as.data.frame(t(exp))
rownames(expt)[rownames(expt) == "T"] <- "TBXT" #hPSC data uses T instead of TBXT
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
load("WNTd.only.RData")
CH <- NormalizeData(CH, verbose = FALSE)
CH <- FindVariableFeatures(CH, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ty <- NormalizeData(ty, verbose = FALSE)
ty <- FindVariableFeatures(ty, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ty$exp <- "ty"
CH$exp <- "hPSC"
int.list = list(ty,CH)
aligned.anchors <- FindIntegrationAnchors(object.list = int.list, dims = 1:30)
tymeso <- IntegrateData(anchorset = aligned.anchors, dims = 1:30)
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
ID <- tymeso@meta.data$cluster_id
ID.replaced <- replace_na(ID, "WNTd")
tymeso <- AddMetaData(tymeso, metadata = ID.replaced, col.name = "cluster_id")


# Supplementary Fig. 5A
pdf()
DimPlot(tymeso, reduction = "umap", group.by = "exp", pt.size = 2.25, order = "ty", cols = c("#3c3c3c","#bbbbbb")) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20))
DimPlot(tymeso, reduction = "umap", group.by = "exp", pt.size = 2.25, order = "ty", cols = c("#3c3c3c","#bbbbbb")) +
theme(axis.title = element_text(size = 18), axis.text = element_text(size = 20)) + NoLegend()
dev.off()

# Supplementary Fig. 5Bi
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
table(ty.query$predicted.id)
pt1 <- as.character(ty@meta.data$orig.ident)
pt2 <- as.character(ty.query@meta.data$predicted.id)
transfer <- c(pt1,pt2)
tymeso <- AddMetaData(tymeso, metadata = transfer, col.name = "transfer")
pt3 <- as.character(ty@meta.data$cluster_id)
transfer2 <- c(pt3,pt2)
tymeso <- AddMetaData(tymeso, metadata = transfer2, col.name = "transfer2")
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

# Supplementary Fig. 5C
DefaultAssay(tymeso) <- "RNA"
Idents(tymeso) <- "exp"
embryo <- WhichCells(tymeso, idents = "ty")
vitro <- WhichCells(tymeso, idents = "hPSC")
CDX_features <- list(c("CDX1","CDX2","CDX4"))
tymeso <- AddModuleScore(tymeso, features = CDX_features, name = 'CDX_features')
pdf()
FeaturePlot(tymeso, cells = embryo, features = c("CDX_features1"), split.by = "exp", reduction = "umap", pt.size = 2.25, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = vitro, features = c("CDX_features1"), split.by = "exp", reduction = "umap", pt.size = 2.25, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = embryo, features = "ALDH1A2", split.by = "exp", reduction = "umap", pt.size = 2.25, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = vitro, features = "ALDH1A2", split.by = "exp", reduction = "umap", pt.size = 2.25, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = embryo, features = "CXCR4", split.by = "exp", reduction = "umap", pt.size = 2.25, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(tymeso, cells = vitro, features = "CXCR4", split.by = "exp", reduction = "umap", pt.size = 2.25, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()


# Supplementary Table 4A
embryo <- SubsetData(tymeso, subset.name = "orig.ident", accept.value = "tyser")
Idents(embryo) <- "cluster_id"
embryo.all.markers <- FindAllMarkers(embryo, logfc = 0.176)
write.table(embryo.all.markers, file="embryo.all.markers.txt", sep="\t")


# Supplementary Fig. 5D
embryo.sub <- SubsetData(embryo, subset.name = "cluster_id", accept.value = c("Nascent Mesoderm","Primitive Streak"))
DefaultAssay(embryo.sub) <- "RNA"
embryo.sub <- FindVariableFeatures(embryo.sub, selection.method = "vst", nfeatures = 2000)
embryo.sub <- ScaleData(embryo.sub, features = rownames(embryo.sub))
embryo.sub <- RunPCA(embryo.sub, features = VariableFeatures(embryo.sub), dims = 50)
embryo.sub <- JackStraw(embryo.sub, num.replicate = 100, dims = 50)
embryo.sub <- ScoreJackStraw(embryo.sub, dims = 1:50)
pdf()
JackStrawPlot(embryo.sub, dims = 1:50)
ElbowPlot(embryo.sub, ndims = 50)
#used a combination of the two plots to decide
dev.off()
embryo.sub <- FindNeighbors(embryo.sub, dims = 1:10)
embryo.sub <- RunUMAP(embryo.sub, dims = 1:10)
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

# Supplementary Fig. 5Di
Idents(embryo.sub) <- "RNA_snn_res.1.2"
embryo.sub <- RenameIdents(embryo.sub,'0'='6','1'='3','2'='4','3'='2','4'='7','5'='5','6'='1')
embryo.sub[["newclusters"]] <- Idents(embryo.sub)
Idents(embryo.sub) <- "newclusters"
my_levels <- c("1","2","3","4","5","6","7")
levels(embryo.sub) <- my_levels
pdf()
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "newclusters", order = c("7","6","5","4","3","2","1")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "newclusters", order = c("7","6","5","4","3","2","1")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
dev.off()

# Supplementary Fig. 5Dii
pdf()
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "cluster_id", cols = c("#00b4f0","#7cae00")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(embryo.sub, reduction = "umap", pt.size = 4, group.by = "cluster_id", cols = c("#00b4f0","#7cae00")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
dev.off()

# Supplementary Fig. 5E
pdf(width = 2, height = 3.5)
VlnPlot(embryo.sub, features = c("ALDH1A2","CXCR4"), log = TRUE, pt.size = 0.25, ncol = 1)
dev.off()

# Supplementary Table 4B
Idents(embryo.sub) <- "newclusters"
embryo.sub.markers <- FindAllMarkers(embryo.sub, logfc = 0.176)
write.table(embryo.sub.markers, file="embryo.sub.markers.txt", sep="\t")

# Supplementary Fig. 5F
CA_features <- list(c("ALDH1A2","CXCR4"))
embryo.sub <- AddModuleScore(embryo.sub, features = CA_features, name = 'CA_features')
CDX_features <- list(c("CDX1","CDX2","CDX4"))
embryo.sub <- AddModuleScore(embryo.sub, features = CDX_features, name = 'CDX_features')
pdf()
FeaturePlot(embryo.sub, features = c("CA_features1"), reduction = "umap", pt.size = 3, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(embryo.sub, features = c("CA_features1"), reduction = "umap", pt.size = 3, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank()) + NoLegend()
FeaturePlot(embryo.sub, features = c("CDX_features1"), reduction = "umap", pt.size = 3, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(embryo.sub, features = c("CDX_features1"), reduction = "umap", pt.size = 3, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank()) + NoLegend()
dev.off()


# Supplementary Fig. 5G
Idents(embryo.sub) <- "newclusters"
markers.all <- FindAllMarkers(embryo.sub, logfc = 0, return.thresh = 1)
write.table(markers.all, file="markers.all.txt", sep="\t")

cluster1 <- filter(markers.all, cluster == "1")
  genes.1 <- as.matrix(cluster1[,7])
  FC.1 <- as.matrix(cluster1[,2])
  genes.FC.1 <- matrix(c(genes.1,FC.1), ncol = 2)
  colnames(genes.FC.1) <- c("Genes","FC")
write.table(genes.FC.1, file = "genes.FC.1.csv", sep = ",", row.names = FALSE)

cluster2 <- filter(markers.all, cluster == "2")
  genes.2 <- as.matrix(cluster2[,7])
  FC.2 <- as.matrix(cluster2[,2])
  genes.FC.2 <- matrix(c(genes.2,FC.2), ncol = 2)
  colnames(genes.FC.2) <- c("Genes","FC")
write.table(genes.FC.2, file = "genes.FC.2.csv", sep = ",", row.names = FALSE)

cluster3 <- filter(markers.all, cluster == "3")
  genes.3 <- as.matrix(cluster3[,7])
  FC.3 <- as.matrix(cluster3[,2])
  genes.FC.3 <- matrix(c(genes.3,FC.3), ncol = 2)
  colnames(genes.FC.3) <- c("Genes","FC")
write.table(genes.FC.3, file = "genes.FC.3.csv", sep = ",", row.names = FALSE)

cluster4 <- filter(markers.all, cluster == "4")
  genes.4 <- as.matrix(cluster4[,7])
  FC.4 <- as.matrix(cluster4[,2])
  genes.FC.4 <- matrix(c(genes.4,FC.4), ncol = 2)
  colnames(genes.FC.4) <- c("Genes","FC")
write.table(genes.FC.4, file = "genes.FC.4.csv", sep = ",", row.names = FALSE)

cluster5 <- filter(markers.all, cluster == "5")
  genes.5 <- as.matrix(cluster5[,7])
  FC.5 <- as.matrix(cluster5[,2])
  genes.FC.5 <- matrix(c(genes.5,FC.5), ncol = 2)
  colnames(genes.FC.5) <- c("Genes","FC")
write.table(genes.FC.5, file = "genes.FC.5.csv", sep = ",", row.names = FALSE)

cluster6 <- filter(markers.all, cluster == "6")
  genes.6 <- as.matrix(cluster6[,7])
  FC.6 <- as.matrix(cluster6[,2])
  genes.FC.6 <- matrix(c(genes.6,FC.6), ncol = 2)
  colnames(genes.FC.6) <- c("Genes","FC")
write.table(genes.FC.6, file = "genes.FC.6.csv", sep = ",", row.names = FALSE)

cluster7 <- filter(markers.all, cluster == "7")
  genes.7 <- as.matrix(cluster7[,7])
  FC.7 <- as.matrix(cluster7[,2])
  genes.FC.7 <- matrix(c(genes.7,FC.7), ncol = 2)
  colnames(genes.FC.7) <- c("Genes","FC")
write.table(genes.FC.7, file = "genes.FC.7.csv", sep = ",", row.names = FALSE)

save.image(file = "tyWNTd-transfer-integrated.RData")
