library(Seurat)

## create embryo (Tyser) dataset based on dataset downloaded from http://www.human-gastrula.net/ on 22 July 2020
exp <- readRDS(file = "expression_values.rds") #obtained from "download gene expression matrix" button
expt <- as.data.frame(t(exp))
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
pdf()
ElbowPlot(ty)
dev.off()
ty <- FindNeighbors(ty, dims = 1:7)
ty <- FindClusters(ty, resolution = 0.5)
ty <- RunUMAP(ty, dims = 1:7)



## Regress cell cycle from Tyser dataset
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ty <- CellCycleScoring(ty, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ty <- RunPCA(ty, features = c(s.genes, g2m.genes))
pdf()
DimPlot(ty, reduction = "pca")
dev.off()
ty <- ScaleData(ty, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ty))
ty <- RunPCA(ty, features = VariableFeatures(ty), nfeatures.print = 10)
ty <- RunPCA(ty, features = c(s.genes, g2m.genes))
pdf()
DimPlot(ty, reduction = "pca")
dev.off()
save(ty, file = "ty-cc.RData")



## Regress cell cycle from WNTd-only dataset
load("WNTd.RData") # this dataset was generated from the original dataset prior to integration with WNTi/IWP2 and cell cycle regression
CH <- RunPCA(CH, features = VariableFeatures(CH), ndims.print = 1:10, nfeatures.print = 10)
CH <- CellCycleScoring(CH, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
CH <- RunPCA(CH, features = c(s.genes, g2m.genes))
CH <- ScaleData(CH, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CH))
CH <- RunPCA(CH, features = VariableFeatures(CH), nfeatures.print = 10)
CH <- RunPCA(CH, features = c(s.genes, g2m.genes))
save(CH, file = "temp-CH-cc.RData")

pdf()
DimPlot(CH, reduction = "pca")
DimPlot(CH, reduction = "umap", pt.size = 0.75, group.by = "Phase")
dev.off()

# Finishing clustering full regression
CH <- RunPCA(CH)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:6)
CH <- FindNeighbors(CH, dims = 1:6)
CH <- FindClusters(CH, resolution = 0.5)
pdf()
DimPlot(CH, reduction = "umap", pt.size = 0.75)
DimPlot(CH, reduction = "umap", pt.size = 0.75, group.by = "Phase")
DimPlot(CH, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID")
dev.off()
save(CH, file = "WNTd-cc.RData")



## Integrate hPSC WNTd with Tyser dataset
ty$exp <- "TY"
CH$exp <- "CHIRSB"
aligned.anchors <- FindIntegrationAnchors(object.list = list(ty, as))
double.aligned <- IntegrateData(anchorset = aligned.anchors, dims = 1:30)
DefaultAssay(double.aligned) <- "integrated"
double.aligned <- ScaleData(double.aligned, verbose = FALSE)
double.aligned <- RunPCA(double.aligned, npcs = 10, verbose = FALSE)
double.aligned <- FindNeighbors(double.aligned, reduction = "pca", dims = 1:8)
double.aligned <- FindClusters(double.aligned, resolution = 0.7)
double.aligned <- RunUMAP(double.aligned, reduction = "pca", dims = 1:8)

# Fig. S5A & S5Bi
pdf()
DimPlot(double.aligned, reduction = "umap", group.by = "exp", pt.size = 1, order = "ty")
DimPlot(double.aligned, reduction = "umap", group.by = "cluster_id", pt.size = 1, order = "ty")
dev.off()


# Fig. S5Bii
load("WNTd.RData")
write.table(CH@meta.data, file="CHmeta.txt", sep="\t")
write.table(double.aligned@meta.data, file="alignedmeta.txt", sep="\t")
# Using Excel, using Vlookup to obtain label for each CHIRSB cell in aligned metadata
# from the CH metadata column labeled "types"
newtypes <- read.table("tychir-newtypes.txt", header = TRUE, sep ="\t")
double.aligned2 = double.aligned@meta.data
double.aligned2["newtypes"] <- newtypes
trim <- subset(double.aligned2, select = c("newtypes"))
double.aligned <- AddMetaData(double.aligned, trim)
pdf()
DimPlot(double.aligned, reduction = "umap", pt.size = 2, group.by = "newtypes", order = c("Ectoderm","Endoderm","Pluripotent","Mesoderm","TY"), cols = c("#c7c7c7","#eb0006","#e5c321","#18b796","#6b8ddd")) + NoLegend()
dev.off()


# Fig. S5C
DefaultAssay(double.aligned) <- "integrated"
pdf()
FeaturePlot(double.aligned, features = c("KDR"), split.by = "exp", reduction = "umap", min.cutoff = 0, order = TRUE, cols = c("#e7e7e7","#eb0006")) + NoLegend()
FeaturePlot(double.aligned, features = c("CXCR4"), split.by = "exp", reduction = "umap", min.cutoff = 0, order = TRUE, cols = c("#e7e7e7","#eb0006")) + NoLegend()
FeaturePlot(double.aligned, features = c("ALDH1A2"), split.by = "exp", reduction = "umap", min.cutoff = 0, order = TRUE, cols = c("#e7e7e7","#eb0006")) + NoLegend()
FeaturePlot(double.aligned, features = c("CYP26A1"), split.by = "exp", reduction = "umap", min.cutoff = 0, order = TRUE, cols = c("#e7e7e7","#eb0006")) + NoLegend()
DimPlot(double.aligned, reduction = "umap", group.by = "integrated_snn_res.0.7", pt.size = 1)
dev.off()


# Fig. S5D
embryo <- OldWhichCells(double.aligned, subset.name = "exp", accept.value = "TY")
Idents(double.aligned) <- "integrated_snn_res.0.7"
A <- OldWhichCells(double.aligned, subset.name = "integrated_snn_res.0.7", accept.value = "0")
H <- OldWhichCells(double.aligned, subset.name = "integrated_snn_res.0.7", accept.value = "7")
I <- OldWhichCells(double.aligned, subset.name = "integrated_snn_res.0.7", accept.value = "8")
J <- OldWhichCells(double.aligned, subset.name = "integrated_snn_res.0.7", accept.value = "9")
L <- OldWhichCells(double.aligned, subset.name = "integrated_snn_res.0.7", accept.value = "11")
embryo.A <- Reduce(intersect, list(embryo,A))
embryo.H <- Reduce(intersect, list(embryo,H))
embryo.I <- Reduce(intersect, list(embryo,I))
embryo.J <- Reduce(intersect, list(embryo,J))
embryo.L <- Reduce(intersect, list(embryo,L))
double.aligned <- SetIdent(double.aligned, cells = embryo.I, value = "embryo.I")
double.aligned <- SetIdent(double.aligned, cells = embryo.A, value = "embryo.A")
double.aligned <- SetIdent(double.aligned, cells = embryo.H, value = "embryo.H")
double.aligned <- SetIdent(double.aligned, cells = embryo.J, value = "embryo.J")
double.aligned <- SetIdent(double.aligned, cells = embryo.L, value = "embryo.L")
DefaultAssay(double.aligned) <- "RNA"
my_levels <- c("embryo.A","embryo.H","embryo.I","embryo.J","embryo.L","0","1","2","3","4","5","6","7","8","9","10","11","12","13")
levels(double.aligned) <- my_levels
pdf()
VlnPlot(double.aligned, features = c("ALDH1A2","CXCR4","CYP26A1","TEK"), log = TRUE, cols = c("#109a54","#e5c321","#4ab4d5","#ed8f08","#eb0006","#fb6d92","#6b8ddd"), pt.size = 1, idents = c("embryo.A","embryo.H","embryo.I","embryo.J","embryo.L"))
dev.off()


# Fig. S5F
I.markers <- FindMarkers(double.aligned, logfc.threshold = 0.176, ident.1 = embryo.I, ident.2 = embryo.A)
write.table(I.markers, file="I-markers.txt", sep="\t")
A.markers <- FindMarkers(double.aligned, logfc.threshold = 0.176, ident.1 = embryo.A, ident.2 = embryo.I)
write.table(A.markers, file="A-markers.txt", sep="\t")


# Fig. S5E
load("WNTd.mesoderm.RData")
write.table(meso@meta.data, file="mesometa.txt", sep="\t")
write.table(double.aligned@meta.data, file="alignedmeta.txt", sep="\t")
# Using Excel, using Vlookup to obtain label for each mesoderm cell in aligned metadata
# from the mesoderm metadata column labeled "RNA_snn_res.0.5"
Fig1C.clusters <- read.table("Fig1C-clusters.txt", header = TRUE, sep ="\t") #overlay
double.aligned2 = double.aligned@meta.data
double.aligned2["Fig1C.clusters"] <- Fig1C.clusters
trim <- subset(double.aligned2, select = c("Fig1C.clusters"))
double.aligned <- AddMetaData(double.aligned, trim)
pdf()
DimPlot(double.aligned, reduction = "umap", pt.size = 2, group.by = "Fig1C.clusters")
dev.off()
