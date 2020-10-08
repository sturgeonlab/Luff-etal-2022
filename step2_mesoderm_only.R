library(Seurat)
library(viridis)

## Setting up WNTd-only dataset
load("CHIR-IWP-before-integration.RData")
Idents(agg) <- "Cell.ID"
CH <- SubsetData(agg, subset.name = "Cell.ID", accept.value = "CHIRSB")
CH <- ScaleData(CH)
CH <- RunPCA(CH, npcs = 30)
pdf()
ElbowPlot(CH)
dev.off()
CH <- FindNeighbors(CH, reduction = "pca", dims = 1:5)
CH <- FindClusters(CH, resolution = 1)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:6)

# Setting identities so that every cell is identified
# come back to here following cell cycle regression
endo2 <- WhichCells(object = CH, expression = APOA1 > 0.1)
ecto4 <- WhichCells(object = CH, expression = DLGAP5 > 0.1)
ecto3 <- WhichCells(object = CH, expression = EPCAM > 0.1)
ecto2 <- WhichCells(object = CH, expression = NPY > 0.1)
ecto <- WhichCells(object = CH, expression = TFAP2A > 0.1)
endo <- WhichCells(object = CH, expression = FOXA2 > 0.1)
pluri <- WhichCells(object = CH, expression = SOX2 > 0.1)
meso7 <- WhichCells(object = CH, expression = TBX6 > 0.1) #paraxial
meso4 <- WhichCells(object = CH, expression = PRRX1 > 0.1) #muscle
meso3 <- WhichCells(object = CH, expression = TMEM88 > 0.1) #cardio
meso2 <- WhichCells(object = CH, expression = MESP1 > 0.1) #cardio
meso6 <- WhichCells(object = CH, expression = MEST > 0.1) #generic meso
meso <- WhichCells(object = CH, expression = KDR > 0.1)
CH <- SetIdent(object = CH, cells = ecto4, value = 'Ectoderm')
CH <- SetIdent(object = CH, cells = ecto3, value = 'Ectoderm')
CH <- SetIdent(object = CH, cells = ecto2, value = 'Ectoderm')
CH <- SetIdent(object = CH, cells = endo2, value = 'Endoderm')
CH <- SetIdent(object = CH, cells = meso7, value = 'Mesoderm')
CH <- SetIdent(object = CH, cells = meso4, value = 'Mesoderm')
CH <- SetIdent(object = CH, cells = meso3, value = 'Mesoderm')
CH <- SetIdent(object = CH, cells = endo, value = 'Endoderm')
CH <- SetIdent(object = CH, cells = ecto, value = 'Ectoderm')
CH <- SetIdent(object = CH, cells = pluri, value = 'Pluripotent')
CH <- SetIdent(object = CH, cells = meso6, value = 'Mesoderm')
CH <- SetIdent(object = CH, cells = meso2, value = 'Mesoderm')
CH <- SetIdent(object = CH, cells = meso, value = 'Mesoderm')
CH[["types"]] <- Idents(CH)
save(CH, file = "WNTd.RData")

# Fig. S2Fi,ii
Idents(CH) <- "types"
my_levels <- c("Mesoderm", "Endoderm", "Ectoderm", "Pluripotent")
levels(CH) <- my_levels
pdf()
FeaturePlot(CH, features = c("KDR"), pt.size = 1.25, reduction = "umap", min.cutoff = 0, order = TRUE) + scale_color_viridis()
DimPlot(CH, reduction = "umap", pt.size = 1.5, order = c("Ectoderm","Endoderm","Pluripotent","Mesoderm"), cols = c("#eb0006","#e5c321","#fb6d92","#6b8ddd"))
dev.off()

# Fig. S2Fiii
features = c("SOX2","POU5F1","NANOG","TFAP2A","DLX5","KRT7",
"FOXA2","SOX17","HNF1B","MEST","KDR","MESP1")
pdf(width = 8, height = 3.5)
DotPlot(CH, features = features, cols = c("#f0f0f0", "#f02207"), col.min = 0) + RotatedAxis()
dev.off()



## Setting up mesoderm-only WNTd dataset
# input can be non-regressed or cell cycle regressed (it doesn't affect this downstream analysis)
meso <- SubsetData(CH, subset.name = "KDR", low.threshold = 0.1)
meso <- NormalizeData(meso)
meso <- FindVariableFeatures(meso, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(meso)
meso <- ScaleData(meso, features = all.genes)
meso <- RunPCA(meso, features = VariableFeatures(object = meso))
pdf()
ElbowPlot(meso)
dev.off()
meso <- FindNeighbors(meso, dims = 1:5)
meso <- FindClusters(meso, resolution = 0.6)
meso <- RunUMAP(meso, dims = 1:5)
save(meso, file = "WNTd.mesoderm.RData")

# Fig. 1C
pdf()
DimPlot(meso, reduction = "umap", pt.size = 2, group.by = "RNA_snn_res.0.6") + NoLegend()
FeaturePlot(meso, features = c("ALDH1A2"), pt.size = 2, reduction = "umap", min.cutoff = 0, order = TRUE, cols = c("#e7e7e7","#eb0006"))
FeaturePlot(meso, features = c("CXCR4"), pt.size = 2, reduction = "umap", min.cutoff = 0, order = TRUE, cols = c("#e7e7e7","#eb0006"))
dev.off()

# Identifying markers of ALDH1A2+ mesodermal cells
aldh.pos <- WhichCells(meso, expression = ALDH1A2 > 0.1)
aldh.neg <- WhichCells(meso, expression = ALDH1A2 < 0.1)
meso <- SetIdent(meso, cells = aldh.neg, value = "aldh.neg")
meso <- SetIdent(meso, cells = aldh.pos, value = "aldh.pos")
aldh.markers <- FindAllMarkers(meso, logfc.threshold = 0.176) #0.176 logfc = 1.5 linear fc
write.table(aldh.markers, file="aldhmesomarkers.txt", sep="\t")










## cell cycle regressions not used in publication
load("CHIR-IWP-before-integration.RData")
Idents(agg) <- "Cell.ID"
CH <- SubsetData(agg, subset.name = "Cell.ID", accept.value = "CHIRSB")
CH <- ScaleData(CH)
CH <- RunPCA(CH, npcs = 30)
CH <- FindNeighbors(CH, reduction = "pca", dims = 1:5)
CH <- FindClusters(CH, resolution = 1)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:6)
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
CH <- RunPCA(CH, features = VariableFeatures(CH), ndims.print = 1:10, nfeatures.print = 10)
CH <- CellCycleScoring(CH, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
CH <- RunPCA(CH, features = c(s.genes, g2m.genes))
pdf()
DimPlot(CH, reduction = "pca")
dev.off()

## partial
CH$CC.Difference <- CH$S.Score - CH$G2M.Score
CH <- ScaleData(CH, vars.to.regress = "CC.Difference", features = rownames(CH))
CH <- RunPCA(CH, features = VariableFeatures(CH), nfeatures.print = 10)
CH <- RunPCA(CH, features = c(s.genes, g2m.genes))
pdf()
DimPlot(CH, reduction = "pca")
dev.off()
CH <- RunPCA(CH)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:6)
CH <- FindNeighbors(CH, dims = 1:6)
CH <- FindClusters(CH, resolution = 0.5)
pdf()
DimPlot(CH, reduction = "umap", pt.size = 0.75)
DimPlot(CH, reduction = "umap", pt.size = 0.75, group.by = "Phase")
dev.off()
# go up to set identities
save(CH, file = "WNTd-partialcc.RData")
 
## full
CH <- ScaleData(CH, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CH))
CH <- RunPCA(CH, features = VariableFeatures(CH), nfeatures.print = 10)
CH <- RunPCA(CH, features = c(s.genes, g2m.genes))
pdf()
DimPlot(CH, reduction = "pca")
dev.off()
CH <- RunPCA(CH)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:6)
CH <- FindNeighbors(CH, dims = 1:6)
CH <- FindClusters(CH, resolution = 0.5)
pdf()
DimPlot(CH, reduction = "umap", pt.size = 0.75)
DimPlot(CH, reduction = "umap", pt.size = 0.75, group.by = "Phase")
DimPlot(CH, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID")
dev.off()
# go up to set identities
save(CH, file = "WNTd-fullcc.RData")
