## Setting up WNTd dataset
load("agg.RData")
Idents(agg) <- "Cell.ID"
CH <- SubsetData(agg, subset.name = "Cell.ID", accept.value = "CHIRSB")
CH <- ScaleData(CH)
CH <- RunPCA(CH, npcs = 30)
pdf()
ElbowPlot(meso)
dev.off()
CH <- FindNeighbors(CH, reduction = "pca", dims = 1:6)
CH <- FindClusters(CH, resolution = 1)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:6)

# Setting identities so that every cell is identified
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

# Fig. S2F
Idents(CH) <- "types"
features = c("MESP1", "KDR", "MEST", "HNF1B", "SOX17", "FOXA2",
  "KRT7", "DLX5", "TFAP2A", "NANOG", "POU5F1", "SOX2")
my_levels <- c("Mesoderm", "Endoderm", "Ectoderm", "Pluripotent")
levels(CH) <- my_levels
pdf()
DimPlot(CH, reduction = "umap", pt.size = 2, order = c("Ectoderm","Endoderm","Pluripotent","Mesoderm"), cols = c("#eb0006","#e5c321","#fb6d92","#6b8ddd"))
DotPlot(CH, features = features, cols = c("#f0f0f0", "#f02207"), col.min = 0) + RotatedAxis()
dev.off()



## Setting up mesoderm-only WNTd dataset
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
DimPlot(meso, reduction = "umap", pt.size = 2, group.by = "RNA_snn_res.0.5") + NoLegend()
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
