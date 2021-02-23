library(Seurat)
library(viridis)
library(ggplot2)
library(clustree)

library(devtools)
install_github("sturgeonlab/Luff-etal-2021/Seurat/averagesc3")
library(averagesc3)

## Setting up WNTd-only dataset
load("CHIRSB-IWP2-before-integration.RData")
CH <- SubsetData(CHIRSB.IWP2, subset.name = "orig.ident", accept.value = "CHIRSB")
CH <- ScaleData(CH, features = rownames(CH))
CH <- RunPCA(CH, npcs = 50)
CH <- JackStraw(CH, num.replicate = 100, dims = 50)
CH <- ScoreJackStraw(CH, dims = 1:50)
pdf()
JackStrawPlot(CH, dims = 1:50)
dev.off()
CH <- FindNeighbors(CH, reduction = "pca", dims = 1:43)
CH <- FindClusters(CH, resolution = 0.7)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:43)
pdf()
DimPlot(CH, reduction = "umap", pt.size = 1.5)
dev.off()

# determine clustering resolution
CH <- FindClusters(CH, resolution = 0.0)
  CH <- FindClusters(CH, resolution = 0.1)
  CH <- FindClusters(CH, resolution = 0.2)
  CH <- FindClusters(CH, resolution = 0.3)
  CH <- FindClusters(CH, resolution = 0.4)
  CH <- FindClusters(CH, resolution = 0.5)
  CH <- FindClusters(CH, resolution = 0.6)
  CH <- FindClusters(CH, resolution = 0.7)
  CH <- FindClusters(CH, resolution = 0.8)
  CH <- FindClusters(CH, resolution = 0.9)
  CH <- FindClusters(CH, resolution = 1)
  CH <- FindClusters(CH, resolution = 1.1)
  CH <- FindClusters(CH, resolution = 1.2)
  CH <- FindClusters(CH, resolution = 1.3)
  CH <- FindClusters(CH, resolution = 1.4)
  pdf(width = 8, height = 14)
  clustree(CH, prefix = "RNA_snn_res.", layout = "sugiyama", use_core_edges = FALSE)
  clustree(CH, prefix = "RNA_snn_res.", layout = "sugiyama", use_core_edges = FALSE, node_colour = "sc3_stability") + scale_color_viridis(option="plasma")
  dev.off()

# sometimes it's challenging to decide, calculating highest average sc3_stability can help
res0 <- CH@meta.data$RNA_snn_res.0
  res0.1 <- CH@meta.data$RNA_snn_res.0.1
  res0.2 <- CH@meta.data$RNA_snn_res.0.2
  res0.3 <- CH@meta.data$RNA_snn_res.0.3
  res0.4 <- CH@meta.data$RNA_snn_res.0.4
  res0.5 <- CH@meta.data$RNA_snn_res.0.5
  res0.6 <- CH@meta.data$RNA_snn_res.0.6
  res0.7 <- CH@meta.data$RNA_snn_res.0.7
  res0.8 <- CH@meta.data$RNA_snn_res.0.8
  res0.9 <- CH@meta.data$RNA_snn_res.0.9
  res1 <- CH@meta.data$RNA_snn_res.1
  res1.1 <- CH@meta.data$RNA_snn_res.1.1
  res1.2 <- CH@meta.data$RNA_snn_res.1.2
  res1.3 <- CH@meta.data$RNA_snn_res.1.3
  res1.4 <- CH@meta.data$RNA_snn_res.1.3
#nrow is the number of cells in dataset
mat <- matrix(c(res0,res0.1,res0.2,res0.3,res0.4,res0.5,res0.6,res0.7,res0.8,
res0.9,res1,res1.1,res1.2,res1.3,res1.4), nrow = 6858, dimnames = list(c(rownames(CH@meta.data)),
c("res0","res0.1","res0.2","res0.3","res0.4","res0.5","res0.6","res0.7","res0.8","res0.9","res1",
"res1.1","res1.2","res1.3","res1.4")))
nodes <- get_tree_nodes(clusterings = mat, prefix = "res", node_aes_list = FALSE, metadata = FALSE)
write.table(nodes, file = "nodes.txt", sep = "\t")
#res0.2 highest average sc3_stability

# Fig. S2Ei
pdf()
DimPlot(CH, reduction = "umap", pt.size = 1.5, group.by = "RNA_snn_res.0.2") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(CH, reduction = "umap", pt.size = 1.5, group.by = "RNA_snn_res.0.2") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend()
dev.off()

## Setting identities so that every cell is labeled
# (round 1) first determine the threshold for key germ layer markers and find unlabeled cells
pdf()
  plot(density(CH@assays$RNA@data['KDR',])) #mesoderm
  abline(v=0.25)
  plot(density(CH@assays$RNA@data['FOXA2',])) #endoderm
  abline(v=0.2)
  plot(density(CH@assays$RNA@data['TFAP2A',])) #ectoderm
  abline(v=0.25)
  plot(density(CH@assays$RNA@data['SOX2',])) #pluripotent
  abline(v=0.25)
  dev.off()
Idents(CH) <- "RNA_snn_res.0.2"
  pluri <- WhichCells(object = CH, expression = SOX2 > 0.25)
  ecto <- WhichCells(object = CH, expression = TFAP2A > 0.25)
  endo <- WhichCells(object = CH, expression = FOXA2 > 0.2)
  meso <- WhichCells(object = CH, expression = KDR > 0.25)
  CH <- SetIdent(CH, cells = pluri, value = 'Pluripotent')
  CH <- SetIdent(CH, cells = ecto, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = endo, value = 'Endoderm')
  CH <- SetIdent(CH, cells = meso, value = 'Mesoderm')
  CH[["types"]] <- Idents(CH)
  pdf()
  DimPlot(CH, reduction = "umap", pt.size = 1.5)
  dev.off()

# run DGE to determine markers for unidentified cells
markers <- FindAllMarkers(CH, only.pos = TRUE) #only.pos = TRUE so that we can identify positive selection markers
write.table(markers, file="markers.txt", sep="\t")

# (round 2) select high scoring markers in cells from unidentified clusters (all clusters)
pdf()
  plot(density(CH@assays$RNA@data['MEST',])) #cluster0, broadly mesoderm
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['MESP1',])) #cluster0, broadly mesoderm
  abline(v=0.25)
  plot(density(CH@assays$RNA@data['KRT19',])) #cluster1, ectoderm
  abline(v=0.35)
  plot(density(CH@assays$RNA@data['DLX5',])) #cluster1, ectoderm
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['POU5F1',])) #cluster2, primitive stream
  abline(v=0.35)
  plot(density(CH@assays$RNA@data['TBXT',])) #cluster2, primitive stream
  abline(v=0.2)
  plot(density(CH@assays$RNA@data['TEK',])) #cluster3, endothelial
  abline(v=0.2)
  plot(density(CH@assays$RNA@data['FLT1',])) #cluster3, endothelial
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['APOA2',])) #cluster4, endoderm
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['APOA1',])) #cluster4, endoderm
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['GATA3',])) #cluster5, ectoderm
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['NANOS3',])) #cluster7, PGC
  abline(v=0.1)
  plot(density(CH@assays$RNA@data['DND1',])) #cluster7, PGC
  abline(v=0.1)
  dev.off()
Idents(CH) <- "RNA_snn_res.0.2"
  pluri2 <- WhichCells(CH, expression = NANOS3 > 0.1)
  pluri3 <- WhichCells(CH, expression = DND1 > 0.1)
  pluri4 <- WhichCells(CH, expression = POU5F1 > 0.35)
  pluri5 <- WhichCells(CH, expression = TBXT > 0.2)
  endo2 <- WhichCells(CH, expression = APOA2 > 0.3)
  endo3 <- WhichCells(CH, expression = APOA1 > 0.3)
  meso2 <- WhichCells(CH, expression = MEST > 0.3)
  meso3 <- WhichCells(CH, expression = MESP1 > 0.25)
  meso4 <- WhichCells(CH, expression = TEK > 0.2)
  meso5 <- WhichCells(CH, expression = FLT1 > 0.3)
  ecto2 <- WhichCells(CH, expression = DLX5 > 0.3)
  ecto3 <- WhichCells(CH, expression = GATA3 > 0.3)
  CH <- SetIdent(CH, cells = c(pluri2,pluri3,pluri4,pluri5), value = 'Pluripotent')
  CH <- SetIdent(CH, cells = c(endo2,endo3), value = 'Endoderm')
  CH <- SetIdent(CH, cells = c(ecto2,ecto3), value = 'Ectoderm')
  CH <- SetIdent(CH, cells = c(meso2,meso3,meso4,meso5), value = 'Mesoderm')
  #add round 1 labels last to pull cells back into labels using best known germ markers
  CH <- SetIdent(CH, cells = pluri, value = 'Pluripotent')
  CH <- SetIdent(CH, cells = ecto, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = endo, value = 'Endoderm')
  CH <- SetIdent(CH, cells = meso, value = 'Mesoderm')
CH[["types"]] <- Idents(CH) #replaces previous

# run DGE to determine markers for remaining unidentified cells
markers <- FindAllMarkers(CH, only.pos = TRUE)
write.table(markers, file="markers.txt", sep="\t")

# few remaining cells labeled based on what cluster they were in and
# the predominate germ layer in that cluster
# SetIdent in this order so that strongest markers are prioritized
Idents(CH) <- "RNA_snn_res.0.2"
  #round2
  CH <- SetIdent(CH, cells = c(ecto2,ecto3), value = 'Ectoderm')
  CH <- SetIdent(CH, cells = c(pluri2,pluri3,pluri4,pluri5), value = 'Pluripotent')
  CH <- SetIdent(CH, cells = c(endo2,endo3), value = 'Endoderm')
  CH <- SetIdent(CH, cells = c(meso2,meso3,meso4,meso5), value = 'Mesoderm')
  #round1
  CH <- SetIdent(CH, cells = ecto, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = pluri, value = 'Pluripotent')
  CH <- SetIdent(CH, cells = endo, value = 'Endoderm')
  CH <- SetIdent(CH, cells = meso, value = 'Mesoderm')
  #leftovers
  meso6 <- WhichCells(CH, idents = c(0,3,6))
  ecto5 <- WhichCells(CH, idents = 1)
  endo4 <- WhichCells(CH, idents = c(4,5))
  pluri6 <- WhichCells(CH, idents = 2)
  CH <- SetIdent(CH, cells = ecto5, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = pluri6, value = 'Pluripotent')
  CH <- SetIdent(CH, cells = endo4, value = 'Endoderm')
  CH <- SetIdent(CH, cells = meso6, value = 'Mesoderm')
  CH[["types"]] <- Idents(CH)


# Fig. S2Eii
pdf()
FeaturePlot(CH, features = c("KDR"), pt.size = 1.5, reduction = "umap", min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(CH, features = c("KDR"), pt.size = 1.5, reduction = "umap", min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank()) + NoLegend()
dev.off()

# Fig. S2Fi
pdf()
DimPlot(CH, reduction = "umap", pt.size = 1, order = c("Mesoderm","Pluripotent","Ectoderm","Endoderm"), cols = c("#00D142","#6D75F2","#EDC834","#E3312D")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(CH, reduction = "umap", pt.size = 1, order = c("Mesoderm","Pluripotent","Ectoderm","Endoderm"), cols = c("#00D142","#6D75F2","#EDC834","#E3312D")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend()
dev.off()

# Fig. S2Fii
Idents(CH) <- "types"
my_levels <- c("Mesoderm","Endoderm","Ectoderm","Pluripotent")
levels(CH) <- my_levels
features = c("SOX2","POU5F1","NANOG","TFAP2A","DLX5","KRT7",
"FOXA2","SOX17","HNF1B","MEST","KDR","MESP1")
pdf(width = 7, height = 2.75) #make it a little taller to get full scale bar
DotPlot(CH, features = features) + RotatedAxis() + scale_colour_gradient2(low="#003cd5", mid="white", high="#f0260a")
dev.off()

# Fig. S2Fiii
Idents(CH) <- "types"
my_levels <- c("Pluripotent","Ectoderm","Endoderm","Mesoderm")
levels(CH) <- my_levels
pdf()
VlnPlot(CH, features = "KDR", cols = c("#EDC834","#6D75F2","#00D142","#E3312D"), pt.size = 0.25)
dev.off()

save(CH, file = "WNTd.only.RData")


## Setting up mesoderm-only WNTd dataset
meso <- SubsetData(CH, subset.name = "KDR", low.threshold = 0.25)
meso <- NormalizeData(meso)
meso <- FindVariableFeatures(meso, selection.method = "vst", nfeatures = 2000)
meso <- ScaleData(meso, features = rownames(meso))
meso <- RunPCA(meso, features = VariableFeatures(object = meso))
meso <- JackStraw(meso, num.replicate = 100, dims = 50)
meso <- ScoreJackStraw(meso, dims = 1:50)
pdf()
JackStrawPlot(meso, dims = 1:50)
dev.off()
meso <- FindNeighbors(meso, dims = 1:44)
meso <- FindClusters(meso)
meso <- RunUMAP(meso, dims = 1:44)

# determine clustering resolution
meso <- FindClusters(meso, resolution = 0.0)
  meso <- FindClusters(meso, resolution = 0.1)
  meso <- FindClusters(meso, resolution = 0.2)
  meso <- FindClusters(meso, resolution = 0.3)
  meso <- FindClusters(meso, resolution = 0.4)
  meso <- FindClusters(meso, resolution = 0.5)
  meso <- FindClusters(meso, resolution = 0.6)
  meso <- FindClusters(meso, resolution = 0.7)
  meso <- FindClusters(meso, resolution = 0.8)
  meso <- FindClusters(meso, resolution = 0.9)
  meso <- FindClusters(meso, resolution = 1)
  meso <- FindClusters(meso, resolution = 1.1)
  meso <- FindClusters(meso, resolution = 1.2)
  meso <- FindClusters(meso, resolution = 1.3)
  meso <- FindClusters(meso, resolution = 1.4)
  pdf(width = 8, height = 14)
  clustree(meso, prefix = "RNA_snn_res.", layout = "sugiyama", use_core_edges = FALSE)
  clustree(meso, prefix = "RNA_snn_res.", layout = "sugiyama", use_core_edges = FALSE, node_colour = "sc3_stability") + scale_color_viridis(option="plasma")
  dev.off()

#calculate highest sc3_stability
res0 <- meso@meta.data$RNA_snn_res.0
  res0.1 <- meso@meta.data$RNA_snn_res.0.1
  res0.2 <- meso@meta.data$RNA_snn_res.0.2
  res0.3 <- meso@meta.data$RNA_snn_res.0.3
  res0.4 <- meso@meta.data$RNA_snn_res.0.4
  res0.5 <- meso@meta.data$RNA_snn_res.0.5
  res0.6 <- meso@meta.data$RNA_snn_res.0.6
  res0.7 <- meso@meta.data$RNA_snn_res.0.7
  res0.8 <- meso@meta.data$RNA_snn_res.0.8
  res0.9 <- meso@meta.data$RNA_snn_res.0.9
  res1 <- meso@meta.data$RNA_snn_res.1
  res1.1 <- meso@meta.data$RNA_snn_res.1.1
  res1.2 <- meso@meta.data$RNA_snn_res.1.2
  res1.3 <- meso@meta.data$RNA_snn_res.1.3
  res1.4 <- meso@meta.data$RNA_snn_res.1.3
#nrow is the number of cells in dataset
mat <- matrix(c(res0,res0.1,res0.2,res0.3,res0.4,res0.5,res0.6,res0.7,res0.8,
res0.9,res1,res1.1,res1.2,res1.3,res1.4), nrow = 2771, dimnames = list(c(rownames(meso@meta.data)),
c("res0","res0.1","res0.2","res0.3","res0.4","res0.5","res0.6","res0.7","res0.8","res0.9","res1",
"res1.1","res1.2","res1.3","res1.4")))
nodes <- get_tree_nodes(clusterings = mat, prefix = "res", node_aes_list = FALSE, metadata = FALSE)
write.table(nodes, file = "nodes.txt", sep = "\t")
#res0.2, res0.3, and res0.4 all have very high res, but lower resolutions fail to resolve mesoderm, so res0.4 chosen

Idents(meso) <- "RNA_snn_res.0.4"
meso <- RenameIdents(meso,'0'='4','1'='5','2'='1','3'='3','4'='2','5'='8','6'='6','7'='7','8'='9')
meso[["newclusters"]] <- Idents(meso)
Idents(meso) <- "newclusters"


# Fig. 1Di
pdf()
DimPlot(meso, reduction = "umap", pt.size = 2, group.by = "newclusters", order = c("9","8","7","6","5","4","3","2","1"),
cols = c("#f8766d","#db72fb","#00b9e3","#00ba38","#d39200","#ff61c3","#619cff","#00c19f","#93aa00")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(meso, reduction = "umap", pt.size = 2, group.by = "newclusters", order = c("9","8","7","6","5","4","3","2","1"),
cols = c("#f8766d","#db72fb","#00b9e3","#00ba38","#d39200","#ff61c3","#619cff","#00c19f","#93aa00")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend()
dev.off()


# Fig. 1Dii,iii
DefaultAssay(meso) <- "RNA"
pdf()
FeaturePlot(meso, features = c("ALDH1A2"), reduction = "umap", order = TRUE, pt.size = 2) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(meso, features = c("ALDH1A2"), reduction = "umap", order = TRUE, pt.size = 2) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend() + ggtitle(element_blank())
FeaturePlot(meso, features = c("CXCR4"), reduction = "umap", order = TRUE, pt.size = 2) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(meso, features = c("CXCR4"), reduction = "umap", order = TRUE, pt.size = 2) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend() + ggtitle(element_blank())
dev.off()

# Supplementary Table 2A
avg <- AverageExpression(meso, assay = "RNA")
write.table(avg, file="table2a.txt", sep="\t")

# Supplementary Table 2B
meso.markers <- FindAllMarkers(meso, only.pos = TRUE, logfc = 0.176)
write.table(meso.markers, file="table2b.txt", sep="\t")

# Supplementary Table 2C
aldh.pos <- WhichCells(meso, expression = ALDH1A2 >= 0.15)
aldh.neg <- WhichCells(meso, expression = ALDH1A2 < 0.15)
meso <- SetIdent(meso, cells = aldh.neg, value = "aldh.neg")
meso <- SetIdent(meso, cells = aldh.pos, value = "aldh.pos")
aldh.markers <- FindAllMarkers(meso, only.pos = TRUE, logfc = 0.176)
write.table(aldh.markers, file="table2c.txt", sep="\t")

save(meso, file = "WNTd.mesoderm.only.RData")
