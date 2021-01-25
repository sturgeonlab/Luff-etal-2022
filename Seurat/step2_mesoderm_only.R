library(Seurat)
library(viridis)
library(ggplot2)

## Setting up WNTd-only dataset
load("CHIR-IWP-before-integration.RData")
Idents(agg) <- "Cell.ID"
CH <- SubsetData(agg, subset.name = "Cell.ID", accept.value = "CHIRSB")
CH <- ScaleData(CH, features = rownames(CH))
CH <- RunPCA(CH, npcs = 50)
CH <- JackStraw(CH, num.replicate = 100, dims = 50)
CH <- ScoreJackStraw(CH, dims = 1:50)
pdf()
JackStrawPlot(CH, dims = 1:50)
dev.off()
CH <- FindNeighbors(CH, reduction = "pca", dims = 1:50)
CH <- FindClusters(CH, resolution = 1)
CH <- RunUMAP(CH, reduction = "pca", dims = 1:50)
pdf()
DimPlot(CH, reduction = "umap", pt.size = 1.5)
dev.off()


## Setting identities so that every cell is labeled
# (round 1) first determine the threshold for key germ layer markers and find unlabeled cells
pdf()
  plot(density(CH@assays$RNA@data['KDR',])) #mesoderm
  abline(v=0.25)
  plot(density(CH@assays$RNA@data['FOXA2',])) #endoderm
  abline(v=0.2)
  plot(density(CH@assays$RNA@data['TFAP2A',])) #ectoderm
  abline(v=0.15)
  plot(density(CH@assays$RNA@data['SOX2',])) #pluripotent
  abline(v=0.25)
  dev.off()
Idents(CH) <- "RNA_snn_res.1"
  pluri <- WhichCells(object = CH, expression = SOX2 > 0.25)
  ecto <- WhichCells(object = CH, expression = TFAP2A > 0.15)
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
markers <- FindAllMarkers(CH, only.pos = TRUE)
write.table(markers, file="markers.txt", sep="\t")


# (round 2) select high scoring markers that are also known to be expressed within a germ layer
pdf()
  plot(density(CH@assays$RNA@data['PITX1',])) #cluster 0, posterior mesoderm
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['MEST',])) #cluster1, broadly mesoderm
  abline(v=0.25)
  plot(density(CH@assays$RNA@data['EPAS1',])) #cluster2, endothelial (mesoderm)
  abline(v=0.25)
  plot(density(CH@assays$RNA@data['IHH',])) #cluster3, endoderm
  abline(v=0.15)
  plot(density(CH@assays$RNA@data['POU5F1',])) #cluster4, primitive streak
  abline(v=0.3)
  plot(density(CH@assays$RNA@data['NODAL',])) #cluster4/6, primitive streak
  abline(v=0.2)
  plot(density(CH@assays$RNA@data['DNAH2',])) #cluster9, ectoderm
  abline(v=0.1)
  dev.off()
Idents(CH) <- "RNA_snn_res.1"
  pluri2 <- WhichCells(CH, expression = POU5F1 > 0.3)
  pluri3 <- WhichCells(CH, expression = NODAL > 0.2)
  ecto2 <- WhichCells(CH, expression = DNAH2 > 0.1)
  endo2 <- WhichCells(CH, expression = IHH > 0.15)
  meso2 <- WhichCells(CH, expression = PITX1 > 0.3)
  meso3 <- WhichCells(CH, expression = MEST > 0.25)
  meso4 <- WhichCells(CH, expression = EPAS1 > 0.25)
  CH <- SetIdent(CH, cells = c(pluri2,pluri3), value = 'Pluripotent')
  CH <- SetIdent(CH, cells = ecto2, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = endo2, value = 'Endoderm')
  CH <- SetIdent(CH, cells = c(meso2,meso3,meso4), value = 'Mesoderm')
  #add round 1 labels last to pull cells back into labels using best known germ markers
  CH[["types"]] <- Idents(CH) #replaces previous


# run DGE to determine markers for unidentified cells
markers <- FindAllMarkers(CH, only.pos = TRUE)
write.table(markers, file="markers.txt", sep="\t")


# (round3) select high scoring markers that are also known to be expressed within a germ layer
pdf()
  plot(density(CH@assays$RNA@data['SSPN',])) #cluster 0, cardiogenic
  abline(v=0.03)
  plot(density(CH@assays$RNA@data['SENCR',])) #cluster1, muscle/endo
  abline(v=0.04)
  plot(density(CH@assays$RNA@data['HNF1B',])) #cluster3, endoderm
  abline(v=0.1)
  plot(density(CH@assays$RNA@data['MYO6',])) #cluster7, cardiogenic
  abline(v=0.2)
  plot(density(CH@assays$RNA@data['NELL1',])) #cluster9, ectoderm
  abline(v=0.02)
  dev.off()
Idents(CH) <- "RNA_snn_res.1"
  meso5 <- WhichCells(CH, expression = SSPN > 0.03)
  meso6 <- WhichCells(CH, expression = SENCR > 0.04)
  endo3 <- WhichCells(CH, expression = HNF1B > 0.1)
  meso7 <- WhichCells(CH, expression = MYO6 > 0.2)
  ecto3 <- WhichCells(CH, expression = NELL1 > 0.02)
  CH <- SetIdent(CH, cells = ecto3, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = endo3, value = 'Endoderm')
  CH <- SetIdent(CH, cells = c(meso5,meso6,meso7), value = 'Mesoderm')
  #add round 2 next to pull cells back into labels using more significant markers
  #add round 1 labels last to pull cells back into labels using best known germ markers
  CH[["types"]] <- Idents(CH) #replaces previous


# few remaining cells did not have enriched markers; therefore they were labeled based on what
# cluster they were in and the predominate germ layer in that cluster
# reverse ordering for rounds so that strongest markers are used more frequently

Idents(CH) <- "RNA_snn_res.1"
  #round3
  CH <- SetIdent(CH, cells = ecto3, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = endo3, value = 'Endoderm')
  CH <- SetIdent(CH, cells = c(pluri2,pluri3), value = 'Pluripotent')
  CH <- SetIdent(CH, cells = c(meso5,meso6,meso7), value = 'Mesoderm')
  #round2
  CH <- SetIdent(CH, cells = ecto2, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = endo2, value = 'Endoderm')
  CH <- SetIdent(CH, cells = c(meso2,meso3,meso4), value = 'Mesoderm')
  #round1
  CH <- SetIdent(CH, cells = ecto, value = 'Ectoderm')
  CH <- SetIdent(CH, cells = pluri, value = 'Pluripotent')
  CH <- SetIdent(CH, cells = endo, value = 'Endoderm')
  CH <- SetIdent(CH, cells = meso, value = 'Mesoderm')
  #leftovers
  meso8 <- WhichCells(CH, idents = c(0,3,5,8,9,12,13))
  pluri4 <- WhichCells(CH, idents = 1)
  CH <- SetIdent(CH, cells = pluri4, value = 'Pluripotent')
  CH <- SetIdent(CH, cells = meso8, value = 'Mesoderm')
  CH[["types"]] <- Idents(CH)


# Fig. S2Ei,ii
Idents(CH) <- "types"
my_levels <- c("Mesoderm", "Endoderm", "Ectoderm", "Pluripotent")
levels(CH) <- my_levels
pdf()
FeaturePlot(CH, features = c("KDR"), pt.size = 1.5, reduction = "umap", min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(CH, features = c("KDR"), pt.size = 1.5, reduction = "umap", min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank()) + NoLegend()
DimPlot(CH, reduction = "umap", pt.size = 1.5, order = c("Ectoderm","Endoderm","Pluripotent","Mesoderm"), cols = c("#eb0006","#e5c321","#fb6d92","#6b8ddd")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(CH, reduction = "umap", pt.size = 1.5, order = c("Ectoderm","Endoderm","Pluripotent","Mesoderm"), cols = c("#eb0006","#e5c321","#fb6d92","#6b8ddd")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend()
dev.off()

save(CH, file = "WNTd.RData")


# Fig. S2Eiii
features = c("SOX2","POU5F1","NANOG","TFAP2A","DLX5","KRT7",
"FOXA2","SOX17","HNF1B","MEST","KDR","MESP1")
pdf(width = 7, height = 2.75)
DotPlot(CH, features = features) + RotatedAxis() + scale_colour_gradient2(low="#003cd5", mid="white", high="#f0260a")
dev.off()


####################################################################


## Setting up mesoderm-only WNTd dataset
# input can be non-regressed or cell cycle regressed
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
meso <- FindNeighbors(meso, dims = 1:40)
meso <- FindClusters(meso)
meso <- RunUMAP(meso, dims = 1:40)
save(meso, file = "WNTd.mesoderm.RData")

# determine clustering resolution
library('clustree')
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
# between 0.8 and 1 gives nearly identical favorable results

Idents(meso) <- "RNA_snn_res.0.8"
meso <- RenameIdents(meso,'0'='9','1'='5','2'='4','3'='2','4'='6','5'='3',
'6'='1','7'='12','8'='7','9'='8','10'='10','11'='11')
meso[["newclusters"]] <- Idents(meso)
Idents(meso) <- "newclusters"


# Fig. 1Ci
pdf()
DimPlot(meso, reduction = "umap", pt.size = 2, group.by = "newclusters", order = c("12","11","10","9","8","7",
"6","5","4","3","2","1"), cols = c("#00bfc4","#7cae00","#00b4f0","#de8c00","#00ba38",
"#b79f00","#619cff","#c77cff","#f8766d","#00c08b","#f564e3","#ff64b0")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(meso, reduction = "umap", pt.size = 2, group.by = "newclusters", order = c("12","11","10","9","8","7",
"6","5","4","3","2","1"), cols = c("#00bfc4","#7cae00","#00b4f0","#de8c00","#00ba38",
"#b79f00","#619cff","#c77cff","#f8766d","#00c08b","#f564e3","#ff64b0")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
dev.off()


# Fig. 1Cii,iv
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

# Fig. 1Ciii
Idents(meso) <- "newclusters"
my_levels <- c("1","2","3","4","5","6","7","8","9","10","11","12")
levels(meso) <- my_levels
pdf(width = 2.5, height = 6)
VlnPlot(meso, idents = c("4","6","7","10","11"), features = c("ALDH1A2","TEK"), assay = "RNA", log = TRUE, pt.size = 0.25, ncol = 1)
dev.off()

# Supplementary Table 2A
avg <- AverageExpression(meso, assay = "RNA")
write.table(avg, file="newclusters-avg.txt", sep="\t")

# Supplementary Table 2B
meso.markers <- FindAllMarkers(meso, logfc.threshold = 0.176)
write.table(meso.markers, file="newclusters-mesomarkers.txt", sep="\t")

# Supplementary Table 2C
pdf()
plot(density(meso@assays$RNA@data['ALDH1A2',]))
abline(v=0.1)
dev.off()
aldh.pos <- WhichCells(meso, expression = ALDH1A2 > 0.1)
aldh.neg <- WhichCells(meso, expression = ALDH1A2 < 0.1)
meso <- SetIdent(meso, cells = aldh.neg, value = "aldh.neg")
meso <- SetIdent(meso, cells = aldh.pos, value = "aldh.pos")
aldh.markers <- FindAllMarkers(meso, logfc.threshold = 0.176) #0.176 logfc = 1.5 linear fc
write.table(aldh.markers, file="aldhmesomarkers.txt", sep="\t")

# Supplementary Table 2D
sub <- SubsetData(meso, subset.name = "newclusters", accept.value = c("1","2","3","4","5","6","10","11","12"))
pdf()
plot(density(sub@assays$RNA@data['ALDH1A2',]))
abline(v=0.15)
dev.off()
aldh.pos <- WhichCells(sub, expression = ALDH1A2 > 0.15)
aldh.neg <- WhichCells(sub, expression = ALDH1A2 < 0.15)
sub <- SetIdent(sub, cells = aldh.neg, value = "aldh.neg")
sub <- SetIdent(sub, cells = aldh.pos, value = "aldh.pos")
aldh.markers <- FindAllMarkers(sub, only.pos = TRUE, logfc.threshold = 0)
write.table(aldh.markers, file="aldh.sub.markers.txt", sep="\t")
