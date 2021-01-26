library(Seurat)
library(viridis)
library(cowplot)
library(ggplot2)

# Workflow using V3.2.2
agg.data <- Read10X(data.dir = "./hg19/")
agg <- CreateSeuratObject(agg.data, assay = "RNA", min.cells = 3, min.features = 200, project = "CHIR-IWP")
agg[["percent.mito"]] <- PercentageFeatureSet(agg, pattern = "^MT-")
agg <- subset(agg, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 10)
agg <- NormalizeData(agg, normalization.method = "LogNormalize", scale.factor = 10000)
agg <- FindVariableFeatures(agg, selection.method = "vst")
all.genes <- rownames(agg)
agg <- ScaleData(agg, features = all.genes, vars.to.regress = c("percent.mito","nCount_RNA"))
agg <- RunPCA(agg, features = VariableFeatures(agg), dims = 50)
agg <- JackStraw(agg, num.replicate = 100, dims = 50)
agg <- ScoreJackStraw(agg, dims = 1:50)
pdf()
JackStrawPlot(agg, dims = 1:50)
dev.off()

# add WNTd (CHIRSB) or WNTi (IWP2) metadata
ID <- read.table("barcode_IDs.txt", header = TRUE, sep ="\t")
IDonly <- ID$Cell.ID
agg <- AddMetaData(agg, metadata = IDonly, col.name = "Cell.ID")

# Fig. S2A
pdf()
VlnPlot(agg, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0, cols = "#20a288")
dev.off()

# Fig. S2B, left panel
agg <- RunUMAP(agg, dims = 1:50)
pdf()
DimPlot(agg, reduction = "umap", pt.size = 1, group.by = "Cell.ID", cols = c("#2a59b8","#e22e2f")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(agg, reduction = "umap", pt.size = 1, group.by = "Cell.ID", cols = c("#2a59b8","#e22e2f")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
dev.off()
save(agg, file = "CHIR-IWP-before-integration.RData")


####################################################################

## Integrate CHIRSB & IWP2 samples
# 22732 features chosen so cell cycle genes will be included for downstream regression
agg.list <- SplitObject(agg, split.by = "Cell.ID")
for (i in 1:length(agg.list)) {
    agg.list[[i]] <- NormalizeData(agg.list[[i]], verbose = FALSE)
    agg.list[[i]] <- FindVariableFeatures(agg.list[[i]], selection.method = "vst",
        nfeatures = 22732, verbose = FALSE)
}

reference.list <- agg.list[c("CHIRSB", "IWP2")]
agg.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
all_genes <- Reduce(intersect, lapply(reference.list, rownames))
agg.integrated <- IntegrateData(anchorset = agg.anchors, dims = 1:50, features.to.integrate = all_genes)
DefaultAssay(agg.integrated) <- "integrated"
agg.integrated <- ScaleData(agg.integrated, verbose = FALSE)
agg.integrated <- RunPCA(agg.integrated, npcs = 50, verbose = FALSE)
agg.integrated <- JackStraw(agg.integrated, num.replicate = 100, dims = 50)
agg.integrated <- ScoreJackStraw(agg.integrated, dims = 1:50)
pdf()
JackStrawPlot(agg.integrated, dims = 1:50)
dev.off()
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:50)
agg.integrated <- FindNeighbors(agg.integrated, reduction = "pca", dims = 1:50)
agg.integrated <- FindClusters(agg.integrated, resolution = 0.5)


# Fig. 1A (also Fig. S2B, right panel)
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 1, group.by = "Cell.ID", cols = c("#e22e2f","#2a59b8"), order = "CHIRSB") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(agg.integrated, reduction = "umap", pt.size = 1, group.by = "Cell.ID", cols = c("#e22e2f","#2a59b8"), order = "CHIRSB") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend()
dev.off()

# Save file
save(agg.integrated, file = "CHIR-IWP-integrated.RData")


#######################################################################

# Fig. S2C
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 1, group.by = "nFeature_RNA") + NoLegend() + scale_color_viridis(discrete=TRUE, direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
DimPlot(agg.integrated, reduction = "umap", pt.size = 1, group.by = "nCount_RNA") + NoLegend() + scale_color_viridis(discrete=TRUE, direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
DimPlot(agg.integrated, reduction = "umap", pt.size = 1, group.by = "percent.mito") + NoLegend() + scale_color_viridis(discrete=TRUE, direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()

# Fig. 1B
Idents(agg.integrated) <- "Cell.ID"
WNTi <- WhichCells(agg.integrated, idents = "IWP2")
WNTd <- WhichCells(agg.integrated, idents = "CHIRSB")
DefaultAssay(agg.integrated) <- "RNA"
pdf()
FeaturePlot(agg.integrated, cells = WNTi, features = c("KDR"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(agg.integrated, cells = WNTd, features = c("KDR"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())

FeaturePlot(agg.integrated, cells = WNTi, features = c("GYPA"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(agg.integrated, cells = WNTd, features = c("GYPA"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())

FeaturePlot(agg.integrated, cells = WNTi, features = c("CDX4"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(agg.integrated, cells = WNTd, features = c("CDX4"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()

# Supplementary Table 1
agg.integrated <- SetIdent(agg.integrated, cells = WNTi, value = 'WNTi')
agg.integrated <- SetIdent(agg.integrated, cells = WNTd, value = 'WNTd')
pdf()
plot(density(agg.integrated@assays$RNA@data['CDX4',]))
abline(v=0.2)
plot(density(agg.integrated@assays$RNA@data['GYPA',]))
abline(v=0.1)
plot(density(agg.integrated@assays$RNA@data['KDR',]))
abline(v=0.25)
dev.off()
WNTd.CDX4 <- WhichCells(agg.integrated, expression = KDR > 0.25 & CDX4 > 0.2, idents = "WNTd")
WNTd.CDX4neg <- WhichCells(agg.integrated, expression = KDR > 0.25 & CDX4 < 0.2, idents = "WNTd")
WNTi.GYPA <- WhichCells(agg.integrated, expression = KDR > 0.25 & GYPA > 0.1, idents = "WNTi")
WNTi.GYPAneg <- WhichCells(agg.integrated, expression = KDR > 0.25 & GYPA < 0.1, idents = "WNTi")
agg.integrated <- SetIdent(agg.integrated, cells = WNTd.CDX4, value = 'WNTd.CDX4')
agg.integrated <- SetIdent(agg.integrated, cells = WNTd.CDX4neg, value = 'WNTd.CDX4neg')
agg.integrated <- SetIdent(agg.integrated, cells = WNTi.GYPA, value = 'WNTi.GYPA')
agg.integrated <- SetIdent(agg.integrated, cells = WNTi.GYPAneg, value = 'WNTi.GYPAneg')
agg.integrated[["subsets"]] <- Idents(agg.integrated)
CDX.markers <- FindMarkers(agg.integrated, ident.1 = "WNTd.CDX4", ident.2 = "WNTd.CDX4neg", only.pos = TRUE, logfc = 0.176)
CDXneg.markers <- FindMarkers(agg.integrated, ident.1 = "WNTd.CDX4neg", ident.2 = "WNTd.CDX4", only.pos = TRUE, logfc = 0.176)
GYP.markers <- FindMarkers(agg.integrated, ident.1 = "WNTi.GYPA", ident.2 = "WNTi.GYPAneg", only.pos = TRUE, logfc = 0.176)
GYPneg.markers <- FindMarkers(agg.integrated, ident.1 = "WNTi.GYPAneg", ident.2 = "WNTi.GYPA", only.pos = TRUE, logfc = 0.176)
write.table(CDX.markers, file="cdx-markers.txt", sep="\t") #Table 1A
write.table(CDXneg.markers, file="cdxneg-markers.txt", sep="\t") #Table 1A
write.table(GYP.markers, file="gyp-markers.txt", sep="\t") #Table 1B
write.table(GYPneg.markers, file="gypneg-markers.txt", sep="\t") #Table 1B


# Fig. S2D
DefaultAssay(agg.integrated) <- "RNA"
Idents(agg.integrated) <- "Cell.ID"
WNTd2 <- WhichCells(agg.integrated, expression = KDR > 0.1 & GYPA < 1, idents = "CHIRSB")
WNTi2 <- WhichCells(agg.integrated, expression = KDR > 0.1 & GYPA > 1, idents = "IWP2")
agg.integrated <- SetIdent(agg.integrated, cells = WNTi2, value = 'WNTi GYPA+')
agg.integrated <- SetIdent(agg.integrated, cells = WNTd2, value = 'WNTd GYPA-')
p1 <- VlnPlot(agg.integrated, features = c("CYP26A1"), idents = "WNTi", pt.size = 0.5, log = TRUE, cols = c("#e12e2e")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p2 <- VlnPlot(agg.integrated, features = c("ALDH1A2"), idents = "WNTi", pt.size = 0.5, log = TRUE, cols = c("#e12e2e")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p3 <- VlnPlot(agg.integrated, features = c("CYP26A1"), idents = "WNTd", pt.size = 0.5, log = TRUE, cols = c("#2959b9")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p4 <- VlnPlot(agg.integrated, features = c("ALDH1A2"), idents = "WNTd", pt.size = 0.5, log = TRUE, cols = c("#2959b9")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
pdf(width = 5, height = 2)
plot_grid(
p1, p2, p3, p4,
nrow = 1)
dev.off()
