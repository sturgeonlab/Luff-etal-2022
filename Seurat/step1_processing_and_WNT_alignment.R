library(Seurat)
library(viridis)
library(cowplot)
library(ggplot2)
library(clustree)

# load WNTd/WNTi data from GEO
CHIRSB.data <- Read10X(data.dir = "./CHIRSB/")
CHIRSB <- CreateSeuratObject(CHIRSB.data, assay = "RNA", min.cells = 3, min.features = 200, project = "CHIRSB")
IWP2.data <- Read10X(data.dir = "./IWP2/")
IWP2 <- CreateSeuratObject(IWP2.data, assay = "RNA", min.cells = 3, min.features = 200, project = "IWP2")
CHIRSB.IWP2 <- merge(CHIRSB,IWP2, project = "CHIRSB.IWP2")
CHIRSB.IWP2[["percent.mito"]] <- PercentageFeatureSet(CHIRSB.IWP2, pattern = "^MT-")
CHIRSB.IWP2[["merged"]] <- "CHIRSB.IWP2"

# Fig. S2A
pdf()
VlnPlot(CHIRSB.IWP2, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1, cols = c("#20a288","#20a288","#20a288"))
VlnPlot(CHIRSB.IWP2, group.by = "merged", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1, cols = c("#20a288","#20a288","#20a288"))
dev.off()

CHIRSB.IWP2 <- subset(CHIRSB.IWP2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 10)
CHIRSB.IWP2 <- NormalizeData(CHIRSB.IWP2, normalization.method = "LogNormalize", scale.factor = 10000)
CHIRSB.IWP2 <- FindVariableFeatures(CHIRSB.IWP2, selection.method = "vst")
CHIRSB.IWP2 <- ScaleData(CHIRSB.IWP2, features = rownames(CHIRSB.IWP2), vars.to.regress = c("percent.mito","nCount_RNA"))
#saved temp
CHIRSB.IWP2 <- RunPCA(CHIRSB.IWP2, features = VariableFeatures(CHIRSB.IWP2), dims = 50)
CHIRSB.IWP2 <- JackStraw(CHIRSB.IWP2, num.replicate = 100, dims = 50)
CHIRSB.IWP2 <- ScoreJackStraw(CHIRSB.IWP2, dims = 1:50)
pdf()
JackStrawPlot(CHIRSB.IWP2, dims = 1:50)
dev.off()


# Fig. S2B, left panel
CHIRSB.IWP2 <- RunUMAP(CHIRSB.IWP2, dims = 1:50)
pdf()
DimPlot(CHIRSB.IWP2, reduction = "umap", pt.size = 1, group.by = "orig.ident", cols = c("#2a59b8","#e22e2f")) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(CHIRSB.IWP2, reduction = "umap", pt.size = 1, group.by = "orig.ident", cols = c("#2a59b8","#e22e2f")) + NoLegend() +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
dev.off()
save(CHIRSB.IWP2, file = "CHIRSB-IWP2-before-integration.RData")


## Integrate CHIRSB & IWP2 samples
CHIRSB.IWP2.list <- SplitObject(CHIRSB.IWP2, split.by = "orig.ident")
for (i in 1:length(CHIRSB.IWP2.list)) {
    CHIRSB.IWP2.list[[i]] <- NormalizeData(CHIRSB.IWP2.list[[i]], verbose = FALSE)
    CHIRSB.IWP2.list[[i]] <- FindVariableFeatures(CHIRSB.IWP2.list[[i]], selection.method = "vst",
        nfeatures = 22401, verbose = FALSE)
}

reference.list <- CHIRSB.IWP2.list[c("CHIRSB", "IWP2")]
CHIRSB.IWP2.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
all_genes <- Reduce(intersect, lapply(reference.list, rownames))
CHIRSB.IWP2.integrated <- IntegrateData(anchorset = CHIRSB.IWP2.anchors, dims = 1:50, features.to.integrate = all_genes)
DefaultAssay(CHIRSB.IWP2.integrated) <- "integrated"
CHIRSB.IWP2.integrated <- ScaleData(CHIRSB.IWP2.integrated, verbose = FALSE)
CHIRSB.IWP2.integrated <- RunPCA(CHIRSB.IWP2.integrated, npcs = 50, verbose = FALSE)
CHIRSB.IWP2.integrated <- JackStraw(CHIRSB.IWP2.integrated, num.replicate = 100, dims = 50)
CHIRSB.IWP2.integrated <- ScoreJackStraw(CHIRSB.IWP2.integrated, dims = 1:50)
pdf()
JackStrawPlot(CHIRSB.IWP2.integrated, dims = 1:50)
dev.off()
CHIRSB.IWP2.integrated <- RunUMAP(CHIRSB.IWP2.integrated, reduction = "pca", dims = 1:50)
CHIRSB.IWP2.integrated <- FindNeighbors(CHIRSB.IWP2.integrated, reduction = "pca", dims = 1:50)
CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.5)

# Fig. 1A (also Fig. S2B, right panel)
pdf()
DimPlot(CHIRSB.IWP2.integrated, reduction = "umap", pt.size = 1, group.by = "orig.ident", cols = c("#e22e2f","#2a59b8"), order = "CHIRSB") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20))
DimPlot(CHIRSB.IWP2.integrated, reduction = "umap", pt.size = 1, group.by = "orig.ident", cols = c("#e22e2f","#2a59b8"), order = "CHIRSB") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + NoLegend()
dev.off()

# Fig. S2C
pdf()
DimPlot(CHIRSB.IWP2.integrated, reduction = "umap", pt.size = 1, group.by = "nFeature_RNA") + NoLegend() + scale_color_viridis(discrete=TRUE, direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
DimPlot(CHIRSB.IWP2.integrated, reduction = "umap", pt.size = 1, group.by = "nCount_RNA") + NoLegend() + scale_color_viridis(discrete=TRUE, direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
DimPlot(CHIRSB.IWP2.integrated, reduction = "umap", pt.size = 1, group.by = "percent.mito") + NoLegend() + scale_color_viridis(discrete=TRUE, direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()

# Fig. 1B
Idents(CHIRSB.IWP2.integrated) <- "orig.ident"
WNTi <- WhichCells(CHIRSB.IWP2.integrated, idents = "IWP2")
WNTd <- WhichCells(CHIRSB.IWP2.integrated, idents = "CHIRSB")
DefaultAssay(CHIRSB.IWP2.integrated) <- "RNA"
cdx_features <- list(c("CDX1","CDX2","CDX4"))
CHIRSB.IWP2.integrated <- AddModuleScore(CHIRSB.IWP2.integrated, features = cdx_features, name = 'cdx_features')
pdf()
FeaturePlot(CHIRSB.IWP2.integrated, cells = WNTi, features = c("KDR"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(CHIRSB.IWP2.integrated, cells = WNTd, features = c("KDR"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(CHIRSB.IWP2.integrated, cells = WNTi, features = c("GYPA"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(CHIRSB.IWP2.integrated, cells = WNTd, features = c("GYPA"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(CHIRSB.IWP2.integrated, cells = WNTi, features = c("cdx_features1"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
FeaturePlot(CHIRSB.IWP2.integrated, cells = WNTd, features = c("cdx_features1"), reduction = "umap", pt.size = 1.5, min.cutoff = 0, order = TRUE) + scale_color_viridis(direction = -1) +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
dev.off()

# Fig. S2D
DefaultAssay(CHIRSB.IWP2.integrated) <- "integrated"
CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.0)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.1)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.2)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.3)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.4)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.5)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.6)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.7)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.8)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 0.9)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 1)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 1.1)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 1.2)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 1.3)
  CHIRSB.IWP2.integrated <- FindClusters(CHIRSB.IWP2.integrated, resolution = 1.4)
pdf(width = 8, height = 14)
clustree(CHIRSB.IWP2.integrated, prefix = "integrated_snn_res.", layout = "sugiyama", use_core_edges = FALSE)
clustree(CHIRSB.IWP2.integrated, prefix = "integrated_snn_res.", layout = "sugiyama", use_core_edges = FALSE, node_colour = "sc3_stability") + scale_color_viridis(option="plasma")
dev.off()
pdf()
DimPlot(CHIRSB.IWP2.integrated, reduction = "umap", pt.size = 1, group.by = "integrated_snn_res.1") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank())
DimPlot(CHIRSB.IWP2.integrated, reduction = "umap", pt.size = 1, group.by = "integrated_snn_res.1") +
theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 20)) + ggtitle(element_blank()) + NoLegend()
dev.off()

# Supplementary table 1
DefaultAssay(CHIRSB.IWP2.integrated) <- "RNA"
Idents(CHIRSB.IWP2.integrated) <- "integrated_snn_res.1"
markers <- FindAllMarkers(CHIRSB.IWP2.integrated, logfc = 0.176)
write.table(markers, file="all-markers.txt", sep="\t")

# Fig. S2E
pdf()
plot(density(CHIRSB.IWP2.integrated@assays$RNA@data['CDX4',]))
abline(v=0.2)
plot(density(CHIRSB.IWP2.integrated@assays$RNA@data['GYPA',]))
abline(v=0.1)
plot(density(CHIRSB.IWP2.integrated@assays$RNA@data['KDR',]))
abline(v=0.25)
dev.off()
WNTd <- WhichCells(CHIRSB.IWP2.integrated, expression = KDR > 0.25 & GYPA < 0.1, idents = "CHIRSB")
WNTi <- WhichCells(CHIRSB.IWP2.integrated, expression = KDR > 0.25 & GYPA > 0.1, idents = "IWP2")
CHIRSB.IWP2.integrated <- SetIdent(CHIRSB.IWP2.integrated, cells = WNTi, value = 'WNTi GYPA+')
CHIRSB.IWP2.integrated <- SetIdent(CHIRSB.IWP2.integrated, cells = WNTd, value = 'WNTd GYPA-')
p1 <- VlnPlot(CHIRSB.IWP2.integrated, features = c("CYP26A1"), idents = "WNTi GYPA+", pt.size = 0.5, log = TRUE, cols = c("#e12e2e")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p2 <- VlnPlot(CHIRSB.IWP2.integrated, features = c("ALDH1A2"), idents = "WNTi GYPA+", pt.size = 0.5, log = TRUE, cols = c("#e12e2e")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p3 <- VlnPlot(CHIRSB.IWP2.integrated, features = c("CYP26A1"), idents = "WNTd GYPA-", pt.size = 0.5, log = TRUE, cols = c("#2959b9")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p4 <- VlnPlot(CHIRSB.IWP2.integrated, features = c("ALDH1A2"), idents = "WNTd GYPA-", pt.size = 0.5, log = TRUE, cols = c("#2959b9")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
pdf(width = 5, height = 2)
plot_grid(
p1, p2, p3, p4,
nrow = 1)
dev.off()

save(CHIRSB.IWP2.integrated, markers, file = "CHIRSB-IWP2-after-integration.RData")
