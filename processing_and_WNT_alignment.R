library(Seurat)
library(dplyr)
library(Matrix)

## First steps of data processing were done with Seurat V2, dataset was update prior
## to alignment and further processing, as shown below
agg.data <- Read10X(data.dir = "hg19/")
agg <- CreateSeuratObject(agg.data, assay = "RNA", min.cells = 3, min.features = 200, project = "CHIR-IWP")
agg[["percent.mt"]] <- PercentageFeatureSet(agg, pattern = "^MT-")
agg <- subset(agg, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normalize data and find variable genes
agg <- NormalizeData(object = agg, normalization.method = "LogNormalize", scale.factor = 10000)
agg <- FindVariableGenes(object = agg, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# Regression
agg <- ScaleData(object = agg, vars.to.regress = c("nUMI", "percent.mito"))

# Find significant PCs
agg <- RunPCA(object = agg, pc.genes = agg@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = agg)

# Dimension reduction
agg <- FindClusters(object = agg, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
agg <- RunTSNE(object = agg, dims.use = 1:20, do.fast = TRUE)



## All steps going forward are using Seurat V3
agg <- UpdateSeuratObject(agg)
save(agg, file = "agg.RData")


# Fig. S2A
library(viridis)
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.5, group.by = "nFeature_RNA") + theme(legend.position = "none") + scale_color_viridis(discrete=TRUE)
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.5, group.by = "nCount_RNA") + theme(legend.position = "none") + scale_color_viridis(discrete=TRUE)
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.5, group.by = "percent.mito") + theme(legend.position = "none") + scale_color_viridis(discrete=TRUE)
dev.off()

# Fig. S2B, left panel
agg <- RunUMAP(agg, dims = 1:6)
pdf()
DimPlot(agg, reduction = "umap", pt.size = 1.5, group.by ="Cell.ID")
dev.off()

# Split combined dataset for integration
agg.list <- SplitObject(agg, split.by = "Cell.ID")
for (i in 1:length(agg.list)) {
    agg.list[[i]] <- NormalizeData(agg.list[[i]], verbose = FALSE)
    agg.list[[i]] <- FindVariableFeatures(agg.list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
}

# Integrate CHIRSB & IWP2 samples
reference.list <- agg.list[c("CHIRSB", "IWP2")]
agg.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:26)
agg.integrated <- IntegrateData(anchorset = agg.anchors, dims = 1:26)

# Process integrated dataset
DefaultAssay(agg.integrated) <- "integrated"
agg.integrated <- ScaleData(agg.integrated, verbose = FALSE)
agg.integrated <- RunPCA(agg.integrated, npcs = 30, verbose = FALSE)
pdf()
ElbowPlot(meso)
dev.off()
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:26)
agg.integrated <- FindNeighbors(agg.integrated, reduction = "pca", dims = 1:26)
agg.integrated <- FindClusters(agg.integrated, resolution = 0.5)

# Fig. S2B, right panel
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75)
dev.off()

# Save file
save(agg.integrated, file = "agg-integrated-26dim.RData")



## Cell Cycle Regression Setup
# Assign cell cycle scores
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
agg.integrated <- CellCycleScoring(agg.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))

# Fig. S2C, no regression
pdf()
DimPlot(agg.integrated, reduction = "pca")
dev.off()



## G2M/S (Partial) Cell Cycle Regression
agg.integrated <- CellCycleScoring(agg.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))
agg.integrated$CC.Difference <- agg.integrated$S.Score - agg.integrated$G2M.Score
agg.integrated <- ScaleData(agg.integrated, vars.to.regress = "CC.Difference", features = rownames(agg.integrated))
agg.integrated <- RunPCA(agg.integrated, features = VariableFeatures(agg.integrated), nfeatures.print = 10)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))

# Fig. S2C, G2M/S regression, PCA plot
pdf()
DimPlot(agg.integrated, reduction = "pca")
dev.off()

# Finishing clustering partial regression
---agg.integrated <- RunPCA(agg.integrated)
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:26)
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75)
dev.off()
save(agg.integrated, file = "agg-integrated-cclite.RData")



## Full Regression; reload dataset to undo partial regression
load("agg-integrated-26dim.RData")
agg.integrated <- ScaleData(agg.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(agg.integrated))
agg.integrated <- RunPCA(agg.integrated, features = VariableFeatures(agg.integrated), nfeatures.print = 10)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))

# Fig. S2C, full regression
pdf()
DimPlot(agg.integrated, reduction = "pca")
dev.off()

# Finishing clustering full regression
---agg.integrated <- RunPCA(agg.integrated)
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:26)
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75)
dev.off()
save(agg.integrated, file = "agg-integrated-cc-clusters.RData")

# Fig. S2D
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "percent.mito")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "nUMI")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "nGene")
dev.off()

# Fig. 1B
pdf()
FeaturePlot(agg.integrated, features = c("KDR"), split.by = "Cell.ID", reduction = "umap", pt.size = 0.5, min.cutoff = 0, cols = c("#f0f0f0", "#f02207"), order = TRUE) + theme(legend.position = "none")
FeaturePlot(agg.integrated, features = c("GYPA"), split.by = "Cell.ID", reduction = "umap", pt.size = 0.5, min.cutoff = 0, cols = c("#f0f0f0", "#f02207"), order = TRUE) + theme(legend.position = "none")
FeaturePlot(agg.integrated, features = c("CDX4"), split.by = "Cell.ID", reduction = "umap", pt.size = 0.5, min.cutoff = 0, cols = c("#f0f0f0", "#f02207"), order = TRUE) + theme(legend.position = "none")
dev.off()

# Fig. S2E
DefaultAssay(agg.integrated) <- "RNA"
Idents(agg.integrated) <- "Cell.ID"
wntd <- WhichCells(object = agg.integrated, expression = KDR > 0.1 & GYPA < 1, idents = "CHIRSB")
wnti <- WhichCells(object = agg.integrated, expression = KDR > 0.1 & GYPA > 1, idents = "IWP2")
agg.integrated <- SetIdent(object = agg.integrated, cells = wnti, value = 'WNTi')
agg.integrated <- SetIdent(object = agg.integrated, cells = wntd, value = 'WNTd')
p1 <- VlnPlot(agg.integrated, features = c("CYP26A1"), idents = "WNTi", pt.size = 0.5, log = TRUE, cols = c("#e12e2e")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p2 <- VlnPlot(agg.integrated, features = c("ALDH1A2"), idents = "WNTi", pt.size = 0.5, log = TRUE, cols = c("#e12e2e")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p3 <- VlnPlot(agg.integrated, features = c("CYP26A1"), idents = "WNTd", pt.size = 0.5, log = TRUE, cols = c("#2959b9")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
p4 <- VlnPlot(agg.integrated, features = c("ALDH1A2"), idents = "WNTd", pt.size = 0.5, log = TRUE, cols = c("#2959b9")) + NoLegend() + FontSize(x.title = 10, y.title = 6) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 0)) + scale_y_continuous(trans="log2")
pdf(width = 5, height = 2)
plot_grid(
p1, p2, p3, p4,
nrow = 1)
dev.off()
