library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)

# Workflow using V3.2.2
agg.data <- Read10X(data.dir = "./hg19/")
agg <- CreateSeuratObject(agg.data, assay = "RNA", min.cells = 3, min.features = 200, project = "CHIR-IWP")
agg[["percent.mito"]] <- PercentageFeatureSet(agg, pattern = "^MT-")
agg <- subset(agg, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 10)
agg <- NormalizeData(agg, normalization.method = "LogNormalize", scale.factor = 10000)
agg <- FindVariableFeatures(agg, selection.method = "vst")
all.genes <- rownames(agg)
agg <- ScaleData(agg, features = all.genes)
agg <- RunPCA(agg, features = VariableFeatures(object = agg))
pdf()
ElbowPlot(agg)
dev.off()
agg <- FindNeighbors(agg, dims = 1:6)
agg <- FindClusters(agg, resolution = 0.5)
ID <- read.table("barcode_IDs.txt", header = TRUE, sep ="\t")
IDonly <- ID$Cell.ID
agg <- AddMetaData(agg, metadata = IDonly, col.name = "Cell.ID")

# Fig. S2A
library(viridis)
pdf()
VlnPlot(agg, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)
dev.off()

# Fig. S2B, left panel
agg <- RunUMAP(agg, dims = 1:6)
pdf()
DimPlot(agg, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID", cols = c("#2a59b8","#e22e2f"))
DimPlot(agg, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID", cols = c("#2a59b8","#e22e2f")) + NoLegend()
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
agg.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
all_genes <- Reduce(intersect, lapply(reference.list, rownames))
agg.integrated <- IntegrateData(anchorset = agg.anchors, dims = 1:30, features.to.integrate = all_genes)
DefaultAssay(agg.integrated) <- "integrated"
agg.integrated <- ScaleData(agg.integrated, verbose = FALSE)
agg.integrated <- RunPCA(agg.integrated, npcs = 30, verbose = FALSE)
pdf()
ElbowPlot(agg.integrated)
dev.off()
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:6)
agg.integrated <- FindNeighbors(agg.integrated, reduction = "pca", dims = 1:6)
agg.integrated <- FindClusters(agg.integrated, resolution = 0.5)


# Fig. S2B, right panel
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID", cols = c("#e22e2f","#2a59b8"), order = "CHIRSB")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID", cols = c("#e22e2f","#2a59b8"), order = "CHIRSB") + NoLegend()
dev.off()

# Save file
save(agg.integrated, file = "CHIR-IWP-integrated.RData")


#######################################################################
## Cell Cycle Regression Setup
# Assign cell cycle scores
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
agg.integrated <- RunPCA(agg.integrated, features = VariableFeatures(agg.integrated), ndims.print = 1:10, nfeatures.print = 10)
agg.integrated <- CellCycleScoring(agg.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))

# Fig. S2Ci, no regression
pdf()
DimPlot(agg.integrated, reduction = "pca")
dev.off()


## G2M/S (Partial) Cell Cycle Regression
agg.integrated$CC.Difference <- agg.integrated$S.Score - agg.integrated$G2M.Score
agg.integrated <- ScaleData(agg.integrated, vars.to.regress = "CC.Difference", features = rownames(agg.integrated))
agg.integrated <- RunPCA(agg.integrated, features = VariableFeatures(agg.integrated), nfeatures.print = 10)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))

# Fig. S2Ci, G2M/S regression, PCA plot
pdf()
DimPlot(agg, reduction = "pca")
dev.off()

# Finishing clustering partial regression
agg.integrated <- RunPCA(agg.integrated)
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:6)
agg.integrated <- FindNeighbors(agg.integrated, dims = 1:6)
agg.integrated <- FindClusters(agg.integrated, resolution = 0.5)
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75)
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "Phase")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID")
dev.off()
save(agg.integrated, file = "CHIR-IWP-integrated-cclite.RData")


## Full Regression; reload dataset to undo partial regression
load("CHIR-IWP-integrated.RData")
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
agg.integrated <- RunPCA(agg.integrated, features = VariableFeatures(agg.integrated), ndims.print = 1:10, nfeatures.print = 10)
agg.integrated <- CellCycleScoring(agg.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))
agg.integrated <- ScaleData(agg.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(agg.integrated))
agg.integrated <- RunPCA(agg.integrated, features = VariableFeatures(agg.integrated), nfeatures.print = 10)
agg.integrated <- RunPCA(agg.integrated, features = c(s.genes, g2m.genes))

# Fig. S2C, full regression
pdf()
DimPlot(agg.integrated, reduction = "pca")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "Phase")
dev.off()

# Finishing clustering full regression
agg.integrated <- RunPCA(agg.integrated)
agg.integrated <- RunUMAP(agg.integrated, reduction = "pca", dims = 1:5)
agg.integrated <- FindNeighbors(agg.integrated, dims = 1:6)
agg.integrated <- FindClusters(agg.integrated, resolution = 0.5)
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75)
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "Phase")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "Cell.ID")
dev.off()
save(agg.integrated, file = "CHIR-IWP-integrated-ccfull.RData")

# Fig. S2D
pdf()
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "percent.mito")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "nCount_RNA")
DimPlot(agg.integrated, reduction = "umap", pt.size = 0.75, group.by = "nFeature_RNA")
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
