library(Seurat)
library(viridis)
library(cowplot)
library(ggplot2)

# load WNTd/WNTi data from GEO
CHIRSB.data <- Read10X(data.dir = "./CHIRSB/")
CHIRSB <- CreateSeuratObject(CHIRSB.data, assay = "RNA", min.cells = 3, min.features = 200, project = "CHIRSB")
IWP2.data <- Read10X(data.dir = "./IWP2/")
IWP2 <- CreateSeuratObject(IWP2.data, assay = "RNA", min.cells = 3, min.features = 200, project = "IWP2")
CHIRSB.IWP2 <- merge(CHIRSB,IWP2, project = "CHIRSB.IWP2")
CHIRSB.IWP2[["percent.mito"]] <- PercentageFeatureSet(CHIRSB.IWP2, pattern = "^MT-")
CHIRSB.IWP2[["merged"]] <- "CHIRSB.IWP2"

# Fig. S2A, left
pdf()
VlnPlot(CHIRSB.IWP2, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1, cols = c("#20a288","#20a288","#20a288"))
VlnPlot(CHIRSB.IWP2, group.by = "merged", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1, cols = c("#20a288","#20a288","#20a288"))
dev.off()

CHIRSB.IWP2 <- subset(CHIRSB.IWP2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 10)

# Fig. S2A, right
pdf()
VlnPlot(CHIRSB.IWP2, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1, cols = c("#20a288","#20a288","#20a288"))
VlnPlot(CHIRSB.IWP2, group.by = "merged", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.1, cols = c("#20a288","#20a288","#20a288"))
dev.off()

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


####################################################################

## Integrate CHIRSB & IWP2 samples
# 22732 features chosen so cell cycle genes will be included for downstream regression
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

# Save file
save(CHIRSB.IWP2.integrated, file = "CHIRSB-IWP2-after-integration.RData")
