library(Seurat)
library(dplyr)
library(Matrix)

# set up Zeng dataset
#obtain matrices from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135202
cs10 <-read.table("GSM3993420_CS10_Body_rawdata.txt", header = T, row.names=1, as.is=T)
cs11 <-read.table("GSM3993421_CS11_CH_rawdata.txt", header = T, row.names=1, as.is=T)
cs13 <-read.table("GSM3993422_CS13_DA_rawdata.txt", header = T, row.names=1, as.is=T)
cs14 <-read.table("GSM3993425_CS14_DA_UMI_raw.txt", header = T, row.names=1, as.is=T)
cs15 <-read.table("GSM3993426_CS15_DA_UMI_raw.txt", header = T, row.names=1, as.is=T)
data <- CreateSeuratObject(cs10, project = "zeng", min.cells = 1)
data2 <- CreateSeuratObject(cs11, project = "zeng", min.cells = 1)
data3 <- CreateSeuratObject(cs13, project = "zeng", min.cells = 1)
data4 <- CreateSeuratObject(cs14, project = "zeng", min.cells = 1)
data5 <- CreateSeuratObject(cs15, project = "zeng", min.cells = 1)
zeng <- merge(data, y = c(data2,data3,data4,data5), add.cell.ids = c("CS10", "CS11", "CS13", "CS14", "CS15"), project = "zeng")
zeng <- NormalizeData(zeng, normalization.method = "LogNormalize", scale.factor = 10000)
zeng <- FindVariableFeatures(zeng, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(zeng)
zeng <- ScaleData(zeng, features = all.genes)

save(zeng, file="zeng.RData")


# identify HE and arterial endothelial cells for SingleR
pdf()
  plot(density(zeng@assays$RNA@data['CDH5',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['CXCR4',]))
  abline(v=0.3)
  plot(density(zeng@assays$RNA@data['GJA5',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['DLL4',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['PTPRC',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['SPN',]))
  abline(v=0.05)
  plot(density(zeng@assays$RNA@data['HEY2',]))
  abline(v=0.15)
  plot(density(zeng@assays$RNA@data['RUNX1',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['HOXA1',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['HOXA2',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['HOXA3',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['HOXA4',]))
  abline(v=0.1)
  plot(density(zeng@assays$RNA@data['HOXA5',]))
  abline(v=0.15)
  plot(density(zeng@assays$RNA@data['HOXA6',]))
  abline(v=0.05)
  plot(density(zeng@assays$RNA@data['HOXA7',]))
  abline(v=0.15)
  plot(density(zeng@assays$RNA@data['HOXA9',]))
  abline(v=0.05)
  plot(density(zeng@assays$RNA@data['HOXA10',]))
  abline(v=0.2)
  plot(density(zeng@assays$RNA@data['ITGA2B',]))
  abline(v=0.1)
  dev.off()
aec <- WhichCells(zeng, expression = CDH5 > 0.2 & CXCR4 > 0.3 & GJA5 > 0.2 & DLL4 > 0.1 & PTPRC < 0.1 & SPN < 0.05 & HEY2 > 0.15)
he1 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA1 > 0.2 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he2 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA2 > 0.1 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he3 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA3 > 0.2 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he4 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA4 > 0.1 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he5 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA5 > 0.15 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he6 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA6 > 0.05 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he7 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA7 > 0.15 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he9 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA9 > 0.05 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
he10 <- WhichCells(zeng, expression =  RUNX1 > 0.1 & CDH5 > 0.2 & HOXA10 > 0.2 & PTPRC < 0.1 & ITGA2B < 0.1 & SPN < 0.05)
#there are no HOXA11/13+ cells in the dataset
allhe <- union(he1,c(he2,he3,he4,he5,he6,he7,he9,he10))
Idents(zeng) <- "exp"
zeng <- SetIdent(zeng, cells = aec, value = "AEC")
zeng <- SetIdent(zeng, cells = allhe, value = "HE")
zeng[["newlabels"]] <- Idents(zeng)
sub <- SubsetData(zeng, subset.name = "newlabels", accept.value = c("HE", "AEC"))


##### SingleR portion, go to bottom and assign functions instead of loading package, as some things behave drastically different between versions
library(matrixStats)
library(outliers)
library(R.utils)
library(pheatmap)

# set up reference dataset
name = "HE_Endo"
bulk <- read.table("bulkseq-input.txt", header=T, sep="\t", row.names = 1)
bulk = as.matrix(bulk)
types = list("IWP","RAi","RAd","CD184")
types = as.character(types)
main_types = list("HE","HE","HE","Endo")
main_types = as.character(main_types)
myref = list(data = bulk, main_types=main_types, types=types)
myref[["de.genes"]] = CreateVariableGeneSet(bulk,types,200)
myref[["de.genes.main"]] = CreateVariableGeneSet(bulk,main_types,300)

# assign categories to cells in cs dataset for "clusters" when plotted later
stage = read.table("stage.txt", sep = "\t", header = TRUE)
stage <- as.list(stage)
labels <- sub@meta.data$newlabels

# create singlecell dataset for comparison to bulk RNAseq
data <- GetAssayData(sub, assay = "RNA", slot = "data")
data <- as.matrix(data)

# run main SingleR function & add names for cluster
obj <- SingleR.CreateObject(data, myref, clusters = NULL, variable.genes = "de", fine.tune = F)
obj[["names"]] <- stage
obj[["names2"]] <- labels

# plot results by hierarchial clustering and "names"
pdf()
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names$stage, order.by.clusters = TRUE)
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names2, order.by.clusters = TRUE)
dev.off()


# export relative correlation scores and names for each cell in cs dataset
scores = obj$SingleR.single.main$r
  datascores <- as.matrix(scores)
  mmax = rowMaxs(datascores)
  mmin = rowMins(datascores)
  datascoresrel = (datascores-mmin)/(mmax-mmin)
  datascoresrel = datascoresrel^3
  write.table(datascoresrel, file="datascores.txt", sep="\t")
  write.table(obj$names2, file="scoresnames.txt", sep="\t")
  write.table(obj$names$stage, file="scoresnames2.txt", sep="\t")


# order cells based on scores
order <- read.table("order.txt", header = TRUE, sep ="\t")
orderv <- as.vector(order$Barcode)

pdf(width = 12, height = 3)
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names2, order.by.clusters = FALSE, cells_order = orderv)
SingleR.DrawHeatmap(obj, top.n = Inf, clusters = obj$names$stage, order.by.clusters = FALSE, cells_order = orderv)
dev.off()




####### functions to be used after loading packages ######


CreateSinglerObject = function(counts,annot=NULL,project.name,
                                 min.genes=0,technology='10X',
                                 species='Human',citation='',
                                 ref.list=list(),normalize.gene.length=F,
                                 variable.genes='de',fine.tune=T,
                                 do.signatures=F,clusters=NULL,
                                 do.main.types=T,reduce.file.size=T,
                                 temp.dir=NULL,numCores = SingleR.numCores) {

    sc.data = ReadSingleCellData(counts2,annot)

    print(paste0('Dimensions of counts data: ',
                 nrow(sc.data$counts),'x',ncol(sc.data$counts)))

    singler = list()


    N = colSums(counts2>0)
    orig.ident = sc.data$orig.ident[N>=min.genes]
    counts = sc.data$counts[,N>=min.genes]

    if (normalize.gene.length == F) {
      sc.data.gl = counts
      rownames(sc.data.gl) = tolower(rownames(sc.data.gl))
    } else {
      if (species == 'Human') {
        sc.data.gl = TPM(counts,human_lengths)
      } else if (species == 'Mouse') {
        sc.data.gl = TPM(counts,mouse_lengths)
      }
    }

    if (length(ref.list)==0) {
      if (species == 'Mouse') {
        #if (!exists('immgen'))
        #  data('Immgen')
        #if (!exists('mouse.rnaseq'))
        #  data('Mouse-RNAseq')
        res = list(SingleR.CreateObject(sc.data.gl,immgen,clusters,species,
                                        citation,technology,
                                        do.main.types=do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores),
                   SingleR.CreateObject(sc.data.gl,mouse.rnaseq,clusters,
                                        species,citation,technology,
                                        do.main.types=do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores)
        )
      } else if (species == 'Human') {
        #if(!exists('hpca'))
        #  data ('HPCA')
        #if (!exists('blueprint_encode'))
        #  data('Blueprint_Encode')
        res = list(SingleR.CreateObject(sc.data.gl,hpca,clusters,species,
                                        citation,technology,
                                        do.main.types = do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores),
                   SingleR.CreateObject(sc.data.gl,blueprint_encode,
                                        clusters,species,citation,technology,
                                        do.main.types = do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores))
      }
    } else {
      res = lapply(ref.list, FUN=function(x) {
        SingleR.CreateObject(sc.data.gl,x,clusters,species,citation,technology,
                             do.main.types=do.main.types,
                             variable.genes=variable.genes,fine.tune=fine.tune,
                             numCores = numCores)
      })
    }

    singler$singler = res

    if (do.signatures==TRUE) {
      signatures = calculateSingScores(sc.data.gl,species=species)
      singler$signatures = signatures

    }

    if (species == 'Human') {
      kang = SingleR.CreateKangAnnotations(sc.data.gl)
      singler$other = kang$kang_annotation
    }

    singler$meta.data = list(project.name=project.name,orig.ident=orig.ident)

    if (reduce.file.size==T) {
      singler = remove.Unnecessary.Data.single(singler)
    }

    singler

  }
ReadSingleCellData = function(counts,annot) {
    if (typeof(counts) == 'character') {
      if (file.info(counts)$isdir==T) {
        counts = as.matrix(Read10X(counts))
      } else if (file.info(counts)$isdir==F) {
        counts <- as.matrix(read.table(counts, header=TRUE, sep="\t",
                                       row.names=1, as.is=TRUE,comment.char='!'))
      } else {
        stop('Cannot find file or directory.')
      }
    }

    A = tolower(rownames(counts))
    dupA = duplicated(A)
    if (sum(dupA)>0) {
      counts = counts[!dupA,]
    }

    #  if (species=='Mouse') {
    #    rownames(counts) = paste(toupper(substr(rownames(counts), 1, 1)), tolower(substr(rownames(counts), 2, nchar(rownames(counts)))), sep="")
    #  }

    if (!is.null(annot)) {
      if (length(annot) == 1) {
        types <- read.table(annot, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
        orig.ident = types[,1]
        names(orig.ident) = rownames(types)
      } else {
        orig.ident = annot
        if (is.null(names(orig.ident))) {
          names(orig.ident) = colnames(counts)
        }
      }
    } else {
      orig.ident = rep('NA',ncol(counts))
      names(orig.ident)=colnames(counts)
    }

    colnames(counts) = make.unique(colnames(counts))
    names(orig.ident) = colnames(counts)

    list(counts=counts,orig.ident=orig.ident)
  }
SingleR.CreateObject <- function(sc.data,ref,clusters=NULL,species='Human',
                                   citation='-',technology='-',variable.genes='de',
                                   fine.tune=T,do.main.types=T,
                                   numCores = SingleR.numCores) {
    types = ref[["types"]]

    print(paste0('Annotating data with ',ref[["name"]],'...'))

    print(paste('Variable genes method:',variable.genes))

    if (variable.genes=='de') {
      if (!is.null(ref[["de.genes"]])) {
        variable.genes = ref[["de.genes"]]
        variable.genes.main = ref[["de.genes.main"]]
      } else {
        variable.genes = CreateVariableGeneSet(ref[["data"]],ref[["types"]],200)
        variable.genes.main = CreateVariableGeneSet(ref[["data"]],ref[["main_types"]],300)
      }
    } else {
      variable.genes.main = variable.genes
    }

    SingleR.single = SingleR("single",sc.data,ref[["data"]],types=types,
                             sd.thres = ref[["sd.thres"]],genes = variable.genes,
                             fine.tune = fine.tune,numCores = numCores)

    if (is.null(clusters)) {
      SingleR.single[["clusters"]] = SingleR.Cluster(SingleR.single,10)
      clusters = SingleR.single[["clusters"]][["cl"]]
    }

    SingleR.clusters = SingleR("cluster",sc.data,ref[["data"]],types=types,
                               clusters = factor(clusters),
                               sd.thres = ref$sd.thres,
                               genes = variable.genes,
                               fine.tune = fine.tune,numCores = numCores)

    about = list(Organism = capitalize(species),Citation=citation,
                 Technology = technology,RefData=ref[["name"]])


    singler = list(SingleR.single = SingleR.single,
                   SingleR.clusters = SingleR.clusters,about=about)

    if (do.main.types==T) {
      print(paste0('Annotating data with ',ref[["name"]],' (Main types)...'))
      types = ref[["main_types"]]
      singler[["SingleR.single.main"]] = SingleR("single",sc.data,ref[["data"]],
                                            types=types,sd.thres = ref[["sd.thres"]],
                                            quantile.use = 0.8,
                                            genes = variable.genes.main,
                                            fine.tune = fine.tune,
                                            numCores = numCores)
      if (is.null(clusters)) {
        singler[["SingleR.single.main"]][["clusters"]] =
          SingleR.Cluster(singler[["SingleR.single.main"]],10)
      }
      singler[["SingleR.clusters.main"]] =
        SingleR("cluster",sc.data,ref[["data"]],types=types,
                clusters=factor(clusters),sd.thres = ref[["sd.thres"]],
                quantile.use = 0.8,genes = variable.genes.main,
                fine.tune = fine.tune,numCores = numCores)
    }


      singler[["about"]][["reference"]] = ref


    singler
  }
CreateSinglerSeuratObject = function(counts,annot=NULL,project.name,
                                     min.genes=200,technology='10X',
                                     species='Human',citation='',
                                     ref.list=list(),normalize.gene.length=F,
                                     variable.genes='de',fine.tune=T,
                                     reduce.file.size=T,do.signatures=F,
                                     min.cells=2,npca=10,regress.out='nUMI',
                                     do.main.types=T,reduce.seurat.object=T,
                                     temp.dir=NULL, numCores = SingleR.numCores) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat needed for this function to work. Please install it.",
         call. = FALSE)
  }}


  medianMatrix = function(mat,groups) {
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if(sum(A)==1) {
      mat[,A]
    } else {
      rowMedians(mat[,A],na.rm=T)
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}

CreateVariableGeneSet = function(ref_data,types,n) {
  mat = medianMatrix(ref_data,types)
  genes = lapply(1:ncol(mat), function(j) {
    lapply(1:ncol(mat), function(i) {
      s=sort(mat[,j]-mat[,i],decreasing=T)
      s=s[s>0]
      tolower(names(s)[1:min(n,length(s))])
    })
  } )
  names(genes) = colnames(mat)
  for (i in 1:length(genes)) {
    names(genes[[i]]) = colnames(mat)
  }
  genes
}
medianMatrix = function(mat,groups) {
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if(sum(A)==1) {
      mat[,A]
    } else {
      rowMedians(mat[,A],na.rm=T)
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}
SingleR <- function(method = "single", sc_data, ref_data, types,
                    clusters = NULL, genes = "de", quantile.use = 0.8,
                    p.threshold = 0.05, fine.tune = TRUE,
                    fine.tune.thres = 0.05,sd.thres=1, do.pvals = T,
                    numCores = SingleR.numCores, ...) {
  rownames(ref_data) = tolower(rownames(ref_data))
  rownames(sc_data) = tolower(rownames(sc_data))
  A = intersect(rownames(ref_data),rownames(sc_data))
  sc_data = as.matrix(sc_data[A,])
  ref_data = ref_data[A,]
  if (ncol(sc_data)>1) {
    not.use = rowSums(is.na(ref_data))>0 | rowSums(is.na(sc_data))>0 |
      rowSums(ref_data)==0
    ref_data = ref_data[!not.use,]
    sc_data = sc_data[!not.use,]
  }

  mat = medianMatrix(ref_data,types)

  if (typeof(genes)=='list') {
    utypes = unique(types)
    n = round(500*(2/3)^(log2(c(ncol(mat)))))
    genes.filtered = unique(unlist(unlist(lapply(utypes,function(j)
      lapply(utypes, function(i) genes[[i]][[j]][1:n])))))
    genes.filtered = intersect(genes.filtered,rownames(mat))
    print(paste0("Number of DE genes:", length(genes.filtered)))
  } else if (genes[1] == "de") {
    n = round(500*(2/3)^(log2(c(ncol(mat)))))
    genes.filtered = unique(unlist(unlist(lapply(1:ncol(mat), function(j) {
      lapply(1:ncol(mat), function(i) {
        s=sort(mat[,j]-mat[,i],decreasing=T);
        s=s[s>0];
        names(s)[1:min(n,length(s))]
      })}))))[-1]
    print(paste0("Number of DE genes:", length(genes.filtered)))
  } else if (genes[1] == "sd") {
    sd =  rowSds(as.matrix(mat))
    genes.filtered=intersect(rownames(mat)[sd>sd.thres],rownames(sc_data))
    print(paste0("Number of genes with SD>",sd.thres,": ",
                 length(genes.filtered)))
  } else {
    genes.filtered=intersect(tolower(genes),intersect(tolower(rownames(sc_data)),
                                                      tolower(rownames(ref_data))))
    print(paste("Number of genes using in analysis:",length(genes.filtered)))

  }

  cell.names = colnames(sc_data)

  if (method == "single") {
    if (dim(sc_data)[2]>2) {
      print(paste("Number of cells:",dim(sc_data)[2]))
    }
  } else if (method == "cluster") {
    n = length(levels(clusters))
    print(paste("Number of clusters:",n))
    data = matrix(nrow=dim(sc_data)[1],ncol=n)
    for (i in 1:n) {
      data[,i] = rowSums(as.matrix(sc_data[,is.element(clusters,
                                                       levels(clusters)[i])]))
    }
    colnames(data) = levels(clusters)
    rownames(data) = rownames(sc_data)
    sc_data = data
  } else {
    print ("Error: method must be 'single' or 'cluster'")
    return(0)
  }
  output = SingleR.ScoreData(sc_data,ref_data,genes.filtered,types,quantile.use,numCores=numCores,...)
  if (do.pvals == T) {
    output$pval = SingleR.ConfidenceTest(output$scores)
  }

  # second round with top labels
  if (fine.tune==TRUE & length(unique(types)) > 2) {
    labels = SingleR.FineTune(sc_data,ref_data,types,output$scores,
                              quantile.use,fine.tune.thres,genes = genes,
                              sd.thres,mat,numCores = numCores)
    output$labels1 = as.matrix(output$labels)
    output$labels = as.matrix(labels)
    output$labels1.thres = c(output$labels)
    if (do.pvals == T)
      output$labels1.thres[output$pval>p.threshold] = "X"
  } else {
    labels = as.matrix(output$labels)
    output$labels.thres = c(output$labels)
    if (do.pvals == T)
      output$labels.thres[output$pval>p.threshold] = "X"
  }
  output$cell.names = cell.names
  output$quantile.use = quantile.use
  output$types = types
  output$method = method

  return (output)
}
SingleR.ScoreData <- function(sc_data,ref_data,genes,types,quantile.use,numCores=1,step=10000) {
  sc_data = as.matrix(sc_data[genes,])
  ref_data = as.matrix(ref_data[genes,])

  if (ncol(sc_data)>step) {
    n = ncol(sc_data)
    s = seq(step+1,n,by=step)
    if (FALSE) {
      cl <- makeCluster(numCores)
      doParallel::registerDoParallel(cl)
      tmpr = foreach (i = 0:length(s)) %dopar% {
        if(i == 0){
          res = data.table::data.table(cor.stable(sc_data[,1:step],ref_data,method='spearman'))
        } else {
          A = seq(s[i],min(s[i]+step-1,n))
          # r=rbind(r,cor(sc_data[,A],ref_data,method='spearman'))
          res = data.table::data.table(cor.stable(sc_data[,A],ref_data,method='spearman'))
        }
        res
      }
      r = data.table::rbindlist(tmpr, use.names = F)
      r = as.matrix(r)
      rownames(r) = colnames(sc_data)
      on.exit(stopCluster(cl))
    } else {
      s = seq(step+1,n,by=step)
      r=cor.stable(sc_data[,1:step],ref_data,method='spearman')
      for (i in 1:length(s)) {
        A = seq(s[i],min(s[i]+step-1,n))
        r = rbind(r,cor.stable(sc_data[,A],ref_data,method='spearman'))
      }
    }

  } else {
    r=cor.stable(sc_data,ref_data,method='spearman')
  }
  agg_scores = quantileMatrix(r,types,quantile.use);
  #agg_scores = aggregate(t(r)~types,FUN = quantile, probs  = quantile.use)
  labels = colnames(agg_scores)[max.col(agg_scores)]
  output = list()


  if (dim(sc_data)[2]>1) {
    names(labels)=t(colnames(sc_data))
  }

  output$scores = as.matrix(t(agg_scores))

  output$labels = as.matrix(labels)
  output$r = r
  output$scores = t(output$scores)

  return(output)
}
fineTuningRound = function(topLabels,types,ref_data,genes,mat,sd.thres,
                           sc_data,quantile.use,fine.tune.thres) {
  labels.use = is.element(types,topLabels)
  ref_data.filtered = as.matrix(ref_data[,labels.use])
  types.filtered = types[labels.use]
  if (typeof(genes)=='list') {
    n = round(1000*(2/3)^(log2(c(ncol(mat)))))
    utypes = colnames(mat)
    genes.filtered = unique(unlist(unlist(lapply(utypes,function(j)
      lapply(utypes, function(i) genes[[i]][[j]][1:n])))))
    genes.filtered = intersect(genes.filtered,rownames(mat))
  } else if (genes[1] == "de") {
    n = round(500*(2/3)^(log2(c(ncol(mat)))))
    genes.filtered = unique(unlist(unlist(
      lapply(1:ncol(mat), function(j) {
        lapply(1:ncol(mat), function(i) {
          s=sort(mat[,j]-mat[,i],decreasing=T);
          s=s[s>0];names(s)[1:min(n,length(s))]})
      }))))[-1]
  } else if (genes[1] == "sd") {
    sd =  rowSds(mat)
    thres = min(sort(sd,decreasing = TRUE)[500],sd.thres)
    genes.filtered = intersect(rownames(ref_data)[sd>=thres],names(sc_data))
  } else {
    genes.filtered=intersect(genes,intersect(rownames(sc_data),
                                             (rownames(ref_data))))
  }

  ref_data.filtered = ref_data.filtered[genes.filtered,]
  sc_data.filtered = as.matrix(sc_data[genes.filtered])
  if (length(genes.filtered)<20) {
    return (topLabels[1])
  }
  if (sd(sc_data.filtered)>0) {
    r=cor.stable(sc_data.filtered,ref_data.filtered,method='spearman')
    agg_scores = quantileMatrix(r,types.filtered,quantile.use)[1,];
    agg_scores = sort(agg_scores,decreasing = T)
    agg_scores = agg_scores[-length(agg_scores)]
    topLabels = names(agg_scores)[agg_scores>=agg_scores[1]-fine.tune.thres]
  } else {
    topLabels = topLabels[1]
  }
  topLabels
}
cor.stable <- function (x, y, method="pearson", ...) {
  omit1 <- which(apply(x, 2, sd) == 0)
  omit2 <- which(apply(y, 2, sd) == 0)
  if (length(omit1) > 0 && length(omit2) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[-omit1,-omit2] = cor(x[,-omit1], y[,-omit2], method=method, ...)
  } else if (length(omit1) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[-omit1,] = cor(x[,-omit1], y, method=method, ...)
  } else if (length(omit2) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[,-omit2] = cor(x, y[,-omit2], method=method, ...)
  } else {
    r = cor(x, y, method=method, ...)
  }
}
quantileMatrix = function(mat,groups,q) {
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if (nrow(mat)==1) {
      quantile(mat[A],na.rm=T,probs=q)
    } else {
      if(sum(A)==1) {
        mat[,A]
      } else {
        rowQuantiles(mat[,A],na.rm=T,probs=q)
      }
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}
SingleR.ConfidenceTest = function(scores) {
  apply(scores, 1,function(x) chisq.out.test(x)$p.value)
}
SingleR.Cluster = function(SingleR,num.clusts=10,normalize_rows=F,
                           normalize_cols=T) {
  if (normalize_rows==T) {
    SingleR$scores = scale(SingleR$scores)
  }
  if (normalize_cols==T) {
    scores = t(scale(t(SingleR$scores^3)))
  } else {
    scores = SingleR$scores
  }
  #r <- cor(t(scores), method="pearson")
  #d <- as.dist(1-r)
  #hc = hclust(d,method='ward.D2')
  hc = hclust(dist(scores,method='euclidean'),method='ward.D2')

  cl = cutree(hc,k=num.clusts)
  list(hc=hc,cl=factor(cl))
}

SingleR.DrawHeatmap = function(SingleR,cells.use = NULL, types.use = NULL,
                               clusters=NULL,top.n=40,normalize=T,
                               order.by.clusters=F,cells_order=NULL,silent=F,
                               fontsize_row=9,...) {
  #scores = SingleR$scores
  #scores = SingleR$SingleR.single.main$scores
  scores = SingleR$SingleR.single.main$r
  if (!is.null(cells.use)) {
    scores = scores[cells.use,]
  }
  if (!is.null(types.use)) {
    scores = scores[,types.use]
  }

  m = apply(t(scale(t(scores))),2,max)

  thres = sort(m,decreasing=TRUE)[min(top.n,length(m))]

  data = as.matrix(scores)

  if (normalize==T) {
    mmax = rowMaxs(data)
    mmin = rowMins(data)
    data = (data-mmin)/(mmax-mmin)
    data = data^3
  }
  data = data[,m>(thres-1e-6)]


  data = t(data)

  if (!is.null(clusters)) {
    clusters = as.data.frame(clusters)
    colnames(clusters) = 'Clusters'
    rownames(clusters) = colnames(data)

  }
  additional_params = list(...)
  if (is.null(additional_params$annotation_colors)) {
    annotation_colors = NA
  } else {
    annotation_colors = additional_params$annotation_colors
  }
  clustering_method = 'ward.D2'
  if (order.by.clusters==T) {
    data = data[,order(clusters$Clusters)]
    clusters = clusters[order(clusters$Clusters),,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,
             clustering_method=clustering_method,fontsize_row=fontsize_row,
             annotation_col = clusters,cluster_cols = F,silent=silent,
             annotation_colors=annotation_colors)
  } else if (!is.null(cells_order)) {
    data = data[,cells_order]
    clusters = clusters[cells_order,,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,
             clustering_method=clustering_method,fontsize_row=fontsize_row,
             annotation_col = clusters,cluster_cols = F,silent=silent,
             annotation_colors=annotation_colors)
  } else {
    if (!is.null(clusters)) {
      pheatmap(data,border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=fontsize_row,
               annotation_col = clusters,silent=silent,
               annotation_colors=annotation_colors)
    } else {
      pheatmap(data[,sample(ncol(data))],border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=fontsize_row,
               silent=silent, annotation_colors=annotation_colors)

    }
  }
}


#### if you want no legend
SingleR.DrawHeatmap = function(SingleR,cells.use = NULL, types.use = NULL,
                               clusters=NULL,top.n=40,normalize=T,
                               order.by.clusters=F,cells_order=NULL,silent=F,
                               fontsize_row=9,...) {
  #scores = SingleR$scores
  #scores = SingleR$SingleR.single.main$scores
  scores = SingleR$SingleR.single.main$r
  if (!is.null(cells.use)) {
    scores = scores[cells.use,]
  }
  if (!is.null(types.use)) {
    scores = scores[,types.use]
  }

  m = apply(t(scale(t(scores))),2,max)

  thres = sort(m,decreasing=TRUE)[min(top.n,length(m))]

  data = as.matrix(scores)

  if (normalize==T) {
    mmax = rowMaxs(data)
    mmin = rowMins(data)
    data = (data-mmin)/(mmax-mmin)
    data = data^3
  }
  data = data[,m>(thres-1e-6)]


  data = t(data)

  if (!is.null(clusters)) {
    clusters = as.data.frame(clusters)
    colnames(clusters) = 'Clusters'
    rownames(clusters) = colnames(data)

  }
  additional_params = list(...)
  if (is.null(additional_params$annotation_colors)) {
    annotation_colors = NA
  } else {
    annotation_colors = additional_params$annotation_colors
  }
  clustering_method = 'ward.D2'
  if (order.by.clusters==T) {
    data = data[,order(clusters$Clusters)]
    clusters = clusters[order(clusters$Clusters),,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,
             clustering_method=clustering_method,fontsize_row=fontsize_row,
             annotation_col = clusters,cluster_cols = F,silent=silent,
             annotation_colors=annotation_colors, annotation_legend = FALSE)
  } else if (!is.null(cells_order)) {
    data = data[,cells_order]
    clusters = clusters[cells_order,,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,
             clustering_method=clustering_method,fontsize_row=fontsize_row,
             annotation_col = clusters,cluster_cols = F,silent=silent,
             annotation_colors=annotation_colors, annotation_legend = FALSE)
  } else {
    if (!is.null(clusters)) {
      pheatmap(data,border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=fontsize_row,
               annotation_col = clusters,silent=silent,
               annotation_colors=annotation_colors, annotation_legend = FALSE)
    } else {
      pheatmap(data[,sample(ncol(data))],border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=fontsize_row,
               silent=silent, annotation_colors=annotation_colors, annotation_legend = FALSE)

    }
  }
}

#### extra, not used
library(rlang)
DoMultiBarHeatmap <- function (object,
                               features = NULL,
                               cells = NULL,
                               group.by = "ident",
                               additional.group.by = NULL,
                               group.bar = TRUE,
                               disp.min = -2.5,
                               disp.max = NULL,
                               slot = "scale.data",
                               assay = NULL,
                               label = TRUE,
                               size = 5.5,
                               hjust = 0,
                               angle = 45,
                               raster = TRUE,
                               draw.lines = TRUE,
                               lines.width = NULL,
                               group.bar.height = 0.02,
                               combine = TRUE)
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot,
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                                                             collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object,
                                                             slot = slot)[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = StashIdent(object = object,
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by)]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    group.use <- groups.use[, c(i, additional.group.by), drop = FALSE]

    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }
    }

    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) *
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(foo=rep(x = levels(x = group.use[[i]]), times = lines.width))
      placeholder.groups[additional.group.by] = NA
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells

      group.levels <- levels(x = group.use[[i]])

      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells

      group.use <- rbind(group.use, placeholder.groups)

      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }

    #group.use = group.use[order(group.use[[i]]), , drop=F]
    group.use <- group.use[with(group.use, eval(parse(text=paste('order(', paste(c(i, additional.group.by), collapse=', '), ')', sep='')))), , drop=F]

    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster,
                            disp.min = disp.min, disp.max = disp.max, feature.order = features,
                            cell.order = rownames(x = group.use), group.by = group.use[[i]])

    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))){
        if (colname == group.by){
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }

        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))), "#FFFFFF")
        } else {
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])


        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)

        plot <- suppressMessages(plot +
                                annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) +
                                annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = grid::gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                coord_cartesian(ylim = c(0, y.max), clip = "off"))

        #temp <- as.data.frame(cols[[colname]][levels(group.use[[colname]])])
        #colnames(temp) <- 'color'
        #temp$x <- temp$y <- 1
        #temp[['name']] <- as.factor(rownames(temp))

        #temp <- ggplot(temp, aes(x=x, y=y, fill=name)) + geom_point(shape=21, size=5) + labs(fill=colname) + theme(legend.position = "bottom")
        #legend <- get_legend(temp)
        #multiplot(plot, legend, heights=3,1)

        if ((colname == group.by) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos),
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity",
                                   data = label.x.pos, aes_string(label = "group",
                                                                  x = "label.x.pos"), y = y.max + y.max *
                                     0.03 * 0.5, angle = angle, hjust = hjust,
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) *
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}
