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
