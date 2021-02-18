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
