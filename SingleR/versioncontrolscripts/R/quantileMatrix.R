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
