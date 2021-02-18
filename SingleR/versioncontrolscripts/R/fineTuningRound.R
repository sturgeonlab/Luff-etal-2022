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
