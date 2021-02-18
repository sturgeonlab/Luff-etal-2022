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
