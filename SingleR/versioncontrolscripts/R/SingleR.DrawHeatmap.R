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
