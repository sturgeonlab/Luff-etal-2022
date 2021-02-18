SingleR.FineTune <- function(sc_data,ref_data,types,scores,quantile.use,
                             fine.tune.thres,genes,sd.thres,mean_mat,
                             numCores = 1) {
  N = dim(sc_data)[2]
  print(paste("Fine-tuning round on top cell types (using", numCores,
              "CPU cores):"))
  labels = pbmclapply(1:N,FUN=function(i){
    max_score = max(scores[i,])
    topLabels = names(scores[i,scores[i,]>=max_score-fine.tune.thres])
    if (length(topLabels)==0) {
      return (names(which.max(scores[i,])))
    } else {
      k=1
      while(length(topLabels)>1) {
        topLabels = fineTuningRound(topLabels,types,ref_data,genes,
                                    mean_mat[,topLabels],sd.thres,
                                    sc_data[,i],quantile.use,fine.tune.thres)
        k=k+1
      }
      return (topLabels)
    }
  },mc.cores=numCores)
  labels = as.matrix(unlist(labels))

  if (dim(sc_data)[2]>1) {
    rownames(labels)=t(colnames(sc_data))
  }
  return(labels)

}
