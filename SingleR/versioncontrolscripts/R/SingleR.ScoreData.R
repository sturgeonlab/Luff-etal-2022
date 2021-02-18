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
