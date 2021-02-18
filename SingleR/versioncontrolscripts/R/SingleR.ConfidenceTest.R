SingleR.ConfidenceTest = function(scores) {
  apply(scores, 1,function(x) chisq.out.test(x)$p.value)
}
