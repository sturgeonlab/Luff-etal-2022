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
