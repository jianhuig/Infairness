#' Nonparametric local constant regression
#' @param S_t phenotyping score S in the labeled dataset
#' @param Y_t outcome Y in the labeled dataset
#' @param S_v phenotyping score S in the unlabeled dataset
#' @param bw bandwidth for kernel smoothing
#' @param Wt optional weights
#' @param kern.mat optional kernel matrix
#' @return vector of predicted value of S_v
#' @noRd

npreg <- function(St, Yt, Sv, bw, Wt = NULL, kern.mat = NULL) {
  nv <- length(Sv)
  nt <- length(St)
  if (is.null(Wt)) {
    Wt <- rep(1, nt)
  }
  if (is.null(kern.mat)) {
    kern.mat <- sapply(1:nt, function(kk) dnorm(Sv - rep(St[kk], nv), sd = bw))
  }
  nw.est <- kern.mat %*% (Yt * Wt) * (1 / (kern.mat %*% Wt))
  return(nw.est)
}
