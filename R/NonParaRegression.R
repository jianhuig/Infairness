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

kernel_bandwidth <- function(S_values, order = 0.45) {
  sqrt(stats::var(S_values)) / (length(S_values)^order)
}

kernel_metric_scale <- function(S, A, basis, threshold, use_ecdf = FALSE) {
  nclass <- sort(unique(A))
  if (length(basis) == 1) {
    basis <- rep(basis, length(nclass))
  }

  S_metric <- S
  if (length(threshold) == 1) {
    threshold_metric <- rep(threshold, length(S))
  } else if (length(threshold) == length(S)) {
    threshold_metric <- threshold
  } else {
    stop("`threshold` must have length 1 or length equal to `S`.")
  }

  if (!use_ecdf) {
    return(list(S = S_metric, threshold = threshold_metric))
  }

  for (a in nclass) {
    if (basis[which(nclass == a)] != "kernel") {
      next
    }

    group_idx <- A == a
    ecdf_a <- stats::ecdf(S[group_idx])
    S_metric[group_idx] <- as.numeric(ecdf_a(S[group_idx]))
    threshold_metric[group_idx] <- as.numeric(ecdf_a(threshold_metric[group_idx]))
  }

  list(S = S_metric, threshold = threshold_metric)
}
