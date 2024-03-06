#' Purpose: Perturbation Resampling
#' Updated: 2022-08-30
#'
#' @param Y Outcome in labeled dataset.
#' @param S Model score in labeled dataset.
#' @param A Group indicator in labeled dataset.
#' @param method Group fairness estimation method. Choose from "supervised" or "spline" or "ks".
#' @param threshold Threshold for classification based on the model score.
#' Default value is 0.5.
#' @export

resample <- function(Y,
                     S,
                     A,
                     method,
                     threshold = 0.5,
                     n_boot = 500,
                     ...) {
  if (method == "supervised") {
    lapply(1:n_boot, function(i) {
      W <- 4 * rbeta(length(Y), 1 / 2, 3 / 2)
      SupervisedFairness(
        Y,
        S,
        A,
        threshold = threshold,
        W = W
      )
    })
  } else {
    lapply(1:n_boot, function(i) {
      W <- 4 * rbeta(sum(is.na(Y)), 1 / 2, 3 / 2)
      Infairness(Y,
        S,
        A,
        W = W,
        method = method,
        ...
      )$metric
    })
  }
}
