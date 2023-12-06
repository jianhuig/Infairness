# Purpose: Estimation of fairness metrics in supervised setting.
# Updated: 2023-07-17

#' Supervised fairness estimation.
#'
#' @param Y Outcome in the observed dataset.
#' @param S Model score in observed dataset.
#' @param A Group indicator in observed  dataset.
#' @param W Weight. Default is Null.
#' @param threshold Threshold for classification based on the model score.
#' Default value is 0.5.
#' @export

SupervisedFairness <- function(Y,
                               S,
                               A,
                               threshold = 0.5,
                               W = NULL) {
  results <- get_metric(Y, S, A, threshold, W = W)

  return(results)
}
