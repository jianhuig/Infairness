# Purpose: Fairness Auditing of a Machine Learning Model

#' Wrapper function for fairness auditing.
#'
#' @param Y Outcome variable, may contain missing values.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Threshold for classification based on the model score.
#' Default is 0.5.
#' @param method Fairness estimation method. Options are "supervised" or
#' "semi-supervised". Default is "semi-supervised".
#' @param X Optional covariates matrix for adjustment in the semi-supervised
#' setting. Default is NULL.
#' @param basis Basis expansion method for score augmentation. Options include
#' "Poly(S)", "Poly(S) + X", "Spline(S)", "Spline(S) + X",
#' "Spline Interaction", "Interaction", "Beta", and "kernel". When `basis` is
#' `NULL`, `Audit_Fairness()` defaults to "Poly(S)" if `X` is `NULL` and to
#' "Poly(S) + X" otherwise.
#' @param ... Additional parameters passed to `SSFairness()`.
#' @return A list containing the fairness audit results.
#'
#' @export
#'

Audit_Fairness <- function(Y,
                           S,
                           A,
                           threshold = 0.5,
                           method = "semi-supervised",
                           X = NULL,
                           basis = NULL,
                           k = 10,
                           ...) {
  if (method == "supervised") {
    # Get labeled data
    labeled_ind <- which(!is.na(Y))
    Y_labeled <- Y[labeled_ind]
    S_labeled <- S[labeled_ind]
    A_labeled <- A[labeled_ind]
    return(SupervisedFairness(Y_labeled, S_labeled, A_labeled, threshold, ...))
  } else if (method == "semi-supervised") {
    if (is.null(basis)) {
      basis <- if (is.null(X)) "Poly(S)" else "Poly(S) + X"
    }
    return(SSFairness(Y, S, A, threshold, X, basis, k = k, ...))
  }
}
