# Purpose: Calibration helpers used by semi-supervised estimators.

#' Parametric calibration models.
#'
#' @param Y_labeled Outcome in labeled dataset.
#' @param S_labeled Model score in labeled dataset.
#' @param A_labeled Group indicator in labeled dataset.
#' @param S_unlabeled Model score in unlabeled dataset.
#' @param method Method to use; Platt scaling or Beta calibration.
#' Default is Platt scaling.
#' @export
FitParametricCalibration <- function(Y_labeled,
                                     S_labeled,
                                     A_labeled,
                                     A_val,
                                     S_unlabeled,
                                     method = "Platt",
                                     W = NULL) {
  S_labeled[S_labeled == 1] <- S_labeled[S_labeled == 1] - 1e-6
  S_labeled[S_labeled == 0] <- S_labeled[S_labeled == 0] + 1e-6

  S_unlabeled[S_unlabeled == 1] <- S_unlabeled[S_unlabeled == 1] - 1e-6
  S_unlabeled[S_unlabeled == 0] <- S_unlabeled[S_unlabeled == 0] + 1e-6

  if (method == "Platt") {
    param_model <- glm(
      Y_labeled[A_labeled == A_val] ~ S_labeled[A_labeled == A_val],
      family = binomial,
      weights = W
    )$coeff

    imp_unlabeled <- expit(cbind(1, S_unlabeled) %*% param_model)
  } else {
    param_model <- betacal::beta_calibration(
      p = S_labeled[A_labeled == A_val],
      y = Y_labeled[A_labeled == A_val]
    )
    imp_unlabeled <- betacal::beta_predict(
      p = S_unlabeled,
      calib = param_model
    )
  }

  return(imp_unlabeled)
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}
