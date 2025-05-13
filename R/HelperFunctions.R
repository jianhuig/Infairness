# Purpose: Helper functions,
# Updated: 2025-05-12

#' Parametric calibration models.
#'
#' @param Y_labeled Outcome in labeled dataset.
#' @param S_labeled Model score in labeled dataset.
#' @param A_labeled Group indicator in labeled dataset.
#' @param S_unlabeled Model score in unlabeled dataset.
#' @param method Method to use; Platt scaling or Beta calibration.
#' Default is Platt scaling.
#' @export
#'
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
    param_model <- glm(Y_labeled[A_labeled == A_val] ~
      S_labeled[A_labeled == A_val], family = binomial, weights = W)$coeff

    imp_unlabeled <- expit(cbind(1, S_unlabeled) %*%
      param_model)
  } else {
    param_model <- glm(
      Y_labeled[A_labeled == A_val] ~
        log(S_labeled[A_labeled == A_val]) +
        log(1 - S_labeled[A_labeled == A_val]),
      family = binomial, weights = W
    )$coeff

    imp_unlabeled <- expit(cbind(1, log(S_unlabeled), log(1 - S_unlabeled)) %*%
      param_model)
  }

  return(imp_unlabeled)
}

# Compute all the metrics
get_metric <- function(Y, S, A, threshold = 0.5) {
  class <- sort(unique(A))
  out <- c()
  D <- 1 * (S > threshold)
  for (i in class) {
    mu_Y <- mean(Y[A == i])
    mu_D <- mean(D[A == i])
    mu_DY <- mean(D[A == i] * Y[A == i])
    mu_SY <- mean(S[A == i] * Y[A == i])
    mu_S2 <- mean(S[A == i]**2)
    
    print(paste0("mu_Y: ", mu_Y))
    print(paste0("mu_D: ", mu_D))
    print(paste0("mu_DY: ", mu_DY))
    tpr <- mu_DY / mu_Y
    print(paste0("tpr: ", tpr))
    fpr <- (mu_D - mu_DY) / (1 - mu_Y)
    ppv <- mu_DY / mu_D
    npv <- (1 - mu_D - mu_Y + mu_DY) / (1 - mu_D)
    f1 <- 2 * mu_DY / (mu_D + mu_Y)
    acc <- 1 - mu_Y - mu_D + 2 * mu_DY
    bs <- mu_S2 - 2 * mu_SY + mu_Y

    out <- c(out, tpr, fpr, ppv, npv, f1, acc, bs)
    print(out)
    
  }
  if (length(class) == 2) {
    out <- cbind(matrix(out, ncol = 2, byrow = FALSE), NA)
    out[, 3] <- out[, 1] - out[, 2]
    colnames(out) <- c(paste0("Group", class), "Delta")
  } else {
    return(error("More than 2 groups not supported"))
  }
  rownames(out) <- c("TPR", "FPR", "PPV", "NPV", "F1", "ACC", "BS")
  tibble::rownames_to_column(as.data.frame(out), "Metric")
}



# Indicator function in R
Indicator <- function(x) {
  ifelse(I(x), 1, 0)
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}
