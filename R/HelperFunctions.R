# Purpose: Helper functions,
# Updated: 2022-08-31

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
    # cat("leng1",length(Y_labeled[A_labeled == A_val]))
    # cat("leng2",length(log(S_labeled[A_labeled == A_val])))
    # cat("leng3",length(log(1 - S_labeled[A_labeled == A_val]))

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

get_metric <- function(Y, S, A, threshold = 0.5, W = NULL) {
  if (is.null(W)) {
    W <- rep(1, length(Y))
  }
  class <- sort(unique(A))
  out <- c()
  for (i in class) {
    C <- 1 * (S[A == i] > threshold)
    P <- sum(Y[A == i] * W[A == i])
    N <- sum((1 - Y[A == i]) * W[A == i])
    TP <- sum(Y[A == i] * C * W[A == i])
    FP <- sum((1 - Y[A == i]) * C * W[A == i])
    TN <- N - FP
    FN <- P - TP
    PP <- TP + FP
    PN <- FN + TN

    tpr <- TP / P
    tnr <- TN / (TN + FP)
    fpr <- FP / N
    fnr <- FN / (TP + FN)
    npv <- TN / PN
    ppv <- TP / PP
    acc <- (TP + TN) / (P + N)
    f1 <- (2 * TP) / (2 * TP + FP + FN)
    bs <- mean(Y[A == i] * W[A == i] +
      (S[A == i])**2 * W[A == i] -
      2 * Y[A == i] * S[A == i] * W[A == i]) / mean(W[A == i])
    out <- c(out, tpr, tnr, fpr, fnr, npv, ppv, acc, f1, bs)
  }
  if (length(class) == 2) {
    out <- cbind(matrix(out, ncol = 2, byrow = FALSE), NA)
    out[, 3] <- out[, 1] - out[, 2]
    colnames(out) <- c(paste0("Group", class), "Delta")
  } else {
    out <- matrix(out, ncol = length(class), byrow = FALSE)
    mad <- as.numeric(apply(out, 1, mad))
    var <- as.numeric(apply(out, 1, var))
    gei <- as.numeric(apply(out, 1, entropy))
    out <- cbind(out, mad, var, gei)
    colnames(out) <- c(paste0("Group", class), "mad", "var", "gei")
  }
  rownames(out) <- c("TPR", "TNR", "FPR", "FNR", "NPV", "PPV", "ACC", "F1", "BS")
  tibble::rownames_to_column(as.data.frame(out), "Metric")
}

# Indicator function in R
Indicator <- function(x) {
  ifelse(I(x), 1, 0)
}


# computes the AUC
AUC.FUN <- function(data) {
  dd <- data[, 1]
  xx <- data[, 2]
  n0 <- sum(1 - dd)
  n1 <- sum(dd)
  x0 <- xx[dd == 0]
  x1 <- xx[dd == 1]
  sum((sum.I(x0, "<=", x1) + sum.I(x0, "<", x1)) / 2) / (n0 * n1)
}


# Computes sums efficiently based on ranks
sum.I <- function(yy, FUN, Yi, Vi = NULL, ties.method = "first") {
  if (FUN == "<" | FUN == ">=") {
    yy <- -yy
    Yi <- -Yi
  }
  pos <- rank(c(yy, Yi), ties.method = ties.method)[1:length(yy)] - rank(yy, ties.method = ties.method)
  if (substring(FUN, 2, 2) == "=") pos <- length(Yi) - pos
  if (!is.null(Vi)) {
    if (substring(FUN, 2, 2) == "=") tmpind <- order(-Yi) else tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind, , drop = F], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  } else {
    return(pos)
  }
}

# Computes the natural spline basis
ns.basis <- function(X, nk) {
  ## X: covariate matrix
  ## nk: number of knots

  X <- as.matrix(X)
  basis <- NULL

  for (i in 1:ncol(X)) {
    X_i <- X[, i]

    # quantiles to determine appropriate knots
    knots <- quantile(X_i, seq(0, 1, length = nk))

    # changes quantiles if there aren't enough unique values
    j <- 0
    while (length(unique(knots)) != nk) {
      j <- j + 1
      knots <- unique(quantile(X_i, seq(0, 1, length = nk + j)))
    }

    # compute the natural spline basis
    d_k <- (trunc.cub(X_i, knots[nk - 1])) / (knots[nk] - knots[nk - 1])
    evals <- sapply(1:(nk - 2), function(ii) {
      d_i <- (trunc.cub(X_i, knots[ii])) / (knots[nk] - knots[ii])
      basis.new <- d_i - d_k
    })

    # bind together and include original variable
    basis <- cbind(basis, cbind(X_i, evals))

    # cbind(X_i, evals))
  }
  return(basis)
}

############################
### Natural Spline Basis ###
############################

# Computes the truncated cubic
trunc.cub <- function(X, x) {
  ## X: variable of interest
  ## x: knot location

  ((X > x) * (X - x))^3
}
