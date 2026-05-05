# Purpose: Shared fairness-metric helpers.

#' Compute fairness metrics by group.
#'
#' @param Y Outcome or imputed outcome.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Classification threshold. Can be a scalar or a vector
#' with one value per observation.
#' @param W Optional observation weights.
#'
#' @export
get_metric <- function(Y, S, A, threshold = 0.5, W = NULL) {
  if (is.null(W)) {
    W <- rep(1, length(Y))
  }

  if (length(threshold) == 1) {
    threshold <- rep(threshold, length(S))
  } else if (length(threshold) != length(S)) {
    stop("`threshold` must have length 1 or length equal to `S`.")
  }

  class <- sort(unique(A))
  out <- c()
  for (i in class) {
    C <- 1 * (S[A == i] > threshold[A == i])
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
    bs <- mean(
      Y[A == i] * W[A == i] +
        (S[A == i])^2 * W[A == i] -
        2 * Y[A == i] * S[A == i] * W[A == i]
    ) / mean(W[A == i])
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

Indicator <- function(x) {
  ifelse(I(x), 1, 0)
}

AUC.FUN <- function(data) {
  dd <- data[, 1]
  xx <- data[, 2]
  n0 <- sum(1 - dd)
  n1 <- sum(dd)
  x0 <- xx[dd == 0]
  x1 <- xx[dd == 1]
  sum((sum.I(x0, "<=", x1) + sum.I(x0, "<", x1)) / 2) / (n0 * n1)
}

#' Efficiently count comparison relations based on ranks.
#'
#' @param yy Vector of values to compare.
#' @param FUN Comparison operator such as `"<"` or `"<="`.
#' @param Yi Reference vector.
#' @param Vi Optional matrix of values to cumulatively aggregate.
#' @param ties.method Tie handling method passed to `rank()`.
#'
#' @export
sum.I <- function(yy, FUN, Yi, Vi = NULL, ties.method = "first") {
  if (FUN == "<" | FUN == ">=") {
    yy <- -yy
    Yi <- -Yi
  }
  pos <- rank(c(yy, Yi), ties.method = ties.method)[1:length(yy)] -
    rank(yy, ties.method = ties.method)
  if (substring(FUN, 2, 2) == "=") {
    pos <- length(Yi) - pos
  }
  if (!is.null(Vi)) {
    if (substring(FUN, 2, 2) == "=") {
      tmpind <- order(-Yi)
    } else {
      tmpind <- order(Yi)
    }
    Vi <- apply(as.matrix(Vi)[tmpind, , drop = FALSE], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  }
  return(pos)
}
