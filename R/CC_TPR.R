#' Control-corrected TPR estimation.
#'
#' @param Y Outcome, with missing values for unlabeled observations.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Threshold for classification based on the model score.
#' @param X Optional covariates used through spline-by-X interactions.
#' @param nknots Number of knots for the natural spline basis.
#' @param D Logical; include `I(S > threshold)` in the imputation basis.
#'
#' @return A list with TPR-only `est`, `var`, and `cv_weights` tables.
#' @export
CC_TPR <- function(Y,
                   S,
                   A,
                   threshold = 0.5,
                   X = NULL,
                   nknots = 3,
                   D = FALSE) {
  labeled_ind <- which(!is.na(Y))
  Y_labeled <- Y[labeled_ind]
  S_labeled <- S[labeled_ind]
  A_labeled <- A[labeled_ind]
  threshold_all <- rep(threshold, length.out = length(S))
  threshold_labeled <- threshold_all[labeled_ind]
  X_all <- if (is.null(X)) NULL else as.matrix(X)
  X_labeled <- if (is.null(X)) NULL else X_all[labeled_ind, , drop = FALSE]

  basis_labeled <- NaturalSplineBasis(as.matrix(S_labeled), nknots, return_knots = TRUE)
  spline_knots <- attr(basis_labeled, "knots")
  basis_all <- NaturalSplineBasis(as.matrix(S), nknots, knots = spline_knots)

  class <- sort(unique(A))
  est <- data.frame(Metric = "TPR")
  var <- data.frame(Metric = "TPR")
  cv_weights <- data.frame(Metric = "TPR")

  for (a in class) {
    idx_labeled <- A_labeled == a
    idx_all <- A == a

    X_fit <- basis_labeled[idx_labeled, , drop = FALSE]
    X_pred <- basis_all[idx_all, , drop = FALSE]

    if (D) {
      X_fit <- cbind(X_fit, as.numeric(S_labeled[idx_labeled] > threshold_labeled[idx_labeled]))
      X_pred <- cbind(X_pred, as.numeric(S[idx_all] > threshold_all[idx_all]))
    }

    if (!is.null(X_all)) {
      X_fit_interaction <- do.call(
        cbind,
        lapply(seq_len(ncol(X_labeled)), function(j) {
          basis_labeled[idx_labeled, , drop = FALSE] * X_labeled[idx_labeled, j]
        })
      )
      X_pred_interaction <- do.call(
        cbind,
        lapply(seq_len(ncol(X_all)), function(j) {
          basis_all[idx_all, , drop = FALSE] * X_all[idx_all, j]
        })
      )
      X_fit <- cbind(X_fit, X_labeled[idx_labeled, , drop = FALSE], X_fit_interaction)
      X_pred <- cbind(X_pred, X_all[idx_all, , drop = FALSE], X_pred_interaction)
    }

    fit <- fit_spline_model(X = X_fit, y = Y_labeled[idx_labeled])
    m_labeled <- as.numeric(predict_spline_model(fit, X_fit))
    m_all <- as.numeric(predict_spline_model(fit, X_pred))

    y_labeled <- Y_labeled[idx_labeled]
    A_labeled_bin <- as.numeric(S_labeled[idx_labeled] > threshold_labeled[idx_labeled])
    A_all_bin <- as.numeric(S[idx_all] > threshold_all[idx_all])

    eta <- mean(y_labeled)
    mu <- mean(m_all)
    tpr_y <- mean(A_labeled_bin * y_labeled) / eta
    tpr_m_labeled <- mean(A_labeled_bin * m_labeled) / mean(m_labeled)
    tpr_m_all <- mean(A_all_bin * m_all) / mu

    phi_y_for_w <- (A_labeled_bin - tpr_y) * y_labeled / eta
    phi_m <- (A_labeled_bin - tpr_m_all) * m_labeled / mu
    w <- stats::cov(phi_y_for_w, phi_m) / stats::var(phi_m)

    tpr <- tpr_y - w * (tpr_m_labeled - tpr_m_all)
    phi_y <- (A_labeled_bin - tpr) * y_labeled / eta
    phi <- phi_y - w * phi_m

    group_col <- paste0("Group", a)
    est[[group_col]] <- tpr
    var[[group_col]] <- stats::var(phi) / length(y_labeled)
    cv_weights[[group_col]] <- w
  }

  if (length(class) == 2) {
    group_cols <- paste0("Group", class)
    est$Delta <- est[[group_cols[1]]] - est[[group_cols[2]]]
    var$Delta <- var[[group_cols[1]]] + var[[group_cols[2]]]
  }

  list(est = est, var = var, cv_weights = cv_weights)
}
