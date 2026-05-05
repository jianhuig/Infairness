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
  S_labeled <- clamp_probabilities(S_labeled)
  S_unlabeled <- clamp_probabilities(S_unlabeled)

  if (method == "Platt") {
    param_model <- glm(
      Y_labeled[A_labeled == A_val] ~ S_labeled[A_labeled == A_val],
      family = binomial,
      weights = W
    )$coeff

    imp_unlabeled <- expit(cbind(1, S_unlabeled) %*% param_model)
  } else {
    S_fit <- S_labeled[A_labeled == A_val]
    Y_fit <- Y_labeled[A_labeled == A_val]
    W_fit <- if (is.null(W)) NULL else W[A_labeled == A_val]

    param_model <- glm(
      Y_fit ~ log(S_fit) + log(1 - S_fit),
      family = binomial,
      weights = W_fit
    )$coefficients

    imp_unlabeled <- expit(
      cbind(1, log(S_unlabeled), log(1 - S_unlabeled)) %*% param_model
    )
  }

  return(imp_unlabeled)
}

beta_design_matrix <- function(S) {
  S_safe <- clamp_probabilities(S)
  cbind(1, log(S_safe), log(1 - S_safe))
}

safe_inverse <- function(mat) {
  tryCatch(solve(mat), error = function(e) MASS::ginv(mat))
}

beta_if_variance_update <- function(var_est,
                                    est,
                                    Y_labeled,
                                    S_labeled,
                                    A_labeled,
                                    S_unlabeled,
                                    A_unlabeled,
                                    threshold_unlabeled,
                                    basis) {
  nclass <- sort(unique(c(A_labeled, A_unlabeled)))
  if (length(basis) == 1) {
    basis <- rep(basis, length(nclass))
  }

  if (length(threshold_unlabeled) == 1) {
    threshold_unlabeled <- rep(threshold_unlabeled, length(S_unlabeled))
  }

  out <- var_est

  for (a in nclass) {
    if (basis[which(nclass == a)] != "Beta") {
      next
    }

    idx_lab <- A_labeled == a
    idx_unlab <- A_unlabeled == a
    n_lab <- sum(idx_lab)
    n_unlab <- sum(idx_unlab)

    if (n_lab < 2 || n_unlab < 1) {
      next
    }

    Y_a <- Y_labeled[idx_lab]
    S_lab_a <- S_labeled[idx_lab]
    S_unlab_a <- S_unlabeled[idx_unlab]
    D_unlab_a <- as.numeric(S_unlab_a > threshold_unlabeled[idx_unlab])

    W_n <- beta_design_matrix(S_lab_a)
    W_N <- beta_design_matrix(S_unlab_a)

    fit <- tryCatch(
      {
        stats::glm.fit(
          x = W_n,
          y = Y_a,
          family = stats::binomial()
        )
      },
      error = function(e) NULL
    )

    if (is.null(fit)) {
      next
    }

    beta_hat <- fit$coefficients
    beta_hat[is.na(beta_hat)] <- 0

    m_n <- drop(boot::inv.logit(W_n %*% beta_hat))
    m_N <- drop(boot::inv.logit(W_N %*% beta_hat))
    w_N <- m_N * (1 - m_N)

    mu_hat <- mean(m_N)
    nu_hat <- mean(1 - m_N)
    p_hat <- mean(D_unlab_a)

    if (!is.finite(mu_hat) || !is.finite(nu_hat) ||
      !is.finite(p_hat) || mu_hat <= 0 || nu_hat <= 0 ||
      p_hat <= 0 || p_hat >= 1 || (p_hat + mu_hat) <= 0) {
      next
    }

    theta_hat <- est[est$Metric == "TPR", paste0("Group", a)]
    gamma_hat <- est[est$Metric == "FPR", paste0("Group", a)]
    f1_hat <- est[est$Metric == "F1", paste0("Group", a)]

    I_beta <- crossprod(W_N, W_N * w_N) / n_unlab
    I_beta_inv <- safe_inverse(I_beta)

    H_TPR <- colMeans(W_N * ((D_unlab_a - theta_hat) * w_N)) / mu_hat
    H_FPR <- -colMeans(W_N * ((D_unlab_a - gamma_hat) * w_N)) / nu_hat
    H_PPV <- colMeans(W_N * (D_unlab_a * w_N)) / p_hat
    H_NPV <- -colMeans(W_N * ((1 - D_unlab_a) * w_N)) / (1 - p_hat)
    H_F1 <- colMeans(W_N * ((2 * D_unlab_a - f1_hat) * w_N)) / (p_hat + mu_hat)
    H_ACC <- colMeans(W_N * ((2 * D_unlab_a - 1) * w_N))
    H_BS <- colMeans(W_N * ((1 - 2 * S_unlab_a) * w_N))

    score_n <- W_n * (Y_a - m_n)
    IF_beta <- t(I_beta_inv %*% t(score_n))

    IF_TPR <- drop(IF_beta %*% H_TPR)
    IF_FPR <- drop(IF_beta %*% H_FPR)
    IF_PPV <- drop(IF_beta %*% H_PPV)
    IF_NPV <- drop(IF_beta %*% H_NPV)
    IF_F1 <- drop(IF_beta %*% H_F1)
    IF_ACC <- drop(IF_beta %*% H_ACC)
    IF_BS <- drop(IF_beta %*% H_BS)

    metric_var <- c(
      TPR = stats::var(IF_TPR) / n_lab,
      FPR = stats::var(IF_FPR) / n_lab,
      PPV = stats::var(IF_PPV) / n_lab,
      NPV = stats::var(IF_NPV) / n_lab,
      F1 = stats::var(IF_F1) / n_lab,
      ACC = stats::var(IF_ACC) / n_lab,
      BS = stats::var(IF_BS) / n_lab
    )

    group_col <- paste0("Group", a)
    out[out$Metric %in% c("TPR", "FNR"), group_col] <- metric_var["TPR"]
    out[out$Metric %in% c("FPR", "TNR"), group_col] <- metric_var["FPR"]
    out[out$Metric == "PPV", group_col] <- metric_var["PPV"]
    out[out$Metric == "NPV", group_col] <- metric_var["NPV"]
    out[out$Metric == "F1", group_col] <- metric_var["F1"]
    out[out$Metric == "ACC", group_col] <- metric_var["ACC"]
    out[out$Metric == "BS", group_col] <- metric_var["BS"]
  }

  if (length(nclass) == 2) {
    out[, "Delta"] <- rowSums(out[, paste0("Group", nclass), drop = FALSE])
  } else if (length(nclass) > 2) {
    k <- length(nclass)
    pest_groups <- est[, paste0("Group", nclass), drop = FALSE]
    out_groups <- out[, paste0("Group", nclass), drop = FALSE]
    pest_mean <- rowMeans(pest_groups)
    out[, "var"] <- 4 / (k - 1)^2 *
      rowSums((pest_groups - pest_mean)^2 * out_groups)
    out[, "gei"] <- 1 / k^2 * rowSums(
      ((pest_groups - pest_mean) / (pest_mean)^2 -
        2 * est[, "gei"] / pest_mean)^2 * out_groups
    )
  }

  out
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}
