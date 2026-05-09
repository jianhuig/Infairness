# Purpose: Standalone control-corrected semi-supervised fairness estimator.

#' Control-corrected semi-supervised fairness estimation.
#'
#' `CC()` is a standalone control-variate estimator. It fits an imputation
#' model for `E(Y | S, X, A)` within each group, then uses labeled-data
#' supervised metrics corrected by the labeled-vs-full imputation contrast.
#'
#' @param Y Outcome, with missing values for unlabeled observations.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Threshold for classification based on the model score.
#' Default value is 0.5. Can be a scalar or one value per observation.
#' @param X Optional covariates matrix used by `"Spline(S) + X"` and
#' `"Spline Interaction"`.
#' @param basis Character vector giving the imputation basis within each group.
#' If a single value is supplied, it is reused for every group. Supported
#' values are `"Spline(S)"`, `"Spline(S) + X"`, `"Spline Interaction"`,
#' `"Beta"`, and `"kernel"`.
#' @param nknots Number of knots for spline-based branches. Default is 3.
#' @param D Logical; if `TRUE`, include the threshold indicator
#' `I(S > threshold)` in spline design matrices. Default is `FALSE`.
#' @param kernel_order Bandwidth exponent used by the kernel smoother.
#' Default is 0.45.
#' @param kernel_use_ecdf Logical; if `TRUE`, transform `S` to the within-group
#' empirical CDF scale before fitting the kernel smoother.
#' @param return_imputation Logical; if `TRUE`, return fitted imputations.
#' @param ... Reserved for future compatibility.
#'
#' @return A list containing `est`, `var`, `alpha`, and `cv_weights`. If
#' `return_imputation = TRUE`, the list also contains `m_labeled` and `m_all`.
#' @export
CC <- function(Y,
               S,
               A,
               threshold = 0.5,
               X = NULL,
               basis = "Spline Interaction",
               nknots = 3,
               D = FALSE,
               kernel_order = 0.45,
               kernel_use_ecdf = FALSE,
               return_imputation = FALSE,
               ...) {
  dots <- list(...)
  if (length(dots) > 0) {
    stop("Unused argument(s): ", paste(names(dots), collapse = ", "))
  }

  fit <- cc_fit_imputation(
    Y = Y,
    S = S,
    A = A,
    threshold = threshold,
    X = X,
    basis = basis,
    nknots = nknots,
    include_D = isTRUE(D),
    kernel_order = kernel_order,
    kernel_use_ecdf = kernel_use_ecdf
  )

  labeled_ind <- fit$labeled_ind
  Y_labeled <- Y[labeled_ind]
  S_labeled <- S[labeled_ind]
  A_labeled <- A[labeled_ind]
  threshold_labeled <- fit$threshold_raw[labeled_ind]

  est_plugin <- get_metric(
    Y = fit$m_all,
    S = S,
    A = A,
    threshold = fit$threshold_raw
  )

  var_plugin <- Influence_curve(
    pest = est_plugin,
    Y = Y_labeled,
    S = S_labeled,
    A = A_labeled,
    m = fit$m_labeled,
    threshold = threshold_labeled,
    method = "semi-supervised"
  )

  cc_update <- cc_apply_control_variate_update(
    est_plugin = est_plugin,
    var_plugin = var_plugin,
    Y_labeled = Y_labeled,
    S_labeled = S_labeled,
    A_labeled = A_labeled,
    m_labeled = fit$m_labeled,
    threshold_labeled = threshold_labeled,
    S_all = S,
    A_all = A,
    m_all = fit$m_all,
    threshold_all = fit$threshold_raw
  )

  result <- list(
    est = cc_update$est,
    var = cc_update$var,
    alpha = fit$alpha,
    cv_weights = cc_update$cv_weights
  )

  if (return_imputation) {
    result$m_labeled <- fit$m_labeled
    result$m_all <- fit$m_all
  }

  result
}

cc_fit_imputation <- function(Y,
                              S,
                              A,
                              threshold,
                              X,
                              basis,
                              nknots,
                              include_D,
                              kernel_order,
                              kernel_use_ecdf) {
  if (!(length(Y) == length(S) && length(S) == length(A))) {
    stop("`Y`, `S`, and `A` must have the same length.")
  }

  labeled_ind <- which(!is.na(Y))
  if (length(labeled_ind) == 0) {
    stop("`Y` must contain at least one labeled observation.")
  }

  nclass <- sort(unique(A))
  if (length(basis) == 1) {
    basis <- rep(basis, length(nclass))
  }
  if (length(basis) != length(nclass)) {
    stop("`basis` must have length 1 or one entry per group.")
  }

  supported_basis <- c("Spline(S)", "Spline(S) + X", "Spline Interaction", "Beta", "kernel")
  if (any(!basis %in% supported_basis)) {
    stop("Unsupported `basis`. Use Spline, Beta, or kernel bases.")
  }

  if (any(basis %in% c("Spline(S) + X", "Spline Interaction")) && is.null(X)) {
    stop("`X` is required for basis = \"Spline(S) + X\" or \"Spline Interaction\".")
  }
  if (!is.null(X)) {
    X <- as.matrix(X)
    if (nrow(X) != length(S)) {
      stop("`X` must have the same number of rows as `S`.")
    }
  }

  if (length(threshold) == 1) {
    threshold_raw <- rep(threshold, length(S))
  } else if (length(threshold) == length(S)) {
    threshold_raw <- threshold
  } else {
    stop("`threshold` must have length 1 or length equal to `S`.")
  }

  S_labeled <- S[labeled_ind]
  A_labeled <- A[labeled_ind]
  Y_labeled <- Y[labeled_ind]
  threshold_labeled <- threshold_raw[labeled_ind]

  if (any(tabulate(match(A_labeled, nclass), nbins = length(nclass)) == 0)) {
    stop("Each group must have at least one labeled observation.")
  }

  kernel_scale <- kernel_metric_scale(
    S = S,
    A = A,
    basis = basis,
    threshold = threshold_raw,
    use_ecdf = kernel_use_ecdf
  )

  C_all <- as.numeric(S > threshold_raw)
  C_labeled <- C_all[labeled_ind]
  D_all <- cc_spline_indicator(C_all, include_D)
  D_labeled <- cc_spline_indicator(C_labeled, include_D)

  m_all <- S
  alpha <- c()

  for (a in nclass) {
    basis_a <- basis[which(nclass == a)]
    idx_all <- A == a
    idx_lab <- A_labeled == a

    if (basis_a == "Spline(S)") {
      alpha <- c(alpha, nknots)

      basis_labeled <- NaturalSplineBasis(as.matrix(S_labeled), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_all <- NaturalSplineBasis(as.matrix(S), nknots, knots = spline_knots)

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[idx_lab, , drop = FALSE],
        D = cc_subset_matrix(D_labeled, idx_lab)
      )
      predictors_all <- build_spline_predictors(
        basis_matrix = basis_all[idx_all, , drop = FALSE],
        D = cc_subset_matrix(D_all, idx_all)
      )

      model <- fit_spline_model(
        X = predictors_labeled,
        y = Y_labeled[idx_lab],
        penalty_factor = build_spline_penalty_factor(
          basis_matrix = basis_labeled[idx_lab, , drop = FALSE],
          D = cc_subset_matrix(D_labeled, idx_lab)
        )
      )

      m_all[idx_all] <- predict_spline_model(model, predictors_all)
    } else if (basis_a == "Spline(S) + X") {
      alpha <- c(alpha, nknots)

      basis_labeled <- NaturalSplineBasis(as.matrix(S_labeled), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_all <- NaturalSplineBasis(as.matrix(S), nknots, knots = spline_knots)
      X_labeled <- X[labeled_ind, , drop = FALSE]

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[idx_lab, , drop = FALSE],
        D = cc_subset_matrix(D_labeled, idx_lab),
        X = X_labeled[idx_lab, , drop = FALSE]
      )
      predictors_all <- build_spline_predictors(
        basis_matrix = basis_all[idx_all, , drop = FALSE],
        D = cc_subset_matrix(D_all, idx_all),
        X = X[idx_all, , drop = FALSE]
      )

      model <- fit_spline_model(
        X = predictors_labeled,
        y = Y_labeled[idx_lab],
        penalty_factor = build_spline_penalty_factor(
          basis_matrix = basis_labeled[idx_lab, , drop = FALSE],
          D = cc_subset_matrix(D_labeled, idx_lab),
          X = X_labeled[idx_lab, , drop = FALSE]
        )
      )

      m_all[idx_all] <- predict_spline_model(model, predictors_all)
    } else if (basis_a == "Spline Interaction") {
      alpha <- c(alpha, nknots)

      basis_labeled <- NaturalSplineBasis(as.matrix(S_labeled), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_all <- NaturalSplineBasis(as.matrix(S), nknots, knots = spline_knots)
      X_labeled <- X[labeled_ind, , drop = FALSE]

      X_labeled_spline_int <- do.call(
        cbind,
        lapply(seq_len(ncol(X_labeled)), function(j) basis_labeled * X_labeled[, j])
      )
      X_all_spline_int <- do.call(
        cbind,
        lapply(seq_len(ncol(X)), function(j) basis_all * X[, j])
      )

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[idx_lab, , drop = FALSE],
        D = cc_subset_matrix(D_labeled, idx_lab),
        X = X_labeled[idx_lab, , drop = FALSE],
        X_interaction = X_labeled_spline_int[idx_lab, , drop = FALSE]
      )
      predictors_all <- build_spline_predictors(
        basis_matrix = basis_all[idx_all, , drop = FALSE],
        D = cc_subset_matrix(D_all, idx_all),
        X = X[idx_all, , drop = FALSE],
        X_interaction = X_all_spline_int[idx_all, , drop = FALSE]
      )

      model <- fit_spline_model(
        X = predictors_labeled,
        y = Y_labeled[idx_lab],
        penalty_factor = build_spline_penalty_factor(
          basis_matrix = basis_labeled[idx_lab, , drop = FALSE],
          D = cc_subset_matrix(D_labeled, idx_lab),
          X = X_labeled[idx_lab, , drop = FALSE],
          X_interaction = X_labeled_spline_int[idx_lab, , drop = FALSE]
        )
      )

      m_all[idx_all] <- predict_spline_model(model, predictors_all)
    } else if (basis_a == "Beta") {
      calibrated <- FitParametricCalibration(
        Y_labeled = Y_labeled,
        S_labeled = S_labeled,
        A_labeled = A_labeled,
        A_val = a,
        S_unlabeled = S,
        method = "Beta"
      )
      m_all[idx_all] <- calibrated[idx_all]
    } else if (basis_a == "kernel") {
      bandwidth <- kernel_bandwidth(
        kernel_scale$S[labeled_ind][idx_lab],
        order = kernel_order
      )
      m_all[idx_all] <- npreg(
        kernel_scale$S[labeled_ind][idx_lab],
        Y_labeled[idx_lab],
        kernel_scale$S[idx_all],
        bandwidth
      )
    }
  }

  list(
    m_all = as.numeric(m_all),
    m_labeled = as.numeric(m_all[labeled_ind]),
    labeled_ind = labeled_ind,
    threshold_raw = threshold_raw,
    alpha = alpha
  )
}

cc_spline_indicator <- function(C, include_D) {
  if (!isTRUE(include_D)) {
    return(NULL)
  }
  as.matrix(C)
}

cc_subset_matrix <- function(x, idx) {
  if (is.null(x)) {
    return(NULL)
  }
  as.matrix(x)[idx, , drop = FALSE]
}

cc_metric_values <- function(metric_table, group_col, metrics) {
  out <- metric_table[match(metrics, metric_table$Metric), group_col]
  names(out) <- metrics
  out
}

cc_recompute_estimate_summaries <- function(est_table, class) {
  group_cols <- paste0("Group", class)

  if (length(class) == 2) {
    est_table[, "Delta"] <- est_table[, group_cols[1]] - est_table[, group_cols[2]]
  } else if (length(class) > 2) {
    est_groups <- est_table[, group_cols, drop = FALSE]
    est_table[, "mad"] <- as.numeric(apply(est_groups, 1, stats::mad))
    est_table[, "var"] <- as.numeric(apply(est_groups, 1, stats::var))
    est_table[, "gei"] <- as.numeric(apply(est_groups, 1, entropy))
  }

  est_table
}

cc_recompute_variance_summaries <- function(var_table, est_table, class) {
  group_cols <- paste0("Group", class)

  if (length(class) == 2) {
    var_table[, "Delta"] <- rowSums(var_table[, group_cols, drop = FALSE])
  } else if (length(class) > 2) {
    k <- length(class)
    pest_groups <- est_table[, group_cols, drop = FALSE]
    var_groups <- var_table[, group_cols, drop = FALSE]
    pest_mean <- rowMeans(pest_groups)
    var_table[, "var"] <- 4 / (k - 1)^2 *
      rowSums((pest_groups - pest_mean)^2 * var_groups)
    var_table[, "gei"] <- 1 / k^2 * rowSums(
      ((pest_groups - pest_mean) / (pest_mean)^2 -
        2 * est_table[, "gei"] / pest_mean)^2 * var_groups
    )
  }

  var_table
}

cc_optimal_weight <- function(phi_y, phi_m) {
  denom <- stats::var(phi_m)
  num <- stats::cov(phi_y, phi_m)

  if (!is.finite(denom) || denom <= .Machine$double.eps || !is.finite(num)) {
    return(1)
  }

  num / denom
}

cc_safe_div <- function(num, denom) {
  if (!is.finite(denom) || denom <= .Machine$double.eps || any(!is.finite(num))) {
    return(rep(NA_real_, length(num)))
  }

  num / denom
}

cc_pack_update <- function(theta_y,
                           theta_m_lab,
                           theta_m_all,
                           phi_y_for_w,
                           phi_m,
                           phi_y_final_fun,
                           n) {
  inputs <- c(theta_y, theta_m_lab, theta_m_all)
  if (length(phi_y_for_w) == 0 || length(phi_m) == 0 || any(!is.finite(inputs))) {
    return(list(theta = NA_real_, var = NA_real_, w = 1))
  }

  w_hat <- cc_optimal_weight(phi_y_for_w, phi_m)
  theta <- theta_y - w_hat * (theta_m_lab - theta_m_all)
  phi_y <- phi_y_final_fun(theta)
  phi <- phi_y - w_hat * phi_m
  var_hat <- stats::var(phi) / n

  if (!is.finite(theta) || !is.finite(var_hat)) {
    return(list(theta = NA_real_, var = NA_real_, w = w_hat))
  }

  list(theta = theta, var = var_hat, w = w_hat)
}

cc_apply_control_variate_update <- function(est_plugin,
                                            var_plugin,
                                            Y_labeled,
                                            S_labeled,
                                            A_labeled,
                                            m_labeled,
                                            threshold_labeled,
                                            S_all,
                                            A_all,
                                            m_all,
                                            threshold_all) {
  class <- sort(unique(A_all))
  metrics <- c("TPR", "FPR", "PPV", "NPV", "F1", "ACC", "BS")
  group_cols <- paste0("Group", class)

  est_out <- est_plugin
  var_out <- var_plugin
  cv_weights <- data.frame(
    Metric = est_out$Metric,
    matrix(NA_real_, nrow = nrow(est_out), ncol = length(group_cols))
  )
  names(cv_weights)[-1] <- group_cols

  for (a in class) {
    idx_lab <- A_labeled == a
    idx_all <- A_all == a

    y_lab <- Y_labeled[idx_lab]
    s_lab <- S_labeled[idx_lab]
    m_lab <- m_labeled[idx_lab]
    cutoff_lab <- threshold_labeled[idx_lab]

    s_all <- S_all[idx_all]
    m_all_a <- m_all[idx_all]
    cutoff_all <- threshold_all[idx_all]

    n <- length(y_lab)
    group_col <- paste0("Group", a)
    corrected_vals <- cc_metric_values(est_out, group_col, metrics)
    corrected_var <- cc_metric_values(var_out, group_col, metrics)

    A_lab <- as.numeric(s_lab > cutoff_lab)
    A_pop <- as.numeric(s_all > cutoff_all)
    B_lab <- 1 - A_lab
    B_pop <- 1 - A_pop
    y0_lab <- 1 - y_lab
    q_lab <- 1 - m_lab
    q_pop <- 1 - m_all_a

    pi1 <- mean(A_pop)
    pi0 <- mean(B_pop)
    pi1_lab <- mean(A_lab)
    pi0_lab <- mean(B_lab)
    eta1 <- mean(y_lab)
    eta0 <- mean(y0_lab)
    mu1 <- mean(m_all_a)
    mu0 <- mean(q_pop)
    kappa <- mean(s_all^2)

    tpr_y <- cc_safe_div(mean(A_lab * y_lab), eta1)
    tpr_m_all <- cc_safe_div(mean(A_pop * m_all_a), mu1)

    fpr_y <- cc_safe_div(mean(A_lab * y0_lab), eta0)
    fpr_m_all <- cc_safe_div(mean(A_pop * q_pop), mu0)

    ppv_y <- cc_safe_div(mean(A_lab * y_lab), pi1_lab)
    ppv_m_all <- cc_safe_div(mean(A_pop * m_all_a), pi1)

    npv_y <- cc_safe_div(mean(B_lab * y0_lab), pi0_lab)
    npv_m_all <- cc_safe_div(mean(B_pop * q_pop), pi0)

    f1_y <- cc_safe_div(2 * mean(A_lab * y_lab), pi1_lab + eta1)
    f1_m_all <- cc_safe_div(2 * mean(A_pop * m_all_a), pi1 + mu1)

    acc_y <- pi0 + mean((2 * A_lab - 1) * y_lab)
    acc_m_all <- pi0 + mean((2 * A_pop - 1) * m_all_a)

    brier_y <- kappa + mean((1 - 2 * s_lab) * y_lab)
    brier_m_all <- kappa + mean((1 - 2 * s_all) * m_all_a)

    metric_updates <- list(
      TPR = cc_pack_update(
        theta_y = tpr_y,
        theta_m_lab = cc_safe_div(mean(A_lab * m_lab), mean(m_lab)),
        theta_m_all = tpr_m_all,
        phi_y_for_w = cc_safe_div((A_lab - tpr_y) * y_lab, eta1),
        phi_m = cc_safe_div((A_lab - tpr_m_all) * m_lab, mu1),
        phi_y_final_fun = function(theta) {
          cc_safe_div((A_lab - theta) * y_lab, eta1)
        },
        n = n
      ),
      FPR = cc_pack_update(
        theta_y = fpr_y,
        theta_m_lab = cc_safe_div(mean(A_lab * q_lab), mean(q_lab)),
        theta_m_all = fpr_m_all,
        phi_y_for_w = cc_safe_div((A_lab - fpr_y) * y0_lab, eta0),
        phi_m = cc_safe_div((A_lab - fpr_m_all) * q_lab, mu0),
        phi_y_final_fun = function(theta) {
          cc_safe_div((A_lab - theta) * y0_lab, eta0)
        },
        n = n
      ),
      PPV = cc_pack_update(
        theta_y = ppv_y,
        theta_m_lab = cc_safe_div(mean(A_lab * m_lab), pi1_lab),
        theta_m_all = ppv_m_all,
        phi_y_for_w = cc_safe_div(A_lab * (y_lab - ppv_y), pi1_lab),
        phi_m = cc_safe_div(A_lab * (m_lab - ppv_m_all), pi1),
        phi_y_final_fun = function(theta) {
          cc_safe_div(A_lab * (y_lab - theta), pi1_lab)
        },
        n = n
      ),
      NPV = cc_pack_update(
        theta_y = npv_y,
        theta_m_lab = cc_safe_div(mean(B_lab * q_lab), pi0_lab),
        theta_m_all = npv_m_all,
        phi_y_for_w = cc_safe_div(B_lab * (y0_lab - npv_y), pi0),
        phi_m = cc_safe_div(B_lab * (q_lab - npv_m_all), pi0),
        phi_y_final_fun = function(theta) {
          cc_safe_div(B_lab * (y0_lab - theta), pi0)
        },
        n = n
      ),
      ACC = cc_pack_update(
        theta_y = acc_y,
        theta_m_lab = pi0 + mean((2 * A_lab - 1) * m_lab),
        theta_m_all = acc_m_all,
        phi_y_for_w = (2 * A_lab - 1) * y_lab - (acc_y - pi0),
        phi_m = (2 * A_lab - 1) * m_lab - (acc_m_all - pi0),
        phi_y_final_fun = function(theta) {
          (2 * A_lab - 1) * y_lab - (theta - pi0)
        },
        n = n
      ),
      F1 = cc_pack_update(
        theta_y = f1_y,
        theta_m_lab = cc_safe_div(2 * mean(A_lab * m_lab), pi1_lab + mean(m_lab)),
        theta_m_all = f1_m_all,
        phi_y_for_w = cc_safe_div(
          2 * A_lab * y_lab - f1_y * (A_lab + y_lab),
          pi1 + eta1
        ),
        phi_m = cc_safe_div(
          2 * A_lab * m_lab - f1_m_all * (A_lab + m_lab),
          pi1 + mu1
        ),
        phi_y_final_fun = function(theta) {
          cc_safe_div(2 * A_lab * y_lab - theta * (A_lab + y_lab), pi1 + eta1)
        },
        n = n
      ),
      BS = cc_pack_update(
        theta_y = brier_y,
        theta_m_lab = kappa + mean((1 - 2 * s_lab) * m_lab),
        theta_m_all = brier_m_all,
        phi_y_for_w = (1 - 2 * s_lab) * y_lab - (brier_y - kappa),
        phi_m = (1 - 2 * s_lab) * m_lab - (brier_m_all - kappa),
        phi_y_final_fun = function(theta) {
          (1 - 2 * s_lab) * y_lab - (theta - kappa)
        },
        n = n
      )
    )

    for (metric in metrics) {
      update <- metric_updates[[metric]]
      if (!is.null(update) && is.finite(update$theta) && is.finite(update$var)) {
        corrected_vals[metric] <- update$theta
        corrected_var[metric] <- update$var
        cv_weights[cv_weights$Metric == metric, group_col] <- update$w
      }
    }

    est_out[est_out$Metric == "TPR", group_col] <- corrected_vals["TPR"]
    est_out[est_out$Metric == "FNR", group_col] <- 1 - corrected_vals["TPR"]
    est_out[est_out$Metric == "FPR", group_col] <- corrected_vals["FPR"]
    est_out[est_out$Metric == "TNR", group_col] <- 1 - corrected_vals["FPR"]
    est_out[est_out$Metric == "PPV", group_col] <- corrected_vals["PPV"]
    est_out[est_out$Metric == "NPV", group_col] <- corrected_vals["NPV"]
    est_out[est_out$Metric == "F1", group_col] <- corrected_vals["F1"]
    est_out[est_out$Metric == "ACC", group_col] <- corrected_vals["ACC"]
    est_out[est_out$Metric == "BS", group_col] <- corrected_vals["BS"]

    var_out[var_out$Metric %in% c("TPR", "FNR"), group_col] <- corrected_var["TPR"]
    var_out[var_out$Metric %in% c("FPR", "TNR"), group_col] <- corrected_var["FPR"]
    var_out[var_out$Metric == "PPV", group_col] <- corrected_var["PPV"]
    var_out[var_out$Metric == "NPV", group_col] <- corrected_var["NPV"]
    var_out[var_out$Metric == "F1", group_col] <- corrected_var["F1"]
    var_out[var_out$Metric == "ACC", group_col] <- corrected_var["ACC"]
    var_out[var_out$Metric == "BS", group_col] <- corrected_var["BS"]

    cv_weights[cv_weights$Metric == "FNR", group_col] <- cv_weights[cv_weights$Metric == "TPR", group_col]
    cv_weights[cv_weights$Metric == "TNR", group_col] <- cv_weights[cv_weights$Metric == "FPR", group_col]
  }

  est_out <- cc_recompute_estimate_summaries(est_out, class)
  var_out <- cc_recompute_variance_summaries(var_out, est_out, class)

  list(est = est_out, var = var_out, cv_weights = cv_weights)
}
