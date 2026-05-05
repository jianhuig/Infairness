# Purpose: Control-variate style bias correction for semi-supervised estimators.

cv_metric_values <- function(metric_table, group_col, metrics) {
  out <- metric_table[match(metrics, metric_table$Metric), group_col]
  names(out) <- metrics
  out
}

cv_recompute_estimate_summaries <- function(est_table, class) {
  group_cols <- paste0("Group", class)

  if (length(class) == 2) {
    est_table[, "Delta"] <- est_table[, group_cols[1]] - est_table[, group_cols[2]]
  } else if (length(class) > 2) {
    group_values <- est_table[, group_cols, drop = FALSE]
    est_table[, "mad"] <- as.numeric(apply(group_values, 1, mad))
    est_table[, "var"] <- as.numeric(apply(group_values, 1, var))
    if (exists("entropy", mode = "function")) {
      est_table[, "gei"] <- as.numeric(apply(group_values, 1, entropy))
    } else {
      est_table[, "gei"] <- NA_real_
    }
  }

  est_table
}

cv_recompute_variance_summaries <- function(var_table, est_table, class) {
  group_cols <- paste0("Group", class)

  if (length(class) == 2) {
    var_table[, "Delta"] <- rowSums(var_table[, group_cols, drop = FALSE])
  } else if (length(class) > 2) {
    k <- length(class)
    pest_groups <- est_table[, group_cols, drop = FALSE]
    out_groups <- var_table[, group_cols, drop = FALSE]
    pest_mean <- rowMeans(pest_groups)
    var_table[, "var"] <- 4 / (k - 1)^2 *
      rowSums((pest_groups - pest_mean)^2 * out_groups)
    var_table[, "gei"] <- 1 / k^2 * rowSums(
      ((pest_groups - pest_mean) / (pest_mean)^2 -
        2 * est_table[, "gei"] / pest_mean)^2 * out_groups
    )
  }

  var_table
}

cv_optimal_weight <- function(psi_sup, psi_cv) {
  denom <- stats::var(psi_cv)
  num <- stats::cov(psi_sup, psi_cv)

  if (!is.finite(denom) || denom <= .Machine$double.eps || !is.finite(num)) {
    return(1)
  }

  num / denom
}

cv_safe_div <- function(num, denom) {
  if (!is.finite(denom) || denom <= .Machine$double.eps || any(!is.finite(num))) {
    return(rep(NA_real_, length(num)))
  }

  num / denom
}

cv_pack_update <- function(theta_y,
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

  w_hat <- cv_optimal_weight(phi_y_for_w, phi_m)
  theta_cv <- theta_y - w_hat * (theta_m_lab - theta_m_all)
  phi_y <- phi_y_final_fun(theta_cv)
  phi <- phi_y - w_hat * phi_m
  var_hat <- stats::var(phi) / n

  if (!is.finite(theta_cv) || !is.finite(var_hat)) {
    return(list(theta = NA_real_, var = NA_real_, w = w_hat))
  }

  list(theta = theta_cv, var = var_hat, w = w_hat)
}

apply_control_variate_update <- function(est_plugin,
                                         var_plugin,
                                         Y_labeled,
                                         S_labeled,
                                         A_labeled,
                                         m_labeled,
                                         threshold_labeled,
                                         S_unlabeled,
                                         A_unlabeled,
                                         m_unlabeled,
                                         threshold_unlabeled) {
  class <- sort(unique(A_labeled))
  metrics <- c("TPR", "FPR", "PPV", "NPV", "F1", "ACC", "BS")
  group_cols <- paste0("Group", class)

  est_sup <- get_metric(
    Y = Y_labeled,
    S = S_labeled,
    A = A_labeled,
    threshold = threshold_labeled
  )
  est_lab <- get_metric(
    Y = m_labeled,
    S = S_labeled,
    A = A_labeled,
    threshold = threshold_labeled
  )

  est_out <- est_plugin
  var_out <- var_plugin
  cv_weights <- data.frame(
    Metric = est_out$Metric,
    matrix(NA_real_, nrow = nrow(est_out), ncol = length(group_cols))
  )
  names(cv_weights)[-1] <- group_cols

  if (length(threshold_labeled) == 1) {
    threshold_labeled <- rep(threshold_labeled, length(S_labeled))
  }
  if (length(threshold_unlabeled) == 1) {
    threshold_unlabeled <- rep(threshold_unlabeled, length(S_unlabeled))
  }

  for (a in class) {
    idx_a <- A_labeled == a
    idx_unlab_a <- A_unlabeled == a
    y_a <- Y_labeled[idx_a]
    m_a <- m_labeled[idx_a]
    s_a <- S_labeled[idx_a]
    m_unlab_a <- m_unlabeled[idx_unlab_a]
    s_unlab_a <- S_unlabeled[idx_unlab_a]
    m_all_a <- c(m_a, m_unlab_a)
    s_all_a <- c(s_a, s_unlab_a)
    threshold_all_a <- c(threshold_labeled[idx_a], threshold_unlabeled[idx_unlab_a])
    n_a <- length(y_a)
    group_col <- paste0("Group", a)
    plugin_vals <- cv_metric_values(est_plugin, group_col, metrics)
    plugin_var_vals <- cv_metric_values(var_plugin, group_col, metrics)

    corrected_vals <- plugin_vals
    corrected_var <- plugin_var_vals
    A_lab <- as.numeric(s_a > threshold_labeled[idx_a])
    A_all <- as.numeric(s_all_a > threshold_all_a)
    y0_lab <- 1 - y_a
    q_lab <- 1 - m_a
    q_all <- 1 - m_all_a

    pi1 <- mean(A_all)
    pi1_lab <- mean(A_lab)
    pi0 <- mean(1 - A_all)
    pi0_lab <- mean(1 - A_lab)
    eta1 <- mean(y_a)
    eta0 <- mean(y0_lab)
    mu1 <- mean(m_all_a)
    mu0 <- mean(q_all)
    kappa <- mean(s_all_a^2)

    metric_updates <- list(
      TPR = cv_pack_update(
        theta_y = cv_safe_div(mean(A_lab * y_a), eta1),
        theta_m_lab = cv_safe_div(mean(A_lab * m_a), mean(m_a)),
        theta_m_all = cv_safe_div(mean(A_all * m_all_a), mu1),
        phi_y_for_w = cv_safe_div((A_lab - cv_safe_div(mean(A_lab * y_a), eta1)) * y_a, eta1),
        phi_m = cv_safe_div((A_lab - cv_safe_div(mean(A_all * m_all_a), mu1)) * m_a, mu1),
        phi_y_final_fun = function(theta) {
          cv_safe_div((A_lab - theta) * y_a, eta1)
        },
        n = n_a
      ),
      FPR = cv_pack_update(
        theta_y = cv_safe_div(mean(A_lab * y0_lab), eta0),
        theta_m_lab = cv_safe_div(mean(A_lab * q_lab), mean(q_lab)),
        theta_m_all = cv_safe_div(mean(A_all * q_all), mu0),
        phi_y_for_w = cv_safe_div((A_lab - cv_safe_div(mean(A_lab * y0_lab), eta0)) * y0_lab, eta0),
        phi_m = cv_safe_div((A_lab - cv_safe_div(mean(A_all * q_all), mu0)) * q_lab, mu0),
        phi_y_final_fun = function(theta) {
          cv_safe_div((A_lab - theta) * y0_lab, eta0)
        },
        n = n_a
      ),
      PPV = cv_pack_update(
        theta_y = cv_safe_div(mean(A_lab * y_a), pi1_lab),
        theta_m_lab = cv_safe_div(mean(A_lab * m_a), pi1_lab),
        theta_m_all = cv_safe_div(mean(A_all * m_all_a), pi1),
        phi_y_for_w = cv_safe_div(
          A_lab * (y_a - cv_safe_div(mean(A_lab * y_a), pi1_lab)),
          pi1_lab
        ),
        phi_m = cv_safe_div(
          A_lab * (m_a - cv_safe_div(mean(A_all * m_all_a), pi1)),
          pi1
        ),
        phi_y_final_fun = function(theta) {
          cv_safe_div(A_lab * (y_a - theta), pi1_lab)
        },
        n = n_a
      ),
      NPV = cv_pack_update(
        theta_y = cv_safe_div(mean((1 - A_lab) * y0_lab), pi0_lab),
        theta_m_lab = cv_safe_div(mean((1 - A_lab) * q_lab), pi0_lab),
        theta_m_all = cv_safe_div(mean((1 - A_all) * q_all), pi0),
        phi_y_for_w = cv_safe_div(
          (1 - A_lab) * (y0_lab - cv_safe_div(mean((1 - A_lab) * y0_lab), pi0_lab)),
          pi0_lab
        ),
        phi_m = cv_safe_div(
          (1 - A_lab) * (q_lab - cv_safe_div(mean((1 - A_all) * q_all), pi0)),
          pi0
        ),
        phi_y_final_fun = function(theta) {
          cv_safe_div((1 - A_lab) * (y0_lab - theta), pi0_lab)
        },
        n = n_a
      ),
      ACC = cv_pack_update(
        theta_y = pi0 + mean((2 * A_lab - 1) * y_a),
        theta_m_lab = pi0 + mean((2 * A_lab - 1) * m_a),
        theta_m_all = pi0 + mean((2 * A_all - 1) * m_all_a),
        phi_y_for_w = (2 * A_lab - 1) * y_a - ((pi0 + mean((2 * A_lab - 1) * y_a)) - pi0),
        phi_m = (2 * A_lab - 1) * m_a - ((pi0 + mean((2 * A_all - 1) * m_all_a)) - pi0),
        phi_y_final_fun = function(theta) {
          (2 * A_lab - 1) * y_a - (theta - pi0)
        },
        n = n_a
      ),
      F1 = cv_pack_update(
        theta_y = cv_safe_div(2 * mean(A_lab * y_a), pi1_lab + eta1),
        theta_m_lab = cv_safe_div(2 * mean(A_lab * m_a), pi1_lab + mean(m_a)),
        theta_m_all = cv_safe_div(2 * mean(A_all * m_all_a), pi1 + mu1),
        phi_y_for_w = cv_safe_div(
          2 * A_lab * y_a -
            cv_safe_div(2 * mean(A_lab * y_a), pi1_lab + eta1) * (A_lab + y_a),
          pi1_lab + eta1
        ),
        phi_m = cv_safe_div(
          2 * A_lab * m_a -
            cv_safe_div(2 * mean(A_all * m_all_a), pi1 + mu1) * (A_lab + m_a),
          pi1 + mu1
        ),
        phi_y_final_fun = function(theta) {
          cv_safe_div(2 * A_lab * y_a - theta * (A_lab + y_a), pi1_lab + eta1)
        },
        n = n_a
      ),
      BS = cv_pack_update(
        theta_y = kappa + mean((1 - 2 * s_a) * y_a),
        theta_m_lab = kappa + mean((1 - 2 * s_a) * m_a),
        theta_m_all = kappa + mean((1 - 2 * s_all_a) * m_all_a),
        phi_y_for_w = (1 - 2 * s_a) * y_a - ((kappa + mean((1 - 2 * s_a) * y_a)) - kappa),
        phi_m = (1 - 2 * s_a) * m_a - ((kappa + mean((1 - 2 * s_all_a) * m_all_a)) - kappa),
        phi_y_final_fun = function(theta) {
          (1 - 2 * s_a) * y_a - (theta - kappa)
        },
        n = n_a
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

  est_out <- cv_recompute_estimate_summaries(est_out, class)
  var_out <- cv_recompute_variance_summaries(var_out, est_out, class)

  list(est = est_out, var = var_out, cv_weights = cv_weights)
}
