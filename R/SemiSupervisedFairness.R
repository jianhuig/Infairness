# Purpose: Estimation of fairness metrics in semi-supervised setting.

#' Semi-supervised fairness estimation.
#'
#' @param Y Outcome, can contain missing values.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Threshold for classification based on the model score.
#' Default value is 0.5.
#' @param X Optional covariates matrix to be adjusted in the semi-supervised setting.
#' Default is NULL.
#' @param basis Character vector giving the basis strategy used within each
#' group. If a single value is supplied, it is reused for every group.
#' Supported values in the current implementation include `"Spline(S)"`,
#' `"Spline(S) + X"`, `"Spline Interaction"`, `"Beta"`, and `"kernel"`.
#' `"Spline(S) + X"` uses a shared spline in `S` plus additive covariate
#' effects, while `"Spline Interaction"` augments that model with spline-by-`X`
#' interaction terms so the shape in `S` can vary with `X`.
#' @param nknots Number of knots used by spline-based branches
#' (`"Spline(S)"`, `"Spline(S) + X"`, and `"Spline Interaction"`). Default is
#' 3. Ignored otherwise.
#' @param W Optional observation weights. Currently reserved for compatibility.
#' @param k Number of folds used for labeled-data cross-fitting when
#' `cross_fit_variance = TRUE`. Default is 10.
#' @param folds Optional named list of precomputed fold assignments, keyed by
#' group value. Supplying the same `folds` object across candidate models makes
#' imputation-quality comparisons use the same labeled-data splits.
#' @param cross_fit_variance Logical; if `TRUE`, use `k`-fold cross-fitted
#' labeled imputations when estimating the variance. If `FALSE`, labeled
#' imputations are computed in-sample.
#' @param return_imputation_quality Logical; if `TRUE`, attach imputation
#' quality metrics and the labeled/unlabeled imputations to the result.
#' @param kernel_order Bandwidth exponent used by the kernel smoother.
#' Default is 0.45.
#' @param kernel_use_ecdf Logical; if `TRUE`, transform `S` to the within-group
#' empirical CDF scale before fitting the kernel smoother. Reported metrics,
#' including Brier score, remain on the original score scale.
#' @param ... Additional arguments. Removed options are rejected explicitly.
#' @return List of estimated fairness metrics and their variances.
#' The returned list always contains `est`, `var`, and `alpha`. When
#' `return_imputation_quality = TRUE`, it also contains `imp_quality`,
#' `m_labeled`, and `m_unlabeled`.
#' @export
#'

SSFairness <- function(
  Y,
  S,
  A,
  threshold = 0.5,
  X = NULL,
  basis = c("Spline(S)", "Spline(S)"),
  nknots = 3,
  W = NULL,
  k = 10,
  folds = NULL,
  cross_fit_variance = FALSE,
  return_imputation_quality = FALSE,
  kernel_order = 0.45,
  kernel_use_ecdf = FALSE,
  ...
) {
  dots <- list(...)
  removed_args <- intersect(
    names(dots),
    c(
      "control_variate",
      "variance_method",
      "bootstrap_reps",
      "bootstrap_resample",
      "bootstrap_seed",
      "bootstrap_cores",
      "spline_ridge",
      "spline_gam",
      "ridge_unpenalize_binary_X",
      "spline_include_D",
      "log_include_D"
    )
  )
  if (length(removed_args) > 0) {
    stop(
      "Removed argument(s): ",
      paste(removed_args, collapse = ", "),
      ". These tuning options are no longer supported."
    )
  }
  removed_basis <- c(
    "Poly(S)", "Poly(S) + X", "Interaction",
    "Log(S)", "Log(S) + X", "Log Interaction"
  )
  if (any(basis %in% removed_basis)) {
    stop(
      "This basis option has been removed. Use Spline, Beta, or kernel ",
      "bases instead."
    )
  }
  spline_ridge <- TRUE
  spline_gam <- FALSE
  ridge_unpenalize_binary_X <- FALSE
  spline_include_D <- TRUE

  if (cross_fit_variance || return_imputation_quality) {
    quality_fit <- compute_imputation_quality(
      Y = Y,
      S = S,
      A = A,
      threshold = threshold,
      X = X,
      basis = basis,
      nknots = nknots,
      k = k,
      folds = folds,
      cross_fit = cross_fit_variance,
      kernel_order = kernel_order,
      kernel_use_ecdf = kernel_use_ecdf,
      ...
    )

    result <- list(
      est = quality_fit$est,
      var = quality_fit$var,
      alpha = quality_fit$alpha
    )

    if (return_imputation_quality) {
      result$imp_quality <- quality_fit$imp_quality
      result$m_labeled <- quality_fit$m_labeled
      result$m_unlabeled <- quality_fit$m_unlabeled
    }

    attr(result, "cv_folds") <- quality_fit$cv_folds

    return(result)
  }

  if (!is.null(W)) {
    W_label <- 4 * rbeta(sum(!is.na(Y)), 1 / 2, 3 / 2)
  } else {
    W_label <- NULL
  }

  labeled_ind <- which(!is.na(Y))

  # Create Indicator
  C <- ifelse(S > threshold, 1, 0)

  # Labeled data.
  Y_labeled <- Y[labeled_ind]
  S_labeled <- S[labeled_ind]
  A_labeled <- A[labeled_ind]
  C_labeled <- C[labeled_ind]

  # Unlabeled data.
  S_unlabeled <- S[-labeled_ind]
  A_unlabeled <- A[-labeled_ind]
  C_unlabeled <- C[-labeled_ind]
  if (!is.null(X)) {
    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]
  }
  nclass <- sort(unique(A))
  if (length(basis) == 1) {
    basis <- rep(basis, length(nclass))
  }

  if (length(threshold) == 1) {
    threshold_raw <- rep(threshold, length(S))
  } else if (length(threshold) == length(S)) {
    threshold_raw <- threshold
  } else {
    stop("`threshold` must have length 1 or length equal to `S`.")
  }

  kernel_scale <- kernel_metric_scale(
    S = S,
    A = A,
    basis = basis,
    threshold = threshold,
    use_ecdf = kernel_use_ecdf
  )
  S_metric <- kernel_scale$S
  threshold_metric <- kernel_scale$threshold
  S_labeled_metric <- S_metric[labeled_ind]
  S_unlabeled_metric <- S_metric[-labeled_ind]
  threshold_labeled_metric <- threshold_metric[labeled_ind]
  threshold_unlabeled_metric <- threshold_metric[-labeled_ind]
  threshold_labeled_raw <- threshold_raw[labeled_ind]
  threshold_unlabeled_raw <- threshold_raw[-labeled_ind]

  # Augmentation in each class
  m_unlabeled <- S_unlabeled # imputed values
  m_labeled <- S_labeled # imputed values
  alphas <- c() # store alpha values

  for (a in nclass) {
    basis_a <- basis[which(nclass == a)]

    if (basis_a == "Spline(S)") {
      alphas <- c(alphas, nknots)

      basis_labeled <- NaturalSplineBasis(S_labeled %>% as.matrix(), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_unlabeled <- NaturalSplineBasis(
        S_unlabeled %>% as.matrix(),
        nknots,
        knots = spline_knots
      )

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D)
      )
      predictors_unlabeled <- build_spline_predictors(
        basis_matrix = basis_unlabeled[A_unlabeled == a, , drop = FALSE],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D)
      )

      gamma <- tryCatch(
        {
          fit_spline_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a],
            use_ridge = spline_ridge,
            use_gam = spline_gam,
            S = S_labeled[A_labeled == a],
            D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
            num_knots = nknots,
            penalty_factor = build_spline_penalty_factor(
              basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
              D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Spline model fitting produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- predict_spline_model(
        gamma,
        predictors_unlabeled,
        S = S_unlabeled[A_unlabeled == a],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D)
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- predict_spline_model(
        gamma,
        predictors_labeled,
        S = S_labeled[A_labeled == a],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D)
      )
      m_labeled[A_labeled == a] <- imputed_labeled
    } else if (basis_a == "Spline(S) + X") {
      alphas <- c(alphas, nknots)

      basis_labeled <- NaturalSplineBasis(S_labeled %>% as.matrix(), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_unlabeled <- NaturalSplineBasis(
        S_unlabeled %>% as.matrix(),
        nknots,
        knots = spline_knots
      )

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X = X_labeled[A_labeled == a, , drop = FALSE]
      )
      predictors_unlabeled <- build_spline_predictors(
        basis_matrix = basis_unlabeled[A_unlabeled == a, , drop = FALSE],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X = X_unlabeled[A_unlabeled == a, , drop = FALSE]
      )

      gamma <- tryCatch(
        {
          fit_spline_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a],
            use_ridge = spline_ridge,
            use_gam = spline_gam,
            S = S_labeled[A_labeled == a],
            D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
            X_cov = X_labeled[A_labeled == a, , drop = FALSE],
            num_knots = nknots,
            penalty_factor = build_spline_penalty_factor(
              basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
              D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
              X = X_labeled[A_labeled == a, , drop = FALSE],
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Spline model fitting produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- predict_spline_model(
        gamma,
        predictors_unlabeled,
        S = S_unlabeled[A_unlabeled == a],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X_cov = X_unlabeled[A_unlabeled == a, , drop = FALSE]
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- predict_spline_model(
        gamma,
        predictors_labeled,
        S = S_labeled[A_labeled == a],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X_cov = X_labeled[A_labeled == a, , drop = FALSE]
      )
      m_labeled[A_labeled == a] <- imputed_labeled
    } else if (basis_a == "Spline Interaction") {
      alphas <- c(alphas, nknots)

      basis_labeled <- NaturalSplineBasis(S_labeled %>% as.matrix(), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_unlabeled <- NaturalSplineBasis(
        S_unlabeled %>% as.matrix(),
        nknots,
        knots = spline_knots
      )

      X_labeled_spline_int <- do.call(
        cbind,
        lapply(seq_len(ncol(X_labeled)), function(j) {
          basis_labeled * X_labeled[, j]
        })
      )
      X_unlabeled_spline_int <- do.call(
        cbind,
        lapply(seq_len(ncol(X_unlabeled)), function(j) {
          basis_unlabeled * X_unlabeled[, j]
        })
      )

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X = X_labeled[A_labeled == a, , drop = FALSE],
        X_interaction = X_labeled_spline_int[A_labeled == a, , drop = FALSE]
      )
      predictors_unlabeled <- build_spline_predictors(
        basis_matrix = basis_unlabeled[A_unlabeled == a, , drop = FALSE],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X = X_unlabeled[A_unlabeled == a, , drop = FALSE],
        X_interaction = X_unlabeled_spline_int[A_unlabeled == a, , drop = FALSE]
      )

      gamma <- tryCatch(
        {
          fit_spline_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a],
            use_ridge = spline_ridge,
            use_gam = spline_gam,
            S = S_labeled[A_labeled == a],
            D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
            X_cov = X_labeled[A_labeled == a, , drop = FALSE],
            interaction = TRUE,
            num_knots = nknots,
            penalty_factor = build_spline_penalty_factor(
              basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
              D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
              X = X_labeled[A_labeled == a, , drop = FALSE],
              X_interaction = X_labeled_spline_int[A_labeled == a, , drop = FALSE],
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Spline model fitting produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- predict_spline_model(
        gamma,
        predictors_unlabeled,
        S = S_unlabeled[A_unlabeled == a],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X_cov = X_unlabeled[A_unlabeled == a, , drop = FALSE]
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- predict_spline_model(
        gamma,
        predictors_labeled,
        S = S_labeled[A_labeled == a],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X_cov = X_labeled[A_labeled == a, , drop = FALSE]
      )
      m_labeled[A_labeled == a] <- imputed_labeled
    } else if (basis_a == "Beta") {
      S_calibrated <- FitParametricCalibration(
        Y_labeled,
        S_labeled,
        A_labeled,
        A_val = a,
        S_unlabeled,
        method = "Beta"
      )
      m_unlabeled[A_unlabeled == a] <- S_calibrated[A_unlabeled == a]
      S_calibrated <- FitParametricCalibration(
        Y_labeled,
        S_labeled,
        A_labeled,
        A_val = a,
        S_labeled,
        method = "Beta"
      )
      m_labeled[A_labeled == a] <- S_calibrated[A_labeled == a]
    } else if(basis_a == "kernel"){
      bandwidth <- kernel_bandwidth(
        S_labeled_metric[A_labeled == a],
        order = kernel_order
      )
      mhat <- npreg(
        S_labeled_metric[A_labeled == a],
        Y_labeled[A_labeled == a],
        S_unlabeled_metric,
        bandwidth
      )
      m_unlabeled[A_unlabeled == a] <- mhat[A_unlabeled == a]

      m_labeled[A_labeled == a] <- npreg(
        S_labeled_metric[A_labeled == a],
        Y_labeled[A_labeled == a],
        S_labeled_metric[A_labeled == a],
        bandwidth
      )
    } else {
      stop("Basis not recognized.")
    }
  }

  est <- get_metric(
    Y = m_unlabeled,
    S = S_unlabeled,
    A = A_unlabeled,
    threshold = threshold_unlabeled_raw,
    W = NULL
  )

  var <- Influence_curve(
    est,
    Y_labeled,
    S_labeled,
    A_labeled,
    m_labeled,
    threshold_labeled_raw,
    method = "semi-supervised"
  )

  if (any(basis == "Beta")) {
    var <- beta_if_variance_update(
      var_est = var,
      est = est,
      Y_labeled = Y_labeled,
      S_labeled = S_labeled,
      A_labeled = A_labeled,
      S_unlabeled = S_unlabeled,
      A_unlabeled = A_unlabeled,
      threshold_unlabeled = threshold_unlabeled_raw,
      basis = basis
    )
  }

  result <- list(est = est, var = var, alpha = alphas)

  return(result)
}
