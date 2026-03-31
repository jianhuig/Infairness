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
#' Supported values in the current implementation include `"Poly(S)"`,
#' `"Poly(S) + X"`, `"Spline(S)"`, `"Spline(S) + X"`, `"Spline Interaction"`,
#' `"Interaction"`, `"Beta"`, and `"kernel"`.
#' `"Spline(S) + X"` uses a shared spline in `S` plus additive covariate
#' effects, while `"Spline Interaction"` augments that model with spline-by-`X`
#' interaction terms so the shape in `S` can vary with `X`.
#' @param nknots Number of knots used by spline-based branches
#' (`"Spline(S)"`, `"Spline(S) + X"`, and `"Spline Interaction"`). Default is
#' 3. Ignored otherwise.
#' @param cross_fit_variance Logical; if `TRUE`, use cross-fitted labeled
#' imputations when estimating the variance through the shared imputation-quality
#' path.
#' @param return_imputation_quality Logical; if `TRUE`, attach imputation
#' quality metrics and the labeled/unlabeled imputations to the result.
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
  basis = c("Poly(S)", "Poly(S)"),
  nknots = 3,
  W = NULL,
  k = 10,
  cross_fit_variance = FALSE,
  return_imputation_quality = FALSE,
  ...
) {
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

  # Augmentation in each class
  m_unlabeled <- S_unlabeled # imputed values
  m_labeled <- S_labeled # imputed values
  alphas <- c() # store alpha values

  for (a in nclass) {
    basis_a <- basis[which(nclass == a)]

    if (basis_a == "Poly(S)") {
      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = C_labeled[A_labeled == a] %>% as.matrix()
          )
        },
        error = function(e) {
          1
        }
      )

      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a]
            ) %>%
              as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a]
          )
        ) %*%
          gamma$coefficients
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
      
      imputed_labeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_labeled[A_labeled == a, -1],
            C_labeled[A_labeled == a]
          )
        ) %*%
          gamma$coefficients
      )
      m_labeled[A_labeled == a] <- imputed_labeled

    } else if (basis_a == "Poly(S) + X") {
      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = cbind(
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>%
              as.matrix()
          )
        },
        error = function(e) {
          1
        }
      )

      # print(alpha)
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>%
              as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled[A_unlabeled == a, ]
          )
        ) %*%
          gamma$coefficients
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
      
      imputed_labeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_labeled[A_labeled == a, -1],
            C_labeled[A_labeled == a],
            X_labeled[A_labeled == a, ]
          )
        ) %*%
          gamma$coefficients
      )
      m_labeled[A_labeled == a] <- imputed_labeled
    } else if (basis_a == "Spline(S)") {
      alphas <- c(alphas, nknots)

      basis_labeled <- NaturalSplineBasis(S_labeled %>% as.matrix(), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_unlabeled <- NaturalSplineBasis(
        S_unlabeled %>% as.matrix(),
        nknots,
        knots = spline_knots
      )

      gamma <- tryCatch(
        {
          SplineRidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, ],
              C_labeled[A_labeled == a]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Spline ridge regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, ],
            C_unlabeled[A_unlabeled == a]
          )
        ) %*% gamma$coefficients
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_labeled[A_labeled == a, ],
            C_labeled[A_labeled == a]
          )
        ) %*% gamma$coefficients
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

      gamma <- tryCatch(
        {
          SplineRidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, ],
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Spline ridge regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, ],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled[A_unlabeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      imputed_labeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_labeled[A_labeled == a, ],
            C_labeled[A_labeled == a],
            X_labeled[A_labeled == a, ]
          )
        ) %*% gamma$coefficients
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

      gamma <- tryCatch(
        {
          SplineRidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, ],
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, , drop = FALSE],
              X_labeled_spline_int[A_labeled == a, ]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Spline ridge regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, ],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled[A_unlabeled == a, , drop = FALSE],
            X_unlabeled_spline_int[A_unlabeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      X_a <- X_labeled[A_labeled == a, , drop = FALSE]
      C_a <- C_labeled[A_labeled == a]

      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)
      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)

        basis_train <- NaturalSplineBasis(
          S_a[train_id] %>% as.matrix(),
          nknots,
          return_knots = TRUE
        )
        fold_knots <- attr(basis_train, "knots")
        basis_test <- NaturalSplineBasis(
          S_a[test_id] %>% as.matrix(),
          nknots,
          knots = fold_knots
        )

        X_train_int <- do.call(
          cbind,
          lapply(seq_len(ncol(X_a)), function(j) {
            basis_train * X_a[train_id, j]
          })
        )
        X_test_int <- do.call(
          cbind,
          lapply(seq_len(ncol(X_a)), function(j) {
            basis_test * X_a[test_id, j]
          })
        )

        gamma <- tryCatch(
          {
            SplineRidgeRegression(
              X = cbind(
                basis_train,
                C_a[train_id],
                X_a[train_id, , drop = FALSE],
                X_train_int
              ) %>% as.matrix(),
              y = Y_a[train_id]
            )
          },
          error = function(e) {
            print("Spline ridge regression produced an error")
            print(e)
          }
        )$coefficients

        imputed_test <- boot::inv.logit(
          as.matrix(
            cbind(
              1,
              basis_test,
              C_a[test_id],
              X_a[test_id, , drop = FALSE],
              X_test_int
            )
          ) %*% gamma
        )
        m_labeled[A_labeled == a][test_id] <- imputed_test
      }
    } else if (basis_a == "Interaction") {
      X_int <- S * as.matrix(X)
      X_int <- cbind(X, X_int)
      X_labeled_int <- X_int[labeled_ind, , drop = FALSE]
      X_unlabeled_int <- X_int[-labeled_ind, , drop = FALSE]

      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = cbind(
              C_labeled[A_labeled == a],
              X_labeled_int[A_labeled == a, ]
            ) %>%
              as.matrix()
          )
        },
        error = function(e) {
          1
        }
      )

      # print(alpha)
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a],
              X_labeled_int[A_labeled == a, ]
            ) %>%
              as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      imputed_unlabeled <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled_int[A_unlabeled == a, ]
          )
        ) %*%
          gamma$coefficients
      )

      m_unlabeled[A_unlabeled == a] <- imputed_unlabeled

      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      X_a <- X_labeled_int[A_labeled == a, , drop = FALSE]
      C_a <- C_labeled[A_labeled == a]

      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)
      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)

        alpha_fold <- tryCatch(
          {
            find_alpha_glm(
              Y = Y_a[train_id],
              covariates_matrix = as.matrix(S_a[train_id]),
              additional_matrix = as.matrix(cbind(
                X_a[train_id, ],
                C_a[train_id]
              ))
            )
          },
          error = function(e) {
            1
          }
        )

        basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
        basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

        gamma <- tryCatch(
          {
            RidgeRegression(
              X = cbind(
                basis_train[, -1],
                C_a[train_id],
                X_a[train_id, ]
              ) %>%
                as.matrix(),
              y = Y_a[train_id],
              weights = NULL
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients

        imputed_test <- boot::inv.logit(
          as.matrix(
            cbind(
              1,
              basis_test[, -1],
              C_a[test_id],
              X_a[test_id, ]
            )
          ) %*%
            gamma
        )
        m_labeled[A_labeled == a][test_id] <- imputed_test
      }
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
      bandwidth <- sqrt(var(S_labeled[A_labeled == a])) / (length(Y_labeled[A_labeled == a])^0.45)
      mhat <- npreg(S_labeled[A_labeled == a], Y_labeled[A_labeled == a], S_unlabeled, bandwidth)
      m_unlabeled[A_unlabeled == a] <- mhat[A_unlabeled == a]
      
      Y_a <- Y_labeled[A_labeled == a]
      S_a <- S_labeled[A_labeled == a]
      C_a <- C_labeled[A_labeled == a]
      fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)
      for (i in 1:k) {
        train_id <- which(fold != i)
        test_id <- which(fold == i)
        
        bandwidth_fold <- sqrt(var(S_a[train_id])) / (length(Y_a[train_id])^0.45)
        mhat_fold <- npreg(S_a[train_id], Y_a[train_id], S_a[test_id], bandwidth_fold)
        m_labeled[A_labeled == a][test_id] <- mhat_fold
      }
      } else {
      stop("Basis not recognized.")
    }
  }

  est <- get_metric(
    Y = m_unlabeled,
    S = S_unlabeled,
    A = A_unlabeled,
    threshold = threshold,
    W = NULL
  )

  var <- Influence_curve(
    est,
    Y_labeled,
    S_labeled,
    A_labeled,
    m_labeled,
    threshold,
    method = "semi-supervised"
  )

  return(list(est = est, var = var, alpha = alphas))
}
