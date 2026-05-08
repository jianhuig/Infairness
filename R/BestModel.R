#' Select best model
#'
#' When candidate models were fit with cross-fitted imputation quality, they
#' should share the same labeled-data folds. Pass the same `folds` argument to
#' each `SSFairness(..., return_imputation_quality = TRUE)` call before
#' comparing them here.
#'
#' @param models List of candidate `SSFairness()` fits with
#' `return_imputation_quality = TRUE`.
#' @param Y Outcome vector with missing values for unlabeled observations.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Threshold for classification based on the model score.
#' @param criterion Model-selection criterion computed from `imp_quality`.
#' Default is `"weighted_mse"`, which uses a TPR-weighted cross-fitted
#' squared-error criterion. `Select_Model()` uses this single criterion
#' directly, not a vote across metrics. Use `"bic"` for the plain BIC score
#' or `"mbic"` for a modified BIC with penalty multiplier
#' `min(n^0.1, log(n))`, stored for regression-basis candidates.
#' @return List containing `est`, `var`, and `alpha`.
#' @export

Select_Model <- function(models, Y, S, A, threshold,
                         criterion = "weighted_mse") {
  criterion <- match.arg(
    criterion,
    choices = c("weighted_mse", "bic", "mbic")
  )

  model_folds <- lapply(models, function(x) attr(x, "cv_folds"))
  non_null_folds <- Filter(Negate(is.null), model_folds)

  if (length(non_null_folds) > 1) {
    same_folds <- vapply(
      non_null_folds[-1],
      identical,
      logical(1),
      y = non_null_folds[[1]]
    )
    if (!all(same_folds)) {
      stop(
        "Candidate models were evaluated on different cross-fitting folds. ",
        "Pass the same `folds` argument to each SSFairness() call before ",
        "calling Select_Model()."
      )
    }
  }

  labeled_ind <- which(!is.na(Y))

  # Labeled data.
  Y_labeled <- Y[labeled_ind]
  S_labeled <- S[labeled_ind]
  A_labeled <- A[labeled_ind]

  # Unlabeled data.
  S_unlabeled <- S[-labeled_ind]
  A_unlabeled <- A[-labeled_ind]

  metric_name <- c(
    weighted_mse = "Weighted_MSE",
    bic = "BIC",
    mbic = "MBIC"
  )[[criterion]]

  valid_models <- Filter(function(x) !all(is.na(x)), models)
  if (length(valid_models) == 0) {
    stop("No valid candidate models were supplied to `Select_Model()`.")
  }
  template_model <- valid_models[[1]]
  if (!metric_name %in% names(template_model$imp_quality)) {
    stop(
      "Selection criterion `", criterion,
      "` is not available in `imp_quality`."
    )
  }

  fallback <- Inf
  model_scores <- matrix(
    fallback,
    nrow = length(models),
    ncol = 2,
    dimnames = list(NULL, c("Group1", "Group2"))
  )

  for (j in seq_along(models)) {
    ss <- models[[j]]
    if (!all(is.na(ss))) {
      model_scores[j, 1] <- ss$imp_quality[1, metric_name]
      model_scores[j, 2] <- ss$imp_quality[2, metric_name]
    }
  }

  best_model_idx <- apply(model_scores, 2, which.min)

  m_unlabeled <- models[[best_model_idx[1]]]$m_unlabeled
  m_labeled <- models[[best_model_idx[1]]]$m_labeled

  a <- sort(unique(A))[2]
  m_unlabeled[A_unlabeled == a] <- models[[best_model_idx[[2]]]]$m_unlabeled[
    A_unlabeled == a
  ]
  m_labeled[A_labeled == a] <- models[[best_model_idx[[2]]]]$m_labeled[
    A_labeled == a
  ]

  alphas <- c(
    models[[best_model_idx[1]]]$alpha[1],
    models[[best_model_idx[2]]]$alpha[2]
  )

  est <- get_metric(
    Y = m_unlabeled,
    S = S_unlabeled,
    A = A_unlabeled,
    threshold = threshold
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

  list(est = est, var = var, alpha = alphas)
}
