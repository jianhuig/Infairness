#' Select best model
#'
#' When candidate models were fit with cross-fitted imputation quality, they
#' should share the same labeled-data folds. Pass the same `folds` argument to
#' each `SSFairness(..., return_imputation_quality = TRUE)` call before
#' comparing them here.
#' @export

Select_Model <- function(models, Y, S, A, threshold) {
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

  # map each metric to the correct column in imp_quality
  metric_col <- c(log_loss = 2, brier = 1, ece = 4)
  metric_col <- c(brier = 1)

  winners <- do.call(
    rbind,
    lapply(names(metric_col), function(m) {
      col <- metric_col[[m]]
      v0 <- v1 <- c()
      for (ss in models) {
        if (!all(is.na(ss))) {
          v0 <- c(v0, ss$imp_quality[1, col])
          v1 <- c(v1, ss$imp_quality[2, col])
        } else {
          v0 <- c(v0, Inf)
          v1 <- c(v1, Inf)
        }
      }
      c(which.min(v0), which.min(v1)) # winner idx for outcome 0 and 1
    })
  )

  best_model_idx <- apply(winners, 2, function(x) {
    as.integer(names(which.max(table(x))))
  })

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

  return(list(est = est, var = var, alpha = alphas))
}
