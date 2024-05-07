#' Semi-supervised fairness estimation.
#'
#' @param Y Outcome, can contain missing values.
#' @param S Model score.
#' @param A Group indicator.
#' @param threshold Threshold for classification based on the model score.
#' @export
#'

Infairness_nobasis <- function(Y,
                       S,
                       A,
                       threshold = 0.5,
                       W = NULL) {
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
  
  # Augmentation in each class
  m_unlabeled <- S_unlabeled # imputed values
  m_labeled <- S_labeled # imputed values
  
  nclass <- sort(unique(A))
  for (a in nclass) {
    X <- cbind(S_labeled[A_labeled == a], C_labeled[A_labeled == a])
    gamma <- glm(Y_labeled[A_labeled == a] ~ X, family = binomial)$coefficients
    index_vector <- A_unlabeled == a
    imputed_unlabeled <- boot::inv.logit(as.matrix(cbind(1, S_unlabeled[index_vector], C_unlabeled[index_vector])) %*% gamma)
    m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
    imputed_labeled <- boot::inv.logit(as.matrix(cbind(1, S_labeled[A_labeled == a], C_labeled[A_labeled == a])) %*% gamma)
    m_labeled[A_labeled == a] <- imputed_labeled
  }
  
  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold, W = W
  )
  var <- Influence_curve(est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold, method = "semi-supervised")
  
  return(list(est = est, var = var, m_unlabeled = m_unlabeled, m_labeled = m_labeled))
}
