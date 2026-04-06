#' # Purpose: Estimation of fairness metrics in supervised setting.
#' # Updated: 2023-07-17
#' 

#' Supervised fairness estimation with Cross-Fitting
#'
#' @param Y Outcome in the observed dataset.
#' @param S Model score in observed dataset.
#' @param A Group indicator in observed dataset.
#' @param threshold Threshold for classification. Default is 0.5.
#' @param W Weight. Default is Null.
#' @param V Number of folds for cross-fitting. Default is 5.
#' @export

SupervisedFairness <- function(Y,
                                  S,
                                  A,
                                  threshold = 0.5,
                                  W = NULL,
                                  V = 10) {
  
  lab <- which(!is.na(Y))
  Y <- Y[lab]; S <- S[lab]; A <- A[lab]
  n <- length(Y)
  if (is.null(W)) W <- rep(1, n)
  folds <- integer(n)
  for (a in sort(unique(A))) {
    idx_a <- which(A == a)
    folds[idx_a] <- sample(rep(1:V, length.out = length(idx_a)))
  }
  
  # Accumulators
  accumulated_sigma <- NULL
  fold_counts <- NULL # Track how many valid folds we have per metric/column
  
  for (v in 1:V) {
    idx_test <- which(folds == v)
    idx_train <- which(folds != v)
    
    # 1. Estimation
    est_v <- get_metric(Y[idx_train], S[idx_train], A[idx_train], threshold)
    
    # 2. Influence Curve (might return NaNs for PPV)
    res_table <- Influence_curve(pest = est_v, Y = Y[idx_test], S = S[idx_test], A = A[idx_test], 
                                 threshold = threshold, method = "supervised")
    
    # 3. Extract Numeric Variances
    vars_v <- as.matrix(res_table[, c("Group0", "Group1", "Delta")])
    
    # 4. Un-scale to get Sum of Squares
    n_test <- length(idx_test)
    sum_sq_v <- vars_v * (n_test^2)
    
    # Initialize if first iteration
    if (is.null(accumulated_sigma)) {
      accumulated_sigma <- matrix(0, nrow=nrow(vars_v), ncol=ncol(vars_v))
      fold_counts <- matrix(0, nrow=nrow(vars_v), ncol=ncol(vars_v))
    }
    
    # 5. Robust Addition: Only add if not NaN
    # Create a mask of valid (non-NA) entries in this fold
    valid_mask <- !is.na(sum_sq_v)
    
    # Add values where valid
    accumulated_sigma[valid_mask] <- accumulated_sigma[valid_mask] + sum_sq_v[valid_mask]
    
    # Increment counters where valid
    fold_counts[valid_mask] <- fold_counts[valid_mask] + 1
  }
  
  # Global Point Estimate
  final_est <- get_metric(Y, S, A, threshold)
  
  # 6. Final Calculation with Dynamic Scaling
  # If a metric failed in all folds, fold_counts is 0 -> Avoid division by zero
  fold_counts[fold_counts == 0] <- 1 
  
  # Scale: We summed 'k' valid folds. We want to estimate the total sum over 'V' folds.
  # Projected Total = (Sum of Valid / k) * V
  # Then Divide by N^2
  
  projected_total_sigma <- (accumulated_sigma / fold_counts) * V
  final_vars_numeric <- projected_total_sigma / (n^2)
  colnames(final_vars_numeric) <- c("Group0", "Group1", "Delta")
  
  # Re-attach names
  final_var_df <- data.frame(Metric = final_est$Metric, final_vars_numeric)
  
  return(list(est = final_est, var = final_var_df))
}
