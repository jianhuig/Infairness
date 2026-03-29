#' @export
ImputeQuality_GIC <- function(Y, S, A, threshold = 0.5, X = NULL,
                             max_alpha = 10, patience = 2) {
  labeled_ind <- which(!is.na(Y))
  unlabeled_ind <- which(is.na(Y))

  C_ind <- ifelse(S > threshold, 1, 0)

  Y_lab <- Y[labeled_ind]
  S_lab <- S[labeled_ind]
  A_lab <- A[labeled_ind]
  C_lab <- C_ind[labeled_ind]
  X_lab <- if (!is.null(X)) X[labeled_ind, , drop = FALSE] else NULL

  # Ensure X has names to prevent lookup errors later
  if (!is.null(X)) {
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("Var", 1:ncol(X))
    }
    # Update X_lab with names if they were missing
    if (!is.null(X_lab)) colnames(X_lab) <- colnames(X)
  }

  final_models <- list()
  nclass <- sort(unique(A_lab))

  # Global m placeholders
  m_global <- rep(NA, length(Y_lab)) # Vector size of Labeled
  m_unlabeled <- rep(NA, length(unlabeled_ind)) # Vector size of Unlabeled
  basis <- list()
  
  # --- Iterate Per Group ---
  for (a in nclass) {
    cat(sprintf("\n=== Processing Group %s ===\n", a))

    idx_a <- which(A_lab == a)
    Y_a <- Y_lab[idx_a]
    S_a <- S_lab[idx_a]
    C_a <- C_lab[idx_a]
    X_a <- if (!is.null(X_lab)) X_lab[idx_a, , drop = FALSE] else NULL

    X_names <- if (!is.null(X_a)) colnames(X_a) else NULL

    folds_a <- caret::createFolds(factor(Y_a), k = 10, list = FALSE)

    alpha_results <- list()
    best_global_loss <- Inf
    no_improve_counter <- 0

    # --- Outer Loop: Iterate Alpha ---
    for (alpha in 1:max_alpha) {
      cat(sprintf("-> Evaluating Alpha %d\n", alpha))
      
      # A. Base Model Setup
      S_poly <- poly(S_a, degree = alpha, raw = TRUE)
      colnames(S_poly) <- paste0("S_deg", 1:alpha)

      base_basis <- cbind(1, S_poly, C_a)
      colnames(base_basis) <- c("Intercept", colnames(S_poly), "Threshold_C")
      
      pX <- ncol(X)
      P_alpha <- (2 + max_alpha) + 2 * pX   # intercept + S_deg(1:alpha) + C + main + interactions
      current_loss <- compute_gbic(Y_a, base_basis, P_total = P_alpha, gamma_ebic = 1)
      current_basis <- base_basis
      selected_vars_indices <- c()

      # B. Forward Selection for X
      if (!is.null(X_a)) {
        candidates <- 1:ncol(X_a)

        repeat {
          best_cand <- NULL
          best_loss_step <- current_loss

          for (j in candidates) {
            test_basis <- cbind(current_basis, X_a[, j, drop = FALSE])
            score <- compute_gbic(Y_a, test_basis, P_total = P_alpha, gamma_ebic = 1)

            if (score< best_loss_step) {
              best_loss_step <- score
              best_cand <- j
            }
          }

          if (!is.null(best_cand)) {
            current_loss <- best_loss_step
            selected_vars_indices <- c(selected_vars_indices, best_cand)

            X_col <- X_a[, best_cand, drop = FALSE]
            colnames(X_col) <- X_names[best_cand]

            current_basis <- cbind(current_basis, X_col)
            candidates <- setdiff(candidates, best_cand)
          } else {
            break
          }
        }
      }

      # C. Interactions Selection (S * X)
      if (length(selected_vars_indices) > 0) {
        candidates_inter <- selected_vars_indices

        repeat{
          best_cand <- NULL
          best_loss_step <- current_loss

          for (j in candidates_inter) {
            interaction_term <- S_a * X_a[, j]
            test_basis <- cbind(current_basis, interaction_term)
            score <- compute_gbic(Y_a, test_basis, P_total = P_alpha, gamma_ebic = 1)

            if (score < best_loss_step) {
              best_loss_step <- score
              best_cand <- j
            }
          }

          if (!is.null(best_cand)) {
            current_loss <- best_loss_step
            interaction_term <- S_a * X_a[, best_cand]

            var_name <- X_names[best_cand]
            new_col_name <- paste0("S_x_", var_name)

            test_basis_named <- cbind(current_basis, interaction_term)
            colnames(test_basis_named)[ncol(test_basis_named)] <- new_col_name

            current_basis <- test_basis_named
            candidates_inter <- setdiff(candidates_inter, best_cand)
          } else {
            break
          }
        }
      }

      full_model_vars <- colnames(current_basis)

      alpha_results[[alpha]] <- list(
        degree = alpha,
        final_loss = current_loss,
        all_vars = full_model_vars
      )

      cat(sprintf(
        "  Alpha %d | Loss: %.5f | Model: [%s]\n",
        alpha, current_loss, paste(full_model_vars, collapse = ", ")
      ))

      if (current_loss < best_global_loss) {
        best_global_loss <- current_loss
        no_improve_counter <- 0
      } else {
        no_improve_counter <- no_improve_counter + 1
      }

      if (no_improve_counter >= patience) {
        cat(sprintf("  >> Stopping Early: No improvement for %d alphas.\n", patience))
        break
      }
    }

    # Winner Selection
    valid_results <- alpha_results[!sapply(alpha_results, is.null)]
    all_losses <- sapply(valid_results, function(x) x$final_loss)
    winner_idx <- which.min(all_losses)
    winner_model <- valid_results[[winner_idx]]

    final_models[[as.character(a)]] <- winner_model
    
    S_poly_all <- poly(S, degree = winner_model$degree, raw = TRUE)
    colnames(S_poly_all) <- paste0("S_deg", 1:winner_model$degree)
    basis_temp <- cbind(Intercept=1, S_poly_all, Threshold_C=C_ind)
    existing_names <- winner_model$all_vars
    extra_vars <- setdiff(existing_names, c("Intercept", colnames(S_poly_all), "Threshold_C"))

    
    for (var_name in extra_vars) {
      if (startsWith(var_name, "S_x_")) {
        orig_x_name <- sub("S_x_", "", var_name)
        # Ensure X_unlab has names
        term <- S * X[, orig_x_name]
      } else {
        term <- X[, var_name]
      }
      basis_temp <- cbind(basis_temp, term)
    }
    
    basis[[as.character(a)]] <- basis_temp
    
}


  return(StepII(Y, S, A, threshold, basis[[as.character(nclass[1])]],
                basis[[as.character(nclass[2])]]))
}


StepII <- function(Y,
                   S,
                   A,
                   threshold = 0.5,
                   basis1 = NULL,
                   basis2 = NULL,
                   k = 10) {
  labeled_ind <- which(!is.na(Y))
  
  # Labeled data.
  Y_labeled <- Y[labeled_ind]
  A_labeled <- A[labeled_ind]
  S_labeled <- S[labeled_ind]
  
  # Unlabeled data.
  A_unlabeled <- A[-labeled_ind]
  S_unlabeled <- S[-labeled_ind]
  
  nclass <- sort(unique(A))
  
  # Augmentation in each class
  m_unlabeled <- rep(NULL, length(A_unlabeled))
  m_labeled <- rep(NULL, length(A_labeled))
  
  g1 <- nclass[1]; g2 <- nclass[2]
  
  for (a in nclass) {
    if (a == g1) X <- as.matrix(basis1) else X <- as.matrix(basis2)
    
    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]
    
    
    gamma <- tryCatch(
      {
        RidgeRegression(
          X = X_labeled[A_labeled == a, -1],
          y = Y_labeled[A_labeled == a],
          coef = log(ncol(X_labeled)),
          weights = NULL
        )
      },
      error = function(e) {
        print("Ridge Regression produced an error")
        print(e)
      }
    )
    
    imputed_unlabeled <- boot::inv.logit(as.matrix(X_unlabeled[A_unlabeled == a, ]
    ) %*% gamma$coefficients)
    
    m_unlabeled[A_unlabeled == a] <- imputed_unlabeled
  }
  
  for (a in nclass) {
    if (a == g1) X <- as.matrix(basis1) else X <- as.matrix(basis2)
    
    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]
    
    Y_a <- Y_labeled[A_labeled == a]
    X_labeled_a <- X_labeled[A_labeled == a, , drop = FALSE]
    X_unlabeled_a <- X_unlabeled[A_labeled == a, , drop = FALSE]
    
    fold <- caret::createFolds(factor(Y_a), k = k, list = FALSE)
    
    for (i in 1:k) {
      train_id <- which(fold != i)
      test_id <- which(fold == i)
      
      gamma <- tryCatch(
        {
          RidgeRegression(
            X = X_labeled_a[train_id, -1],
            y = Y_a[train_id],
            coef = log(ncol(X_labeled_a)),
            weights = NULL
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )$coefficients
      
      
      imputed_test <- boot::inv.logit(
        as.matrix(X_labeled_a[test_id, , drop = FALSE]
        )
        %*% gamma)
      
      m_labeled[which(A_labeled == a)[test_id] ] <- imputed_test
    }
  }
  
  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold, W = NULL
  )
  
  var <- Influence_curve(
    est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold,
    method = "semi-supervised"
  )
  
  return(list(est = est, var = var))
}