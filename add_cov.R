library(doParallel)
library(dplyr)
library(SSFairnessAudit)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(dplyr)
  library(SSFairnessAudit)
  library(glmnet)
})

#args <- commandArgs(TRUE)
args = c("balanced", "misspecified 3")
if (args[1] == "balanced") {
  prev <- 0.5 # Prevalence of the protected attribute.
} else {
  prev <- 0.4 # Prevalence of the protected attribute.
}
print(prev)

nclass <- 2

n <- 500 * nclass # Labeled data size.
N <- 1e4 * nclass # Unlabeled data size.
rho <- 0.4
threshold <- 0.5
model <- args[2]
print(model)

set.seed(1234)
indep <- DataGeneration(
  n_labeled = 3000,
  N_unlabeled = 0,
  prot_att_prevalence = prev,
  model = model,
  rho = rho
)
mean(indep$Y[indep$A == 0])
mean(indep$Y[indep$A == 1])

model_0 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 0) %>% select(-Y_miss, -A, -X_11:-X_15))
model_1 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 1) %>% select(-Y_miss, -A, -X_11:-X_15))

clusterExport(cl = cl, varlist = list(
  "indep", "nclass", "n", "N",
  "p", "rho", "threshold", "model", "prev", "model_0", "model_1"
))
clusterEvalQ(cl, source("TwoStep.R"))

nsim <- 1e4
result <- pbapply::pblapply(cl = cl, X = 1:nsim, FUN = function(j) {
  # generate the main dataset
  dat <- DataGeneration(
    n_labeled = n,
    N_unlabeled = N,
    prot_att_prevalence = prev,
    model = model,
    p = 15,
    rho = rho
  )
  
  # using indepdent data to train S
  dat$S <- rep(NA, nrow(dat))
  # for (a in unique(indep$A)) {
  #   model_S <- glm(Y ~ .,
  #                  family = binomial(),
  #                  data = indep %>% filter(A == a) %>% select(-Y_miss, -A))
  #
  #   # Predict S for the corresponding group in dat
  #   dat$S[dat$A == a] <- predict(model_S, newdata = dat %>% filter(A == a), type = "response")
  # }
  dat$S[dat$A == 0] <- predict(model_0, newdata = dat %>% filter(A == 0), type = "response")
  dat$S[dat$A == 1] <- predict(model_1, newdata = dat %>% filter(A == 1), type = "response")
  
  # prepare main data
  dat$C <- ifelse(dat$S > threshold, 1, 0)
  labeled <- dat %>% filter(!is.na(Y_miss))
  unlabeled <- dat %>% filter(is.na(Y_miss))
  
  # check AUC
  # pROC::auc(dat$Y, dat$S)
  # pROC::auc(dat$Y[dat$A==0], dat$S[dat$A==0])
  # pROC::auc(dat$Y[dat$A==1], dat$S[dat$A==1])
  
  # oracle
  oracle <- Audit_Fairness(Y = dat$Y,
                           S = dat$S,
                           A = dat$A,
                           threshold = threshold,
                           method = "supervised")
  
  # sup
  sup <- Audit_Fairness(Y = labeled$Y,
                        S = labeled$S,
                        A = labeled$A,
                        threshold = threshold,
                        method = "supervised")
  
  # naive
  # naive <- SupervisedFairness(
  #   Y = dat$S,
  #   S = dat$S,
  #   A = dat$A,
  #   threshold = threshold
  # )
  
  
  # GLM(S)
  # ss_s <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = NULL,
  #   basis = "glm(S)"
  # )$est
  
  # ss_s_all <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = dat %>% select(starts_with("X")) %>% as.matrix(),
  #   basis = "glm(S + X)"
  # )$est
  
  ss_ridge_s <- Audit_Fairness(Y = dat$Y_miss,
                               S = dat$S,
                               A = dat$A,
                               threshold = threshold,
                               method = "semi-supervised",
                               basis = "ridge(S)")
  
  ss_ridge_sw <- Audit_Fairness(Y = dat$Y_miss,
                                S = dat$S,
                                A = dat$A,
                                threshold = threshold,
                                method = "semi-supervised",
                                X = dat %>% select(X_11:X_15) %>% as.matrix(),
                                basis = "ridge(S + X)")
  
  # GLM(Spline(S))
  ss_spline_s <- Audit_Fairness(Y = dat$Y_miss,
                             S = dat$S,
                             A = dat$A,
                             threshold = threshold,
                             method = "semi-supervised",
                             basis = "Spline(S)")
    
  ss_spline_s_w <- Audit_Fairness(Y = dat$Y_miss,
                             S = dat$S,
                             A = dat$A,
                             threshold = threshold,
                             method = "semi-supervised",
                             X = dat %>% select(X_11:X_15) %>% as.matrix(),
                             basis = "Spline(S) + X")  
  
  ss_spline_sw <- Audit_Fairness(Y = dat$Y_miss,
                                  S = dat$S,
                                  A = dat$A,
                                  threshold = threshold,
                                  method = "semi-supervised",
                                  X = dat %>% select(X_11:X_15) %>% as.matrix(),
                                  basis = "Spline(S, X)")  
  
  # Poly(S)
  ss_poly_s <-  Audit_Fairness(Y = dat$Y_miss,
                             S = dat$S,
                             A = dat$A,
                             threshold = threshold,
                             method = "semi-supervised",
                             basis = "Poly(S)")

  # Poly(S) + X
  ss_poly_s_w <- Audit_Fairness(Y = dat$Y_miss,
                             S = dat$S,
                             A = dat$A,
                             threshold = threshold,
                             method = "semi-supervised",
                             X = dat %>% select(X_11:X_15) %>% as.matrix(),
                             basis = "Poly(S) + X")
  
  
  # Poly(S) + all
  ss_poly_sw <- Audit_Fairness(Y = dat$Y_miss,
                             S = dat$S,
                             A = dat$A,
                             threshold = threshold,
                             method = "semi-supervised",
                             X = dat %>% select(X_11:X_15) %>% as.matrix(),
                             basis = "Poly(S, X)")
  
  ss_inter <- Audit_Fairness(Y = dat$Y_miss,
                             S = dat$S,
                             A = dat$A,
                             threshold = threshold,
                             method = "semi-supervised",
                             X = dat %>% select(X_11:X_15) %>% as.matrix(),
                             basis = "interaction")
  
  
  # ss_platt <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = NULL,
  #   basis = "Platt"
  # )
  # 
  # # Beta Calibration
  # ss_Beta <- Infairness(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = NULL,
  #   basis = "Beta"
  # )
  
  # ss_twostep_s <- Infairness_twostep(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = NULL,
  #   basis = "S"
  # )
  # 
  # ss_twostep_sx <- Infairness_twostep(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = dat %>% select(X_1:X_10),
  #   basis = "S,X"
  # )
  # 
  # ss_twostep_sall <- Infairness_twostep(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = dat %>% select(starts_with("X")),
  #   basis = "S,X"
  # )
  # 
  # ss_twostep_x <- Infairness_twostep(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = dat %>% select(contains("X")),
  #   basis = "X"
  # )
  # 
  # ss_twostep_xall <- Infairness_twostep(
  #   Y = dat$Y_miss,
  #   S = dat$S,
  #   A = dat$A,
  #   X = dat %>% select(contains("X")),
  #   basis = "X"
  # )
  
  
  list(
    oracle = oracle,
    sup = sup,
    ss_ridge_s = ss_ridge_s,
    ss_ridge_sw = ss_ridge_sw,
    ss_spline_s = ss_spline_s,
    ss_spline_s_w = ss_spline_s_w,
    ss_spline_sw = ss_spline_sw,
    ss_poly_s = ss_poly_s,
    ss_poly_s_w = ss_poly_s_w,
    ss_poly_sw = ss_poly_sw,
    ss_inter = ss_inter
  )
  
})
stopCluster(cl)
saveRDS(result, file = paste0("0324_", args[1], "_", args[2], ".rds"))
