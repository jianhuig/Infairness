library(dplyr)
library(Infairness)
library(glmnet)

args = c("unbalanced", "misspecified 1")
if (args[1] == "balanced") {
  prev <- 0.5 # Prevalence of the protected attribute.
} else {
  prev <- 0.4 # Prevalence of the protected attribute.
}

nclass <- 2

n <- 500 * nclass # Labeled data size.
N <- 1e4 * nclass # Unlabeled data size.
rho <- 0.4
threshold <- 0.5
model <- args[2]

set.seed(1234)
indep <- DataGeneration(
  n_labeled = 3000,
  N_unlabeled = 0,
  prot_att_prevalence = prev,
  model = model,
  rho = rho
)

model_0 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 0)
               %>% select(Y, contains("X")))
model_1 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 1) 
               %>% select(Y, contains("X")))

nsim <- 1e4
results <- list()
result <- for (i in 1:nsim) {
  # generate the main dataset
  dat <- DataGeneration(
    n_labeled = n,
    N_unlabeled = N,
    prot_att_prevalence = prev,
    model = model,
    rho = rho
  )
  
  # using indepdent data to train S
  dat$S <- rep(NA, nrow(dat))
  dat$S[dat$A == 0] <- predict(model_0, newdata = dat %>% filter(A == 0), 
                               type = "response")
  dat$S[dat$A == 1] <- predict(model_1, newdata = dat %>% filter(A == 1), 
                               type = "response")
  
  # prepare main data
  dat$C <- ifelse(dat$S > threshold, 1, 0)
  labeled <- dat %>% filter(!is.na(Y_miss))
  unlabeled <- dat %>% filter(is.na(Y_miss))
  
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
  
  # Infairness (S)
  ss_s <-  Audit_Fairness(Y = dat$Y_miss,
                               S = dat$S,
                               A = dat$A,
                               threshold = threshold,
                               method = "Infairness",
                               basis = "Poly(S)")
  
  # Infairness (S, W)
  ss_sw <- Audit_Fairness(Y = dat$Y_miss,
                       S = dat$S,
                       A = dat$A,
                       threshold = threshold,
                       method = "semi-supervised",
                       X = dat %>% select(contains("W")) %>% as.matrix(),
                       basis = "Poly(S, X)")
  
  results[[i]] <- list(
    oracle = oracle,
    sup = sup,
    ss_s = ss_s,
    ss_sw = ss_sw
  )
}
