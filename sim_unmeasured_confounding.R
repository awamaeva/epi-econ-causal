# ======================================================================================
# epi-econ-causal: Simulation study
#   - ATT estimators: Matching, IPTW, G-comp, DR, DiD
#   - LATE estimator: IV
#   - Local effect @ cutoff: RDD (sharp)
# epi-econ-causal  
# Replication and simulation code for the working paper:  
# “A Unified Causal Inference Framework Bridging Epidemiology and Econometrics”  
# Authors: Awa Diop, Ababacar S. Gueye, Cheikh I. Nokho, Kossi D. Abalo, 
# Janick Weberpals,Nima Hejazi, Denis Talbot
# Note:  The data-generating process includes an unmeasured confounder U
# that affects both treatment and outcome. Since U is omitted from all estimation models, 
# the simulation intentionally evaluates estimator performance under unmeasured confounding.
# Moreover, the treatment effect is constant across individuals. 
# As a result, the LATE identified by the IV estimator is numerically close to the ATT. 
# This equivalence arises from the data-generating process and does not hold in general.
# DiD, IV and RDD are generated from separate designs.
# ========================================================================================


suppressPackageStartupMessages({
  library(MASS)
  library(MatchIt)
  library(ivreg)
  library(tidyverse)
  library(sandwich)
  library(lmtest)
})

set.seed(4423)

# ------------------------------
# Global parameters
# ------------------------------
params <- list(
  R = 100,
  n = 5000,
  true_ATT = 2,
  Bg = 200,
  gamma_U = 0.8,   # U -> treatment
  delta_U = 1.0    # U -> outcome level
)

# ------------------------------
# Helpers
# ------------------------------
ci_cover <- function(est, se, true, alpha = 0.05) {
  z <- qnorm(1 - alpha/2)
  ci <- est + c(-1, 1) * z * se
  true >= ci[1] && true <= ci[2]
}

hc_se <- function(fit, term, type = "HC1") {
  sqrt(vcovHC(fit, type = type)[term, term])
}

# ============================================================
# 1) Data Generators
# ============================================================

gen_att_data <- function(n, true_ATT, gamma_U, delta_U) {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2)
  L <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  colnames(L) <- c("L1", "L2")
  
  U <- rnorm(n)
  
  ps_true <- plogis(0.5 * L[,1] - 0.25 * L[,2] + gamma_U * U)
  A <- rbinom(n, 1, ps_true)
  
  Y <- true_ATT * A + 0.8 * L[,1] - 0.5 * L[,2] + delta_U * U + rnorm(n)
  
  tibble(Y = Y, A = A, L1 = L[,1], L2 = L[,2], U = U)
}

gen_did_data <- function(df, true_ATT, delta_U) {
  n <- nrow(df)
  Y0 <- 5 + 0.8*df$L1 - 0.5*df$L2 + delta_U*df$U + rnorm(n)
  Y1 <- Y0 + true_ATT*df$A + rnorm(n)
  tibble(A = df$A, dY = Y1 - Y0)
}

gen_iv_data <- function(n, true_ATT, gamma_U, delta_U) {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2)
  L <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  colnames(L) <- c("L1", "L2")
  U <- rnorm(n)
  
  Z <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(-0.5 + Z + 0.5*L[,1] + gamma_U*U))
  Y <- true_ATT*A + 0.8*L[,1] - 0.5*L[,2] + delta_U*U + rnorm(n)
  
  tibble(Y = Y, A = A, Z = Z, L1 = L[,1], L2 = L[,2])
}

gen_rdd_data <- function(n, true_ATT, delta_U) {
  U <- rnorm(n)
  S <- runif(n, 0, 100)
  A <- as.numeric(S >= 50)
  Y <- true_ATT*A + 0.05*S + delta_U*U + rnorm(n)
  tibble(Y = Y, A = A, S = S)
}

# ============================================================
# 2) Estimators
# ============================================================

est_matching_att <- function(df) {
  m <- matchit(A ~ L1 + L2, data = df, method = "nearest", replace = TRUE)
  md <- match.data(m)
  
  fit <- lm(Y ~ A + L1 + L2, data = md, weights = weights)
  list(est = coef(fit)[["A"]], se = hc_se(fit, "A"))
}

est_iptw_att <- function(df) {
  ps_hat <- predict(glm(A ~ L1 + L2, family = binomial, data = df), type = "response")
  w <- ifelse(df$A == 1, 1, ps_hat / (1 - ps_hat))
  
  fit <- lm(Y ~ A, data = df, weights = w)
  list(est = coef(fit)[["A"]], se = hc_se(fit, "A"))
}

est_gcomp_att <- function(df, Bg = 200) {
  fit <- lm(Y ~ A + L1 + L2, data = df)
  mu1 <- predict(fit, newdata = transform(df, A = 1))
  mu0 <- predict(fit, newdata = transform(df, A = 0))
  est <- mean(mu1[df$A == 1] - mu0[df$A == 1])
  
  boot <- replicate(Bg, {
    idx <- sample.int(nrow(df), replace = TRUE)
    dfb <- df[idx, ]
    fb <- lm(Y ~ A + L1 + L2, data = dfb)
    m1b <- predict(fb, newdata = transform(dfb, A = 1))
    m0b <- predict(fb, newdata = transform(dfb, A = 0))
    mean(m1b[dfb$A == 1] - m0b[dfb$A == 1])
  })
  
  list(est = est, se = sd(boot))
}

est_dr_aiptw_att <- function(df) {
  ps_hat <- predict(glm(A ~ L1 + L2, family = binomial, data = df), type = "response")
  pi1 <- mean(df$A)
  
  m1_hat <- predict(lm(Y ~ L1 + L2, data = df[df$A == 1, ]), newdata = df)
  m0_hat <- predict(lm(Y ~ L1 + L2, data = df[df$A == 0, ]), newdata = df)
  
  est <- mean(
    df$A*(df$Y - m1_hat)/pi1 -
      (1-df$A)*(ps_hat/(1-ps_hat))*(df$Y - m0_hat)/pi1 +
      df$A*(m1_hat - m0_hat)/pi1
  )
  
  phi <- df$A*(df$Y - m1_hat)/pi1 -
    (1-df$A)*(ps_hat/(1-ps_hat))*(df$Y - m0_hat)/pi1 +
    df$A*(m1_hat - m0_hat - est)/pi1
  
  list(est = est, se = sd(phi)/sqrt(nrow(df)))
}

est_did_att <- function(did_df) {
  est <- mean(did_df$dY[did_df$A == 1]) - mean(did_df$dY[did_df$A == 0])
  n1 <- sum(did_df$A == 1); n0 <- sum(did_df$A == 0)
  se <- sqrt(var(did_df$dY[did_df$A == 1]) / n1 + var(did_df$dY[did_df$A == 0]) / n0)
  list(est = est, se = se)
}

est_iv_late <- function(df_iv) {
  fit <- ivreg(Y ~ A + L1 + L2 | Z + L1 + L2, data = df_iv)
  list(est = coef(fit)[["A"]], se = sqrt(vcovHC(fit, type = "HC1")["A","A"]))
}

est_rdd <- function(df_rdd) {
  fit <- lm(Y ~ A + S, data = df_rdd)
  list(est = coef(fit)[["A"]], se = hc_se(fit, "A"))
}

# ============================================================
# 3) One simulation run
# ============================================================

one_sim <- function(seed, p) {
  set.seed(seed)
  
  # ATT world (confounding + unmeasured U)
  df <- gen_att_data(p$n, p$true_ATT, p$gamma_U, p$delta_U)
  
  out <- list(
    Matching = est_matching_att(df),
    IPTW     = est_iptw_att(df),
    `G-comp` = est_gcomp_att(df, Bg = p$Bg),
    DR       = est_dr_aiptw_att(df),
    DiD      = est_did_att(gen_did_data(df, p$true_ATT, p$delta_U)),
    IV       = est_iv_late(gen_iv_data(p$n, p$true_ATT, p$gamma_U, p$delta_U)),
    RDD      = est_rdd(gen_rdd_data(p$n, p$true_ATT, p$delta_U))
  )
  
  tibble(
    method = names(out),
    est    = map_dbl(out, "est"),
    se     = map_dbl(out, "se"),
    truth  = p$true_ATT,                 # see note below
    cover  = map2_lgl(map_dbl(out, "est"), map_dbl(out, "se"), ~ ci_cover(.x, .y, p$true_ATT))
  )
}

# ============================================================
# 4) Run + Summarize
# ============================================================

sim_results <- map_dfr(1:params$R, one_sim, p = params)

summary_tbl <- sim_results %>%
  group_by(method) %>%
  summarise(
    Bias = mean(est - truth),
    RMSE = sqrt(mean((est - truth)^2)),
    Coverage = mean(cover),
    .groups = "drop"
  )

print(summary_tbl)

# ============================================================
# 5) Plot
# ============================================================

sim_results %>%
  mutate(Bias = est - truth) %>%
  ggplot(aes(x = method, y = Bias, fill = method)) +
  geom_boxplot(alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = paste0("Bias Across Estimators (R = ", params$R, ", n = ", params$n, ")"),
    y = "Bias (Estimate - Truth)", x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
# ============================================================
# 5) end
# ============================================================
