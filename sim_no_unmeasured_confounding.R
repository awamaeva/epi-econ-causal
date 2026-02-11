# ======================================================================================
# epi-econ-causal: Simulation study
#   - ATT estimators: Matching, IPTW, G-comp, DR, DiD
#   - LATE estimator: IV
#   - Local effect @ cutoff: RDD (sharp; bandwidth restriction)
# epi-econ-causal
# Replication and simulation code for the working paper:
# “A Unified Causal Inference Framework Bridging Epidemiology and Econometrics”
# Authors: Awa Diop, Ababacar S. Gueye, Cheikh I. Nokho, Kossi D. Abalo,
# Janick Weberpals, Nima Hejazi, Denis Talbot
# Note: This script uses a DGP with NO unmeasured confounding.
# The treatment effect is constant across individuals (true_ATT = 2).
#  DiD, IV and RDD are generated from separate designs.
# ======================================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(MatchIt)
  library(ivreg)
  library(tidyverse)
  library(sandwich)
  library(lmtest)
})

set.seed(123)

# ------------------------------
# Global parameters
# ------------------------------
params <- list(
  R = 1000,
  n = 5000,
  true_ATT = 2,
  Bg = 200
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

gen_att_data_noU <- function(n, true_ATT) {
  Sigma <- matrix(c(1, 0.5,
                    0.5, 1), 2)
  L <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  colnames(L) <- c("L1", "L2")
  
  ps_true <- plogis(0.5 * L[,1] - 0.25 * L[,2])
  A <- rbinom(n, 1, ps_true)
  
  Y <- true_ATT * A + 0.8 * L[,1] - 0.5 * L[,2] + rnorm(n)
  
  tibble(Y = Y, A = A, L1 = L[,1], L2 = L[,2])
}

gen_did_data_noU <- function(df, true_ATT) {
  n <- nrow(df)
  Y0 <- 1.5 + 0.8 * df$L1 - 0.5 * df$L2 + rnorm(n)
  Y1 <- Y0 + true_ATT * df$A + rnorm(n)
  tibble(A = df$A, dY = Y1 - Y0)
}

gen_iv_data_noU <- function(n, true_ATT) {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2)
  L <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  colnames(L) <- c("L1", "L2")
  
  Z <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(-0.5 + 1 * Z + 0.5 * L[,1]))
  Y <- true_ATT * A + 0.8 * L[,1] - 0.5 * L[,2] + rnorm(n)
  
  tibble(Y = Y, A = A, Z = Z, L1 = L[,1], L2 = L[,2])
}

gen_rdd_data_noU <- function(n, true_ATT) {
  S <- runif(n, 0, 100)
  cutoff <- 50
  A <- as.numeric(S >= cutoff)
  Y <- true_ATT * A + 0.05 * S + rnorm(n)
  tibble(Y = Y, A = A, S = S, cutoff = cutoff)
}

# ============================================================
# 2) Estimators
# ============================================================

est_matching_att <- function(df) {
  m <- matchit(
    A ~ L1 + L2,
    data = df,
    method = "nearest",
    replace = TRUE,
    caliper = 0.2
  )
  md <- match.data(m)
  
  fit <- lm(Y ~ A + L1 + L2, data = md, weights = weights)
  list(est = coef(fit)[["A"]], se = hc_se(fit, "A"))
}

est_iptw_att <- function(df) {
  ps_fit <- glm(A ~ L1 + L2, family = binomial, data = df)
  ps_hat <- predict(ps_fit, type = "response")
  
  w <- ifelse(df$A == 1, 1, ps_hat / (1 - ps_hat))
  fit <- lm(Y ~ A, data = df, weights = w)
  
  list(est = coef(fit)[["A"]], se = hc_se(fit, "A"))
}

est_gcomp_att <- function(df, Bg = 200) {
  fit <- lm(Y ~ A + L1 + L2, data = df)
  mu1 <- predict(fit, newdata = transform(df, A = 1))
  mu0 <- predict(fit, newdata = transform(df, A = 0))
  
  est <- mean((mu1 - mu0)[df$A == 1])
  
  boot <- replicate(Bg, {
    idx <- sample.int(nrow(df), replace = TRUE)
    dfb <- df[idx, ]
    
    fb <- lm(Y ~ A + L1 + L2, data = dfb)
    m1b <- predict(fb, newdata = transform(dfb, A = 1))
    m0b <- predict(fb, newdata = transform(dfb, A = 0))
    
    mean((m1b - m0b)[dfb$A == 1])
  })
  
  list(est = est, se = sd(boot))
}

est_dr_aiptw_att <- function(df) {
  ps_hat <- predict(glm(A ~ L1 + L2, family = binomial, data = df), type = "response")
  pi1 <- mean(df$A)
  
  m1_hat <- predict(lm(Y ~ L1 + L2, data = df[df$A == 1,]), newdata = df)
  m0_hat <- predict(lm(Y ~ L1 + L2, data = df[df$A == 0,]), newdata = df)
  
  est <- mean(
    df$A * (df$Y - m1_hat) / pi1 -
      (1 - df$A) * (ps_hat / (1 - ps_hat)) * (df$Y - m0_hat) / pi1 +
      df$A * (m1_hat - m0_hat) / pi1
  )
  
  phi <- df$A * (df$Y - m1_hat) / pi1 -
    (1 - df$A) * (ps_hat / (1 - ps_hat)) * (df$Y - m0_hat) / pi1 +
    df$A * (m1_hat - m0_hat - est) / pi1
  
  list(est = est, se = sd(phi) / sqrt(nrow(df)))
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

est_rdd <- function(df_rdd, bw = 10) {
  cutoff <- unique(df_rdd$cutoff)
  df_bw <- df_rdd %>% filter(abs(S - cutoff) <= bw)
  
  fit <- lm(Y ~ A + S, data = df_bw)
  list(est = coef(fit)[["A"]], se = hc_se(fit, "A"))
}

# ============================================================
# 3) One simulation run
# ============================================================

one_sim <- function(seed, p) {
  set.seed(seed)
  
  # ATT world (no unmeasured confounding)
  df <- gen_att_data_noU(p$n, p$true_ATT)
  
  out <- list(
    Matching = est_matching_att(df),
    IPTW     = est_iptw_att(df),
    `G-comp` = est_gcomp_att(df, Bg = p$Bg),
    DR       = est_dr_aiptw_att(df),
    DiD      = est_did_att(gen_did_data_noU(df, p$true_ATT)),
    IV       = est_iv_late(gen_iv_data_noU(p$n, p$true_ATT)),
    RDD      = est_rdd(gen_rdd_data_noU(p$n, p$true_ATT), bw = 10)
  )
  
  tibble(
    method = names(out),
    est    = map_dbl(out, "est"),
    se     = map_dbl(out, "se"),
    truth  = p$true_ATT,
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
    Bias     = mean(est - truth),
    RMSE     = sqrt(mean((est - truth)^2)),
    Coverage = mean(cover),
    .groups  = "drop"
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
    title = paste0("Bias Across Estimators (No Unmeasured Confounding) (R = ", params$R, ", n = ", params$n, ")"),
    y = "Bias (Estimate - Truth)",
    x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# ============================================================
# 5) end
# ============================================================
