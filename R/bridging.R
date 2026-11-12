# Companion-diagnostic bridging simulation (Chapter 8). All functions are tested in tests/testthat.

# Simulate a randomised all-comers trial with a clinical trial assay (CTA) measured in everyone
# and a final CDx measured only where a bridging sample is available.
# miss: "none", "MCAR", "MAR" (depends on CTA), "MNAR" (depends on the CDx result itself),
#       "MNAR_resp" (depends on the CDx result and on the trial outcome)
sim_trial <- function(n = 1200, prev = 0.40,
                      cta = c(se = 0.93, sp = 0.90), cdx = c(se = 0.97, sp = 0.95),
                      p_ctrl = 0.20, p_trt_pos = 0.50, p_trt_neg = 0.20,
                      miss = c("MCAR", "MAR", "MNAR", "MNAR_resp", "none"), p_obs = 0.60) {
  miss <- match.arg(miss)
  B <- rbinom(n, 1, prev)
  arm <- rbinom(n, 1, 0.5)
  cta_r <- rbinom(n, 1, ifelse(B == 1, cta[["se"]], 1 - cta[["sp"]]))
  cdx_r <- rbinom(n, 1, ifelse(B == 1, cdx[["se"]], 1 - cdx[["sp"]]))
  p_resp <- ifelse(arm == 0, p_ctrl, ifelse(B == 1, p_trt_pos, p_trt_neg))
  resp <- rbinom(n, 1, p_resp)
  p_o <- switch(miss,
    none = rep(1, n),
    MCAR = rep(p_obs, n),
    MAR  = ifelse(cta_r == 1, 0.75, 0.35),
    MNAR = ifelse(cdx_r == 1, 0.80, 0.40),
    MNAR_resp = ifelse(cdx_r == 1, ifelse(resp == 1, 0.85, 0.50), 0.40))
  data.frame(arm, B, cta = cta_r, cdx = cdx_r, resp, obs = rbinom(n, 1, p_o))
}

# Risk difference (treatment minus control) with Wald variance; optional weights.
rd_est <- function(arm, resp, w = NULL) {
  if (is.null(w)) w <- rep(1, length(arm))
  m <- function(a) { i <- arm == a; sum(w[i] * resp[i]) / sum(w[i]) }
  neff <- function(a) { i <- arm == a; sum(w[i])^2 / sum(w[i]^2) }
  p1 <- m(1); p0 <- m(0)
  c(est = p1 - p0, se = sqrt(p1 * (1 - p1) / neff(1) + p0 * (1 - p0) / neff(0)))
}

# True effect in the population who test CDx+ (analytic)
true_effect_cdx_pos <- function(prev = 0.40, cdx = c(se = 0.97, sp = 0.95),
                                p_ctrl = 0.20, p_trt_pos = 0.50, p_trt_neg = 0.20) {
  pi1 <- cdx[["se"]] * prev / (cdx[["se"]] * prev + (1 - cdx[["sp"]]) * (1 - prev))
  pi1 * (p_trt_pos - p_ctrl) + (1 - pi1) * (p_trt_neg - p_ctrl)
}

# Estimators of the effect among CDx+ subjects --------------------------------
est_full <- function(d) rd_est(d$arm[d$cdx == 1], d$resp[d$cdx == 1])          # if CDx were known for all
est_cc   <- function(d) { k <- d$obs == 1 & d$cdx == 1; rd_est(d$arm[k], d$resp[k]) }
est_ipw  <- function(d) {
  w <- if (all(d$obs == 1)) rep(1, nrow(d)) else 1 / fitted(glm(obs ~ cta + arm, binomial, d))
  k <- d$obs == 1 & d$cdx == 1
  rd_est(d$arm[k], d$resp[k], w[k])
}
est_mi <- function(d, M = 20) {
  o <- d[d$obs == 1, ]
  f <- glm(cdx ~ cta + arm * resp, binomial, o)   # arm x resp: effect differs by biomarker status
  b <- coef(f); V <- vcov(f)
  X <- model.matrix(~ cta + arm * resp, d)
  res <- t(sapply(seq_len(M), function(m) {
    bm <- MASS::mvrnorm(1, b, V)
    p <- plogis(X %*% bm)
    cd <- ifelse(d$obs == 1, d$cdx, rbinom(nrow(d), 1, p))
    k <- cd == 1
    rd_est(d$arm[k], d$resp[k])
  }))
  qbar <- mean(res[, "est"]); ubar <- mean(res[, "se"]^2); bm <- var(res[, "est"])
  c(est = qbar, se = sqrt(ubar + (1 + 1 / M) * bm))
}

# Delta-adjusted tipping-point analysis (pattern mixture). The CDx result of subjects
# without a sample is imputed from the MAR model of est_mi(), with a logit offset
#   delta * arm * (2 * resp - 1)
# added for the missing subjects only. delta < 0 says that, among subjects without a
# sample, treated responders are less likely and treated non-responders more likely to be
# CDx+ than the MAR model predicts, which moves the CDx+ effect towards zero. delta = 0
# reproduces est_mi(). The tipping point is the delta at which the conclusion changes.
est_tipping_delta <- function(d, delta, M = 30) {
  o <- d[d$obs == 1, ]
  f <- glm(cdx ~ cta + arm * resp, binomial, o)
  b <- coef(f); V <- vcov(f)
  X <- model.matrix(~ cta + arm * resp, d)
  shift <- delta * d$arm * (2 * d$resp - 1)
  res <- t(sapply(seq_len(M), function(m) {
    bm <- MASS::mvrnorm(1, b, V)
    p <- plogis(X %*% bm + shift)
    cd <- ifelse(d$obs == 1, d$cdx, rbinom(nrow(d), 1, p))
    k <- cd == 1
    rd_est(d$arm[k], d$resp[k])
  }))
  c(est = mean(res[, "est"]), se = sqrt(mean(res[, "se"]^2) + (1 + 1 / M) * var(res[, "est"])))
}
