# -----------------------------------------------------------------------------
# Program    : R/helpers.R
# Purpose    : Shared statistical helpers sourced by the chapters and the mock SAP:
#              interval methods, EP05 ANOVA, EP17 LoB/LoD, method comparison, ROC,
#              two-cut-off derivation and performance, sample size, sources table
# Inputs     : none (pure functions); show_sources() reads claims_ledger.csv
# Outputs    : none (functions only)
# Author     : Ignacio G. Ribelles
# Created    : 2026-10-16
# Change log : 2026-10-16  IGR  Initial version
#              2026-11-12  IGR  EP17 estimators: pooled SD by sample, small-sample
#                               multiplier, rank rule 0.5 + 0.95 B (quantile type 5),
#                               nonparametric LoD; guard in derive_two_cutoffs();
#                               show_sources() and sim_ep05() added
# Notes      : Every function is unit-tested in tests/testthat/test-helpers.R.
#              deming_hand() takes delta = var(err_y) / var(err_x), the inverse of
#              mcr's error.ratio.
# -----------------------------------------------------------------------------
if (!exists("%||%", baseenv())) `%||%` <- function(x, y) if (is.null(x)) y else x  # R < 4.4 compat (lme4 printing)

# ---- Proportion confidence intervals -------------------------------------
ci_wilson <- function(x, n, conf = 0.95) {
  stopifnot(n > 0, x >= 0, x <= n)
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(est = p, lower = max(0, centre - half), upper = min(1, centre + half))
}

ci_clopper_pearson <- function(x, n, conf = 0.95) {
  stopifnot(n > 0, x >= 0, x <= n)
  a <- 1 - conf
  lo <- if (x == 0) 0 else qbeta(a / 2, x, n - x + 1)
  hi <- if (x == n) 1 else qbeta(1 - a / 2, x + 1, n - x)
  c(est = x / n, lower = lo, upper = hi)
}

# ---- 2x2 diagnostic report ------------------------------------------------
# truth and test are logical/0-1 vectors (1 = condition present / test positive)
report_2x2 <- function(test, truth, conf = 0.95, ci = ci_wilson) {
  test <- as.logical(test); truth <- as.logical(truth)
  tp <- sum(test & truth); fn <- sum(!test & truth)
  fp <- sum(test & !truth); tn <- sum(!test & !truth)
  row <- function(name, x, n) {
    r <- ci(x, n, conf)
    data.frame(measure = name, x = x, n = n, est = r[["est"]],
               lower = r[["lower"]], upper = r[["upper"]])
  }
  rbind(row("Sensitivity", tp, tp + fn),
        row("Specificity", tn, tn + fp),
        row("PPV", tp, tp + fp),
        row("NPV", tn, tn + fn))
}

# PPV/NPV at a chosen prevalence from sensitivity and specificity (Bayes)
predictive_values <- function(sens, spec, prev) {
  ppv <- sens * prev / (sens * prev + (1 - spec) * (1 - prev))
  npv <- spec * (1 - prev) / (spec * (1 - prev) + (1 - sens) * prev)
  c(PPV = ppv, NPV = npv, LRpos = sens / (1 - spec), LRneg = (1 - sens) / spec)
}

# ---- EP05-style nested precision (balanced) -------------------------------
# Balanced design: day (D) > run (R) > replicate (K). Returns variance components
# by the classical ANOVA method of moments (what VCA::anovaVCA also does).
ep05_anova <- function(y, day, run, nrep) {
  d <- data.frame(y = y, day = factor(day), run = factor(paste(day, run)))
  D <- nlevels(d$day); R <- nlevels(d$run) / D; K <- nrep
  grand <- mean(d$y)
  ybar_day <- tapply(d$y, d$day, mean)
  ybar_run <- tapply(d$y, d$run, mean)
  ss_day <- R * K * sum((ybar_day - grand)^2)
  run_day <- sapply(levels(d$run), function(r) sub(" .*", "", r))
  ss_run <- K * sum((ybar_run - ybar_day[run_day])^2)
  ss_err <- sum((d$y - ybar_run[as.character(d$run)])^2)
  ms_day <- ss_day / (D - 1)
  ms_run <- ss_run / (D * (R - 1))
  ms_err <- ss_err / (D * R * (K - 1))
  v_rep <- ms_err
  v_run <- max(0, (ms_run - ms_err) / K)
  v_day <- max(0, (ms_day - ms_run) / (R * K))
  c(day = v_day, run = v_run, repeatability = v_rep,
    within_lab = v_day + v_run + v_rep)
}

# ---- EP17-style detection capability --------------------------------------
lob_param <- function(blank, z = 1.645) mean(blank) + z * sd(blank)

# Pooled SD of low-level results. With `group` (one id per low-level sample) the SD is
# pooled within samples, so differences between sample means do not inflate it.
sd_pooled <- function(y, group = NULL) {
  if (is.null(group)) return(sd(y))
  g <- split(y, group)
  dfs <- lengths(g) - 1
  sqrt(sum(dfs * vapply(g, var, numeric(1))) / sum(dfs))
}

# Parametric LoD = LoB + c * SD_low. `cp = TRUE` uses the small-sample multiplier
# 1.645 / (1 - 1 / (4 * (L - J))), L results in J low-level samples, as described in
# secondary literature on EP17-A2; `cp = FALSE` uses the plain 1.645.
lod_param <- function(lob, low, group = NULL, z = 1.645, cp = FALSE) {
  mult <- z
  if (cp) {
    J <- if (is.null(group)) 1 else length(unique(group))
    mult <- z / (1 - 1 / (4 * (length(low) - J)))
  }
  lob + mult * sd_pooled(low, group)
}

# Rank-based 95th percentile of the blanks: rank position 0.5 + 0.95 * B with linear
# interpolation, which is quantile type 5.
lob_nonparam <- function(blank, alpha = 0.05) {
  unname(quantile(blank, 1 - alpha, type = 5))
}

# Nonparametric LoD: LoB + (median - 5th percentile) of the low-level results, each
# sample centred at its own median before pooling.
lod_nonparam <- function(lob, low, group = NULL, beta = 0.05) {
  centred <- if (is.null(group)) low - median(low) else low - ave(low, group, FUN = median)
  lob + (0 - unname(quantile(centred, beta, type = 5)))
}

# ---- Sample size for sensitivity / specificity (Buderer 1996) --------------
# Total n so that the CI for sens (or spec) has half-width `d`, given prevalence.
n_buderer <- function(p, d, prev, conf = 0.95, what = c("sens", "spec")) {
  what <- match.arg(what)
  z <- qnorm(1 - (1 - conf) / 2)
  n_cases <- z^2 * p * (1 - p) / d^2
  if (what == "sens") ceiling(n_cases / prev) else ceiling(n_cases / (1 - prev))
}

# ---- Two-cut-off classification --------------------------------------------
# Returns factor: "rule-out" (x <= lo), "intermediate" (lo < x <= hi), "rule-in" (x > hi)
classify_two_cut <- function(x, lo, hi) {
  stopifnot(lo <= hi)
  factor(ifelse(x <= lo, "rule-out", ifelse(x > hi, "rule-in", "intermediate")),
         levels = c("rule-out", "intermediate", "rule-in"))
}

# ---- Sources ----------------------------------------------------------------
# Print the sourced claims of one chapter from claims_ledger.csv as a markdown table
# (use inside a chunk with results: asis).
show_sources <- function(chapter, path = "claims_ledger.csv") {
  l <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("claim", "source", "verified_on", "confidence", "chapter") %in% names(l)))
  l <- l[l$chapter == chapter, c("claim", "source", "verified_on", "confidence")]
  if (!nrow(l)) return(invisible(cat("No sourced claims for this chapter.\n")))
  link <- ifelse(grepl("^https?://", l$source), sprintf("<%s>", l$source), l$source)
  cat("| Claim | Source | Checked | Confidence |\n|---|---|---|---|\n")
  cat(sprintf("| %s | %s | %s | %s |\n", gsub("\\|", "/", l$claim), gsub("\\|", "/", link),
              l$verified_on, l$confidence), sep = "")
  invisible(l)
}

# ---- Simulated EP05-style precision study (shared by Chapters 2 and 7) --------
sim_ep05 <- function(seed = 20, D = 20, R = 2, K = 2, mu = 100,
                     sd_day = 1.2, sd_run = 1.5, sd_rep = 2) {
  set.seed(seed)
  d <- expand.grid(rep = 1:K, run = 1:R, day = 1:D)
  d$y <- mu + rnorm(D, 0, sd_day)[d$day] +
    rnorm(D * R, 0, sd_run)[(d$day - 1) * R + d$run] + rnorm(nrow(d), 0, sd_rep)
  d$day <- factor(d$day); d$run <- factor(d$run)
  d
}

# ---- Probit LoD (hit-rate design) ------------------------------------------
# hits/n at each concentration; returns the concentration with `target` detection prob.
lod_probit <- function(conc, hits, n, target = 0.95) {
  fit <- glm(cbind(hits, n - hits) ~ log10(conc), family = binomial(link = "probit"))
  b <- coef(fit)
  10^((qnorm(target) - b[[1]]) / b[[2]])
}

# ---- Linearity: deviation of the linear fit from a polynomial fit ----------
# Returns F-test p-value for the quadratic+cubic terms and %deviation at each level.
linearity_check <- function(x, y, degree = 3) {
  lin <- lm(y ~ x)
  poly_fit <- lm(y ~ poly(x, degree, raw = TRUE))
  a <- anova(lin, poly_fit)
  lv <- sort(unique(x))
  dev <- 100 * (predict(lin, data.frame(x = lv)) - predict(poly_fit, data.frame(x = lv))) /
    predict(poly_fit, data.frame(x = lv))
  list(p_nonlinear = a$`Pr(>F)`[2], F = a$F[2], deviation_pct = setNames(dev, lv))
}

# ---- Method comparison: hand implementations used as cross-checks ----------
# delta = var(error_y) / var(error_x); delta = 1 is orthogonal regression.
deming_hand <- function(x, y, delta = 1) {
  sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
  b <- (syy - delta * sxx + sqrt((syy - delta * sxx)^2 + 4 * delta * sxy^2)) / (2 * sxy)
  c(intercept = mean(y) - b * mean(x), slope = b)
}

# Passing-Bablok point estimates (shifted median of pairwise slopes)
paba_hand <- function(x, y) {
  n <- length(x)
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  dx <- x[idx[, 2]] - x[idx[, 1]]; dy <- y[idx[, 2]] - y[idx[, 1]]
  keep <- !(dx == 0 & dy == 0)
  s <- ifelse(dx == 0, ifelse(dy > 0, Inf, -Inf), dy / dx)[keep]
  s <- s[s != -1]                       # slopes of exactly -1 are discarded
  N <- length(s); K <- sum(s < -1)
  s <- sort(s)
  b <- if (N %% 2 == 1) s[(N + 1) / 2 + K] else mean(s[N / 2 + K + 0:1])
  c(intercept = median(y - b * x), slope = b)
}

# ---- ROC: DeLong AUC and variance by hand (cross-check for pROC) -----------
auc_delong <- function(case, control) {
  psi <- function(a, b) (a > b) + 0.5 * (a == b)
  m <- outer(case, control, psi)
  v10 <- rowMeans(m); v01 <- colMeans(m)
  auc <- mean(v10)
  var_auc <- var(v10) / length(case) + var(v01) / length(control)
  c(auc = auc, se = sqrt(var_auc))
}

# ---- Agreement with a non-reference comparator: PPA / NPA / OPA -------------
# FDA 2007 guidance: do NOT call these sensitivity / specificity.
report_agreement <- function(candidate, comparator, conf = 0.95, ci = ci_wilson) {
  r <- report_2x2(candidate, comparator, conf, ci)
  tp <- sum(candidate & comparator); tn <- sum(!candidate & !comparator)
  n <- length(candidate)
  out <- r[r$measure %in% c("Sensitivity", "Specificity"), ]
  out$measure <- c("PPA", "NPA")
  opa <- ci(tp + tn, n, conf)
  rbind(out, data.frame(measure = "OPA", x = tp + tn, n = n, est = opa[["est"]],
                        lower = opa[["lower"]], upper = opa[["upper"]]))
}

# ---- ROC figure from a list of pROC objects --------------------------------
roc_figure <- function(rocs) {
  stopifnot(is.list(rocs), !is.null(names(rocs)))
  d <- do.call(rbind, lapply(names(rocs), function(nm) {
    r <- rocs[[nm]]
    data.frame(curve = sprintf("%s (AUC %.3f)", nm, as.numeric(pROC::auc(r))),
               fpr = 1 - r$specificities, tpr = r$sensitivities)
  }))
  ggplot2::ggplot(d, ggplot2::aes(fpr, tpr, colour = curve)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_abline(linetype = 3, colour = "grey50") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "1 - specificity", y = "Sensitivity", colour = NULL) +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
}

# ---- Two cut-off derivation on a derivation cohort -------------------------
# Convention (as in the FDA decision summary): negative if x <= lo, positive if x > hi.
# lo: largest cut-off for which sensitivity of "not negative" is >= se_target
# hi: smallest cut-off for which specificity of "not positive" is >= sp_target
derive_two_cutoffs <- function(x, truth, se_target, sp_target) {
  truth <- as.logical(truth)
  cases <- sort(x[truth]); ctrls <- sort(x[!truth])
  ok_lo <- cases[vapply(cases, function(thr) mean(cases > thr) >= se_target, logical(1))]
  ok_hi <- ctrls[vapply(ctrls, function(thr) mean(ctrls <= thr) >= sp_target, logical(1))]
  if (!length(ok_lo) || !length(ok_hi))
    stop("No observed cut-off meets the sensitivity or specificity target")
  lo <- max(ok_lo); hi <- min(ok_hi)
  if (lo > hi) {
    # The rule-out cut-off lies above the rule-in cut-off: both targets cannot hold with an
    # intermediate zone. Return a single cut-off and say so; the targets are then not both met.
    warning("Targets overlap: returning a single cut-off; check achieved Se and Sp")
    m <- (lo + hi) / 2; lo <- hi <- m
  }
  c(lo = lo, hi = hi)
}

# Performance of a locked two-cut-off scheme
two_cut_performance <- function(x, truth, lo, hi, conf = 0.95) {
  truth <- as.logical(truth)
  cat3 <- classify_two_cut(x, lo, hi)
  tab <- table(cat3, truth)[, c("TRUE", "FALSE")]
  colnames(tab) <- c("D+", "D-")
  n <- length(x)
  list(
    table = tab,
    intermediate = c(ci_wilson(sum(cat3 == "intermediate"), n, conf)),
    ruleout_sens = ci_wilson(sum(tab[c("intermediate", "rule-in"), "D+"]), sum(tab[, "D+"]), conf),
    ruleout_npv  = ci_wilson(tab["rule-out", "D-"], sum(tab["rule-out", ]), conf),
    rulein_spec  = ci_wilson(sum(tab[c("rule-out", "intermediate"), "D-"]), sum(tab[, "D-"]), conf),
    rulein_ppv   = ci_wilson(tab["rule-in", "D+"], sum(tab["rule-in", ]), conf),
    # "including indeterminates": intermediates counted as wrong, like PPA/NPA in the FDA summary
    ppa_incl = ci_wilson(tab["rule-in", "D+"], sum(tab[, "D+"]), conf),
    npa_incl = ci_wilson(tab["rule-out", "D-"], sum(tab[, "D-"]), conf)
  )
}

# ---- Sample size: performance goal with an exact binomial test --------------
# Smallest number of subjects n (cases for sensitivity) such that, if the true value
# is p1, the one-sided exact test of H0: p <= p0 rejects with probability >= power.
# Sawtooth in n: returns the smallest n that reaches the power (see also n_exact_goal_stable).
power_exact_goal <- function(n, p1, p0, alpha = 0.025) {
  crit <- qbinom(1 - alpha, n, p0)          # reject if X > crit  (P(X > crit | p0) <= alpha)
  while (pbinom(crit, n, p0, lower.tail = FALSE) > alpha) crit <- crit + 1
  pbinom(crit, n, p1, lower.tail = FALSE)
}
n_exact_goal <- function(p1, p0, alpha = 0.025, power = 0.80, nmax = 5000) {
  for (n in 5:nmax) if (power_exact_goal(n, p1, p0, alpha) >= power) return(n)
  NA_integer_
}
# Because of the sawtooth, require the power to hold for all larger n up to nmax.
n_exact_goal_stable <- function(p1, p0, alpha = 0.025, power = 0.80, nmax = 1500) {
  ok <- sapply(5:nmax, function(n) power_exact_goal(n, p1, p0, alpha) >= power)
  last_bad <- if (any(!ok)) max(which(!ok)) else 0
  (5:nmax)[last_bad + 1]
}

# ---- Sample size: AUC against a null value (Hanley-McNeil variance) ---------
var_auc_hm <- function(auc, n1, n0) {
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  (auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)) / (n1 * n0)
}
# n1 = cases; controls = ratio * cases. Two-sided alpha, one-sample Z test of AUC = auc0.
n_auc_hm <- function(auc1, auc0 = 0.5, alpha = 0.05, power = 0.80, ratio = 1) {
  za <- qnorm(1 - alpha / 2); zb <- qnorm(power)
  f <- function(n1) {
    n0 <- ratio * n1
    (auc1 - auc0) / (za * sqrt(var_auc_hm(auc0, n1, n0)) + zb * sqrt(var_auc_hm(auc1, n1, n0))) - 1
  }
  ceiling(uniroot(f, c(3, 1e5))$root)
}
