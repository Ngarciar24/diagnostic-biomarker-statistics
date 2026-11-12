source(file.path(rprojroot_dir <- normalizePath(file.path("..", "..")), "R", "helpers.R"))

test_that("Wilson and Clopper-Pearson agree with binom package", {
  skip_if_not_installed("binom")
  w <- binom::binom.confint(37, 50, methods = "wilson")
  cp <- binom::binom.confint(37, 50, methods = "exact")
  expect_equal(unname(ci_wilson(37, 50)[2:3]), c(w$lower, w$upper), tolerance = 1e-8)
  expect_equal(unname(ci_clopper_pearson(37, 50)[2:3]), c(cp$lower, cp$upper), tolerance = 1e-8)
})

test_that("CP handles boundaries", {
  expect_equal(ci_clopper_pearson(0, 20)[["lower"]], 0)
  expect_equal(ci_clopper_pearson(20, 20)[["upper"]], 1)
})

test_that("report_2x2 counts are right", {
  truth <- c(rep(1, 10), rep(0, 20)); test <- c(rep(1, 8), rep(0, 2), rep(1, 3), rep(0, 17))
  r <- report_2x2(test, truth)
  expect_equal(r$est[r$measure == "Sensitivity"], 0.8)
  expect_equal(r$est[r$measure == "Specificity"], 17 / 20)
  expect_equal(r$est[r$measure == "PPV"], 8 / 11)
})

test_that("predictive values match hand calculation", {
  pv <- predictive_values(0.9, 0.8, 0.1)
  expect_equal(pv[["PPV"]], 0.09 / (0.09 + 0.18))
  expect_equal(pv[["LRpos"]], 4.5)
})

test_that("ep05_anova recovers truth on a large balanced simulation", {
  set.seed(1)
  D <- 400; R <- 2; K <- 2
  d <- expand.grid(rep = 1:K, run = 1:R, day = 1:D)
  d$y <- 100 + rnorm(D, 0, 3)[d$day] + rnorm(D * R, 0, 2)[(d$day - 1) * R + d$run] + rnorm(nrow(d), 0, 1.5)
  v <- ep05_anova(d$y, d$day, d$run, K)
  expect_equal(unname(v[c("day", "run", "repeatability")]), c(9, 4, 2.25), tolerance = 0.15)
})

test_that("ep05_anova matches VCA::anovaVCA", {
  skip_if_not_installed("VCA")
  set.seed(2)
  D <- 20; R <- 2; K <- 2
  d <- expand.grid(rep = 1:K, run = 1:R, day = 1:D)
  d$y <- 50 + rnorm(D, 0, 2)[d$day] + rnorm(D * R, 0, 1)[(d$day - 1) * R + d$run] + rnorm(nrow(d))
  d$day <- factor(d$day); d$run <- factor(d$run)
  fit <- VCA::anovaVCA(y ~ day/run, d)
  vca <- fit$aov.tab[c("day", "day:run", "error"), "VC"]
  mine <- ep05_anova(d$y, d$day, d$run, K)
  expect_equal(unname(mine[c("day", "run", "repeatability")]), unname(vca), tolerance = 1e-6)
})

test_that("LoB/LoD parametric formulas", {
  b <- c(1, 2, 3, 4, 5); lob <- lob_param(b)
  expect_equal(lob, 3 + 1.645 * sd(b))
  expect_equal(lod_param(lob, c(10, 12, 14)), lob + 1.645 * 2)
})

test_that("Buderer sample size matches hand calc", {
  # sens 0.9, d 0.05, prev 0.2: 1.96^2*.09/.0025 = 138.3 cases -> /0.2 = 691.5 -> 692
  expect_equal(n_buderer(0.9, 0.05, 0.2), 692)
  expect_equal(n_buderer(0.9, 0.05, 0.2, what = "spec"), ceiling(138.2976 / 0.8))
})

test_that("two-cut classification", {
  f <- classify_two_cut(c(1, 5, 9), 2, 8)
  expect_equal(as.character(f), c("rule-out", "intermediate", "rule-in"))
  expect_error(classify_two_cut(1, 5, 2))
})

test_that("lod_probit recovers the 95% point and matches MASS::dose.p", {
  conc <- c(2, 5, 10, 20, 40); b0 <- -2; b1 <- 2   # true model: probit(p) = -2 + 2*log10(conc)
  p <- pnorm(b0 + b1 * log10(conc)); n <- rep(1e6, 5); hits <- round(p * n)
  true_lod <- 10^((qnorm(0.95) - b0) / b1)
  expect_equal(lod_probit(conc, hits, n), true_lod, tolerance = 1e-3)
  skip_if_not_installed("MASS")
  fit <- glm(cbind(hits, n - hits) ~ log10(conc), family = binomial("probit"))
  dp <- MASS::dose.p(fit, p = 0.95)
  expect_equal(10^as.numeric(dp), lod_probit(conc, hits, n), tolerance = 1e-8)
})

test_that("linearity_check flags curvature and not a straight line", {
  x <- rep(1:9, each = 3)
  set.seed(3)
  expect_gt(linearity_check(x, 2 * x + rnorm(27, 0, 0.3))$p_nonlinear, 0.01)
  expect_lt(linearity_check(x, 2 * x - 0.08 * x^2 + rnorm(27, 0, 0.1))$p_nonlinear, 1e-6)
})

test_that("Deming hand estimate matches mcr and recovers truth without noise bias", {
  set.seed(4); t <- runif(200, 10, 100)
  x <- t + rnorm(200, 0, 2); y <- 3 + 1.1 * t + rnorm(200, 0, 2)
  h <- deming_hand(x, y, 1)
  expect_equal(unname(h["slope"]), 1.1, tolerance = 0.05)
  skip_if_not_installed("mcr")
  m <- mcr::mcreg(x, y, method.reg = "Deming", error.ratio = 1, method.ci = "analytical")
  expect_equal(unname(h), unname(m@para[, "EST"]), tolerance = 1e-6)
})

test_that("Passing-Bablok hand estimate matches mcr", {
  set.seed(5); x <- runif(60, 5, 80); y <- 1 + 1.05 * x + rnorm(60, 0, 2)
  skip_if_not_installed("mcr")
  m <- mcr::mcreg(x, y, method.reg = "PaBa", method.ci = "analytical")
  expect_equal(unname(paba_hand(x, y)), unname(m@para[, "EST"]), tolerance = 1e-6)
})

test_that("auc_delong matches pROC (AUC and DeLong SE)", {
  skip_if_not_installed("pROC")
  set.seed(6); case <- rnorm(40, 1); ctrl <- rnorm(60)
  r <- pROC::roc(response = rep(c(1, 0), c(40, 60)), predictor = c(case, ctrl), quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  h <- auc_delong(case, ctrl)
  expect_equal(unname(h["auc"]), as.numeric(pROC::auc(r)), tolerance = 1e-10)
  expect_equal(unname(h["auc"] + 1.96 * h["se"]), min(ci[3], 1) , tolerance = 5e-3)
})

test_that("auc equals Mann-Whitney U / (n1 n0)", {
  set.seed(7); a <- rnorm(30, .8); b <- rnorm(50)
  w <- unname(wilcox.test(a, b)$statistic) / (30 * 50)
  expect_equal(unname(auc_delong(a, b)["auc"]), w)
})

test_that("report_agreement labels PPA/NPA/OPA", {
  r <- report_agreement(c(1, 1, 0, 0, 1), c(1, 0, 0, 0, 1))
  expect_equal(r$measure, c("PPA", "NPA", "OPA"))
  expect_equal(r$est, c(1, 2 / 3, 0.8))
})

test_that("roc_figure returns a ggplot", {
  skip_if_not_installed("pROC")
  set.seed(8); y <- rep(0:1, 50); x <- rnorm(100) + y
  r <- pROC::roc(y, x, quiet = TRUE)
  expect_s3_class(roc_figure(list(test = r)), "ggplot")
})

test_that("derive_two_cutoffs hits the sensitivity/specificity targets on the derivation data", {
  set.seed(9); y <- rbinom(600, 1, .4); x <- rnorm(600, y * 1.5)
  co <- derive_two_cutoffs(x, y, 0.95, 0.95)
  expect_gte(mean(x[y == 1] > co[["lo"]]), 0.95)
  expect_gte(mean(x[y == 0] <= co[["hi"]]), 0.95)
  expect_lt(co[["lo"]], co[["hi"]])
  # tightest: next higher lo value would break the target
  nxt <- min(x[y == 1][x[y == 1] > co[["lo"]]])
  expect_lt(mean(x[y == 1] > nxt), 0.95)
})

test_that("two_cut_performance counts add up", {
  set.seed(10); y <- rbinom(300, 1, .4); x <- rnorm(300, y * 1.5)
  p <- two_cut_performance(x, y, -0.2, 1.2)
  expect_equal(sum(p$table), 300)
  expect_equal(unname(p$intermediate["est"]), sum(x > -0.2 & x <= 1.2) / 300)
})

test_that("exact goal sample size (stable) matches MKmisc::power.diagnostic.test", {
  skip_if_not_installed("MKmisc")
  for (a in c(0.05, 0.025)) {
    r <- MKmisc::power.diagnostic.test(sens = 0.9, delta = 0.1, power = 0.8, sig.level = a, prev = 0.2)
    expect_equal(n_exact_goal_stable(0.9, 0.8, alpha = a, power = 0.8), r$n)
  }
})

test_that("power_exact_goal is monotone-ish and bounded", {
  expect_gt(power_exact_goal(300, 0.9, 0.8), power_exact_goal(50, 0.9, 0.8))
  expect_lte(power_exact_goal(100, 0.8, 0.8), 0.05)
})

test_that("Hanley-McNeil variance equals known SE for AUC 0.85, 50/50", {
  expect_equal(sqrt(var_auc_hm(0.85, 50, 50)), 0.0427, tolerance = 0.01)
})

test_that("n_auc_hm increases as the effect shrinks", {
  expect_gt(n_auc_hm(0.75), n_auc_hm(0.90))
})

test_that("two-cut classification boundaries follow the FDA convention (negative <= lo, positive > hi)", {
  expect_equal(as.character(classify_two_cut(c(2, 8, 8.001), 2, 8)), c("rule-out", "intermediate", "rule-in"))
})

test_that("pooled SD within samples ignores differences between sample means", {
  y <- c(10, 11, 12, 20, 21, 22); g <- rep(1:2, each = 3)
  expect_equal(sd_pooled(y, g), 1)
  expect_equal(sd_pooled(y), sd(y))
  # pooled-SD LoD uses the within-sample SD, not the spread of all results
  expect_equal(lod_param(0, y, group = g), 1.645)
})

test_that("small-sample LoD multiplier", {
  y <- c(10, 11, 12, 20, 21, 22); g <- rep(1:2, each = 3)    # L = 6 results, J = 2 samples
  expect_equal(lod_param(0, y, group = g, cp = TRUE), 1.645 / (1 - 1 / (4 * (6 - 2))))
})

test_that("nonparametric LoB uses rank 0.5 + 0.95 B (quantile type 5)", {
  b <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4)    # B = 20
  r <- 0.5 + 0.95 * 20                                                   # rank 19.5
  s <- sort(b)
  expect_equal(lob_nonparam(b), s[19] + 0.5 * (s[20] - s[19]))
})

test_that("nonparametric LoD adds the median-to-5th-percentile distance", {
  low <- c(0:19, 100:119); g <- rep(1:2, each = 20)    # each sample centred: -9.5 ... 9.5
  centred <- c(0:19 - 9.5, 0:19 - 9.5)
  expect_equal(lod_nonparam(1, low, g), 1 - unname(quantile(centred, 0.05, type = 5)))
})

test_that("Deming: hand delta = var(err_y)/var(err_x) is the inverse of mcr's error.ratio", {
  skip_if_not_installed("mcr")
  set.seed(11); t <- runif(80, 1, 100); x <- t + rnorm(80, 0, 2); y <- 1 + 1.1 * t + rnorm(80, 0, 6)
  for (r in c(0.25, 2)) {
    m <- mcr::mcreg(x, y, method.reg = "Deming", error.ratio = r, method.ci = "analytical")
    expect_equal(unname(deming_hand(x, y, delta = 1 / r)), unname(m@para[, "EST"]), tolerance = 1e-6)
  }
})

test_that("derive_two_cutoffs fails clearly when no cut-off meets a target and warns on overlap", {
  x <- c(1, 2, 3, 4); y <- c(0, 1, 0, 1)
  expect_error(derive_two_cutoffs(x, y, 1.01, 0.5), "No observed cut-off")
  set.seed(12); y2 <- rbinom(300, 1, .5); x2 <- rnorm(300, y2 * 4)   # well separated, modest targets
  expect_warning(co <- derive_two_cutoffs(x2, y2, 0.5, 0.5), "overlap")
  expect_equal(co[["lo"]], co[["hi"]])
})

test_that("sim_ep05 is reproducible and has the balanced design", {
  a <- sim_ep05(seed = 20); b <- sim_ep05(seed = 20)
  expect_identical(a, b)
  expect_equal(nrow(a), 20 * 2 * 2)
  expect_equal(nlevels(a$day), 20)
})

test_that("show_sources prints a markdown table for one chapter", {
  f <- tempfile(fileext = ".csv")
  write.csv(data.frame(claim = c("A | b", "C"), source = c("https://x.org", "book"),
                       verified_on = "2026-01-01", confidence = "high", chapter = c(1, 2)), f, row.names = FALSE)
  out <- capture.output(show_sources(1, f))
  expect_equal(out[1], "| Claim | Source | Checked | Confidence |")
  expect_match(out[3], "A / b \\| <https://x.org>", fixed = FALSE)
  expect_length(out, 3)
})
