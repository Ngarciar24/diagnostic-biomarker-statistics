root <- normalizePath(file.path("..", ".."))
source(file.path(root, "R", "helpers.R")); source(file.path(root, "R", "bridging.R"))

test_that("sim_trial returns the right columns and missingness", {
  set.seed(1); d <- sim_trial(2000, miss = "MCAR", p_obs = 0.6)
  expect_named(d, c("arm", "B", "cta", "cdx", "resp", "obs"))
  expect_equal(mean(d$obs), 0.6, tolerance = 0.05)
  expect_equal(mean(sim_trial(500, miss = "none")$obs), 1)
})

test_that("with no missing data, complete-case, IPW and full estimators coincide", {
  set.seed(2); d <- sim_trial(800, miss = "none")
  expect_equal(est_cc(d), est_full(d))
  expect_equal(unname(est_ipw(d)["est"]), unname(est_full(d)["est"]), tolerance = 1e-8)
})

test_that("rd_est matches a hand calculation", {
  r <- rd_est(c(1, 1, 0, 0, 0), c(1, 0, 0, 0, 1))
  expect_equal(unname(r["est"]), 0.5 - 1 / 3)
})

test_that("full-information estimator is unbiased for the analytic truth with a large trial", {
  set.seed(3); d <- sim_trial(400000, miss = "none")
  expect_equal(unname(est_full(d)["est"]), true_effect_cdx_pos(), tolerance = 0.01)
})

test_that("analytic truth: perfect CDx gives the biomarker-positive effect", {
  expect_equal(true_effect_cdx_pos(cdx = c(se = 1, sp = 1)), 0.30)
})

test_that("multiple imputation is unbiased under MAR, complete case may not be", {
  set.seed(4); truth <- true_effect_cdx_pos()
  r <- replicate(120, { d <- sim_trial(1200, miss = "MAR"); c(cc = est_cc(d)[["est"]], mi = est_mi(d, 10)[["est"]]) })
  expect_lt(abs(mean(r["mi", ]) - truth), 0.02)
})

test_that("complete case is unbiased when missingness depends only on the CDx stratum, biased when it also depends on outcome", {
  set.seed(5); truth <- true_effect_cdx_pos()
  m1 <- mean(replicate(150, est_cc(sim_trial(1200, miss = "MNAR"))[["est"]]))
  m2 <- mean(replicate(150, est_cc(sim_trial(1200, miss = "MNAR_resp"))[["est"]]))
  expect_lt(abs(m1 - truth), 0.015)
  expect_gt(abs(m2 - truth), 0.02)
})

test_that("delta-adjusted tipping point: delta = 0 matches MI and negative delta moves the effect towards zero", {
  set.seed(6); d <- sim_trial(3000, miss = "MNAR")
  set.seed(7); e0 <- est_tipping_delta(d, 0, M = 40)[["est"]]
  set.seed(7); mi <- est_mi(d, 40)[["est"]]
  expect_equal(e0, mi, tolerance = 0.01)
  set.seed(8); eneg <- est_tipping_delta(d, -4, M = 40)[["est"]]
  set.seed(8); epos <- est_tipping_delta(d, 2, M = 40)[["est"]]
  expect_lt(eneg, e0); expect_gt(epos, e0)
})
