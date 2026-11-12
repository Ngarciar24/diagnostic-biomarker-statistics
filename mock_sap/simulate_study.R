# Simulated two-cut-off pTau217-like clinical performance study ("SIM-217-001").
# Run from the repository root:  Rscript mock_sap/simulate_study.R
# All data are simulated. No real assay or patient data are used.
source("R/helpers.R")

marker <- function(y, egfr, shift = TRUE) {
  lx <- ifelse(y == 1, rnorm(length(y), -0.20, 0.30), rnorm(length(y), -0.80, 0.25))
  if (shift) lx <- lx + ifelse(egfr >= 60, 0, ifelse(egfr >= 30, 0.10, 0.30))   # assumed eGFR effect
  10^lx
}

# 1. Derivation cohort (separate from the performance study) -> locked cut-offs
set.seed(9001)
n0 <- 700; y0 <- rbinom(n0, 1, 0.40); egfr0 <- round(pmin(120, pmax(8, rnorm(n0, 78, 22))))
x0 <- marker(y0, egfr0)
co <- derive_two_cutoffs(x0, y0, se_target = 0.94, sp_target = 0.93)
locked <- data.frame(cutoff = c("rule_out_le", "rule_in_gt"), pg_per_ml = round(unname(co), 3))
write.csv(locked, "mock_sap/data/locked_cutoffs.csv", row.names = FALSE)

# 2. Performance study
set.seed(9002)
N <- 1100
setting <- sample(c("primary", "secondary"), N, replace = TRUE, prob = c(0.25, 0.75))
age <- round(pmin(85, pmax(55, rnorm(N, 71, 7))))
sex <- sample(c("F", "M"), N, replace = TRUE, prob = c(0.55, 0.45))
mmse <- pmin(30, pmax(15, round(rnorm(N, 26, 3))))
egfr <- round(pmin(120, pmax(8, rnorm(N, 78, 22))))
prev <- ifelse(setting == "primary", 0.22, 0.45) + ifelse(age > 70, 0.05, -0.05)
pet <- rbinom(N, 1, prev)
ptau <- marker(pet, egfr)
d <- data.frame(
  USUBJID = sprintf("SIM-%04d", seq_len(N)),
  SETTING = setting, AGE = age, SEX = sex, MMSE = mmse, EGFR = egfr,
  PET_VISUAL = pet,                                  # reference standard (simulated truth)
  PTAU217 = round(ptau, 3),
  PET_INTERVAL_D = round(abs(rnorm(N, 38, 29))),
  stringsAsFactors = FALSE)
d$AGEGR <- cut(d$AGE, c(54, 70, 85), labels = c("55-70", "71-85"))
d$EGFRGR <- cut(d$EGFR, c(0, 30, 60, Inf), right = FALSE, labels = c("<30", "30-59", ">=60"))
d$MMSEGR <- cut(d$MMSE, c(14, 24, 27, 30), labels = c("15-24", "25-27", "28-30"))
# invalid assay results (2%) and unreadable PET (3%) -> missing
d$PTAU217[sample(N, round(0.02 * N))] <- NA
d$PET_VISUAL[sample(N, round(0.03 * N))] <- NA
d$PPFL <- ifelse(d$PET_INTERVAL_D <= 180 & runif(N) > 0.02, "Y", "N")      # no major protocol deviation
write.csv(d, "mock_sap/data/sim_ptau217.csv", row.names = FALSE)
cat("Locked cut-offs:", paste(locked$cutoff, locked$pg_per_ml, collapse = "; "), "\n")
cat("Study rows:", nrow(d), " PET+ prevalence:", round(mean(d$PET_VISUAL, na.rm = TRUE), 3), "\n")
