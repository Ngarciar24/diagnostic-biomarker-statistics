# -----------------------------------------------------------------------------
# Program    : run_all.R
# Purpose    : Rebuild and check everything in one command: R unit tests, the
#              book (every chapter stops on a failed cross-check), and the
#              Python/SQL re-derivation of the mock SAP endpoints
# Inputs     : the repository; R packages from renv.lock; Quarto CLI; python3
# Outputs    : _book/ (rendered book), outputs/sap_primary_endpoints.csv,
#              outputs/qc/edit_checks.csv, outputs/qc/python_crosscheck.csv
# Author     : Ignacio G. Ribelles
# Created    : 2026-11-12
# Change log : 2026-11-12  IGR  Initial version
#              2026-11-12  IGR  Stop if pytest is missing instead of skipping it
# Notes      : Run from the repository root:  Rscript run_all.R
#              Stops at the first failing step with a non-zero exit status.
# -----------------------------------------------------------------------------

step <- function(label, ok) {
  cat(sprintf("\n== %s\n", label))
  if (!isTRUE(ok)) stop(label, " failed", call. = FALSE)
}

run <- function(cmd, args) identical(system2(cmd, args), 0L)

step("1/3 R unit tests", {
  res <- testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
  !any(as.data.frame(res)$failed > 0) && !any(as.data.frame(res)$error)
})

step("2/3 Render the book", run("quarto", "render"))

if (!run("python3", c("-c", shQuote("import pytest"))))
  stop("pytest is not installed: run  pip install pytest  and try again", call. = FALSE)
step("3/3 Python / SQL re-derivation of the SAP endpoints",
     run("python3", "python/check_sap.py") && run("python3", c("-m", "pytest", "-q", "python")))

cat("\nAll checks passed. Book: _book/index.html\n")
