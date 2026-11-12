"""Unit tests for check_sap.py. Run from the repository root:  python3 -m pytest python"""

import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_sap as cs  # noqa: E402


def test_wilson_matches_r_helper():
    # values from ci_wilson(37, 50) in R/helpers.R
    est, lo, hi = cs.wilson(37, 50)
    assert est == pytest.approx(0.74, abs=1e-15)
    assert lo == pytest.approx(0.604468427351512, abs=1e-12)
    assert hi == pytest.approx(0.841284725064476, abs=1e-12)


def test_wilson_boundaries():
    assert cs.wilson(0, 20)[1] == pytest.approx(0.0, abs=1e-12)
    assert cs.wilson(0, 20)[2] == pytest.approx(0.161125158052819, abs=1e-12)
    assert cs.wilson(20, 20)[2] == 1.0
    with pytest.raises(ValueError):
        cs.wilson(5, 0)


def _db(rows):
    con = sqlite3.connect(":memory:")
    cols = ", ".join(f"{c} {t}" for c, t in cs.STUDY_COLUMNS.items())
    con.execute(f"CREATE TABLE study ({cols})")
    con.executemany(f"INSERT INTO study VALUES ({', '.join('?' * len(cs.STUDY_COLUMNS))})", rows)
    con.execute("CREATE TABLE cutoffs (cutoff TEXT PRIMARY KEY, pg_per_ml REAL)")
    con.executemany("INSERT INTO cutoffs VALUES (?, ?)", [("rule_out_le", 0.2), ("rule_in_gt", 0.4)])
    return con


GOOD = ("SIM-0001", "primary", 60, "F", 28, 80, 1, 0.5, 30, "55-70", ">=60", "28-30", "Y")


def test_edit_checks_pass_on_clean_record():
    con = _db([GOOD])
    assert all(c["failures"] == 0 for c in cs.run_edit_checks(con))


def test_edit_checks_catch_injected_errors():
    dup = GOOD
    bad_age = ("SIM-0002", "primary", 40, "F", 28, 80, 1, 0.5, 30, "55-70", ">=60", "28-30", "Y")
    bad_group = ("SIM-0003", "primary", 75, "F", 28, 45, 1, 0.5, 30, "55-70", ">=60", "28-30", "Y")
    bad_ppfl = ("SIM-0004", "primary", 60, "F", 28, 80, 1, 0.5, 200, "55-70", ">=60", "28-30", "Y")
    con = _db([GOOD, dup, bad_age, bad_group, bad_ppfl])
    failed = {c["check"]: c["failures"] for c in cs.run_edit_checks(con) if c["failures"]}
    assert failed["duplicate subject IDs"] == 1
    assert failed["age outside 55-85"] == 1
    assert failed["age group inconsistent with age"] == 1      # SIM-0002 (age 40) and SIM-0003 (age 75)
    assert failed["eGFR group inconsistent with eGFR"] == 1
    assert failed["per-protocol flag with PET interval > 180 days"] == 1


def test_category_boundaries_follow_the_sap_convention():
    # equal to the lower cut-off -> rule-out; equal to the upper -> intermediate; above -> rule-in
    rows = [(f"SIM-{i:04d}", "primary", 60, "F", 28, 80, pet, x, 30, "55-70", ">=60", "28-30", "Y")
            for i, (x, pet) in enumerate([(0.2, 0), (0.4, 1), (0.41, 1), (None, 1)])]
    con = _db(rows)
    con.execute(cs.CATEGORY_SQL)
    cats = dict(con.execute("SELECT USUBJID, cat FROM evaluable"))
    assert cats == {"SIM-0000": "rule-out", "SIM-0001": "intermediate", "SIM-0002": "rule-in"}


def test_full_check_agrees_with_r_outputs():
    assert cs.main() == 0
