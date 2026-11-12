"""Cross-language check of the mock SAP (SIM-217-001), independent of the R code.

1. Load the raw study CSV and the locked cut-offs into an in-memory SQLite database.
2. Run data-management edit checks in SQL (duplicates, ranges, codelists, derived groups).
3. Derive the result categories and the primary and secondary endpoints in SQL,
   compute Wilson 95% intervals in plain Python.
4. Compare with the endpoints exported by the R analysis (outputs/sap_primary_endpoints.csv)
   and exit non-zero on any edit-check failure or any difference above the tolerance.

Standard library only. Run from the repository root:  python3 python/check_sap.py
"""

from __future__ import annotations

import csv
import math
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "mock_sap" / "data"
OUT = ROOT / "outputs"
TOL = 1e-10

STUDY_COLUMNS = {
    "USUBJID": "TEXT", "SETTING": "TEXT", "AGE": "INTEGER", "SEX": "TEXT", "MMSE": "INTEGER",
    "EGFR": "INTEGER", "PET_VISUAL": "INTEGER", "PTAU217": "REAL", "PET_INTERVAL_D": "INTEGER",
    "AGEGR": "TEXT", "EGFRGR": "TEXT", "MMSEGR": "TEXT", "PPFL": "TEXT",
}

# Each edit check is a query returning the offending rows; an empty result is a pass.
EDIT_CHECKS = {
    "duplicate subject IDs":
        "SELECT USUBJID FROM study GROUP BY USUBJID HAVING COUNT(*) > 1",
    "age outside 55-85":
        "SELECT USUBJID FROM study WHERE AGE IS NULL OR AGE NOT BETWEEN 55 AND 85",
    "MMSE outside 0-30":
        "SELECT USUBJID FROM study WHERE MMSE IS NULL OR MMSE NOT BETWEEN 0 AND 30",
    "eGFR missing or not positive":
        "SELECT USUBJID FROM study WHERE EGFR IS NULL OR EGFR <= 0",
    "marker not positive":
        "SELECT USUBJID FROM study WHERE PTAU217 IS NOT NULL AND PTAU217 <= 0",
    "PET read not 0/1":
        "SELECT USUBJID FROM study WHERE PET_VISUAL IS NOT NULL AND PET_VISUAL NOT IN (0, 1)",
    "sex not in codelist":
        "SELECT USUBJID FROM study WHERE SEX NOT IN ('F', 'M')",
    "setting not in codelist":
        "SELECT USUBJID FROM study WHERE SETTING NOT IN ('primary', 'secondary')",
    "PET interval negative":
        "SELECT USUBJID FROM study WHERE PET_INTERVAL_D < 0",
    "per-protocol flag not Y/N":
        "SELECT USUBJID FROM study WHERE PPFL NOT IN ('Y', 'N')",
    "per-protocol flag with PET interval > 180 days":
        "SELECT USUBJID FROM study WHERE PPFL = 'Y' AND PET_INTERVAL_D > 180",
    "age group inconsistent with age":
        """SELECT USUBJID FROM study
           WHERE AGEGR <> CASE WHEN AGE <= 70 THEN '55-70' ELSE '71-85' END""",
    "eGFR group inconsistent with eGFR":
        """SELECT USUBJID FROM study
           WHERE EGFRGR <> CASE WHEN EGFR < 30 THEN '<30' WHEN EGFR < 60 THEN '30-59' ELSE '>=60' END""",
    "MMSE group inconsistent with MMSE":
        """SELECT USUBJID FROM study
           WHERE MMSEGR <> CASE WHEN MMSE <= 24 THEN '15-24' WHEN MMSE <= 27 THEN '25-27' ELSE '28-30' END""",
}


def wilson(x: int, n: int, z: float = 1.959963984540054) -> tuple[float, float, float]:
    """Wilson score interval, the same formula as ci_wilson() in R/helpers.R."""
    if n <= 0 or not 0 <= x <= n:
        raise ValueError("need 0 <= x <= n and n > 0")
    p = x / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return p, max(0.0, centre - half), min(1.0, centre + half)


def _value(raw: str, sql_type: str):
    if raw in ("", "NA"):
        return None
    if sql_type == "INTEGER":
        return int(float(raw))
    if sql_type == "REAL":
        return float(raw)
    return raw


def load(con: sqlite3.Connection, study_csv: Path, cutoff_csv: Path) -> None:
    cols = ", ".join(f"{c} {t}" for c, t in STUDY_COLUMNS.items())
    con.execute(f"CREATE TABLE study ({cols})")
    with study_csv.open(newline="") as f:
        rows = [tuple(_value(r[c], t) for c, t in STUDY_COLUMNS.items()) for r in csv.DictReader(f)]
    con.executemany(f"INSERT INTO study VALUES ({', '.join('?' * len(STUDY_COLUMNS))})", rows)
    con.execute("CREATE TABLE cutoffs (cutoff TEXT PRIMARY KEY, pg_per_ml REAL)")
    with cutoff_csv.open(newline="") as f:
        con.executemany("INSERT INTO cutoffs VALUES (?, ?)",
                        [(r["cutoff"], float(r["pg_per_ml"])) for r in csv.DictReader(f)])


def run_edit_checks(con: sqlite3.Connection) -> list[dict]:
    out = []
    for name, query in EDIT_CHECKS.items():
        ids = [r[0] for r in con.execute(query)]
        out.append({"check": name, "failures": len(ids), "example_ids": ";".join(ids[:5])})
    return out


# Category rule as in the SAP and the FDA convention: negative if <= lower, positive if > upper.
CATEGORY_SQL = """
CREATE VIEW evaluable AS
SELECT s.USUBJID, s.PET_VISUAL AS pet,
       CASE WHEN s.PTAU217 <= lo.pg_per_ml THEN 'rule-out'
            WHEN s.PTAU217 >  hi.pg_per_ml THEN 'rule-in'
            ELSE 'intermediate' END AS cat
FROM study s
CROSS JOIN (SELECT pg_per_ml FROM cutoffs WHERE cutoff = 'rule_out_le') lo
CROSS JOIN (SELECT pg_per_ml FROM cutoffs WHERE cutoff = 'rule_in_gt')  hi
WHERE s.PTAU217 IS NOT NULL AND s.PET_VISUAL IS NOT NULL
"""

ENDPOINT_SQL = {
    # endpoint: (numerator, denominator)
    "ppv_positive": ("SELECT COUNT(*) FROM evaluable WHERE cat = 'rule-in' AND pet = 1",
                     "SELECT COUNT(*) FROM evaluable WHERE cat = 'rule-in'"),
    "npv_negative": ("SELECT COUNT(*) FROM evaluable WHERE cat = 'rule-out' AND pet = 0",
                     "SELECT COUNT(*) FROM evaluable WHERE cat = 'rule-out'"),
    "ppa_incl": ("SELECT COUNT(*) FROM evaluable WHERE cat = 'rule-in' AND pet = 1",
                 "SELECT COUNT(*) FROM evaluable WHERE pet = 1"),
    "npa_incl": ("SELECT COUNT(*) FROM evaluable WHERE cat = 'rule-out' AND pet = 0",
                 "SELECT COUNT(*) FROM evaluable WHERE pet = 0"),
    "intermediate_fraction": ("SELECT COUNT(*) FROM evaluable WHERE cat = 'intermediate'",
                              "SELECT COUNT(*) FROM evaluable"),
}


def derive_endpoints(con: sqlite3.Connection) -> dict[str, dict]:
    con.execute(CATEGORY_SQL)
    res = {}
    for name, (num_q, den_q) in ENDPOINT_SQL.items():
        x = con.execute(num_q).fetchone()[0]
        n = con.execute(den_q).fetchone()[0]
        est, lo, hi = wilson(x, n)
        res[name] = {"x": x, "n": n, "est": est, "lower": lo, "upper": hi}
    return res


def compare(py: dict[str, dict], r_csv: Path, tol: float = TOL) -> list[dict]:
    with r_csv.open(newline="") as f:
        r_rows = {r["endpoint"]: r for r in csv.DictReader(f)}
    if set(r_rows) != set(py):
        raise ValueError(f"endpoint sets differ: R {sorted(r_rows)} vs Python {sorted(py)}")
    out = []
    for name, p in py.items():
        r = r_rows[name]
        for field in ("x", "n", "est", "lower", "upper"):
            diff = abs(float(r[field]) - p[field])
            out.append({"endpoint": name, "field": field, "R": float(r[field]), "Python_SQL": p[field],
                        "abs_diff": diff, "ok": diff <= tol})
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    con = sqlite3.connect(":memory:")
    load(con, DATA / "sim_ptau217.csv", DATA / "locked_cutoffs.csv")
    n = con.execute("SELECT COUNT(*) FROM study").fetchone()[0]

    checks = run_edit_checks(con)
    write_csv(OUT / "qc" / "edit_checks.csv", checks)
    failed_checks = [c for c in checks if c["failures"]]
    print(f"Edit checks on {n} records: {len(checks) - len(failed_checks)}/{len(checks)} pass")
    for c in failed_checks:
        print(f"  FAIL {c['check']}: {c['failures']} (e.g. {c['example_ids']})")

    endpoints = derive_endpoints(con)
    cmp = compare(endpoints, OUT / "sap_primary_endpoints.csv")
    write_csv(OUT / "qc" / "python_crosscheck.csv", cmp)
    bad = [c for c in cmp if not c["ok"]]
    for name, e in endpoints.items():
        print(f"  {name:<22} {e['x']:>4}/{e['n']:<4} {100 * e['est']:5.1f}% "
              f"({100 * e['lower']:.1f}-{100 * e['upper']:.1f})")
    print(f"R vs Python/SQL: {len(cmp) - len(bad)}/{len(cmp)} values agree (tolerance {TOL:g})")
    for c in bad:
        print(f"  MISMATCH {c['endpoint']}.{c['field']}: R {c['R']} vs Python {c['Python_SQL']}")
    return 1 if failed_checks or bad else 0


if __name__ == "__main__":
    sys.exit(main())
