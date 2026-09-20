# -*- coding: utf-8 -*-
"""Public reproduction orchestrator. Stages skip when outputs exist."""
import os, sys, subprocess, argparse, json, time

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
os.makedirs(RES, exist_ok=True)

PIPELINE = ["01_download_data.py", "01b_download_estimate.py",
            "02_build_matrices.py",
            "03_nested_compression.py", "04_adipose_ablation.py",
            "05_cross_cohort_axis.py", "06_estimate_baseline.py",
            "07_transfer.py", "08_clinical_subtypes.py"]

def get(d, path):
    cur = d
    for part in path.split("/"):
        cur = cur[part]
    return float(cur)

# verification manifest (path uses / to walk nested dicts)
EXPECTED = {
    "nested_compression.json": {
        "full/acc": 0.989, "panelA/acc": 0.947, "panelB/acc": 0.958,
        "panelC/acc": 0.987, "k1/acc": 0.970, "null/acc": 0.920,
        "randDE/acc": 0.958},
    "adipose_ablation.json": {
        "before/acc": 0.989, "after/acc": 0.909},
    "cross_cohort_axis.json": {
        "tcga_tumor": -1.374, "tcga_adjacent_normal": 0.400,
        "gtex_breast": 0.388, "gtex_adipose": 0.696},
    "transfer.json": {
        "gtex_normal_frac": 0.994, "metabric_frac_tumor": 0.924},
    "clinical_subtypes.json": {
        "histology_bh": 1.34e-04, "er_bh": 1.61e-02,
        "pr_bh": 3.97e-02, "pam50_chi2_without_normal": 2.25e-03},
}
TOL = 0.02

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    for script in PIPELINE:
        base = script[3:-3]
        out = os.path.join(RES, base + ".json")
        if script.startswith(("01", "02")):
            print(f"[run ] {script} (idempotente: baixa apenas o que falta)",
                  flush=True)
        elif os.path.exists(out) and not args.force:
            print(f"[done] {script}")
            continue
        r = subprocess.run([sys.executable, os.path.join(HERE, script)],
                           cwd=HERE)
        if r.returncode != 0:
            raise SystemExit(f"FAILED: {script}")

    print("\n=== verification (headline numbers) ===")
    n_ok = n_bad = 0
    for jf, checks in EXPECTED.items():
        p = os.path.join(RES, jf)
        if not os.path.exists(p):
            print(f"  [skip] {jf}")
            continue
        d = json.load(open(p, encoding="utf-8"))
        for path, exp in checks.items():
            v = get(d, path)
            ok = abs(v - exp) <= TOL
            n_ok += ok
            n_bad += (not ok)
            print(f"  [{'PASS' if ok else 'FAIL'}] {jf}::{path} = "
                  f"{v:.4f} (esperado {exp})")
    print(f"\n{n_ok} PASS / {n_bad} FAIL "
          f"({time.time()-t0:.0f}s)")
    sys.exit(0 if n_bad == 0 else 1)

if __name__ == "__main__":
    main()
