#!/usr/bin/env python3
"""
run_all.py  -  run every tests/test_*.py (Python) and test_*.f90 (Fortran) and
report PASS/FAIL.

Each test prints a final 'PASS'/'FAIL' line and exits 0 (pass) / nonzero (fail).
Python tests run directly. Fortran tests are compiled with gfortran against
their self-contained source dependencies (declared in FORTRAN_DEPS below), then
run. Exit code is the number of failed tests (0 = all pass), CI-friendly.

Usage:  python3 tests/run_all.py
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Fortran tests and the (self-contained) module sources they compile against.
FORTRAN_DEPS = {
    "test_superres_rates.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_eph_cptp.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_screening.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_carrier_carrier.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_vg_basis_nb.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_material_registry.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_auger_cptp.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_eph_interk_cptp.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_eph_detailed_balance.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_mp_kmap.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_ii_interk_cptp.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_auger_interk_cptp.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_rana_2d.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_rana_auger_cptp.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_bath_corr.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_sfsb_kernel.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_colmem_lines.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_dressed_ref.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_rana_saturation.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_colmem_2d.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_dirac_te_fit.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
    "test_vg_sumrule.f90": ["src/ssbe/sbe_superres_ssbe.f90"],
}


def run_python(t):
    r = subprocess.run([sys.executable, t], capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.stderr.strip():
        sys.stdout.write(r.stderr)
    return r.returncode == 0


def run_fortran(t):
    name = os.path.basename(t)
    fc = shutil.which("gfortran")
    if fc is None:
        print("  SKIP: gfortran not found"); return True   # don't fail CI on missing toolchain
    deps = [os.path.join(ROOT, d) for d in FORTRAN_DEPS.get(name, [])]
    with tempfile.TemporaryDirectory() as tmp:
        exe = os.path.join(tmp, "a.out")
        cmd = [fc, "-J", tmp, "-ffree-line-length-none", *deps, t, "-o", exe]
        c = subprocess.run(cmd, capture_output=True, text=True, cwd=tmp)
        if c.returncode != 0:
            sys.stdout.write(c.stdout + c.stderr)
            print("  FAIL: compilation error"); return False
        r = subprocess.run([exe], capture_output=True, text=True)
        sys.stdout.write(r.stdout)
        if r.stderr.strip():
            sys.stdout.write(r.stderr)
        return r.returncode == 0


def main():
    tests = sorted(glob.glob(os.path.join(HERE, "test_*.py")) +
                   glob.glob(os.path.join(HERE, "test_*.f90")))
    if not tests:
        print("no tests found"); return 1
    failed = []
    for t in tests:
        name = os.path.basename(t)
        print(f"\n=== {name} ===")
        ok = run_fortran(t) if t.endswith(".f90") else run_python(t)
        if not ok:
            failed.append(name)
    print("\n" + "=" * 50)
    print(f"  {len(tests) - len(failed)}/{len(tests)} passed")
    if failed:
        print("  FAILED: " + ", ".join(failed))
    return len(failed)


if __name__ == "__main__":
    sys.exit(main())
