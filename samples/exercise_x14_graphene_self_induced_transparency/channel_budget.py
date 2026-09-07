#!/usr/bin/env python3
"""Where the absorbed energy goes: the per-channel ledger of a dissipative sheet run.

The solver writes a CUMULATIVE per-cell ledger to `*_sbe_channels.data`, one pair of
columns per ring channel: dN, the conduction-population change (pairs created > 0), and
dE [Ha], the eigenvalue-weighted energy the electrons GAINED from that channel
(`ring_ledger` in bloch_solver_ssbe.f90 accumulates sum(eval * dpop)/nk, so dE is the
change this channel makes to Tr(rho H); negative means the carriers gave it up).

Read `*_sbe_rt_energy.data` with that in mind. Its Eall is NOT Tr(rho H): realtime_ssbe
accumulates `energy += (E_tot . -J) * volume * dt`, which is the WORK the local field
does on the sheet. So

    dE_all  =  W_field                    (an integration check, not a statement)
    E_electronic  =  W_field + sum_ch dE_ch                                    (1)

and the two are different numbers whenever a channel has a bath on the other side.
Exactly one does: e-ph, whose -dE is what the carriers hand to the phonons. Auger and
impact ionization move energy WITHIN the electron gas, so their dE must come out near
zero while their dN does not -- which is the check that they are doing what they claim.

Every column is reported against a DARK run when one is given, and it should be: with
the ring on, the solver pumps a field-independent current AND a field-independent energy
with the drive switched off (README SS7.11), and on this mesh that offset exceeds the
whole work a 100 kV/cm pulse does. The dark column also says WHICH channel is
responsible, which no other diagnostic here can. Subtracting it makes the remaining
columns physics; it does not rehabilitate the low-field transmissions, because the dark
current is a current and the fields it condemns stay condemned.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transmission import (read_rt, LZ_BOHR, AREA_BOHR2, AU_EV, AU_E_VM,  # noqa: E402
                          TruncatedRun)
from drude_check import run_variable, density_cm2, KB_AU                    # noqa: E402

HA_MEV = AU_EV * 1e3
BOHR_CM = 0.529177210903e-8
# t, then (dN, dE) for e-ph, impact ionization, Auger, Rana
COLS = dict(eph=(1, 2), ii=(3, 4), auger=(5, 6), rana=(7, 8))


def budget(rt_path):
    """(W_field, dE_all, {channel: (dN, dE)}) per cell, in eV, at the end of the run."""
    d = read_rt(rt_path)
    t = d['Time']
    check = len(t)
    ax = max('xyz', key=lambda a: np.max(np.abs(d[f'E_ext_{a}'])))
    nl = int(run_variable(rt_path, 'sbe_sheet_nlayers', 1) or 1)
    J = -d[f'Jm_{ax}'] * LZ_BOHR * nl
    Et = d.get(f'E_tot_{ax}')
    if Et is None:
        Et = d[f'E_ext_{ax}']
    trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    W = float(trapz(J * Et, t)) * AREA_BOHR2 * AU_EV        # eV / cell

    stem = rt_path[:-len('_rt.data')] if rt_path.endswith('_rt.data') else rt_path
    de_all = np.nan
    ep = stem + '_rt_energy.data'
    if os.path.exists(ep):
        de_all = float(np.loadtxt(ep)[-1, 1])               # eV / cell
    ch = {}
    cp = stem + '_channels.data'
    if os.path.exists(cp):
        row = np.loadtxt(cp)[-1]
        for name, (i, j) in COLS.items():
            ch[name] = (float(row[i]), float(row[j]) * AU_EV)
    return W, de_all, ch, check


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+', help='*_sbe_rt.data of the driven runs (globs ok)')
    ap.add_argument('--dark', default=None,
                    help='*_sbe_rt.data of the zero-field control, subtracted from every '
                         'column. Give it whenever the ring is on')
    args = ap.parse_args(argv)

    files = sorted(set(sum((glob.glob(p) for p in args.files), [])))
    dW = dE = 0.0
    dch = {k: (0.0, 0.0) for k in COLS}
    if args.dark:
        dW, dE, dch, _ = budget(args.dark)
        print(f'# DARK control, {args.dark}')
        print(f'#   W_field = {dW * 1e3:+.4f} meV/cell   dE_all = {dE * 1e3:+.4f} meV/cell')
        for k, (n, e) in dch.items():
            if abs(e) > 0 or abs(n) > 0:
                print(f'#   {k:6s} dN = {n:+.3e} /cell   dE = {e * 1e3:+.4f} meV/cell')
        print('#   -- every column below has this subtracted --')

    rows = []
    for f in files:
        try:
            W, ea, ch, _ = budget(f)
        except TruncatedRun as exc:
            print(f'# SKIPPED: {exc}')
            continue
        d = read_rt(f)
        ax = max('xyz', key=lambda a: np.max(np.abs(d[f'E_ext_{a}'])))
        e0 = float(np.max(np.abs(d[f'E_ext_{ax}']))) * AU_E_VM / 1e5        # kV/cm
        rows.append((e0, W - dW, ea - dE,
                     {k: (ch[k][0] - dch[k][0], ch[k][1] - dch[k][1]) for k in ch}, f))
    rows.sort(key=lambda r: r[0])

    ef = float(run_variable(files[0], 'sbe_ef_ev', 0.0) or 0.0) if files else 0.0
    tk = float(run_variable(files[0], 'sbe_temp_init_k', 300.0) or 300.0) if files else 300.0
    n_cell = density_cm2(abs(ef) / AU_EV, KB_AU * tk) * AREA_BOHR2 * BOHR_CM**2

    print(f'{"E0":>7} {"W_field":>10} {"E_elec":>10} {"to lattice":>11} {"":>7} '
          f'{"dN_eph":>11} {"dN_rana":>11} {"dE_rana":>10} {"Eall-W":>9}')
    print(f'{"kV/cm":>7} {"meV/cell":>10} {"meV/cell":>10} {"meV/cell":>11} {"% of W":>7} '
          f'{"/cell":>11} {"/cell":>11} {"meV/cell":>10} {"rel.":>9}')
    for e0, W, ea, ch, f in rows:
        lat = -ch['eph'][1] * 1e3
        e_el = (W + sum(ch[k][1] for k in ch)) * 1e3
        chk = (ea - W) / max(abs(W), 1e-30)
        print(f'{e0:7.1f} {W * 1e3:10.4f} {e_el:10.4f} {lat:11.4f} '
              f'{100 * lat / max(W * 1e3, 1e-30):7.1f} {ch["eph"][0]:11.3e} '
              f'{ch["rana"][0]:11.3e} {ch["rana"][1] * 1e3:10.2e} {chk:+9.1e}')
    if n_cell > 0:
        print(f'# doping: {n_cell:.4e} carriers per cell (E_F = {ef:g} eV, T = {tk:g} K), so a '
              f'dN of 1e-6 is {1e-6 / n_cell * 100:.1e} % of them')
    print('# E_elec = W_field + sum dE_channels, Eq. (1): what is left in the electron gas.')
    print('# "to lattice" = -dE_eph: e-ph is the only channel with a bath on the other side.')
    print('# Eall-W is an integration check (trapezoid here against the code rectangle rule),')
    print('#   not physics: Eall IS the accumulated work, so it can only differ by that.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
