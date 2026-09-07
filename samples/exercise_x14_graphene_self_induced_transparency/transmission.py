#!/usr/bin/env python3
"""Sheet (2D) transmission from the SBE current -- "field before and after".

The single-cell SBE driver (realtime_ssbe) writes the EXTERNAL field it was
driven with and the matter current density Jm of the graphene cell (per cell
VOLUME, vacuum included) into *_sbe_rt.data. Its Ac_tot/E_tot columns are
identical to the external ones: there is NO self-consistent field in the
single-cell driver. For a 2D sheet at normal incidence the transmitted field
follows from the sheet boundary condition (E continuous across the sheet, H
jumps by the sheet current), Hartree atomic units:

    E_t = (2 E_inc - Z0 J_s) / (1 + n_sub),      E_r = E_t - E_inc
    Z0  = 4 pi / c
    J_s = physical charge sheet current = -Jm * L_z

(Jm is the ELECTRON current per volume -- the energy ledger in realtime_ssbe
is dW = -E.Jm V dt, so the charge current is J_phys = -Jm; L_z is the cell
height along the vacuum axis, so Jm*L_z is the current per unit width.)

"field before" = E_inc (column E_ext), "field after" = E_t.

Power (FLUENCE-integrated, Parseval-exact -- the primary numbers):
    T = n_sub int E_t^2 / int E_inc^2,  R = int E_r^2 / int E_inc^2,  A = 1 - T - R
with the pointwise identity (c/4pi)(E_inc^2 - E_t^2 - E_r^2) = E_t J_s
(n_sub = 1; unit-tested), i.e. A is the Joule absorption of the sheet in the
TRANSMITTED (local) field. The SBE, however, was driven by E_inc, and its own
energy ledger gives A_E = int E_inc J_s / F_inc. Exactly:

    A = A_E - S_rr,     S_rr = (Z0/2) int J_s^2 / F_inc     (radiation reaction)

S_rr is the self-consistency correction the single-cell driver lacks. On a
k-mesh that resolves the resonance shell it is O((Z0 sigma)^2/2) ~ 3e-4 for the
universal sigma = e^2/4hbar (negligible against the 2.3 % absorption); on an
UNRESOLVED mesh a few near-resonant k-points produce a large REACTIVE current
and S_rr becomes comparable to A itself -- the printed S_rr is therefore the
reliability flag of the perturbative estimate (the exact routes: a radiation-
reaction term in the single-cell driver, or theory='maxwell_sbe' with hx_m =
L_z; see README). A single-FFT-bin "T at the carrier" is NOT reported: pulse
reshaping moves spectral weight between bins and such a ratio is not bounded
by 1; instead the BAND-integrated values over the incident pulse's FWHM band
and the band-averaged Re sigma / sigma_univ (= 1 for the universal sheet) are
given as secondary diagnostics.

Usage
    python3 transmission.py runs/*/graphene_sit_sbe_rt.data
            [--lz-bohr 37.794523] [--area-bohr2 18.7154] [--n-sub 1.0] [--plot]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

C_AU = 137.035999084          # speed of light [a.u.]
AU_E_VM = 5.14220675e11       # 1 a.u. of electric field in V/m
AU_T_FS = 0.02418884326505    # 1 a.u. of time in fs
AU_EV = 27.211386245988
BOHR_ANG = 0.52917721067
Z0_SI = 376.730313668         # vacuum impedance [Ohm]
SIGMA_UNIV = 0.25             # e^2/(4 hbar) = 1/4 a.u. (universal interband sheet conductance)
# graphene primitive cell (a = 2.46 Ang, 20 Ang vacuum) -- defaults for the x14 inputs
A_BOHR = 4.648726
AREA_BOHR2 = A_BOHR**2 * np.sqrt(3.0) / 2.0     # 18.7154
LZ_BOHR = 37.794523
VF_EPM_AU = 0.439             # 43-PW Ramanujam pi-EPM cone slope (0.96e6 m/s), tests/test_graphene_dirac_levels.py
B_LEN_AU = 4.0 * np.pi / (np.sqrt(3.0) * A_BOHR)   # |b1| = |b2| of the hexagonal mesh


# ---------------------------------------------------------------- physics core
def sheet_fields(E_inc, J_sheet, n_sub=1.0):
    """Transmitted / reflected fields of a current sheet at normal incidence.

    E_inc, J_sheet: arrays (a.u.); J_sheet = physical charge current per unit
    width. n_sub: refractive index of a semi-infinite substrate behind the
    sheet (1 = free-standing). Returns (E_t, E_r)."""
    Z0 = 4.0 * np.pi / C_AU
    E_t = (2.0 * E_inc - Z0 * J_sheet) / (1.0 + n_sub)
    E_r = E_t - E_inc
    return E_t, E_r


def fluence_tra(t, E_inc, E_t, E_r, n_sub=1.0):
    """Fluence-integrated (T, R, A) -- the primary, Parseval-exact numbers."""
    Fi = np.trapezoid(E_inc**2, t)
    Ft = n_sub * np.trapezoid(E_t**2, t)
    Fr = np.trapezoid(E_r**2, t)
    return Ft / Fi, Fr / Fi, 1.0 - (Ft + Fr) / Fi


def radiation_reaction_term(t, E_inc, J_sheet):
    """S_rr = (Z0/2) int J_s^2 dt / F_inc  (= A_E - A exactly for n_sub = 1)."""
    Z0 = 4.0 * np.pi / C_AU
    return (Z0 / 2.0) * np.trapezoid(J_sheet**2, t) / ((1.0 / Z0) * np.trapezoid(E_inc**2, t))


def band_tra(t, E_inc, E_t, E_r, J_sheet, n_sub=1.0, halfwidth=None):
    """Band-integrated (T, R, A, Re sigma/sigma_univ, omega0) over the incident
    pulse's FWHM band around its spectral peak (secondary diagnostic)."""
    dt = t[1] - t[0]
    w = 2.0 * np.pi * np.fft.rfftfreq(len(t), d=dt)
    Ei = np.fft.rfft(E_inc)
    Et = np.fft.rfft(E_t)
    Er = np.fft.rfft(E_r)
    Js = np.fft.rfft(J_sheet)
    P = np.abs(Ei)**2
    i0 = np.argmax(P[1:]) + 1
    w0 = w[i0]
    if halfwidth is None:
        above = np.where(P >= 0.5 * P[i0])[0]
        halfwidth = max(w[above].max() - w0, w0 - w[above].min(), w[1] - w[0])
    m = np.abs(w - w0) <= halfwidth
    T = n_sub * np.sum(np.abs(Et[m])**2) / np.sum(P[m])
    R = np.sum(np.abs(Er[m])**2) / np.sum(P[m])
    re_sig = np.sum((Js[m] * np.conj(Ei[m])).real) / np.sum(P[m])
    # complex band-averaged sheet conductivity referred to the LOCAL (transmitted)
    # field -- the least-squares sigma of J_s = sigma E_t over the band. Its real part
    # is the dissipative conductance, its imaginary part the inductive (Drude
    # reactive) one; both are needed for the sheet impedance z = Z0 sigma.
    den = np.sum(np.abs(Et[m])**2)
    sig_c = (np.sum(Js[m] * np.conj(Et[m])) / den) if den > 0 else 0.0 + 0.0j
    return T, R, 1.0 - T - R, re_sig / SIGMA_UNIV, w0, sig_c / SIGMA_UNIV


def stack_from_one_layer(t, E_inc, E_t, J_sheet, n_layers, n_sub=1.0, floor=1e-8):
    """Fluence transmission of N identical decoupled sheets, from a ONE-layer run.

    The single-frequency estimate T_N = |2/(2+N z)|^2 with one band-averaged z is too
    crude for a single-cycle transient: sigma(omega) varies across a band as wide as
    its own centre, and the fluence ratio is not the ratio at the centroid. Here the
    sheet's own frequency-resolved response is used instead,

        sigma(omega) = J_s(omega) / E_t(omega),     T_N = int |2 E_inc/(2 + N Z0 sigma)|^2 / int |E_inc|^2,

    which reproduces the explicitly propagated two-layer run to ~1 % where the single-
    frequency formula is off by 12 %. Bins with negligible incident weight are dropped
    (`floor`, relative to the spectral peak) so the division does not amplify noise.
    Valid while the response is LINEAR. Each layer of a real stack sees the field
    reduced by both currents, so once drift saturation is under way the prediction
    fails in a definite direction: the weaker local field pushes every layer back
    towards its more conductive small-u response, and the true stack is DARKER than
    predicted. Against explicit two-layer runs at E_F = 0.6 eV, 72^2 (u = A_0/k_F):

        E_0 [kV/cm]      3      10      30     100     300
        u             0.04    0.12    0.37    1.24    3.71
        predicted   0.4292  0.4239  0.3758  0.3868  0.6004
        propagated  0.4362  0.4341  0.4128  0.3778  0.4715

    i.e. ~2 % while u < 0.15 and useless by u ~ 4. Nothing about a stack of doped
    sheets in a strong field can be obtained by scaling one run."""
    Ei = np.fft.rfft(E_inc)
    Et = np.fft.rfft(E_t)
    Js = np.fft.rfft(J_sheet)
    P = np.abs(Ei)**2
    m = P > floor * P.max()
    z = np.zeros_like(Js)
    z[m] = (4.0 * np.pi / C_AU) * Js[m] / Et[m]
    tN = 2.0 / (2.0 + n_layers * z)
    return n_sub * np.sum(np.abs(tN[m] * Ei[m])**2) / np.sum(P[m])


def spectral_tra(t, E_inc, E_t, E_r, n_sub=1.0, w0=None):
    """Single-bin (T, R, A, omega) at the carrier. Kept for the unit test on a
    long quasi-monochromatic pulse; NOT used in the report (see module doc)."""
    dt = t[1] - t[0]
    w = 2.0 * np.pi * np.fft.rfftfreq(len(t), d=dt)
    Ei = np.fft.rfft(E_inc)
    Et = np.fft.rfft(E_t)
    Er = np.fft.rfft(E_r)
    i0 = (np.argmax(np.abs(Ei)[1:]) + 1) if w0 is None else int(np.argmin(np.abs(w - w0)))
    T = n_sub * abs(Et[i0])**2 / abs(Ei[i0])**2
    R = abs(Er[i0])**2 / abs(Ei[i0])**2
    return T, R, 1.0 - T - R, w[i0]


def linear_reference(sigma=SIGMA_UNIV, n_sub=1.0):
    """(T, R, A) of a LINEAR sheet with real conductance sigma [a.u.]."""
    Z0 = 4.0 * np.pi / C_AU
    tt = 2.0 / (1.0 + n_sub + Z0 * sigma)
    T = n_sub * tt**2
    R = (tt - 1.0)**2
    return T, R, 1.0 - T - R


def shell_resolution(hw_au, nk_lin, vf_au=VF_EPM_AU, b_len=B_LEN_AU):
    """Points per resonance-shell RADIUS on an nk_lin x nk_lin hexagonal MP mesh:
    k_res = hw/(2 v_F) divided by the mesh spacing |b|/nk_lin."""
    return (hw_au / (2.0 * vf_au)) / (b_len / nk_lin)


# ---------------------------------------------------------------- file readers
def _header_columns(path, key):
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') and re.search(r'\b1:' + key + r'\[', line):
                return {m.group(2): (int(m.group(1)) - 1, m.group(3))
                        for m in re.finditer(r'(\d+):([\w\-]+)\[([^\]]*)\]', line)}
    raise RuntimeError(f'no column header with "{key}" in {path}')


def _dedup_time(d, path):
    """Drop the seam a checkpoint-restart leaves behind.

    `yn_sbe_checkpoint_restart` appends, and the resume begins at the last checkpoint,
    which is a few steps BEFORE the record got to when the run was killed. The file then
    contains a short block of repeated times, and every integral downstream -- fluence,
    energy ledger, FFT -- assumes a uniform monotone grid. Keep the LAST occurrence of
    each timestamp (the resumed values supersede the pre-kill ones) and sort. A record
    that was never resumed is returned untouched.
    """
    t = d[:, 0]
    if np.all(np.diff(t) > 0):
        return d
    order = np.argsort(t, kind='stable')
    d = d[order]
    keep = np.append(np.diff(d[:, 0]) > 0, True)          # last of each repeated time
    n_dropped = len(d) - int(keep.sum())
    print(f'# {os.path.basename(os.path.dirname(path))}: checkpoint-restart seam, '
          f'{n_dropped} duplicated row(s) dropped')
    return d[keep]


def read_rt(path):
    """*_sbe_rt.data -> dict of columns converted to ATOMIC UNITS."""
    cols = _header_columns(path, 'Time')
    d = _dedup_time(np.loadtxt(path, comments='#'), path)
    out = {}
    for name, (i, unit) in cols.items():
        v = d[:, i].copy()
        if name == 'Time':
            if unit.startswith('fs'):
                v = v / AU_T_FS
        elif name.startswith('E_'):
            if 'V/Angstrom' in unit:
                v = v * 1e10 / AU_E_VM
        elif name.startswith('Ac_'):
            if 'fs' in unit:
                v = v * 1e10 / AU_E_VM / AU_T_FS
        elif name.startswith('Jm_'):
            if 'fs' in unit or 'Angstrom' in unit:
                v = v / ((BOHR_ANG / AU_T_FS) / BOHR_ANG**3)
        out[name] = v
    return out


def read_energy_delta(path):
    """Final Eall-Eall0 [Ha] from *_sbe_rt_energy.data (unit from the header)."""
    cols = _header_columns(path, 'Time')
    d = np.loadtxt(path, comments='#')
    for name, (i, unit) in cols.items():
        if name.startswith('Eall-Eall0'):
            v = d[-1, i]
            return v / AU_EV if unit.lower().startswith('ev') else v
    raise RuntimeError(f'no Eall-Eall0 column in {path}')


def count_kpoints(rt_path):
    """nk from the *_k.data sitting next to the rt file (None if absent)."""
    kfile = rt_path.replace('_sbe_rt.data', '_k.data')
    if not os.path.exists(kfile):
        return None
    n = 0
    with open(kfile) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith('#'):
                n += 1
    return n


# ---------------------------------------------------------------- analysis
def run_variable(rt_path, key, default=None):
    """Value of an input variable from the variables.log next to the rt file (None if absent)."""
    vlog = os.path.join(os.path.dirname(rt_path) or '.', 'variables.log')
    if not os.path.exists(vlog):
        return default
    pat = re.compile(r'^#\s*' + re.escape(key) + r'\s*=\s*(\S+)')
    with open(vlog) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                return m.group(1).strip("'\"")
    return default


class TruncatedRun(RuntimeError):
    """The *_sbe_rt.data record is shorter than the nt the run was asked for.

    Analysing a still-running (or crashed) job silently gives wrong physics: the
    fluence integrals stop mid-pulse, so T is meaningless and A is whatever fraction
    of the drive happened to have arrived. Every entry point raises this instead."""


def check_complete(path, nrows, allow_partial=False):
    """Raise TruncatedRun unless the record holds the nt steps of variables.log."""
    nt = run_variable(path, 'nt', None)
    try:
        nt = int(float(nt))
    except (TypeError, ValueError):
        return None                      # no variables.log next to the file: nothing to check
    if nrows < nt and not allow_partial:
        raise TruncatedRun(
            f'{path}: {nrows} rows against nt = {nt} ({100.0 * nrows / nt:.0f} % of the '
            f'requested propagation) -- the run is still going or died. Wait for it, or '
            f'pass allow_partial/--allow-partial if you really want the partial record.')
    return nt


def dark_fraction(J_s, J_dark):
    """|J_dark|max / |J_s|max -- how much of a driven run is the zero-field artifact.

    A dissipative run started from a Fermi-Dirac occupation carries a current even at
    zero field: the ring's fixed point is not exactly the initial distribution, and on
    a finite mesh its relaxation is not isotropic. That current is FIELD-INDEPENDENT,
    so it matters in proportion to how weak the drive is. Measured at 72^2, E_F = 0.6
    eV (|J_dark|max = 2.42e-8 a.u.): 36 % of the peak current at 3 kV/cm, 14 % at 10,
    3.3 % at 30, 0.6 % at 300.

    The effect on T, from subtracting the dark run and re-applying the sheet BC:
    +0.047 at 3 kV/cm, +0.004 at 10, below 0.001 from 30 kV/cm up. Above ~10 % the
    point is unusable -- and NOT repairable by subtraction either, because the dark
    current of a driven run is not the field-free one (the field changes the
    distribution the ring acts on), so the first-order correction over-shoots.
    """
    p = np.max(np.abs(J_s))
    return float(np.max(np.abs(J_dark)) / p) if p > 0 else np.nan


def analyze(path, lz=LZ_BOHR, area=AREA_BOHR2, n_sub=1.0, n_layers=None, allow_partial=False,
            dark=None):
    d = read_rt(path)
    t = d['Time']
    check_complete(path, len(t), allow_partial)
    ax = max('xyz', key=lambda a: np.max(np.abs(d[f'E_ext_{a}'])))
    E_inc = d[f'E_ext_{ax}']
    Jm = d[f'Jm_{ax}']
    # sbe_sheet_nlayers identical decoupled sheets in the same local field: J_s = -nlayers L_z Jm,
    # and the per-cell energy ledger counts one layer. Auto-read from variables.log unless given.
    if n_layers is None:
        n_layers = int(run_variable(path, 'sbe_sheet_nlayers', 1) or 1)
    J_s = -Jm * lz * n_layers
    E_bc, E_rbc = sheet_fields(E_inc, J_s, n_sub)           # perturbative sheet BC
    E_totc = d.get(f'E_tot_{ax}', E_inc)
    selfc = np.max(np.abs(E_totc - E_inc)) > 1e-9 * max(np.max(np.abs(E_inc)), 1e-300)
    if selfc:
        # yn_sbe_sheet_field run: the driver already propagated in the local field,
        # E_tot IS the transmitted field (free-standing sheet, normal incidence). The
        # BC reconstruction is kept as a consistency check (differs by the explicit lag).
        E_t = E_totc
        E_r = E_t - E_inc
        dEt = np.max(np.abs(E_t - E_bc)) / max(np.max(np.abs(E_inc)), 1e-300)
    else:
        E_t, E_r, dEt = E_bc, E_rbc, 0.0
    T, R, A = fluence_tra(t, E_inc, E_t, E_r, n_sub)
    Tb, Rb, Ab, resig, w0, sig_c = band_tra(t, E_inc, E_t, E_r, J_s, n_sub)
    E0 = float(np.max(np.abs(E_inc)))
    res = dict(file=path, axis=ax, mode='SC' if selfc else 'pert', dEt=dEt, E0_au=E0, E0_kvcm=E0 * AU_E_VM / 1e5, n_layers=n_layers,
               I_wcm2=(E0 * AU_E_VM)**2 / (2.0 * Z0_SI) / 1e4, hw_ev=w0 * AU_EV,
               T=T, R=R, A=A, T_band=Tb, R_band=Rb, A_band=Ab, resig=resig, sig_c=sig_c,
               S_rr=radiation_reaction_term(t, E_inc, J_s), A_energy=np.nan, nk=count_kpoints(path))
    if dark is not None:
        dd = read_rt(dark)
        Jd = -dd[f'Jm_{ax}'] * lz * n_layers
        m = min(len(Jd), len(J_s))
        res['dark_frac'] = dark_fraction(J_s, Jd[:m])
        Etd, Erd = sheet_fields(E_inc[:m], J_s[:m] - Jd[:m], n_sub)
        res['T_dark_sub'] = fluence_tra(t[:m], E_inc[:m], Etd, Erd, n_sub)[0]
    epath = path.replace('_sbe_rt.data', '_sbe_rt_energy.data')
    if os.path.exists(epath):
        dE = read_energy_delta(epath)                                    # Ha per cell
        F_inc = (C_AU / (4.0 * np.pi)) * np.trapezoid(E_inc**2, t)      # Ha per bohr^2
        res['A_energy'] = n_layers * dE / (area * F_inc) if F_inc > 0 else np.nan
    if res['nk']:
        res['shell_pts'] = shell_resolution(w0, int(round(np.sqrt(res['nk']))))
    return res, (t, E_inc, E_t, E_r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+', help='*_sbe_rt.data files (globs ok)')
    ap.add_argument('--lz-bohr', type=float, default=LZ_BOHR, help='cell height along the vacuum axis [bohr]')
    ap.add_argument('--area-bohr2', type=float, default=AREA_BOHR2, help='2D cell area [bohr^2] (energy cross-check)')
    ap.add_argument('--n-sub', type=float, default=1.0, help='substrate refractive index (1 = free-standing)')
    ap.add_argument('--n-layers', type=int, default=None,
                    help='identical decoupled sheets in the same field (default: sbe_sheet_nlayers of the run, else 1); '
                         'pert mode only -- a self-consistent run already carries it in E_tot')
    ap.add_argument('--plot', action='store_true', help='write transmission_scan.png (needs matplotlib)')
    ap.add_argument('--predict-layers', type=int, default=None, metavar='N',
                    help='also print the fluence transmission of N identical decoupled '
                         'sheets predicted from this ONE-layer run via its own '
                         'frequency-resolved sigma(omega) (see stack_from_one_layer)')
    ap.add_argument('--dark', default=None, metavar='DARK_RT',
                    help='the zero-field control run of the same mesh and doping. Reports what '
                         'fraction of each driven current is the field-independent ring artifact '
                         'and what T becomes without it. Above ~10 %% the point is unusable (see '
                         'dark_fraction)')
    ap.add_argument('--allow-partial', action='store_true',
                    help='analyse records shorter than the run\'s own nt (default: refuse -- '
                         'the fluence integrals of a half-finished run are meaningless)')
    args = ap.parse_args(argv)

    files = sorted(sum((glob.glob(f) for f in args.files), []))
    if not files:
        sys.exit('no input files')
    Tl, Rl, Al = linear_reference(SIGMA_UNIV, args.n_sub)
    print(f'# sheet transmission  (L_z = {args.lz_bohr:.4f} bohr, A_2D = {args.area_bohr2:.4f} bohr^2, n_sub = {args.n_sub})')
    print(f'# linear universal-sheet reference (sigma = e^2/4hbar): T = {Tl:.5f}  R = {Rl:.2e}  A = {Al:.5f}'
          f'   [pi*alpha = {np.pi / C_AU:.5f}]')
    print('# mode: SC = yn_sbe_sheet_field run (E_tot IS the transmitted field; dEt = |E_tot - E_BC|/E0 consistency),'
          ' pert = sheet BC applied to Jm of an E_ext-driven run')
    print('# T,R,A = fluence-integrated (primary); A_E = SBE energy ledger; S_rr = (Z0/2) int J_s^2/F'
          ' (pert: A_E - A = S_rr, must be << A for the perturbative estimate to hold)')
    print('# T_band, Re(sigma)/sigma_univ = over the incident FWHM band (secondary); shell = mesh points per resonance-shell radius')
    print(f'{"E0[kV/cm]":>10} {"I[W/cm2]":>10} {"hw[eV]":>7} {"mode":>4} {"dEt":>8} {"T":>9} {"R":>9} {"A":>9} {"A_E":>9} '
          f'{"S_rr":>9} {"T_band":>8} {"Res/s0":>7} {"shell":>6}  file')
    rows = []
    for f in files:
        try:
            r, _ = analyze(f, args.lz_bohr, args.area_bohr2, args.n_sub, args.n_layers,
                           allow_partial=args.allow_partial, dark=args.dark)
        except TruncatedRun as exc:
            print(f'# SKIPPED (incomplete): {exc}')
            continue
        rows.append(r)
        if args.predict_layers:
            dd = read_rt(f)
            axx = r['axis']
            Ei = dd[f'E_ext_{axx}']
            Ett = dd.get(f'E_tot_{axx}', Ei)
            Jss = -dd[f'Jm_{axx}'] * args.lz_bohr * r['n_layers']
            r['T_stack'] = stack_from_one_layer(dd['Time'], Ei, Ett, Jss,
                                                args.predict_layers / r['n_layers'], args.n_sub)
        sp = r.get('shell_pts', np.nan)
        print(f'{r["E0_kvcm"]:10.3f} {r["I_wcm2"]:10.3e} {r["hw_ev"]:7.3f} {r["mode"]:>4} {r["dEt"]:8.1e} {r["T"]:9.6f} '
              f'{r["R"]:9.2e} {r["A"]:9.6f} {r["A_energy"]:9.6f} {r["S_rr"]:9.2e} {r["T_band"]:8.5f} {r["resig"]:7.3f} '
              f'{sp:6.2f}  ' + (f'T{args.predict_layers}={r["T_stack"]:.5f}  ' if 'T_stack' in r else '')
              + (f'dark={100 * r["dark_frac"]:.1f}%{"!" if r["dark_frac"] > 0.10 else ""} '
                 f'T_sub={r["T_dark_sub"]:.5f}  ' if 'dark_frac' in r else '') + f'{f}')
    # In a self-consistent run E_tot already carries the radiation reaction, so the
    # fluence absorption and the SBE energy ledger must agree. They do to ~1 % on a
    # clean run; a large gap means the current is carrying something the field did not
    # put there (the zero-field ring drift is the usual culprit -- see dark_fraction).
    bad = [r for r in rows if r['mode'] == 'SC' and np.isfinite(r['A_energy'])
           and abs(r['A'] - r['A_energy']) > 0.25 * max(abs(r['A']), 1e-9)]
    if bad:
        print('# *** FLUENCE / LEDGER MISMATCH (self-consistent runs, should agree to ~1 %):')
        for r in bad:
            print(f'#     {r["E0_kvcm"]:.0f} kV/cm: A = {r["A"]:.4f} from fluences but '
                  f'A_E = {r["A_energy"]:.4f} from the electron energy '
                  f'({abs(r["A"] / max(r["A_energy"], 1e-12)):.1f}x)')
        print('#     The current holds something the drive did not put there. Run the dark control.')
    neg = [r for r in rows if r['A'] < -1e-6]
    if neg:
        print('# *** NEGATIVE ABSORPTION in ' + ', '.join(f'{r["E0_kvcm"]:.0f} kV/cm (A = {r["A"]:+.3f})'
                                                          for r in neg))
        print('#     A passive sheet cannot emit. This means the electronic state is INVERTED and the'
              '\n#     sheet has gain -- check the ground state for "ZERO partially occupied k-points"'
              '\n#     (a doping the mesh cannot represent) before using any number from that run.')
    if rows and rows[0].get('nk'):
        sp = rows[0].get('shell_pts', np.nan)
        flag = 'RESOLVED' if sp >= 3 else ('marginal' if sp >= 1 else 'UNRESOLVED -- the mesh sees a few discrete'
                                                                   ' near-resonant k-points, not the shell')
        print(f'# resolution advisory: {rows[0]["nk"]} k-points -> {sp:.2f} mesh points per resonance-shell radius '
              f'at hw = {rows[0]["hw_ev"]:.3f} eV: {flag}')
    if args.plot and len(rows) >= 2:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print('# (matplotlib not available -- no plot)')
            return 0
        E = np.array([r['E0_kvcm'] for r in rows])
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
        ax[0].semilogx(E, [r['T'] for r in rows], 'o-', label='T (sheet BC, fluence)')
        ax[0].axhline(Tl, ls='--', c='gray', label='linear universal sheet')
        ax[0].set_xlabel('peak E [kV/cm]'); ax[0].set_ylabel('transmission'); ax[0].legend(fontsize=8)
        ax[1].semilogx(E, [r['A'] for r in rows], 's-', label='A (sheet BC)')
        ax[1].semilogx(E, [r['A_energy'] for r in rows], 'x--', label='A_E (energy ledger)')
        ax[1].semilogx(E, [r['S_rr'] for r in rows], '^:', label='S_rr (radiation reaction)')
        ax[1].axhline(Al, ls='--', c='gray', label=r'$\pi\alpha$ sheet')
        ax[1].set_xlabel('peak E [kV/cm]'); ax[1].set_ylabel('absorption'); ax[1].legend(fontsize=8)
        fig.suptitle('graphene sheet: field before/after -> transmission')
        fig.tight_layout()
        fig.savefig('transmission_scan.png', dpi=140)
        print('# wrote transmission_scan.png')
    return 0


if __name__ == '__main__':
    sys.exit(main())
