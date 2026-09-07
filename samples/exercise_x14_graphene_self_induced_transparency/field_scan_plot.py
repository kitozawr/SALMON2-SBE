#!/usr/bin/env python3
"""The x14 headline figure: transmission of a DOPED and an INTRINSIC graphene sheet
against the peak driving field, on the same mesh and the same pulse, plus the sheet
conductivity of the doped run against the value a measurement implies.

    python3 field_scan_plot.py --doped 'runs_doped/*/graphene_sit_sbe_rt.data' \
                               --intrinsic 'runs_intrinsic/*/graphene_sit_sbe_rt.data' \
                               --t-meas 0.60 0.70 --n-sub 1.65 --out doped_vs_intrinsic.png

Left panel: T(E0) for both occupations. Only the initial occupation differs, so the
two curves isolate what the doping carriers do. The intrinsic sheet darkens
monotonically (Landau-Zener pair creation, wiki/12 sec. 7); the doped sheet is flat
in the linear regime and BRIGHTENS once the vector-potential excursion A_0 exceeds
the Fermi radius k_F -- the drift velocity saturates at v_F and the chord
conductivity falls as k_F/A_0 (wiki/12 sec. 4a.5). The saturation field A_0 = k_F is
drawn from the run's own E_F. Any dip BEFORE that rise is the k-mesh's
representation of the Fermi disc, not physics: --continuum draws the continuum
drift-saturation prediction so the residual mesh artifact is read off directly
(wiki/12 sec. 4a.5.5).

Right panel: Re sigma of the doped sheet over the incident band, with the sheet
conductance(s) implied by the measured transmissions passed in --t-meas (the
substrate's Fresnel loss is divided out with --n-sub). This is the panel that
compares calculation and experiment: both are sheet conductances in units of
sigma_univ = e^2/4hbar, no fitting in between.

Both run sets must be the same mesh, pulse and dt -- otherwise the comparison is
between meshes, not between occupations.
"""
import argparse
import glob
import re
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transmission import analyze, TruncatedRun, AU_EV, VF_EPM_AU  # noqa: E402
from drude_check import sheet_from_transmission, Z0, run_variable  # noqa: E402
from drift_saturation import g_continuum  # noqa: E402

SIGMA_UNIV = 0.25
A0_PER_KVCM = 6.213e-4      # peak |A| [a.u.] per kV/cm of the scaled DAST transient
# Beyond this displacement the "displaced Fermi disc on a cone" picture behind G(u)
# stops being about a cone at all: at u = 12 the excursion is 40 % of |b|, the EPM
# band is strongly warped, and the sheet keeps a conductivity the drift model says
# should have gone. The continuum overlay is therefore drawn only up to here.
CONE_U_MAX = 2.0


def collect(patterns):
    """(E0 [kV/cm], T, R, A, Re sigma/sigma_univ, E_F [eV], sigma_c) per run, by field.

    Returns (array of the first five columns, E_F, list of complex band sigma /
    sigma_univ referred to the LOCAL field -- the anchor of the continuum curve)."""
    files = sorted(sum((glob.glob(p) for p in patterns), []))
    rows = []
    for f in files:
        try:
            r, _ = analyze(f)
        except TruncatedRun as exc:
            print(f'# SKIPPED, this field is missing from the figure: {exc}')
            continue
        ef = float(run_variable(f, 'sbe_ef_ev', 0.0) or 0.0)
        tk = float(run_variable(f, 'sbe_temp_init_k', 300.0) or 300.0)
        rows.append((r['E0_kvcm'], r['T'], r['R'], r['A'], r['resig'], ef, r['sig_c'], tk,
                     r['n_layers'], r['nk']))
    rows.sort(key=lambda r: r[0])
    return (np.array([r[:5] for r in rows]),
            (rows[0][5] if rows else 0.0),
            [r[6] for r in rows],
            (rows[0][7] if rows else 300.0),
            (rows[0][8] if rows else 1),
            (rows[0][9] if rows else None))


def continuum_curve(e0, sig_c0, ef_ev, t_init_k, t0=None, e0_anchor=None):
    """Continuum drift-saturation prediction of T(E_0), anchored on one field.

    The anchor is the lowest field given unless e0_anchor names another one -- which it
    must whenever the lowest field of a scan is not the one being trusted (a field its
    dark control condemned, say), and whenever e0 is a plotting grid rather than the run
    fields themselves. sig_c0 and t0 have to come from that same anchor run.

    The velocity-gauge displacement k -> k + A turns the direction of a Dirac-cone
    velocity but not its modulus, so the chord response of the whole doped disc is
    sigma(E_0) = sigma_lin * G(u)/u with u = A_0/k_F (wiki/12 Eqs. 4a.9-4a.13). Only
    that field dependence is the model. The anchor is the run's OWN lowest-field
    point: the PHASE of z = Z_0 sigma is the band-averaged complex conductivity of
    that run (how resistive vs how inductive the sheet is), and its MODULUS is fixed
    by inverting the sheet BC on that run's own fluence transmission,

        T = |2/(2+z)|^2   ->   |z| = 2[-cos(phi) + sqrt(cos^2(phi) - 1 + 1/T)],

    so the curve starts exactly on the data and everything after it is prediction.
    The gap that opens at higher fields is what the mesh (and, at large u, trigonal
    warping and pair creation) adds to the continuum answer.
    """
    kF = (abs(ef_ev) / AU_EV) / VF_EPM_AU
    if kF <= 0:
        return None
    zb = Z0 * SIGMA_UNIV * sig_c0                     # band-averaged complex impedance
    cphi = (zb.real / abs(zb)) if abs(zb) > 0 else 1.0
    if t0 is None:
        r0 = abs(zb)
    else:
        disc = cphi**2 - 1.0 + 1.0 / max(t0, 1e-12)
        r0 = 2.0 * (-cphi + np.sqrt(max(disc, 0.0)))
    z0 = r0 * (zb / abs(zb)) if abs(zb) > 0 else r0 + 0j
    u = np.asarray(e0, float) * A0_PER_KVCM / kF
    ef_au = abs(ef_ev) / AU_EV
    u_a = (float(e0_anchor) * A0_PER_KVCM / kF) if e0_anchor is not None else float(u[0])
    g0 = g_continuum(u_a, t_init_k, ef_au)
    ratio = np.array([g_continuum(float(x), t_init_k, ef_au) / g0 for x in u])
    return np.abs(2.0 / (2.0 + z0 * ratio))**2


def report_residual(label, e0, tdat, tcont, kF):
    """Per-point residual of the data against the continuum drift curve.

    Split at u = 1 because the two sides miss for different reasons: below the
    saturation field a shortfall is the k-mesh artifact of the Fermi disc
    (README SS7.11), above it the cone stops being a cone (trigonal warping,
    Landau-Zener pairs) and, if the ring is on, momentum relaxation keeps the
    disc from following A(t) rigidly -- which is the one assumption G(u) makes.
    Points past CONE_U_MAX are not reported: the curve is not drawn there either.
    """
    e0 = np.asarray(e0, float)
    u = e0 * A0_PER_KVCM / kF
    dev = 100.0 * (np.asarray(tdat, float) - tcont) / np.maximum(tcont, 1e-12)
    keep = u <= CONE_U_MAX
    print(f'# residual T_data / T_continuum - 1, anchored on this series own lowest '
          f'field E_0 = {e0[0]:g} kV/cm  [{label}]')
    for ui, ei, di in zip(u[keep], e0[keep], dev[keep]):
        print(f'#   E_0 = {ei:8.1f} kV/cm   u = {ui:6.3f}   {di:+7.2f} %')
    lo = keep & (u <= 1.0)
    hi = keep & (u > 1.0)
    if lo.any():
        print(f'#   largest DIP below the curve at u <= 1 (k-mesh artifact of the '
              f'Fermi disc): {np.min(dev[lo]):+.1f} %')
    if hi.any():
        print(f'#   deviation at u > 1 (warping / pair creation / momentum relaxation, '
              f'NOT the mesh): {np.min(dev[hi]):+.1f} %')
    if (~keep).any():
        print(f'#   {int((~keep).sum())} point(s) past u = {CONE_U_MAX:g} not compared: '
              f'the displaced-disc picture has left the cone')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--doped', nargs='+', required=True, help='*_sbe_rt.data of the doped scan (globs ok)')
    ap.add_argument('--intrinsic', nargs='*', default=[], help='*_sbe_rt.data of the undoped control scan')
    ap.add_argument('--t-meas', nargs='*', type=float, default=[], help='measured transmissions, substrate included')
    ap.add_argument('--series', nargs='*', default=[],
                    help='extra doped series to overlay, each "label:glob[:glob...]" (e.g. '
                         '"with e-ph ring:diss/runs/*/..._rt.data"); several colon-separated '
                         'globs are unioned, so a series can be built from a subset of a run '
                         'directory. Plotted on the T and sigma panels beside --doped')
    ap.add_argument('--n-sub', type=float, default=1.65, help='substrate index for --t-meas (PET ~ 1.65)')
    ap.add_argument('--continuum', action='store_true',
                    help='overlay the parameter-free continuum drift-saturation curve '
                         'T(sigma_lin G(u)/u), anchored on the lowest-field run: the gap '
                         'to it is the k-mesh artifact of the Fermi disc')
    ap.add_argument('--excluded', nargs='*', default=[],
                    help='runs to DRAW but never fit: hollow markers, no line, no '
                         'continuum anchor. For fields a dark control condemned '
                         '(dark/peak > 10 %%, README SS7.11) -- they belong in the '
                         'figure as evidence, not in the comparison')
    ap.add_argument('--out', default='doped_vs_intrinsic.png')
    ap.add_argument('--title', default=None)
    ap.add_argument('--doped-label', default=None,
                    help="legend label for the --doped series (default: 'doped, E_F = ... (metal)')")
    args = ap.parse_args(argv)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available -- no figure'); return 0

    dop, ef, dsig, tinit, nlay, nk = collect(args.doped)
    if dop.size == 0:
        print('# no doped runs matched'); return 1
    ins = collect(args.intrinsic)[0] if args.intrinsic else np.empty((0, 5))
    kF = (abs(ef) / AU_EV) / VF_EPM_AU if ef else 0.0
    e_sat = kF / A0_PER_KVCM if kF else None

    npan = 3 if ins.size else 2
    extra = []
    for spec in args.series:
        lab, _, pat = spec.partition(':')
        # the remainder may hold several globs separated by colons or whitespace, so a
        # series can be assembled from a subset of a run directory (e.g. dropping a
        # field whose dark control showed it contaminated)
        arr, _ef, asig, atk, _nl, _nk = collect([g for g in re.split(r'[:\s]+', pat) if g])
        if arr.size:
            extra.append((lab, arr, asig, _ef or ef, atk))
    exc = collect(args.excluded)[0] if args.excluded else np.empty((0, 5))
    fig, ax = plt.subplots(1, npan, figsize=(5.5 * npan, 4.0))
    dlab = args.doped_label or f'doped, $E_F$ = {ef:g} eV (metal)'
    ax[0].semilogx(dop[:, 0], dop[:, 1], 'o-', color='#c0392b', label=dlab)
    for (lab, arr, asig, aef, atk), col, mk in zip(extra, ('#e67e22', '#16a085', '#8e44ad'), ('D-', 'v-', '^-')):
        ax[0].semilogx(arr[:, 0], arr[:, 1], mk, ms=5, color=col, label=lab)
        if args.continuum:
            tc = continuum_curve(arr[:, 0], asig[0], aef, atk, t0=arr[0, 1])
            if tc is not None:
                akF = (abs(aef) / AU_EV) / VF_EPM_AU if aef else kF
                mv = (arr[:, 0] * A0_PER_KVCM / akF <= CONE_U_MAX) if akF else slice(None)
                ax[0].semilogx(arr[mv, 0], tc[mv], ':', lw=1.4, color=col, alpha=0.75)
                if akF:
                    report_residual(lab, arr[:, 0], arr[:, 1], tc, akF)
    uu = dop[:, 0] * A0_PER_KVCM / kF if kF else np.zeros(len(dop))
    if args.continuum:
        tc = continuum_curve(dop[:, 0], dsig[0], ef, tinit, t0=dop[0, 1])
        if tc is not None:
            mv = (uu <= CONE_U_MAX) if kF else slice(None)
            ax[0].semilogx(dop[mv, 0], tc[mv], ':', lw=1.6, color='#c0392b',
                           label=f'continuum drift saturation $G(u)/u$, $u\\leq{CONE_U_MAX:g}$')
            report_residual(dlab, dop[:, 0], dop[:, 1], tc, kF)
    if exc.size:
        ax[0].semilogx(exc[:, 0], exc[:, 1], 'o', ms=7, mfc='none', mew=1.4,
                       color='#7f8c8d', ls='none',
                       label='dark control condemned (not fitted)')
    if ins.size:
        ax[0].semilogx(ins[:, 0], ins[:, 1], 's-', color='#2c3e50', label='intrinsic (semimetal)')
    if e_sat:
        ax[0].axvline(e_sat, ls='--', c='#7f8c8d', lw=1)
        ax[0].axvspan(e_sat, max(dop[:, 0].max(), 1.0), color='#f1c40f', alpha=0.13)
        ax[0].annotate('$A_0 = k_F$\n(current saturation)',
                       xy=(e_sat * 1.15, 0.97), xycoords=('data', 'axes fraction'),
                       fontsize=8, va='top', color='#7f8c8d')
    zs_meas = []
    for tm in args.t_meas:
        zs, t_bare = sheet_from_transmission(tm, args.n_sub)
        zs_meas.append(zs / Z0 / SIGMA_UNIV)
        ax[0].axhline(tm / t_bare, ls=':', c='#27ae60', lw=1.2)
    if args.t_meas:
        ax[0].text(0.995, 0.02, 'green: measurement, substrate divided out',
                   transform=ax[0].transAxes, fontsize=7.5, color='#27ae60',
                   ha='right', va='bottom', zorder=5,
                   bbox=dict(fc='white', ec='none', alpha=0.75, pad=1.5))
    ax[0].set_xlabel('peak field $E_0$ [kV/cm]'); ax[0].set_ylabel('transmission $T$')
    ax[0].set_title(args.title or ('Same mesh, same pulse: only the initial occupation differs'
                                   if ins.size else 'Transmission of a doped Dirac sheet against the peak field'),
                    fontsize=9)
    ax[0].legend(fontsize=8, loc='upper left', framealpha=0.9); ax[0].grid(alpha=0.25)

    # middle panel: the extinction the DOPING carriers add, 1 - T_doped/T_intrinsic.
    # This is the observable that isolates them: the intrinsic Landau-Zener darkening
    # divides out, and what is left is the Drude response of the doped carriers, which
    # collapses once the drift saturates.
    if ins.size:
        tin = {round(r[0], 3): r[1] for r in ins}
        common = np.array([[r[0], 1.0 - r[1] / tin[round(r[0], 3)]] for r in dop
                           if round(r[0], 3) in tin])
        if common.size:
            ax[1].semilogx(common[:, 0], 100.0 * common[:, 1], 'o-', color='#8e44ad')
            ax[1].axhline(0.0, c='k', lw=0.6)
            if e_sat:
                ax[1].axvline(e_sat, ls='--', c='#7f8c8d', lw=1)
                ax[1].axvspan(e_sat, common[:, 0].max(), color='#f1c40f', alpha=0.13)
            ax[1].set_xlabel('peak field $E_0$ [kV/cm]')
            ax[1].set_ylabel(r'extinction added by the doping, $1-T_{\rm doped}/T_{\rm intrinsic}$  [%]')
            ax[1].set_title('The doping carriers alone: their extinction\npeaks at the saturation field and collapses',
                            fontsize=9)
            ax[1].grid(alpha=0.25)
    ia = npan - 1
    ax[ia].semilogx(dop[:, 0], dop[:, 4], 'o-', color='#c0392b', label=f'Re $\\sigma$, {dlab}')
    for (lab, arr, _s, _e, _t), col, mk in zip(extra, ('#e67e22', '#16a085', '#8e44ad'), ('D-', 'v-', '^-')):
        ax[ia].semilogx(arr[:, 0], arr[:, 4], mk, ms=5, color=col, label=f'Re $\\sigma$, {lab}')
    if exc.size:
        ax[ia].semilogx(exc[:, 0], exc[:, 4], 'o', ms=7, mfc='none', mew=1.4,
                        color='#7f8c8d', ls='none', label='not fitted')
    for z, tm in zip(zs_meas, args.t_meas):
        ax[ia].axhline(z, ls=':', c='#27ae60')
        ax[ia].annotate(f'measured {z:.1f} $\\sigma_{{univ}}$ (T = {100 * tm:.0f} %)',
                        xy=(0.985, z), xycoords=('axes fraction', 'data'), fontsize=7.5,
                        color='#27ae60', ha='right', va='bottom', zorder=5,
                        bbox=dict(fc='white', ec='none', alpha=0.75, pad=1.5))
    if e_sat:
        ax[ia].axvline(e_sat, ls='--', c='#7f8c8d', lw=1)
    if len(dop) > 1:
        imax = int(np.argmax(dop[:, 4]))
        drop = 100.0 * (dop[-1, 4] / dop[imax, 4] - 1.0)
        ax[ia].set_title(rf'Sheet conductivity: {drop:+.0f} % from its peak to {dop[-1, 0]:.0f} kV/cm', fontsize=9)
    ax[ia].set_xlabel('peak field $E_0$ [kV/cm]'); ax[ia].set_ylabel(r'Re $\sigma$ / $\sigma_{univ}$')
    ax[ia].legend(fontsize=8, loc='lower left', framealpha=0.9); ax[ia].grid(alpha=0.25)
    # ---- the provenance line the reader needs before comparing with a measurement:
    # what the calculation contains (no substrate, N sheets) and what was done to the
    # measured numbers (the substrate's Fresnel factor DIVIDED OUT, not added).
    layers = 'a single monolayer' if nlay == 1 else f'{nlay} electronically decoupled layers'
    note = ('FILM: the calculation is FREE-STANDING -- vacuum on both sides, no substrate '
            f'in the propagation -- and carries {layers} in one local field '
            '(sbe_sheet_nlayers).')
    if args.t_meas:
        tb = sheet_from_transmission(args.t_meas[0], args.n_sub)[1]
        note += (f'  The measured values have the substrate DIVIDED OUT (not added): '
                 f'T_sheet = T_meas / {tb:.3f}, n_sub = {args.n_sub:g}, so '
                 + ', '.join(f'{t:.2f}->{t / tb:.2f}' for t in args.t_meas) + '.')
    if nk:
        note += f'  Mesh: {nk} k-points.'
    fig.text(0.008, 0.008, note, fontsize=7.0, color='#34495e', wrap=True)
    fig.tight_layout(rect=(0, 0.045, 1, 1)); fig.savefig(args.out, dpi=150)
    print('# ' + note)
    print(f'# wrote {args.out}')

    print(f'{"E0[kV/cm]":>9} {"T doped":>9} {"T intr":>9} {"R doped":>9} {"A doped":>9} {"Re s/s0":>9}')
    tin = {round(r[0], 3): r[1] for r in ins} if ins.size else {}
    for r in dop:
        ti = tin.get(round(r[0], 3))
        print(f'{r[0]:9.2f} {r[1]:9.5f} ' + (f'{ti:9.5f} ' if ti is not None else f'{"--":>9} ')
              + f'{r[2]:9.2e} {r[3]:9.5f} {r[4]:9.3f}')
    if e_sat:
        print(f'# k_F = {kF:.5f} a.u. (E_F = {ef:g} eV) -> current saturation A_0 = k_F at E_0 = {e_sat:.1f} kV/cm')
    return 0


if __name__ == '__main__':
    sys.exit(main())
