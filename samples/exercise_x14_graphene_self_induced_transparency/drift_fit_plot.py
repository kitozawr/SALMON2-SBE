#!/usr/bin/env python3
"""One graphene layer against the analytic drift-saturation law.

The velocity-gauge displacement k -> k + A(t) turns the direction of a Dirac-cone
velocity but never its modulus, so the current of a rigidly displaced Fermi disc is

    J(A) = n e v_F G(u),   u = A_0 / k_F,        (wiki/12 Eqs. 4a.9-4a.13)

and the sheet response follows the CHORD ratio G(u)/u -- exactly 1 - u^2/8 while the
disc still contains the origin, 1/u times (1 - 1/8u^2) once it no longer does. That is
the whole model: no fitted parameter, one number (k_F) taken from the doping.

Two ways to test it, drawn side by side, because they fail differently:

  LEFT   T(E_0).  What an experiment sees, but a compressed view of the model: the
         sheet boundary condition T = |2/(2+z)|^2 is nonlinear in z, so the same
         fractional change in the sheet response shows up differently depending on
         where the sheet sits, and on the PHASE of z (a reactive sheet and a
         resistive one of equal |z| transmit differently).  Anchored on the lowest
         field of each series, so the curve starts on the data by construction.

  RIGHT  D(u)/D(u_0).  The model's actual claim, with the boundary condition and the
         phase divided out.  D is inverted from sigma(w) bin by bin -- no fit, no
         window (drude_check.analyze) -- and its companion tau(u) is drawn beside it:
         if tau moves with the field, a fall in the sheet conductance is scattering,
         not drift saturation, and the left panel cannot tell the two apart.

Only the SHAPE is tested.  The absolute Drude weight is not: a dissipative run can
carry a suppressed D (see the printed D/D_eq column) without that touching how D
depends on the field.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transmission import AU_EV, VF_EPM_AU, TruncatedRun          # noqa: E402
from drude_check import (analyze as drude_analyze, run_variable,  # noqa: E402
                         drude_weight_au, KB_AU, Z0)
from drift_saturation import g_continuum                         # noqa: E402
from field_scan_plot import (A0_PER_KVCM, CONE_U_MAX, SIGMA_UNIV,  # noqa: E402
                             collect, continuum_curve)

# tau larger than this is not a scattering time, it is a collisionless run whose
# sigma(w) has no measurable imaginary-to-real ratio; drawing it would only compress
# the axis of the runs that do relax.
TAU_COLLISIONLESS_FS = 1000.0


def series_rows(patterns):
    """(E0, T, D_spec [eV], tau_spec [fs], D_eq [eV], E_F, T_init) per run, by field."""
    files = sorted(set(sum((glob.glob(p) for p in patterns), [])))
    rows = []
    for f in files:
        try:
            r = drude_analyze(f, 0.6, 300.0)
        except TruncatedRun as exc:
            print(f'# SKIPPED: {exc}')
            continue
        rows.append((r['E0_kvcm'], r['T'], r['D_spec_ev'], r['tau_spec_fs'],
                     r['D_ev'], r['ef_ev'], r['temp_init']))
    rows.sort(key=lambda r: r[0])
    return np.array([r[:5] for r in rows]), (rows[0][5] if rows else 0.0), (rows[0][6] if rows else 300.0)


def parse_series(spec):
    lab, _, pat = spec.partition(':')
    return lab, [g for g in pat.split() if g]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--series', nargs='+', required=True,
                    help='"label:glob [glob ...]" -- one entry per curve. Several '
                         'space-separated globs are unioned so a series can drop the '
                         'fields a dark control condemned')
    ap.add_argument('--excluded', nargs='*', default=[],
                    help='runs to DRAW but never fit (hollow markers): fields whose '
                         'dark control puts the field-independent ring current above '
                         '10 %% of peak (README SS7.11)')
    ap.add_argument('--anchor-kvcm', type=float, default=None,
                    help='field every series is normalised at (default: each series own '
                         'lowest field). Give one when the series start at different '
                         'fields and you want the comparison to be like for like')
    ap.add_argument('--u-max', type=float, default=CONE_U_MAX,
                    help='draw the law out to this u (default %(default)g, where the '
                         'displaced disc leaves the cone). Give a larger value to see '
                         'the whole scan; everything past CONE_U_MAX is shaded, because '
                         'there the law is being shown, not tested')
    ap.add_argument('--predict-ef', nargs='*', type=float, default=[],
                    help='draw the analytic T(E_0) prediction for these dopings [eV] on the '
                         'left panel. No run needed: the linear sheet response is the '
                         'equilibrium Drude weight of (E_F, T), the field dependence is '
                         'G(u)/u with that doping own k_F, and the two go through the sheet '
                         'boundary condition. Only tau is assumed -- see --predict-tau-fs')
    ap.add_argument('--predict-dscale', nargs='*', type=float, default=[1.0],
                    help='scale(s) on the equilibrium Drude weight for --predict-ef. Give '
                         'two and the band between them is shaded: the runs carry only '
                         '0.69 D_eq at 0.6 eV (README SS7.16), an unexplained deficit that '
                         'a prediction from D_eq alone does not know about')
    ap.add_argument('--predict-tau-fs', type=float, default=45.0,
                    help='momentum-relaxation time for --predict-ef (default %(default)g fs, '
                         'the range the 0.6 eV runs measure)')
    ap.add_argument('--e0-max', type=float, default=None,
                    help='draw the LEFT panel (and the law on it) out to this peak field '
                         'in kV/cm, past the last run if you like: the law keeps bending '
                         'after the data stop, and that bend is most of what it says. '
                         'Independent of --u-max, which governs the right panel')
    ap.add_argument('--out', default='drift_fit_1layer.png')
    ap.add_argument('--title', default=None)
    args = ap.parse_args(argv)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available -- no figure')
        return 0

    fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.6))
    colours = ('#c0392b', '#e67e22', '#16a085', '#8e44ad')
    marks = ('o-', 'D-', 'v-', '^-')
    kF_ref = None
    summary = []
    umax = max(args.u_max, 0.1)
    e0_max = 0.0
    earr, _eef, _etk = series_rows(args.excluded) if args.excluded else (np.empty((0, 5)), 0.0, 300.0)
    e0_floor = float(earr[:, 0].min()) if earr.size else 0.0

    for (spec, col, mk) in zip(args.series, colours, marks):
        lab, pats = parse_series(spec)
        arr, ef, tk = series_rows(pats)
        if not arr.size:
            print(f'# no runs matched for "{lab}"')
            continue
        kF = (abs(ef) / AU_EV) / VF_EPM_AU
        kF_ref = kF_ref or kF
        u = arr[:, 0] * A0_PER_KVCM / kF
        e0_max = max(e0_max, float(arr[:, 0].max()))
        ia = 0 if args.anchor_kvcm is None else int(np.argmin(abs(arr[:, 0] - args.anchor_kvcm)))

        # ---- left: the transmission, and the model propagated through the sheet BC.
        # collect() is reused rather than re-derived: it carries the complex band
        # conductivity that fixes the PHASE of z, which the boundary condition needs
        # and the drift law says nothing about.
        ax[0].semilogx(arr[:, 0], arr[:, 1], mk, ms=5, color=col, label=lab)
        carr, cef, csig, ctk, _nl, _nk = collect(pats)
        # Drawn on a dense grid, not on the run fields, and carried DOWN through the
        # fields the dark control threw out: the law is anchored at one trusted field and
        # says something everywhere else, so the reader should be able to see where a
        # discarded point would have fallen.
        e_lo = min(float(arr[0, 0]), e0_floor if e0_floor else float(arr[0, 0]))
        e_hi = args.e0_max or float(arr[-1, 0]) * 1.2
        grid = np.geomspace(e_lo * 0.6, max(e_hi, float(arr[-1, 0]) * 1.05), 500)
        tc = continuum_curve(grid, csig[ia], cef, ctk, t0=carr[ia, 1],
                             e0_anchor=float(arr[ia, 0]))
        if tc is not None:
            ax[0].semilogx(grid, tc, ':', lw=1.6, color=col)
        e0_max = max(e0_max, float(grid[-1]))

        # ---- right: the model's own claim, D(u)/D(u_0) against G(u)/u.
        g = np.array([g_continuum(float(x), tk, abs(ef) / AU_EV) for x in u])
        gr = g / g[ia]
        dr = arr[:, 2] / arr[ia, 2]
        # Only the on-cone points are joined: beyond CONE_U_MAX the displaced disc has
        # left the cone the law is about, so a line through those points would invite
        # exactly the extrapolation the law does not make.
        shown = u <= umax
        tested = u <= CONE_U_MAX
        ax[1].plot(u[shown], dr[shown], mk, ms=6, color=col, label=f'{lab}')
        if (shown & ~tested).any():
            ax[1].plot(u[shown & ~tested], dr[shown & ~tested], mk[0], ms=7, mfc='white',
                       mew=1.6, ls='none', color=col)
        tau = arr[:, 3]
        if np.all(np.isfinite(tau)) and np.median(tau) < TAU_COLLISIONLESS_FS:
            ax[1].plot(u[shown], (tau / tau[ia])[shown], '--', lw=1.3, color='#7f8c8d',
                       marker='x', ms=5, alpha=0.9,
                       label='$\\tau(u)/\\tau(u_0)$ -- flat, so the fall is weight, not scattering')
        resid = 100 * (dr / gr - 1.0)
        fit = tested & (np.arange(len(u)) != ia)
        worst = np.max(np.abs(resid[fit])) if fit.any() else 0.0
        summary.append((lab, col, worst, int(fit.sum())))

        print(f'# ---- {lab}:  E_F = {ef:g} eV, k_F = {kF:.5f} a.u., '
              f'A_0 = k_F at {kF / A0_PER_KVCM:.0f} kV/cm, anchor {arr[ia, 0]:g} kV/cm')
        print(f'#  {"E0":>7} {"u":>6} {"T":>8} {"D[eV]":>8} {"D/D_eq":>7} {"tau[fs]":>8} '
              f'{"D/D_0":>7} {"G(u)/u":>7} {"resid":>8}')
        for i in range(len(arr)):
            flag = '' if u[i] <= CONE_U_MAX else '  (u > %g: off the cone)' % CONE_U_MAX
            tt = '   inf' if not np.isfinite(tau[i]) or tau[i] >= TAU_COLLISIONLESS_FS else f'{tau[i]:8.1f}'
            print(f'#  {arr[i, 0]:7.0f} {u[i]:6.3f} {arr[i, 1]:8.5f} {arr[i, 2]:8.4f} '
                  f'{arr[i, 2] / arr[i, 4]:7.3f} {tt:>8} {dr[i]:7.3f} {gr[i]:7.3f} '
                  f'{100 * (dr[i] / gr[i] - 1):+7.1f} %{flag}')

    if earr.size:
        ax[0].semilogx(earr[:, 0], earr[:, 1], 'o', ms=7, mfc='none', mew=1.4,
                       color='#7f8c8d', ls='none',
                       label='dark control condemned (not fitted)')

    # ---- the law as a PREDICTION at another doping ------------------------------
    # Nothing is fitted here and no run is needed: D_eq(E_F, T) fixes the linear sheet
    # response, G(u)/u with that doping's own k_F gives the field dependence, and the
    # sheet boundary condition turns the pair into a transmission. The saturation field
    # moves as k_F, i.e. linearly in E_F, so a lighter doping saturates earlier -- which
    # is the whole reason to want the curve before spending the runs.
    for ef_p, colp in zip(args.predict_ef, ('#2980b9', '#8e44ad', '#27ae60')):
        kF_p = (abs(ef_p) / AU_EV) / VF_EPM_AU
        e_hi = args.e0_max or (e0_max if e0_max else 1000.0)
        grid = np.geomspace(max(e_hi * 1e-3, 0.3), e_hi, 400)
        u_p = grid * A0_PER_KVCM / kF_p
        kT = KB_AU * 300.0
        D = drude_weight_au(abs(ef_p) / AU_EV, kT)            # Ha
        tau = args.predict_tau_fs / 0.0241888                  # fs -> a.u.
        w = 0.0108 / AU_EV                                     # the drive's band centre
        g0 = g_continuum(1e-6, 300.0, abs(ef_p) / AU_EV)
        gr = np.array([g_continuum(float(x), 300.0, abs(ef_p) / AU_EV) / g0 for x in u_p])
        curves = []
        for sc in (args.predict_dscale or [1.0]):
            sig = (sc * D / np.pi) / (1.0 / tau - 1j * w) * gr     # a.u., sheet
            curves.append(np.abs(2.0 / (2.0 + Z0 * sig))**2)
        lab = (f'PREDICTED, $E_F$ = {ef_p:g} eV ($A_0=k_F$ at '
               f'{kF_p / A0_PER_KVCM:.0f} kV/cm)')
        if len(curves) > 1:
            ax[0].fill_between(grid, np.min(curves, axis=0), np.max(curves, axis=0),
                               color=colp, alpha=0.18, lw=0, label=lab + ', $D$ band')
            for c in curves:
                ax[0].semilogx(grid, c, '-.', lw=1.2, color=colp, alpha=0.85)
        else:
            ax[0].semilogx(grid, curves[0], '-.', lw=1.6, color=colp, label=lab)
        ax[0].axvline(kF_p / A0_PER_KVCM, ls=':', lw=1, color=colp, alpha=0.7)

    # the analytic law itself, once, over the range where a displaced disc on a cone
    # is still what the sheet is doing
    if kF_ref:
        uu = np.linspace(1e-3, umax, 200)
        gg = np.array([g_continuum(float(x), 300.0, kF_ref * VF_EPM_AU) for x in uu])
        u0 = (args.anchor_kvcm or 1e-6) * A0_PER_KVCM / kF_ref
        g0 = g_continuum(max(u0, 1e-6), 300.0, kF_ref * VF_EPM_AU)
        ax[1].plot(uu, gg / g0, '-', lw=2.0, color='#2c3e50', alpha=0.8, zorder=0,
                   label='analytic $G(u)/u$ (no free parameter)')
        if umax > CONE_U_MAX:
            # past here the excursion is no longer small against the zone: the band is
            # warped, pairs are made, and G(u) is an extrapolation of a picture that has
            # stopped applying. Shaded, not hidden -- the reader should see where the
            # curve is a prediction and where it is a quotation.
            e_cone = CONE_U_MAX * kF_ref / A0_PER_KVCM
            ax[0].axvspan(e_cone, max(e_cone, e0_max), color='#95a5a6',
                          alpha=0.13, zorder=0)
            ax[1].axvspan(CONE_U_MAX, umax, color='#95a5a6', alpha=0.13, zorder=0)
            ax[1].annotate('off the cone: warping + LZ pairs\n(law quoted, not tested)',
                           xy=(CONE_U_MAX * 1.03, 0.72), xycoords=('data', 'axes fraction'),
                           fontsize=7.5, color='#7f8c8d')
        ax[0].axvline(kF_ref / A0_PER_KVCM, ls='--', c='#7f8c8d', lw=1)
        ax[0].annotate('$A_0 = k_F$', xy=(kF_ref / A0_PER_KVCM * 1.1, 0.96),
                       xycoords=('data', 'axes fraction'), fontsize=8, va='top', color='#7f8c8d')
    ax[1].axvline(1.0, ls='--', c='#7f8c8d', lw=1)
    ax[1].axhline(1.0, c='k', lw=0.6)

    ax[0].set_xlabel('peak field $E_0$ [kV/cm]')
    ax[0].set_ylabel('transmission $T$')
    ax[0].set_title(args.title or 'One layer: transmission, and the model through the sheet BC',
                    fontsize=9)
    ax[0].legend(fontsize=7.5, loc='upper left', framealpha=0.9)
    ax[0].grid(alpha=0.25)
    ax[1].set_xlabel('$u = A_0 / k_F$')
    ax[1].set_ylabel('$D(u)/D(u_0)$,  $\\tau(u)/\\tau(u_0)$')
    ax[1].set_xlim(0, umax)
    ax[1].set_title('The law itself: Drude weight against $G(u)/u$, with $\\tau$ beside it',
                    fontsize=9)
    ax[1].legend(fontsize=6.8, loc='lower left', framealpha=0.9)
    ax[1].grid(alpha=0.25)
    if summary:
        lines = [f'worst |measured / G(u)/u - 1| where the law applies ($u\\leq{CONE_U_MAX:g}$):']
        for lab, col, w, npt in summary:
            lines.append(f'  {lab}:  {w:.0f} %  ({npt} point{"s" if npt != 1 else ""})')
        ax[1].text(0.985, 0.97, '\n'.join(lines), transform=ax[1].transAxes,
                   fontsize=7.2, ha='right', va='top', color='#34495e',
                   bbox=dict(fc='white', ec='#bdc3c7', alpha=0.9, pad=3))

    note = ('Dotted curves on the left are the analytic law propagated through '
            'T = |2/(2+z)|^2 with the phase of z taken from the run; markers on the right are '
            'D inverted from sigma(w) with no fit. Free-standing monolayer, no substrate. '
            'Only the field DEPENDENCE is compared -- the D/D_eq column above says what each '
            'run does to the absolute weight.')
    fig.text(0.008, 0.008, note, fontsize=7.0, color='#34495e', wrap=True)
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.savefig(args.out, dpi=150)
    print('# wrote ' + args.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
