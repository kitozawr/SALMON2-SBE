#!/usr/bin/env python3
"""THz permittivity of a bulk SBE run: what the record can say, and what it cannot.

Reads a *_sbe_rt.data / *_sbe_nex_nonad.data pair and reports the response over
0-4 THz two ways, because on a record taken WITHOUT the pure-gauge restoration the two
disagree by three orders of magnitude, and the disagreement is the point.

  MEASURED   sigma(w) = J(w)/E(w) straight off the record, zero-padded onto a fine grid.
             A uniform A is a pure gauge, so a COMPLETE basis gives zero current from a
             filled sea at every A. An nb-band basis does not: it leaves eta N_e A/V,
             with eta the oscillator strength the truncation fails to capture (the
             startup banner prints it). This is a BASIS effect, not a mesh effect -- a
             denser k-mesh does not touch it. Measure eta on the LOWEST field available:
             at high field the slope of J against A also carries the real free carriers
             (60 % of it for GaAs at 1043 kV/cm), so reading eta there overstates it. On
             this archive the clean numbers are 0.73 % (Si) and 1.03 % (GaAs), both at
             104 kV/cm, where the carriers contribute nothing to the slope. What is left is reversible, tracking A(t)
             rather than E(t), and at THz it is ruinous: the residue scales as A = E/w
             while the polarisation current it has to be compared with scales as w E, so
             the contamination grows as 1/w^2. eta of a per cent is invisible in the
             optical and fatal at 2 THz. `yn_sbe_vg_sumrule = 'y'` removes it by
             construction (wiki/12 SS6a), subtracting the adiabatic ground-state current
             of the SAME truncated H_k(A) -- parameter-free, exact at every A, one ZHEEV
             per k per current evaluation, and material-agnostic (it reads only eigen,
             p_tm and the reference occupation, so Si and GaAs are on the same footing
             as graphene -- cleaner, in fact, since their sea is undoped).

  MODELLED   Drude from the REAL carrier density in *_sbe_nex_nonad.data (the
             dressed-reference column, which is what the ring dissipators see), with m*
             and eps_inf from the literature. A model, not a measurement -- but built on
             the one number in the archive the residue does not touch, and it answers
             what the measurement cannot: where the plasma edge sits relative to 0-4 THz.

The 0-4 THz grid is finer than the records resolve (1/T = 0.33-0.69 THz for the archive
this was written for), so the padded curve interpolates and adds no information; the
Drude curve is analytic and can be evaluated anywhere. Both are drawn, so the reader can
see which is which.
"""
import argparse
import glob
import os
import sys

import numpy as np

E_SI, EPS0, C_SI, ME = 1.602176634e-19, 8.8541878128e-12, 2.99792458e8, 9.1093837015e-31
J_CONV = 1.602176634e16          # Jm [1/(fs Angstrom^2)] -> A/m^2

# conductivity effective mass and the THz background permittivity
MATERIAL = {'Si':   dict(mstar=0.26, eps_inf=11.68, cell_A3=5.431**3 / 4),
            'GaAs': dict(mstar=0.067, eps_inf=12.9, cell_A3=5.653**3 / 4)}


def load(run):
    rt = np.loadtxt(glob.glob(os.path.join(run, '*_sbe_rt.data'))[0])
    t = rt[:, 0] * 1e-15
    ax = int(np.argmax([np.abs(rt[:, 4 + i]).max() for i in range(3)]))
    return dict(t=t, dt=t[1] - t[0],
                E=rt[:, 4 + ax] * 1e10,            # V/m
                A=rt[:, 7 + ax],                   # shape only, units irrelevant here
                J=-rt[:, 13 + ax] * J_CONV)        # A/m^2, electric current


def real_density(run):
    """Non-adiabatic carrier density [m^-3]: (nex_dref, nex_proj)."""
    f = glob.glob(os.path.join(run, '*_sbe_nex_nonad.data'))
    if not f:
        return np.nan, np.nan
    a = np.loadtxt(f[0])
    return a[-1, 2] * 1e6, a[-1, 1] * 1e6


def measured_eps(d, eps_inf, fmax=4.0, pad=16):
    """eps(w) from the record. numpy's rfft is the exp(+iwt) convention, so conjugate to
    the physics one before using eps = eps_inf + i sigma / (eps0 w)."""
    n = len(d['t'])
    w = np.fft.rfftfreq(pad * n, d=d['dt'])
    Ew = np.conj(np.fft.rfft(d['E'], pad * n))
    Jw = np.conj(np.fft.rfft(d['J'], pad * n))
    m = (w > 0) & (w <= fmax * 1e12)
    sig = np.where(np.abs(Ew[m]) > 0, Jw[m] / np.where(np.abs(Ew[m]) == 0, 1, Ew[m]), 0)
    return w[m] / 1e12, sig, eps_inf + 1j * sig / (EPS0 * 2 * np.pi * w[m])


def drude_eps(f_thz, n_m3, mstar, eps_inf, tau_fs):
    w = 2 * np.pi * np.asarray(f_thz) * 1e12
    wp2 = n_m3 * E_SI**2 / (EPS0 * mstar * ME)
    return eps_inf - wp2 / (w**2 + 1j * w / (tau_fs * 1e-15)), np.sqrt(wp2) / (2 * np.pi) / 1e12


def absorbed_energy(run, cell_A3):
    """Absorbed energy per cell [meV], integrated to the last near-zero of A(t).

    int J.E dt is the one quantity a record dominated by the reversible residue can
    still deliver -- but only if the window is chosen right. The residue J ~ -cA does
    work c[A^2/2] taken between the ends, so it cancels EXACTLY when A starts and
    finishes at zero, and not at all otherwise. These records stop mid-pulse with A at
    20-40 % of its peak, where the leftover is 148, 94 and 2.8 % of the raw integral for
    Si at 104, 1086 and 3258 kV/cm and 0.5, 5.9 and 4.3 % for GaAs at 104, 1043 and 3129.
    Cutting at the last sample where |A| is smallest -- on these records the last zero
    crossing of A -- drops it below 1e-6 of the integral, out of sight. Projecting the residue
    out instead over-subtracts (it takes some of the real dissipative response with it)
    and returns negative absorption, so the window, not the projection, is the repair.
    """
    a = np.loadtxt(glob.glob(os.path.join(run, '*_sbe_rt.data'))[0])
    t = a[:, 0] * 1e-15
    ax = int(np.argmax([np.abs(a[:, 4 + i]).max() for i in range(3)]))
    # A in SI (V s / m): the column is fs*V/Angstrom, so 1e-15 * 1e10. Without this the
    # residue estimate below is out by 1e5 -- it was, and it made a leftover of 1e-7 meV
    # look like 12 % of the absorbed energy.
    A = a[:, 7 + ax] * 1e-5
    Et = a[:, 10 + ax] * 1e10
    Ei = a[:, 4 + ax] * 1e10
    J = -a[:, 13 + ax] * J_CONV
    lo = int(0.85 * len(A))
    i = lo + int(np.argmin(np.abs(A[lo:])))
    conv = cell_A3 * 1e-30 / E_SI * 1e3                    # J -> meV per cell
    W = np.trapezoid(J[:i + 1] * Et[:i + 1], t[:i + 1]) * conv
    c = np.dot(J, A) / np.dot(A, A)
    residue = abs(c * (A[i]**2 - A[0]**2) / 2 * conv)
    F = EPS0 * C_SI * np.trapezoid(Ei[:i + 1]**2, t[:i + 1]) * 1e3 / 1e4   # mJ/cm^2
    alpha = (W * 1e-3 * E_SI / (cell_A3 * 1e-24)) / (F * 1e-3)             # 1/cm
    return dict(W_meV=W, residue_frac=residue / max(abs(W), 1e-30), t_cut_fs=t[i] * 1e15,
                F_mJcm2=F, alpha_cm=alpha, A_end_rel=abs(A[i]) / np.abs(A).max())


def optics(eps, f_thz):
    nt = np.sqrt(np.asarray(eps, dtype=complex))
    nt = np.where(nt.imag < 0, np.conj(nt), nt)
    w = 2 * np.pi * np.asarray(f_thz) * 1e12
    return nt, 2 * w * nt.imag / C_SI, np.abs((1 - nt) / (1 + nt))**2


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='+', help='run directories holding *_sbe_rt.data')
    ap.add_argument('--tau-fs', type=float, default=30.0,
                    help='momentum-relaxation time for the Drude model (default %(default)g fs); '
                         'a record dominated by the gauge residue cannot fix it, so the figure '
                         'brackets it with 10 and 100 fs')
    ap.add_argument('--slab-um', type=float, default=10.0,
                    help='slab thickness for the transmission panel (default %(default)g um), '
                         'single pass with both Fresnel faces, no Fabry-Perot')
    ap.add_argument('--out', default='thz_permittivity.png')
    ap.add_argument('--absorbed-out', default=None,
                    help='also write the absorbed-energy figure here: W(E_0) and the '
                         'effective absorption coefficient, which is where the two '
                         'regimes -- bleaching then avalanche -- are visible')
    args = ap.parse_args(argv)

    fgrid = np.linspace(0.05, 4.0, 400)
    rows = []
    for run in args.runs:
        name = os.path.basename(run.rstrip('/'))
        mat = 'GaAs' if 'GaAs' in name else 'Si'
        p = MATERIAL[mat]
        d = load(run)
        n_dref, n_proj = real_density(run)
        f, sig, eps_meas = measured_eps(d, p['eps_inf'])
        ab = absorbed_energy(run, p['cell_A3'])
        rows.append(dict(name=name, mat=mat, d=d, f=f, sig=sig, eps=eps_meas, ab=ab,
                         n_dref=n_dref, n_proj=n_proj,
                         rA=np.corrcoef(d['J'], d['A'])[0, 1],
                         rE=np.corrcoef(d['J'], d['E'])[0, 1],
                         Epk=np.abs(d['E']).max() / 1e5, **p))

    print(f'{"run":20} {"E_pk[kV/cm]":>11} {"r(J,A)":>7} {"r(J,E)":>7} {"n_real[cm^-3]":>14} '
          f'{"f_p[THz]":>9} {"eps1(1THz)":>10} {"eps2(1THz)":>10} {"alpha[1/cm]":>12} '
          f'{"R":>6} {"eps1 MEASURED":>14}')
    for r in rows:
        eps_d, fp = drude_eps(fgrid, r['n_dref'], r['mstar'], r['eps_inf'], args.tau_fs)
        i1 = int(np.argmin(abs(fgrid - 1.0)))
        j1 = int(np.argmin(abs(r['f'] - 1.0)))
        _, alpha, R = optics(eps_d[i1], 1.0)
        print(f'{r["name"]:20} {r["Epk"]:11.0f} {r["rA"]:+7.3f} {r["rE"]:+7.3f} '
              f'{r["n_dref"]/1e6:14.3e} {fp:9.1f} {eps_d[i1].real:10.2f} {eps_d[i1].imag:10.2f} '
              f'{alpha/100:12.3e} {R:6.3f} {r["eps"][j1].real:14.3e}')
        r['eps_drude'], r['fp'] = eps_d, fp
    print('# r(J,A) near -1 with r(J,E) near 0: the recorded current is REVERSIBLE, tracking')
    print('#   the vector potential. The last column is what that residue implies for eps1 --')
    print('#   compare it with eps_inf = 11.7 (Si) / 12.9 (GaAs).')
    print(f'# Drude columns: tau = {args.tau_fs:g} fs, m* and eps_inf from the literature,')
    print('#   n from nex_dref. alpha and R are at 1 THz, normal incidence, semi-infinite.')
    print()
    print('# ABSORBED ENERGY -- the one thing a residue-dominated record still gives, once the')
    print('# window ends where A does (see absorbed_energy).')
    print(f'{"run":20} {"E_pk[kV/cm]":>11} {"t_cut[fs]":>10} {"residue":>8} {"W[meV/cell]":>12} '
          f'{"F[mJ/cm2]":>10} {"alpha_eff[1/cm]":>16}')
    for r in rows:
        a = r['ab']
        print(f'{r["name"]:20} {r["Epk"]:11.0f} {a["t_cut_fs"]:10.1f} {a["residue_frac"]:7.1%} '
              f'{a["W_meV"]:12.4f} {a["F_mJcm2"]:10.4f} {a["alpha_cm"]:16.3e}')
    for mat in ('Si', 'GaAs'):
        rr = sorted([r for r in rows if r['mat'] == mat], key=lambda x: x['Epk'])
        for i in range(len(rr) - 1):
            e0, e1 = rr[i]['Epk'], rr[i + 1]['Epk']
            w0, w1 = rr[i]['ab']['W_meV'], rr[i + 1]['ab']['W_meV']
            a0, a1 = rr[i]['ab']['alpha_cm'], rr[i + 1]['ab']['alpha_cm']
            print(f'#   {mat}: {e0:.0f} -> {e1:.0f} kV/cm   W ~ E^{np.log(w1 / w0) / np.log(e1 / e0):.2f} '
                  f'(linear absorption = E^2),   alpha x{a1 / a0:.2f}')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# no matplotlib -- table only')
        return 0

    fig, ax = plt.subplots(1, 4, figsize=(21.0, 4.7))
    cols = {'Si': ['#c0392b', '#e67e22', '#d35400'], 'GaAs': ['#2980b9', '#16a085', '#8e44ad']}
    used = {'Si': 0, 'GaAs': 0}
    for r in rows:
        c = cols[r['mat']][used[r['mat']] % 3]; used[r['mat']] += 1
        lab = f'{r["mat"]} {r["Epk"]:.0f} kV/cm'
        ax[0].plot(r['f'], r['eps'].real, '-', lw=1.2, color=c, label=lab)
        ax[1].plot(fgrid, r['eps_drude'].real, '-', lw=1.8, color=c,
                   label=f'{lab},  $f_p$ = {r["fp"]:.0f} THz')
        ax[2].plot(fgrid, r['eps_drude'].imag, '-', lw=1.8, color=c, label=lab)
        for tau, st in ((10.0, ':'), (100.0, '--')):
            e2, _ = drude_eps(fgrid, r['n_dref'], r['mstar'], r['eps_inf'], tau)
            ax[2].plot(fgrid, e2.imag, st, lw=0.9, color=c, alpha=0.55)
        _, al, R = optics(r['eps_drude'], fgrid)
        ax[3].semilogy(fgrid, (1 - R)**2 * np.exp(-al * args.slab_um * 1e-6), '-', lw=1.8,
                       color=c, label=lab)

    ax[0].axhline(11.68, ls='--', c='#7f8c8d', lw=1)
    ax[0].annotate("$\\epsilon_\\infty$ = 11.7 (Si)", xy=(1.6, 11.68), fontsize=8, color='#7f8c8d')
    ax[0].set_yscale('symlog', linthresh=100)
    ax[0].set_title("MEASURED off the record\n$\\epsilon'$ is the gauge residue, not the material",
                    fontsize=9)
    ax[1].axhline(0, c='k', lw=0.6)
    ax[1].set_yscale('symlog', linthresh=10)
    ax[1].set_title("MODELLED from the real carrier density\n$\\epsilon'(\\omega)$, Drude", fontsize=9)
    ax[2].set_yscale('log')
    ax[2].set_title("$\\epsilon''(\\omega)$, Drude\ndotted $\\tau$ = 10 fs, dashed $\\tau$ = 100 fs",
                    fontsize=9)
    # the undamped plasma edge, where eps' would cross zero if tau were long
    for r in rows:
        fe = r['fp'] / np.sqrt(r['eps_inf'])
        if fe <= 4.0:
            ax[1].axvline(fe, ls=':', lw=1.2, color='#7f8c8d')
            ax[1].annotate(f'plasma edge {fe:.1f} THz', xy=(fe * 1.03, 0.06),
                           xycoords=('data', 'axes fraction'), fontsize=7.5, color='#7f8c8d')
    ax[3].axhline(1.0, c='k', lw=0.6)
    ax[3].set_ylim(1e-12, 2)
    ax[3].set_title(f'Transmission of a {args.slab_um:g} um slab\n(both faces, single pass)',
                    fontsize=9)
    for a, yl in zip(ax, ("$\\epsilon'$ from $J/E$", "$\\epsilon'$", "$\\epsilon''$",
                          'transmission')):
        a.set_xlabel('frequency [THz]'); a.set_ylabel(yl)
        a.set_xlim(0, 4); a.grid(alpha=0.25); a.legend(fontsize=7, framealpha=0.9)
    note = ('Left: sigma(w) = J(w)/E(w) from the records, zero-padded x16 -- the records resolve only '
            '0.33-0.69 THz, so the curve interpolates rather than resolves. Middle and right: Drude '
            'from nex_dref with m* and eps_inf from the literature. Both panels describe the SAME runs; '
            'they disagree because the records were taken without yn_sbe_vg_sumrule.')
    fig.text(0.008, 0.008, note, fontsize=7, color='#34495e', wrap=True)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(args.out, dpi=150)
    print(f'# wrote {args.out}')

    if args.absorbed_out:
        f2, ax2 = plt.subplots(1, 2, figsize=(11.6, 4.6))
        style = {'Si': ('#c0392b', 'o'), 'GaAs': ('#2980b9', 's')}
        for mat in ('Si', 'GaAs'):
            rr = sorted([r for r in rows if r['mat'] == mat], key=lambda x: x['Epk'])
            if not rr:
                continue
            c, mk = style[mat]
            E = np.array([r['Epk'] for r in rr])
            W = np.array([r['ab']['W_meV'] for r in rr])
            al = np.array([r['ab']['alpha_cm'] for r in rr])
            ax2[0].loglog(E, W, mk + '-', ms=8, color=c, label=mat)
            ax2[1].loglog(E, al, mk + '-', ms=8, color=c, label=mat)
            for i in range(len(E) - 1):
                pwr = np.log(W[i + 1] / W[i]) / np.log(E[i + 1] / E[i])
                ax2[0].annotate(f'$E^{{{pwr:.1f}}}$',
                                xy=(np.sqrt(E[i] * E[i + 1]), np.sqrt(W[i] * W[i + 1])),
                                fontsize=9, color=c, ha='center', va='bottom')
                ax2[1].annotate(f'$\\times{al[i + 1] / al[i]:.2f}$',
                                xy=(np.sqrt(E[i] * E[i + 1]), np.sqrt(al[i] * al[i + 1])),
                                fontsize=9, color=c, ha='center', va='bottom')
        # what linear absorption would look like: W proportional to the fluence, E^2
        rr = sorted([r for r in rows if r['mat'] == 'Si'], key=lambda x: x['Epk'])
        if rr:
            E0, W0 = rr[0]['Epk'], rr[0]['ab']['W_meV']
            eg = np.geomspace(80, 4000, 50)
            ax2[0].loglog(eg, W0 * (eg / E0)**2, '--', lw=1.2, color='#7f8c8d',
                          label='linear absorption, $W\\propto E_0^2$')
        ax2[0].set_ylabel('absorbed energy [meV per cell]')
        ax2[0].set_title('Absorbed energy against peak field\n'
                         r'$\int J\cdot E\,dt$ to the last zero of $A(t)$', fontsize=9)
        ax2[1].set_ylabel(r'effective $\alpha$ = W / fluence [cm$^{-1}$]')
        ax2[1].set_title('The same, divided by the fluence:\nbleaching, then avalanche',
                         fontsize=9)
        for a in ax2:
            a.set_xlabel('peak field $E_0$ [kV/cm]')
            a.grid(alpha=0.25, which='both'); a.legend(fontsize=8)
        note = ('Window ends at the last near-zero of A(t): the reversible velocity-gauge residue '
                'does work c[A^2/2] taken between the ends, so it cancels exactly there. Residue left: '
                + ', '.join(f'{r["name"].split("_")[0]} {r["Epk"]:.0f} kV/cm {r["ab"]["residue_frac"]:.1%}'
                            for r in rows) + '.')
        f2.text(0.008, 0.008, note, fontsize=6.8, color='#34495e', wrap=True)
        f2.tight_layout(rect=(0, 0.06, 1, 1)); f2.savefig(args.absorbed_out, dpi=150)
        print(f'# wrote {args.absorbed_out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
