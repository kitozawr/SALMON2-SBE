!
!  Nonlocal "super-compute" mode (Part C) -- rate and final-state-search
!  PRIMITIVES. This module contains only PURE functions and cited data tables;
!  it has no SALMON dependencies and no state, so it can be unit-tested on its
!  own (tests/test_superres_rates.f90) and reused by the dynamics integration in
!  later increments. Nothing here is called from the time loop yet -- the whole
!  nonlocal mode is gated behind yn_sbe_superres='n' (default).
!
!  Contents:
!    - nu_saturation        : collision-rate saturation nu(eps) (smooth, no min())
!    - bose_factor          : phonon occupation N_B(hw, T)
!    - gaussian_bin/rect_bin: normalized broadened delta(Delta E) (energy bins, C3)
!    - frohlich_hi_factor   : GaAs polar-optical high-energy asinh factor
!    - ii_rate_general      : general impact-ionization fit P (eps-E_th)^a
!    - bgr_gap_shift_ev     : bandgap-renormalization cube-root law (C7)
!    - Si / GaAs intervalley deformation-potential data tables (C5)
!
!  All energies in a single consistent unit chosen by the caller (the dynamics
!  uses Hartree atomic units); the data tables carry their own stated units.
!
module sbe_superres_ssbe
    implicit none
    private

    real(8), parameter :: PI = 3.14159265358979323846d0
    ! Unit conversions to Hartree atomic units (a.u.). [CODATA 2018]
    real(8), parameter :: HA_EV    = 27.211386245988d0   ! 1 Hartree in eV
    real(8), parameter :: BOHR_CM  = 0.529177210903d-8   ! 1 Bohr in cm
    real(8), parameter :: BOHR_ANG = 0.529177210903d0    ! 1 Bohr in Angstrom
    real(8), parameter :: ME_G     = 9.1093837015d-28    ! electron mass in g

    public :: nu_saturation, bose_factor, gaussian_bin, rect_bin, &
              frohlich_hi_factor, ii_rate_general, bgr_gap_shift_ev, &
              gaussian_shape, amp_damp_channel, &
              mev_to_ha, d_evcm_to_au, d_evang_to_au, rho_gcm3_to_au, &
              golden_rule_prefactor, eph_thermal_split, eph_thermal_split_de, &
              eps_thomas_fermi, tf_kappa2_degenerate, debye_kappa2, &
              lindhard_F, eps_lindhard_static, plasmon_freq2, lopc_branches, &
              eps_cdrb, interk_vq, build_vq_table, build_acscreen_table, t_ring_opts, &
              dirac_n_2d, dirac_mu_2d, dirac_e_2d, dirac_fit_te, rana_qtf, dirac_plasmon_2d, &
              rana_rcccv, rana_auger_dpop, &
              energy_partner_weights, fermi_dirac, fit_fermi_dirac, &
              carrier_carrier_relax, eph_interk_dpop, ii_interk_dpop, &
              auger_interk_dpop, mp_grid_triple, mp_partner_triple, &
              vg_eta_admixture, vg_trunc_shift2, vg_conv_error, vg_ptop_exceeds, &
              bath_t2_high_t, bath_corr_table, sfsb_nc_series, &
              colmem_lines, colmem_response, colmem_pop_filter, colmem_pop_init, &
              dressed_ref_delta, vg_sumrule_eta, vg_gs_current_k, &
              get_material_params

    ! =====================================================================
    ! Silicon intervalley deformation potentials -- Pop "new" set (default).
    ! Order: g-TA, g-LA, g-LO, f-TA, f-LA, f-TO. E [meV], D [1e8 eV/cm].
    ! g-type couple same-<100>-axis valleys; f-type orthogonal-axis valleys.
    ! [Jacoboni & Reggiani, Rev. Mod. Phys. 55, 645 (1983); Pop set compiled in
    !  Jacoboni-Lugli; anchor Canali et al., Phys. Rev. B 15, 3994 (1977)]
    ! =====================================================================
    integer, parameter, public :: SI_N_PHONON = 6
    real(8), parameter, public :: SI_PHONON_E_MEV(SI_N_PHONON) = &
        (/ 10.0d0, 19.0d0, 63.0d0, 19.0d0, 51.0d0, 57.0d0 /)
    real(8), parameter, public :: SI_PHONON_D_1E8EVCM(SI_N_PHONON) = &
        (/  0.3d0,  1.5d0,  6.0d0,  0.5d0,  3.5d0,  1.5d0 /)
    ! Silicon acoustic/transport constants [Jacoboni-Reggiani 1983; Canali 1977]
    real(8), parameter, public :: SI_XI_D_EV     = 9.0d0      ! acoustic def. potential [eV]
    real(8), parameter, public :: SI_RHO_GCM3    = 2.33d0     ! mass density [g/cm^3]
    real(8), parameter, public :: SI_VLA_CMS     = 9.01d5     ! v_LA [cm/s]
    real(8), parameter, public :: SI_VTA_CMS     = 5.23d5     ! v_TA [cm/s]
    real(8), parameter, public :: SI_ALPHA_NP    = 0.5d0      ! non-parabolicity [1/eV]

    ! =====================================================================
    ! GaAs intervalley deformation potentials. Order: Gamma-L, Gamma-X, L-L,
    ! L-X, X-X. E [meV], D [eV/Angstrom]. Plus Frohlich polar-optical constants.
    ! [Fischetti & Laux, Phys. Rev. B 38, 9721 (1988); IEEE TED 38, 634 (1991);
    !  Frohlich: Fawcett-Boardman-Swain, J. Phys. Chem. Solids 31, 1963 (1970)]
    ! =====================================================================
    integer, parameter, public :: GAAS_N_IV = 5
    real(8), parameter, public :: GAAS_IV_E_MEV(GAAS_N_IV) = &
        (/ 27.8d0, 29.9d0, 29.0d0, 29.3d0, 29.9d0 /)
    real(8), parameter, public :: GAAS_IV_D_EVANG(GAAS_N_IV) = &
        (/ 10.0d0, 10.0d0, 10.0d0,  5.0d0,  7.0d0 /)
    real(8), parameter, public :: GAAS_HW_LO_MEV = 36.0d0    ! Frohlich LO phonon [meV]
    real(8), parameter, public :: GAAS_HW_TO_MEV = 33.6d0    ! TO phonon [meV] (LOPC)
    real(8), parameter, public :: GAAS_ALPHA_FR  = 0.068d0   ! Frohlich coupling
    real(8), parameter, public :: GAAS_EPS0      = 12.9d0    ! static dielectric
    real(8), parameter, public :: GAAS_EPS_INF   = 10.89d0   ! high-freq dielectric
    real(8), parameter, public :: GAAS_M_GAMMA   = 0.067d0   ! Gamma effective mass [m_e]
    real(8), parameter, public :: GAAS_M_HH      = 0.51d0    ! heavy-hole mass [m_e]
    real(8), parameter, public :: GAAS_M_LH      = 0.082d0   ! light-hole mass [m_e]

    ! =====================================================================
    ! Silicon screening constants (non-polar: no LO-phonon-plasmon coupling).
    ! [eps: std Si; masses: Si effective-mass tables (m_l/m_t, hh/lh)]
    ! =====================================================================
    real(8), parameter, public :: SI_EPS    = 11.7d0    ! static dielectric (non-polar)
    real(8), parameter, public :: SI_M_L    = 0.98d0    ! longitudinal mass [m_e]
    real(8), parameter, public :: SI_M_T    = 0.19d0    ! transverse mass [m_e]
    real(8), parameter, public :: SI_M_HH   = 0.49d0    ! heavy-hole mass [m_e]
    real(8), parameter, public :: SI_M_LH   = 0.16d0    ! light-hole mass [m_e]

    ! =====================================================================
    ! Wurtzite CdS (P6_3mc) -- polar, NON-centrosymmetric II-VI.
    ! Structure: Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967).
    ! DISSIPATION CHANNELS (cited in the CdS physics-methods spec / md):
    !   * Frohlich polar-optical e-ph -- the PRIMARY room-T channel
    !     (hw_LO = 38 meV [Raman]; alpha ~ 0.5 [cyclotron, arXiv:cond-mat/
    !     0107481; path-integral arXiv:2205.11780]). nu_sat is DERIVED from the
    !     cited coupling: the Frohlich rate scale ~ alpha * omega_LO.
    !   * Coulomb HF (static dielectric eps0 ~ 9.0 isotropic avg, eps_inf ~ 5.3
    !     [md; anisotropic ||/_|_ c]).
    !   * Impact ionization: Keldysh soft threshold E_th ~ 1.5 E_g = 3.6 eV
    !     [md, (3/2)E_g rule]; the PREFACTOR is a fit parameter (CdS-specific
    !     value scarce -- md), so the user must set sbe_ii_prefactor explicitly.
    ! Carrier-carrier (e-e/e-h) and Auger recombination are FORBIDDEN for CdS:
    ! the CdS literature only fixes an e-e *timescale* (sub-100fs thermalization
    ! at n >= 1e18 cm^-3 [Shah et al., IEEE JQE 22, 1728 (1986); Elsaesser PRL 66,
    ! 1757 (1991)]), NOT a carrier-carrier *rate*. There is NO verified CdS Auger
    ! coefficient: the previously-cited "C = 2.0e-30 cm^6/s [Haury, PRB 57, 11513
    ! (1998)]" was a FABRICATED reference -- the real Haury et al. paper is PRL 79,
    ! 511 (1997) on ferromagnetism in CdMnTe quantum wells, unrelated to Auger in
    ! CdS, and the coefficient is unconfirmed anywhere. It has been REMOVED; CdS
    ! Auger is gated off (auger_ok=.false.). A user with a verified rate can still
    ! opt in explicitly via sbe_auger_c_cm6s. (Physically CdS is wide-gap, Eg~2.4
    ! eV, so direct Auger is exponentially suppressed and any future channel would
    ! be the phonon/defect-assisted one -- see wiki/07, the nonlocal-Auger task.)
    ! Piezoelectric acoustic (e33/e31/e15 [Berlincourt PR 129,1009]) and
    ! deformation-potential acoustic are cited but NOT yet SBE Lindblad channels.
    ! =====================================================================
    real(8), parameter, public :: CDS_A_BOHR    = 7.8159d0  ! a = 4.136 Ang  [BC1967 Table I]
    real(8), parameter, public :: CDS_ASQ3_BOHR = 13.5376d0 ! b = a*sqrt(3) (orthorhombic)
    real(8), parameter, public :: CDS_C_BOHR    = 12.6877d0 ! c = 6.714 Ang (c/a=1.623) [BC1967]
    real(8), parameter, public :: CDS_U_INT     = 0.375d0   ! internal parameter u=3/8 [BC1967]
    real(8), parameter, public :: CDS_EG_EV     = 2.58d0    ! direct gap at Gamma (low-T) [BC1967 Table I]
    real(8), parameter, public :: CDS_EPS0      = 8.9d0     ! static dielectric (isotropic avg) [Berlincourt 1963]
    real(8), parameter, public :: CDS_EPS_INF   = 5.3d0     ! high-frequency dielectric [Berlincourt 1963]
    ! (NO verified CdS Auger coefficient exists -- the former CDS_AUGER_C with its
    !  "Haury PRB 57, 11513 (1998)" citation was fabricated and has been removed.)
    real(8), parameter, public :: CDS_EE_ACT_N  = 1.0d18    ! e-e activation density [cm^-3] [Shah JQE 22, 1728]
    real(8), parameter, public :: CDS_HW_LO_MEV = 38.0d0    ! Frohlich LO ~305 cm^-1 [Raman; md]
    real(8), parameter, public :: CDS_ALPHA_FR  = 0.5d0     ! Frohlich coupling alpha [md]
    ! Frohlich e-ph rate scale nu_sat = alpha * omega_LO (omega_LO = hw_LO/hbar):
    ! 0.5 * (38e-3 eV / 6.582e-16 eV.s) = 2.89e13 s^-1 (~35 fs; md: sub-100 fs).
    real(8), parameter, public :: CDS_NU_SAT_SI = 2.89d13   ! derived from cited alpha, hw_LO [md]
    real(8), parameter, public :: CDS_II_ETH_EV = 3.6d0     ! II threshold ~1.5 E_g [md, (3/2)E_g]
    real(8), parameter, public :: CDS_RHO_GCM3  = 4.82d0    ! mass density [g/cm^3]
    real(8), parameter, public :: CDS_E33       = 0.385d0   ! piezo [C/m^2] [Berlincourt PR 129,1009]
    real(8), parameter, public :: CDS_E31       = -0.262d0  ! piezo [C/m^2] [Berlincourt 1963]
    real(8), parameter, public :: CDS_E15       = -0.183d0  ! piezo [C/m^2] [Berlincourt 1963]

    ! =====================================================================
    ! Monolayer graphene -- GAPLESS Dirac semimetal, honeycomb (D6h), V^A=0
    ! (centrosymmetric). a = 2.46 Ang = 4.6511 Bohr.
    ! ELECTRON-PHONON (the cited, enabled channel): the two strongly-coupled
    ! optical modes responsible for the Kohn anomalies and dominant hot-carrier
    ! energy relaxation:
    !   * E2g at Gamma (intra-valley LO/TO, q~0): hw = 196 meV (~1580 cm^-1),
    !     EPC <g^2_Gamma> = 0.0405 eV^2.
    !   * A1' at K   (INTER-valley K<->K', the Kohn-anomaly mode): hw = 160 meV,
    !     EPC <g^2_K> = 0.0994 eV^2 (DFT); enhanced ~x2 by GW.
    ! [Piscanec, Lazzeri, Mauri, Ferrari, Robertson, PRL 93, 185503 (2004);
    !  GW enhancement: Lazzeri, Attaccalite, Wirtz, Mauri, PRB 78, 081406 (2008)].
    ! Both modes RELAX carriers DOWN THE CONE / across valleys, i.e. to a DIFFERENT
    ! k -- so graphene e-ph is physical only through the INTER-K ring path
    ! (yn_sbe_superres='y'); the k-local same-k search is inert on the 2-band cone.
    ! nu_sat = optical-phonon emission rate cap; OP emission is the dominant hot-
    ! carrier cooling channel, time scale ~10-100 fs [hot-carrier cooling lit.].
    ! Conservative cap ~5e13 s^-1 (~20 fs); tunable via sbe_eph_nu_sat.
    real(8), parameter, public :: GRAPH_A_BOHR    = 4.6511d0  ! a = 2.46 Ang
    real(8), parameter, public :: GRAPH_HW_E2G_MEV = 196.0d0  ! Gamma E2g optical [Piscanec 2004]
    real(8), parameter, public :: GRAPH_HW_A1P_MEV = 160.0d0  ! K A1' optical (Kohn anomaly) [Piscanec 2004]
    real(8), parameter, public :: GRAPH_G2_E2G    = 0.0405d0  ! <g^2_Gamma> [eV^2] [Piscanec 2004]
    real(8), parameter, public :: GRAPH_G2_A1P    = 0.0994d0  ! <g^2_K> [eV^2] DFT [Piscanec 2004]
    real(8), parameter, public :: GRAPH_GW_K      = 2.0d0     ! GW enhancement of the K coupling [Lazzeri 2008]
    real(8), parameter, public :: GRAPH_NU_SAT_SI = 5.0d13    ! OP-emission rate cap [~20 fs; hot-carrier cooling]
    ! 2D Rana Auger / carrier-multiplication constants [R07 = F. Rana,
    ! PRB 76, 155431 (2007); journal text verified, wiki/07 sec.6]:
    real(8), parameter, public :: GRAPH_VF_AU     = 1.0d8 / 2.18769126364d8
    !                              ^ Dirac velocity v_F = 1e8 cm/s [R07] in a.u.
    real(8), parameter, public :: GRAPH_RANA_EPS  = 10.0d0
    !                              ^ background eps_r of the R07 Fig.4 lifetime
    !                                benchmarks (Al2O3-like); SUBSTRATE-dependent
    !                                (R07 Fig.5: eps=4 SiO2 doubles the rate) --
    !                                override via sbe_coulomb_epsilon.

    ! =====================================================================
    ! Material registry -- the SINGLE place that maps a material name to all
    ! the per-material constants the SBE dissipation channels need (dielectric,
    ! impact-ionization fit, electron-phonon table). Adding a material is one
    ! `case` block in get_material_params() plus its name in MAT_SUPPORTED;
    ! every channel then auto-selects through the same struct. All numbers are
    ! the cited constants declared above -- the registry only assembles them.
    ! =====================================================================
    integer, parameter, public :: MAT_MAXPH = 8         ! capacity of a phonon table
    character(*), parameter, public :: MAT_SUPPORTED = 'GaAs, Si, Si_cb, CdS, graphene'

    type, public :: s_material_params
        logical       :: found        = .false.
        character(16) :: name         = ''
        ! crystal structure
        real(8)       :: a_lattice_au = 0d0       ! in-plane lattice constant a [Bohr]
        real(8)       :: cell_au(3)   = 0d0       ! expected &system al(1:3) box [Bohr]
        logical       :: is_diamond   = .false.   ! diamond (V^A=0) vs zincblende/wurtzite
        ! -- PROVENANCE GATES --------------------------------------------------
        ! A channel may be enabled for a material ONLY if its constants are
        ! backed by a cited source for THAT material. No source => the flag is
        ! .false. => the SBE init aborts (error stop) if the channel is enabled.
        ! Never fall back to another material's constants.
        logical       :: ii_ok        = .false.   ! impact-ionization fit cited?
        logical       :: eph_ok       = .false.   ! e-ph rate (nu_sat) cited?
        logical       :: eeh_ok       = .false.   ! carrier-carrier rate cited?
        logical       :: coulomb_ok   = .false.   ! dielectric for Coulomb cited?
        ! A7: Sigma^HF uses the 2D sheet kernel V(q) = 2 pi/(eps A (q+kappa))
        ! instead of the 3D 4 pi/(eps Omega q^2) (graphene: 2D slab geometry).
        logical       :: coulomb_2d   = .false.
        logical       :: auger_ok     = .false.   ! Auger coeff C cited?
        ! Auger recombination (Sec 13): R = C n^3, density-gated above n_gate.
        real(8)       :: auger_c_cm6s    = 0d0    ! Auger coefficient C [cm^6/s]
        real(8)       :: auger_n_gate_cm3 = 0d0   ! activation density [cm^-3]
        ! 2D Rana Auger/CM branch (graphene [R07]): auger_ok with this flag set
        ! means the CITED channel is the gapless 2D collinear one -- NOT C*n^3
        ! and NOT the gap-threshold ring kernel. Ring-gated (needs the gather).
        logical       :: auger_2d_rana = .false.
        real(8)       :: rana_vf_au    = 0d0      ! Dirac velocity [a.u.]
        real(8)       :: rana_eps_r    = 0d0      ! background eps_r [R07 Fig.4]
        ! dynamic free-carrier screening lambda^2(n(t)) in the ring II/Auger
        ! |V(q)|^2 (Debye/TF from the Part-G primitives, evaluated on the
        ! gathered carrier density). Si stays .false. -- lambda=0 is CORRECT by
        ! Burt's dynamical argument [L90]: the ~1 eV Auger transition frequency
        ! far exceeds the carrier plasma frequency, so the static free-carrier
        ! screen does not act. GaAs (polar, LOPC-prone) takes the density-
        ! dependent lambda as the cited Part-G refinement.
        logical       :: dyn_lambda_ok = .false.
        ! A2: hole-initiated II strength relative to the electron channel,
        ! Cp/Cn from the source-verified wiki/07 tables (0 = channel off).
        real(8)       :: ii_cpcn = 0d0
        ! A4: quasi-elastic acoustic deformation-potential mode constants
        ! (0 = not cited for this material -> acoustic mode unavailable).
        real(8)       :: eph_ac_xi_ev   = 0d0   ! deformation potential Xi_d [eV]
        real(8)       :: eph_ac_cs_cmps = 0d0   ! LA sound velocity [cm/s]
        ! dielectric (Coulomb HF exchange / carrier screening)
        real(8)       :: eps0         = 1d0       ! static dielectric
        real(8)       :: eps_inf      = 1d0       ! high-frequency dielectric
        ! impact-ionization fit defaults
        character(20) :: ii_form      = 'stobbe_quartic'
        real(8)       :: ii_exponent  = 4d0       ! exponent a in P*(eps-Eth)^a
        real(8)       :: ii_prefactor = 2d12      ! P [s^-1 eV^-a]
        real(8)       :: ii_threshold_ev = 2.1d0  ! Eth above the CBM [eV]
        ! electron-phonon population-relaxing Lindblad table
        real(8)       :: eph_nu_sat_si = 1d14     ! saturation rate [s^-1]
        logical       :: eph_polar    = .false.   ! has a Frohlich polar-LO branch
        integer       :: eph_nph      = 0
        real(8)       :: eph_hw_mev(MAT_MAXPH) = 0d0   ! phonon energies [meV]
        real(8)       :: eph_wraw(MAT_MAXPH)   = 0d0   ! raw (un-normalized) D^2/hw weights
    end type s_material_params

    ! Optional refinements bundle for the ring II/Auger kernels (all fields
    ! default to the inert value -> omitting the argument is bit-identical).
    type :: t_ring_opts   ! exported via the module `public ::` list above
        logical :: use_tab = .false.       ! B1: use the precomputed vq table
        real(8), allocatable :: vq_tab(:)  !     signed-difference table
        real(8) :: vq_floor = 0d0          ! B3: skip quadruples with vq < floor (absolute)
        real(8) :: fk_theta = 0d0          ! A5: FK electro-optic width [Ha]; 0 = hard threshold
        real(8) :: phassist = 0d0          ! A1: phonon-assisted sideband strength; 0 = off
        integer :: nph = 0                 !     phonon table for the sidebands
        real(8), allocatable :: hw(:), nbb(:), wrel(:)
        real(8) :: pref_h = 0d0            ! A2: hole-channel prefactor (pref*Cp/Cn); 0 = off
        real(8) :: evbm = 0d0              !     valence-band maximum [Ha]
    end type t_ring_opts

contains

    ! Look up all per-material SBE constants by name. mp%found = .false. for an
    ! unknown material (callers that need a registry value must check it and stop
    ! with a helpful message). This is the ONLY function to extend when adding a
    ! material -- add a `case` and update MAT_SUPPORTED above.
    pure function get_material_params(name) result(mp)
        character(*), intent(in) :: name
        type(s_material_params)  :: mp
        integer :: p
        mp%name = name
        select case (trim(name))
        case ('GaAs')
            mp%found = .true.
            mp%a_lattice_au = 10.68d0
            mp%cell_au = (/ 10.68d0, 10.68d0, 10.68d0 /)   ! cubic
            mp%is_diamond   = .false.
            mp%ii_ok = .true.; mp%eph_ok = .true.; mp%eeh_ok = .true.; mp%coulomb_ok = .true.
            mp%eps0 = GAAS_EPS0;  mp%eps_inf = GAAS_EPS_INF
            mp%ii_form = 'stobbe_quartic'; mp%ii_exponent = 4d0
            mp%ii_prefactor = 2d12;        mp%ii_threshold_ev = 2.1d0
            ! dynamic free-carrier lambda^2(n(t)) in the ring II/Auger screen:
            ! GaAs is the polar material where the free-carrier (Debye/TF)
            ! screen acts on the collision kernel (Part-G; static Lindhard is
            ! the default screen class for GaAs). Si stays lambda=0 [L90/Burt].
            mp%dyn_lambda_ok = .true.
            ! hole-initiated II/Auger: Cp/Cn = (2.2+3.1)/1.1 ~ 4.8 -- S14's
            ! "hhe ~ 5x eeh" at Eg=1.43 eV [S14, source-verified wiki/07 sec.7]
            mp%ii_cpcn = 4.8d0
            ! acoustic deformation potential Xi_d = 7.0 eV [Fischetti-Laux,
            ! PRB 38, 9721 (1988), maintainer-supplied]; LA sound velocity
            ! 5.24e5 cm/s [same source tables -- verify against the PDF]
            mp%eph_ac_xi_ev = 7.0d0;  mp%eph_ac_cs_cmps = 5.24d5
            mp%eph_nu_sat_si = 1.0d14      ! [Fischetti IEEE TED 38, 634 (1991)]
            mp%eph_polar = .true.
            ! Frohlich polar-LO (mode 1) + 5 intervalley modes
            mp%eph_nph = GAAS_N_IV + 1
            mp%eph_hw_mev(1) = GAAS_HW_LO_MEV
            do p = 1, GAAS_N_IV
                mp%eph_hw_mev(p + 1) = GAAS_IV_E_MEV(p)
                mp%eph_wraw(p + 1)   = GAAS_IV_D_EVANG(p)**2 / GAAS_IV_E_MEV(p)
            end do
            mp%eph_wraw(1) = sum(mp%eph_wraw(2:GAAS_N_IV + 1))   ! polar LO ~ dominant
        case ('Si', 'Si_cb')
            mp%found = .true.
            mp%a_lattice_au = 10.26d0
            mp%cell_au = (/ 10.26d0, 10.26d0, 10.26d0 /)   ! cubic
            mp%is_diamond   = .true.
            mp%ii_ok = .true.; mp%eph_ok = .true.; mp%eeh_ok = .true.; mp%coulomb_ok = .true.
            mp%eps0 = SI_EPS;  mp%eps_inf = SI_EPS     ! non-polar: eps_inf = eps0
            mp%ii_form = 'keldysh_quadratic'; mp%ii_exponent = 2d0
            mp%ii_prefactor = 2d12;           mp%ii_threshold_ev = 1.1d0
            ! hole-initiated II/Auger: Cp/Cn = 0.99/2.8 = 0.354
            ! [Dziewior-Schmid via L90, source-verified wiki/07 sec.7]
            mp%ii_cpcn = 0.99d0 / 2.8d0
            ! acoustic: Xi_d = 9.0 eV, c_LA = 9.04e5 cm/s [Jacoboni-Reggiani,
            ! Rev. Mod. Phys. 55, 645 (1983) -- the already-cited SI_XI_D_EV]
            mp%eph_ac_xi_ev = SI_XI_D_EV;  mp%eph_ac_cs_cmps = 9.04d5
            mp%eph_nu_sat_si = 1.3d14      ! [Meng PRB 91, 075201 (2015)]
            mp%eph_polar = .false.
            mp%eph_nph = SI_N_PHONON       ! 6 intervalley g/f modes
            do p = 1, SI_N_PHONON
                mp%eph_hw_mev(p) = SI_PHONON_E_MEV(p)
                mp%eph_wraw(p)   = SI_PHONON_D_1E8EVCM(p)**2 / SI_PHONON_E_MEV(p)
            end do
        case ('CdS')
            ! Wurtzite (P6_3mc), polar, NON-centrosymmetric. Orthorhombic cell
            ! (a, a*sqrt3, c) via &system al(1:3) -- NOT cubic.
            ! [structure: Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967)]
            !
            ! CITED, ENABLED channels (the constants below are CdS-specific):
            !   * e-ph   : Frohlich polar-optical LO (PRIMARY room-T channel);
            !              hw_LO=38 meV [Raman], nu_sat=alpha*omega_LO=2.89e13
            !              from alpha=0.5 [cyclotron].
            !   * Coulomb: static dielectric eps0=8.9 [Berlincourt PR 129,1009 (1963)].
            !   * impact : Keldysh soft threshold, E_th=1.5*E_g=3.6 eV; the
            !              PREFACTOR is an uncited fit parameter -> sentinel (-1),
            !              so the run aborts unless the user sets sbe_ii_prefactor.
            !
            ! FORBIDDEN channel (no cited CdS *rate* -> eeh_ok = .false.):
            !   * carrier-carrier (yn_sbe_eeh). The implemented channel is an
            !     FD-thermalization parametrised by a single rate scale nu_cc.
            !     CdS has NO cited carrier-carrier rate, so enabling it with the
            !     default would borrow the generic 1e14 s^-1 scale that is cited
            !     ONLY for GaAs/Si (Goodnick-Lugli PRB 37, 2578; Fischetti-Laux
            !     PRB 38, 9721) -- forbidden by the strict provenance rule. The
            !     CdS citations (sub-100fs thermalization @ n>=1e18 [Shah 1986;
            !     Elsaesser 1991]) describe a *timescale*, not a rate. A user who
            !     supplies their own rate can still opt in via sbe_eeh_nu_sat (the
            !     same explicit-input escape hatch as the II prefactor).
            !
            ! FORBIDDEN channel (NO verified CdS Auger coefficient -> auger_ok=.false.):
            !   * Auger recombination (yn_sbe_auger). The former default
            !     "C = 2.0e-30 cm^6/s [Haury, PRB 57, 11513 (1998)]" was a
            !     FABRICATED reference (the real Haury 1997 PRL is on CdMnTe
            !     ferromagnetism, not CdS Auger) and was removed. CdS Auger is
            !     gated off; a user with a verified C may still opt in via
            !     sbe_auger_c_cm6s.
            mp%found = .true.
            mp%a_lattice_au = CDS_A_BOHR
            mp%cell_au = (/ CDS_A_BOHR, CDS_ASQ3_BOHR, CDS_C_BOHR /)  ! orthorhombic
            mp%is_diamond   = .false.                  ! V^A != 0 (broken inversion)
            mp%coulomb_ok = .true.; mp%eph_ok = .true.; mp%ii_ok = .true.
            mp%eeh_ok = .false.   ! no cited CdS carrier-carrier rate (see above)
            ! Auger recombination FORBIDDEN for CdS: no verified coefficient (the
            ! former "Haury PRB 57, 11513" C was fabricated, removed). A density
            ! gate is kept (Shah 1986) only for a user who supplies their own C
            ! via sbe_auger_c_cm6s; with auger_ok=.false. the default C is 0, so a
            ! plain yn_sbe_auger='y' run aborts (provenance) unless C is given.
            mp%auger_ok = .false.
            mp%auger_c_cm6s = 0d0;  mp%auger_n_gate_cm3 = CDS_EE_ACT_N
            mp%eps0 = CDS_EPS0;  mp%eps_inf = CDS_EPS_INF
            ! A4-CdS acoustic deformation potential: E1 = 14.5 eV [D. L. Rode,
            ! PRB 2, 1012 (1970), maintainer-supplied]; LA sound velocity
            ! ~4.25e5 cm/s from the Rode elastic constants (verify vs PDF).
            ! CRITICAL: this anomalously LARGE E1 makes the BARE DP channel
            ! unphysically strong at working densities n >= 1e18 cm^-3 -- the
            ! acoustic mode is therefore ALWAYS applied with the free-carrier
            ! screening factor S(q) = [q/(q+q_TF)]^2 built from the gathered
            ! carrier density (see apply_ring_channels); screening cuts the
            ! small-q divergence and leaves acoustics as the physical fallback
            ! for carriers cooled below hbar*omega_LO.
            mp%eph_ac_xi_ev = 14.5d0;  mp%eph_ac_cs_cmps = 4.25d5
            ! Frohlich polar-optical: a single dominant LO mode at 38 meV; the
            ! rate scale nu_sat = alpha*omega_LO is the cited Frohlich coupling.
            mp%eph_polar = .true.
            mp%eph_nph = 1
            mp%eph_hw_mev(1) = CDS_HW_LO_MEV
            mp%eph_wraw(1)   = 1d0                     ! single mode -> weight 1
            mp%eph_nu_sat_si = CDS_NU_SAT_SI
            ! Impact ionization: Keldysh soft threshold; E_th cited, prefactor is
            ! a fit parameter (left to the user via sbe_ii_prefactor).
            mp%ii_form = 'keldysh_quadratic'; mp%ii_exponent = 2d0
            mp%ii_threshold_ev = CDS_II_ETH_EV
            mp%ii_prefactor = -1d0     ! sentinel: NO cited CdS prefactor; the
            ! init requires the user to set sbe_ii_prefactor explicitly (fit param).
        case ('graphene')
            ! GAPLESS Dirac semimetal (honeycomb D6h, V^A=0, centrosymmetric).
            ! a minimal pi-model (Ramanujam 3-form-factor EPM); the Python ref
            ! emits the GS, so this registry entry only supplies the DISSIPATION
            ! constants. Hexagonal 2-atom primitive (or the rectangular 2-fold
            ! folded cell): cell_au below is informational only.
            !
            ! CITED, ENABLED channel:
            !   * e-ph: the two Kohn-anomaly optical modes E2g(Gamma,196 meV) and
            !     A1'(K,160 meV), EPC <g^2_Gamma>=0.0405, <g^2_K>=0.0994 eV^2 (x2 GW)
            !     [Piscanec PRL 93,185503 (2004); Lazzeri PRB 78,081406 (2008)].
            !     ** Both modes are NON-LOCAL on the cone (carriers relax to a
            !        DIFFERENT k), so graphene e-ph is valid ONLY through the
            !        inter-k RING (yn_sbe_superres='y'); the k-local same-k search
            !        is inert on the 2-band Dirac cone. **
            !
            ! CITED, ENABLED channel (the 2D Rana branch, wiki/07 sec.6):
            !   * Auger / carrier multiplication (auger_ok=.true. + auger_2d_rana):
            !     the gapless 2D collinear CCCV/CVVV recombination and its CVCC
            !     generation partner [R07 = F. Rana, PRB 76, 155431 (2007),
            !     Eqs. (13),(14),(17); journal text verified; lifetime benchmarks
            !     unit-tested in test_rana_2d]. Applied as the net CPTP relaxation
            !     R - G of the Dirac-cone pair density (rana_auger_dpop); at
            !     equilibrium R = G exactly (detailed balance). There is NO single
            !     C [cm^6/s] -- the rate comes from the R07 integrals evaluated on
            !     the instantaneous quasi-Fermi levels. ** Ring-gated
            !     (yn_sbe_superres='y'): needs the global population gather, like
            !     graphene e-ph. ** eps_r default = 10 (the R07 Fig.4 benchmark
            !     substrate); override via sbe_coulomb_epsilon (R07 Fig.5: eps=4
            !     SiO2 roughly doubles the rate).
            !
            ! FORBIDDEN channels (GAPLESS physics -- the gap-based maps don't apply):
            !   * impact ionization (ii_ok=.false.): the Stobbe/Keldysh fit is a
            !     (eps_kin - E_th)^a law with E_th a GAP threshold -- meaningless for
            !     a gapless cone. Graphene carrier multiplication IS the generation
            !     side of the 2D Rana channel above (thresholdless, collinear) --
            !     NOT this gap-threshold channel.
            !   * carrier-carrier (eeh_ok=.false.): no cited graphene FD-thermalization
            !     RATE (the graphene e-e coupling alpha=2.2/eps_eff is a coupling, not
            !     a rate); the gapless e-e/CM physics is the Rana branch.
            !   * Coulomb HF: ENABLED as the 2D SHEET kernel (A7, coulomb_2d):
            !     V_2D(q) = 2 pi e^2/(eps_r A (q + kappa)) with the substrate
            !     eps_r (default = the R07 benchmark 10; override via
            !     sbe_coulomb_epsilon) -- unlocks the Dirac-velocity/excitonic
            !     HF renormalization on the cone.
            !
            ! Also: a Kuhn-Zurek (single-particle wave-packet) dephasing is
            ! UNPHYSICAL for gapless Dirac carriers (coherence loss is many-body) --
            ! the SBE init aborts on graphene + sbe_decoh_* (see bloch_solver init).
            mp%found = .true.
            mp%a_lattice_au = GRAPH_A_BOHR
            mp%cell_au = (/ GRAPH_A_BOHR, GRAPH_A_BOHR*sqrt(3d0), 0d0 /)  ! informational
            mp%is_diamond = .true.    ! V^A = 0 (centrosymmetric)
            mp%eph_ok = .true.
            mp%ii_ok = .false.
            ! 2D Rana Auger/CM [R07] -- cited & enabled (ring-gated; see above)
            mp%auger_ok = .true.;  mp%auger_2d_rana = .true.
            mp%rana_vf_au = GRAPH_VF_AU;  mp%rana_eps_r = GRAPH_RANA_EPS
            mp%eeh_ok = .false.
            ! A7: 2D-sheet Sigma^HF enabled; eps_r = the substrate dielectric
            ! (R07 benchmark 10 by default, sbe_coulomb_epsilon overrides).
            mp%coulomb_ok = .true.;  mp%coulomb_2d = .true.
            mp%eps0 = GRAPH_RANA_EPS;  mp%eps_inf = GRAPH_RANA_EPS
            ! non-polar optical deformation-potential modes (no Frohlich LO branch)
            mp%eph_polar = .false.
            mp%eph_nph = 2
            mp%eph_hw_mev(1) = GRAPH_HW_E2G_MEV
            mp%eph_hw_mev(2) = GRAPH_HW_A1P_MEV
            ! raw per-mode weight = <g^2>/hw (analogue of D^2/hw); the A1' K-mode
            ! carries the x2 GW enhancement and dominates (Kohn anomaly).
            mp%eph_wraw(1) = GRAPH_G2_E2G            / GRAPH_HW_E2G_MEV
            mp%eph_wraw(2) = GRAPH_G2_A1P * GRAPH_GW_K / GRAPH_HW_A1P_MEV
            ! acoustic: D = 16 eV, v_ph = 2.0e6 cm/s [Hwang & Das Sarma,
            ! PRB 77, 195412 (2008), maintainer-supplied]
            mp%eph_ac_xi_ev = 16.0d0;  mp%eph_ac_cs_cmps = 2.0d6
            mp%eph_nu_sat_si = GRAPH_NU_SAT_SI
        case default
            mp%found = .false.
        end select
    end function get_material_params



    ! Collision-rate saturation: nu(eps) = nu_sat [1 - exp(-(eps/eps0)^n)].
    ! Smooth (no hard min() cutoff -- the derivative discontinuity destabilizes
    ! the stiff solver). Rises through ~0.5-1.5 eV and saturates at nu_sat.
    ! [Meng et al., PRB 91, 075201 (2015); Fischetti & Laux, PRB 38, 9721 (1988)]
    pure function nu_saturation(eps, nu_sat, eps0, nexp) result(nu)
        real(8), intent(in) :: eps, nu_sat, eps0, nexp
        real(8) :: nu, x
        if (eps <= 0d0 .or. eps0 <= 0d0) then
            nu = 0d0
        else
            x = (eps / eps0) ** nexp
            nu = nu_sat * (1d0 - exp(-x))
        end if
    end function nu_saturation

    ! Bose-Einstein phonon occupation N_B = 1/(exp(hw/kT) - 1). hw and kT in the
    ! same units. Guards the hw/kT -> 0 (classical) and large-argument limits.
    pure function bose_factor(hw, kT) result(nb)
        real(8), intent(in) :: hw, kT
        real(8) :: nb, x
        if (kT <= 0d0) then
            nb = 0d0                       ! T=0: no absorption
        else
            x = hw / kT
            if (x < 1d-8) then
                nb = kT / max(hw, 1d-300)  ! classical limit N_B ~ kT/hw
            else if (x > 7.0d2) then
                nb = 0d0                   ! exp overflow guard
            else
                nb = 1d0 / (exp(x) - 1d0)
            end if
        end if
    end function bose_factor

    ! Normalized Gaussian energy bin (broadened delta): area 1, width sigma.
    ! delta(Delta E) -> (1/(sqrt(2 pi) sigma)) exp(-dE^2 / (2 sigma^2)).
    pure function gaussian_bin(dE, sigma) result(w)
        real(8), intent(in) :: dE, sigma
        real(8) :: w, z
        if (sigma <= 0d0) then
            w = 0d0
        else
            z = dE / sigma
            w = exp(-0.5d0 * z * z) / (sqrt(2d0 * PI) * sigma)
        end if
    end function gaussian_bin

    ! Normalized unit-area rectangle energy bin: 1/width for |dE| < width/2, else 0.
    pure function rect_bin(dE, width) result(w)
        real(8), intent(in) :: dE, width
        real(8) :: w
        if (width <= 0d0) then
            w = 0d0
        else if (abs(dE) < 0.5d0 * width) then
            w = 1d0 / width
        else
            w = 0d0
        end if
    end function rect_bin

    ! GaAs polar-optical (Frohlich) EMISSION high-energy factor:
    !   f(E) = (1/sqrt(E)) asinh(sqrt(E/hw0 - 1))   for E > hw0, else 0.
    ! Decays slower than 1/sqrt(E); omitting the asinh under-counts ~20-30% at
    ! 2-3 eV. asinh(y) = log(y + sqrt(y^2+1)). [Fawcett-Boardman-Swain 1970]
    pure function frohlich_hi_factor(E, hw0) result(f)
        real(8), intent(in) :: E, hw0
        real(8) :: f, y
        if (E <= hw0 .or. hw0 <= 0d0) then
            f = 0d0
        else
            y = sqrt(E / hw0 - 1d0)
            f = log(y + sqrt(y * y + 1d0)) / sqrt(E)
        end if
    end function frohlich_hi_factor

    ! General impact-ionization fit rate gamma = P (eps_kin - E_th)^a, gated at
    ! threshold (Theta). a=4 GaAs Stobbe (hard), a=2 Si Keldysh (soft), a=4.6 Si
    ! full-band. All consistent units (caller supplies P in 1/(energy^a time)).
    ! [Stobbe-Redmer-Schattke PRB 49, 4494; Keldysh JETP 21, 1135; Kamakura JAP 75, 3500]
    pure function ii_rate_general(eps_kin, E_th, P, a) result(g)
        real(8), intent(in) :: eps_kin, E_th, P, a
        real(8) :: g, d
        d = eps_kin - E_th
        if (d <= 0d0) then
            g = 0d0
        else
            g = P * d ** a
        end if
    end function ii_rate_general

    ! Peak-normalized Gaussian energy-conservation weight (dimensionless, in
    ! [0,1]): exp(-dE^2 / (2 sigma^2)). Used to weight a discrete partner
    ! transition by how well it conserves energy (vs gaussian_bin, which is
    ! area-normalized for converting a state sum into a rate).
    pure function gaussian_shape(dE, sigma) result(w)
        real(8), intent(in) :: dE, sigma
        real(8) :: w, z
        if (sigma <= 0d0) then
            w = merge(1d0, 0d0, dE == 0d0)
        else
            z = dE / sigma
            w = exp(-0.5d0 * z * z)
        end if
    end function gaussian_shape

    ! Exact finite-time amplitude-damping (population-transfer) CPTP map for a
    ! single GKLS jump L = sqrt(gamma) |dst><src| in the basis where src,dst are
    ! basis states (here the Houston/adiabatic basis):
    !   rho_ss -> e^{-g} rho_ss ;  rho_dd -> rho_dd + (1 - e^{-g}) rho_ss
    !   rho_sb -> e^{-g/2} rho_sb  and  rho_bs -> e^{-g/2} rho_bs  (b /= s)
    ! with g = gamma*tau. Trace-preserving and completely positive for any
    ! gamma,tau >= 0 (a genuine GKLS map). This is the building block of the
    ! population-relaxing electron-phonon channel; the impact-ionization channel
    ! uses an identical (separately validated) copy in bloch_solver_ssbe.
    ! [Lindblad, Commun. Math. Phys. 48, 119 (1976)]
    subroutine amp_damp_channel(nba, rho, i_src, i_dst, gamma, tau)
        integer,    intent(in)    :: nba, i_src, i_dst
        complex(8), intent(inout) :: rho(nba, nba)
        real(8),    intent(in)    :: gamma, tau
        real(8) :: g, gh, transfer
        integer :: b
        if (i_src == i_dst) return
        if (gamma * tau < 1d-14) return
        g  = exp(-gamma * tau)
        gh = sqrt(g)
        transfer = (1d0 - g) * real(rho(i_src, i_src))
        do b = 1, nba
            if (b == i_src) cycle
            rho(i_src, b) = gh * rho(i_src, b)
            rho(b, i_src) = gh * rho(b, i_src)
        end do
        rho(i_src, i_src) = g * rho(i_src, i_src)
        rho(i_dst, i_dst) = rho(i_dst, i_dst) + transfer
    end subroutine amp_damp_channel

    ! Electron-hole-plasma bandgap-renormalization cube-root law (C7):
    !   Delta E_gap [eV] = -K * (n [cm^-3])^(1/3),  K ~ 1.9e-8 eV.cm.
    ! Returns a NEGATIVE shift (gap shrinks). [Vashishta & Kalia, PRB 25, 6492 (1982)]
    pure function bgr_gap_shift_ev(n_cm3, K) result(dE)
        real(8), intent(in) :: n_cm3, K
        real(8) :: dE
        if (n_cm3 <= 0d0) then
            dE = 0d0
        else
            dE = -K * n_cm3 ** (1d0 / 3d0)
        end if
    end function bgr_gap_shift_ev

    ! --- unit conversions to a.u. (for golden-rule e-ph rates, C5) -----------
    pure function mev_to_ha(E_mev) result(E)
        real(8), intent(in) :: E_mev
        real(8) :: E
        E = E_mev * 1d-3 / HA_EV
    end function mev_to_ha

    ! Deformation potential D [eV/cm] -> [Ha/Bohr]. (Si tables are in 1e8 eV/cm,
    ! so pass D_table*1e8.)
    pure function d_evcm_to_au(D_evcm) result(D)
        real(8), intent(in) :: D_evcm
        real(8) :: D
        D = D_evcm * BOHR_CM / HA_EV
    end function d_evcm_to_au

    ! Deformation potential D [eV/Angstrom] -> [Ha/Bohr]. (GaAs intervalley tables.)
    pure function d_evang_to_au(D_eva) result(D)
        real(8), intent(in) :: D_eva
        real(8) :: D
        D = D_eva * BOHR_ANG / HA_EV
    end function d_evang_to_au

    ! Mass density rho [g/cm^3] -> [m_e/Bohr^3].
    pure function rho_gcm3_to_au(rho) result(r)
        real(8), intent(in) :: rho
        real(8) :: r
        r = rho * BOHR_CM**3 / ME_G
    end function rho_gcm3_to_au

    ! Golden-rule deformation-potential rate PREFACTOR pi D^2/(rho omega), all
    ! arguments in a.u. (rate per unit final-state density-of-states; the DOS /
    ! energy-bin factor is applied by the caller). [Jacoboni-Reggiani RMP 55, 645]
    pure function golden_rule_prefactor(D_au, rho_au, omega_au) result(g)
        real(8), intent(in) :: D_au, rho_au, omega_au
        real(8) :: g
        if (rho_au <= 0d0 .or. omega_au <= 0d0) then
            g = 0d0
        else
            g = PI * D_au * D_au / (rho_au * omega_au)
        end if
    end function golden_rule_prefactor

    ! Normalized thermal emission/absorption split for one phonon mode:
    !   fe = (N_B+1)/(2 N_B+1),  fa = N_B/(2 N_B+1).
    ! fe+fa = 1 (so the per-mode total is the mode weight) and fe/fa = (N_B+1)/N_B
    ! (detailed balance). At N_B=0: fe=1, fa=0 (spontaneous emission only).
    ! Detailed balance on the ENERGY ACTUALLY EXCHANGED, not on the nominal hw_p.
    !
    ! The energy matching is Gaussian-broadened with width sigma (finite mesh, finite
    ! lifetime), so a source is connected to partners at |dE| anywhere within a few
    ! sigma of hw_p. Weighting all of them by N_B(hw_p) makes the channel enforce a
    ! Boltzmann factor exp(hw_p/kT) on a transition that actually moves delta, and the
    ! stationary state is then NOT the bath's Fermi-Dirac: for delta >> hw_p the upward
    ! rate is too large by exp((delta-hw_p)/kT), and the distribution heats until its
    ! temperature is set by sigma instead of by kT. With hw_ac = 5.4 meV and
    ! sigma = 0.1 eV this ran a doped graphene sheet to T_e ~ 2300 K with no field at
    ! all (x14 README SS7.16).
    !
    ! Evaluating N_B at delta makes the forward/backward ratio exp(delta/kT) for every
    ! pair the broadening admits, so FD(T_bath) is stationary whatever sigma is. It
    ! reduces to the nominal split when delta = hw_p, and to the symmetric 1/2 : 1/2 of
    ! a degenerate pair as delta -> 0.
    pure subroutine eph_thermal_split_de(delta, kT, fe, fa)
        real(8), intent(in)  :: delta, kT
        real(8), intent(out) :: fe, fa
        real(8) :: Nb, denom
        Nb = min(bose_factor(abs(delta), kT), 1d12)
        denom = 2d0 * Nb + 1d0
        fe = (Nb + 1d0) / denom
        fa = Nb / denom
    end subroutine eph_thermal_split_de

    pure subroutine eph_thermal_split(Nb, fe, fa)
        real(8), intent(in)  :: Nb
        real(8), intent(out) :: fe, fa
        real(8) :: denom
        denom = 2d0 * Nb + 1d0
        fe = (Nb + 1d0) / denom
        fa = Nb / denom
    end subroutine eph_thermal_split

    ! =====================================================================
    ! Part G -- dielectric screening models for the carrier-carrier channel.
    ! W(q) = V(q)/eps(q[,omega]). Pure functions; the caller supplies the
    ! carrier density, Fermi wavevector and phonon frequencies in a.u.
    ! =====================================================================

    ! Static dielectric, model (a) Thomas-Fermi/Debye: eps(q) = 1 + kappa^2/q^2.
    ! [Ashcroft & Mermin; long-wavelength limit of Lindhard theory]
    pure function eps_thomas_fermi(q, kappa2) result(eps)
        real(8), intent(in) :: q, kappa2
        real(8) :: eps
        if (q <= 0d0) then
            eps = huge(1d0)            ! q->0: perfect screening
        else
            eps = 1d0 + kappa2 / (q * q)
        end if
    end function eps_thomas_fermi

    ! Degenerate (T->0) Thomas-Fermi screening wavevector squared [a.u.]:
    !   kappa_TF^2 = 4 (3 n / pi)^(1/3) / eps  ( = 4 k_F/(pi eps) ).
    ! [arXiv:2312.13059; Ashcroft & Mermin]
    pure function tf_kappa2_degenerate(n_au, eps_bg) result(k2)
        real(8), intent(in) :: n_au, eps_bg
        real(8) :: k2
        if (n_au <= 0d0) then
            k2 = 0d0
        else
            k2 = 4d0 * (3d0 * n_au / PI) ** (1d0 / 3d0) / max(eps_bg, 1d-12)
        end if
    end function tf_kappa2_degenerate

    ! Nondegenerate (Debye-Huckel) screening wavevector squared [a.u.]:
    !   kappa_D^2 = 4 pi n / (eps k_B T).  [Ashcroft & Mermin]
    pure function debye_kappa2(n_au, eps_bg, kT_au) result(k2)
        real(8), intent(in) :: n_au, eps_bg, kT_au
        real(8) :: k2
        if (n_au <= 0d0 .or. kT_au <= 0d0) then
            k2 = 0d0
        else
            k2 = 4d0 * PI * n_au / (max(eps_bg, 1d-12) * kT_au)
        end if
    end function debye_kappa2

    ! 3D static Lindhard function F(x), x = q/(2 k_F):
    !   F(x) = 1/2 + (1-x^2)/(4x) ln|(1+x)/(1-x)|.
    ! F(0)=1 (recovers Thomas-Fermi), F(1)=1/2 (the 2k_F kink -> Friedel
    ! oscillations). [J. Lindhard, Mat.-Fys. Medd. 28, 8 (1954); Ashcroft & Mermin]
    pure function lindhard_F(x) result(F)
        real(8), intent(in) :: x
        real(8) :: F
        if (abs(x) < 1d-8) then
            F = 1d0
        else if (abs(x - 1d0) < 1d-8) then
            F = 0.5d0
        else
            F = 0.5d0 + (1d0 - x * x) / (4d0 * x) * log(abs((1d0 + x) / (1d0 - x)))
        end if
    end function lindhard_F

    ! Static Lindhard/RPA dielectric, model (b, recommended default):
    !   eps(q,0) = 1 + (kappa_TF^2/q^2) F(q/2k_F).
    ! Correct across the q range that dominates the collision integral (Thomas-
    ! Fermi over-screens for q>2k_F). [Lindhard 1954; arXiv:1206.2003]
    pure function eps_lindhard_static(q, kF, kappa2) result(eps)
        real(8), intent(in) :: q, kF, kappa2
        real(8) :: eps
        if (q <= 0d0) then
            eps = huge(1d0)
        else if (kF <= 0d0) then
            eps = 1d0 + kappa2 / (q * q)        ! no carriers: bare TF form
        else
            eps = 1d0 + kappa2 / (q * q) * lindhard_F(q / (2d0 * kF))
        end if
    end function eps_lindhard_static

    ! =====================================================================
    ! CDRB model dielectric function for the inter-k Coulomb channels
    ! [Cappellini, Del Sole, Reining, Bechstedt, PRB 47, 9892 (1993); used for
    !  Auger by Kioupakis et al. (K15), PRB 92, 035207 (2015), Eq. (8)]:
    !   eps(q) = 1 + [ (eps_inf-1)^-1 + alpha (q/q_TF)^2 + q^4/(4 wp^2) ]^-1
    ! in a.u. (hbar = m = e = 1), alpha = 1.563. Limits: eps(0) = eps_inf,
    ! eps(q->inf) -> 1 -- the large-q transparency that makes umklapp
    ! transitions strong (static eps0 over-screens them by ~10x [L90]).
    ! q_TF and wp are those of the VALENCE electron gas (n = nelec/volume).
    ! =====================================================================
    pure function eps_cdrb(q2, eps_inf, qtf2, wp2) result(eps)
        implicit none
        real(8), intent(in) :: q2, eps_inf, qtf2, wp2
        real(8) :: eps, denom
        real(8), parameter :: CDRB_ALPHA = 1.563d0
        if (eps_inf <= 1d0 + 1d-12 .or. qtf2 <= 0d0 .or. wp2 <= 0d0) then
            eps = 1d0
            return
        end if
        denom = 1d0 / (eps_inf - 1d0) + CDRB_ALPHA * q2 / qtf2 + q2 * q2 / (4d0 * wp2)
        eps = 1d0 + 1d0 / denom
    end function eps_cdrb

    ! =====================================================================
    ! Screened inter-k Coulomb weight WITH the umklapp G-sum, in the CARTESIAN
    ! metric of the actual (possibly non-orthogonal) cell:
    !   vq = sum_G 1 / [ eps(|q+G|) (|q+G|^2 + lambda^2 + q2reg) ]
    ! over the 27 neighbouring reciprocal images G = n1 b1 + n2 b2 + n3 b3,
    ! n_i in {-1,0,1} (bmat rows = b1..b3 [a.u.]; df = fractional coordinate
    ! difference of the two grid points, in (-1,1) per component so the 27
    ! images bracket the minimum image). This implements the two [L90]
    ! "must-not-drop" pieces at the overlap-free level (I(G) -> 1): the
    ! umklapp G-sum and the q-dependent eps(q); the Bloch-overlap factors
    ! I_{13}(G) I_{24}(G'-G) need the plane-wave coefficients the SBE does not
    ! carry and remain a refinement. lambda2 = free-carrier screening (Si: 0
    ! by Burt's dynamical argument [L90]); q2reg = grid-scale q->0 regulariser
    ! of the discrete BZ sum (numerical, refines with the k-grid -- NOT a
    ! physics constant, unlike the old fixed reduced-space kappa2 = 0.05).
    ! =====================================================================
    pure function interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg) result(vq)
        implicit none
        real(8), intent(in) :: df(3), bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
        real(8) :: vq, qc(3), q2
        integer :: g1, g2, g3
        vq = 0d0
        do g3 = -1, 1
            do g2 = -1, 1
                do g1 = -1, 1
                    qc(1:3) = (df(1) + g1) * bmat(1, 1:3) &
                            + (df(2) + g2) * bmat(2, 1:3) &
                            + (df(3) + g3) * bmat(3, 1:3)
                    q2 = dot_product(qc, qc)
                    vq = vq + 1d0 / (eps_cdrb(q2, eps_inf, qtf2, wp2) * (q2 + lambda2 + q2reg))
                end do
            end do
        end do
    end function interk_vq

    ! Unpack the optional t_ring_opts into inert-default locals (shared by the
    ! ring II/Auger kernels; omitted argument -> bit-identical legacy behavior).
    pure subroutine ring_opts_unpack(opts, havetab, vfloor, fk, pa, nphl, prefh, &
                                     evbm, hwl, nbbl, wrell)
        implicit none
        type(t_ring_opts), intent(in), optional :: opts
        logical, intent(out) :: havetab
        real(8), intent(out) :: vfloor, fk, pa, prefh, evbm
        integer, intent(out) :: nphl
        real(8), allocatable, intent(out) :: hwl(:), nbbl(:), wrell(:)
        havetab = .false.; vfloor = 0d0; fk = 0d0; pa = 0d0
        nphl = 0; prefh = 0d0; evbm = 0d0
        if (present(opts)) then
            havetab = opts%use_tab .and. allocated(opts%vq_tab)
            vfloor = opts%vq_floor;  fk = opts%fk_theta
            prefh  = opts%pref_h;    evbm = opts%evbm
            if (opts%phassist > 0d0 .and. opts%nph > 0 .and. allocated(opts%hw)) then
                pa = opts%phassist;  nphl = opts%nph
                hwl = opts%hw(1:nphl); nbbl = opts%nbb(1:nphl); wrell = opts%wrel(1:nphl)
            end if
        end if
        if (nphl == 0) then
            allocate(hwl(1), nbbl(1), wrell(1))
            hwl = 0d0; nbbl = 0d0; wrell = 0d0
        end if
    end subroutine ring_opts_unpack

    ! A1: energy-conservation shape with optional phonon-assisted sidebands.
    ! Forward (II) orientation: phonon EMISSION (N+1) shifts the electronic
    ! surplus to etgt = +hw -> delta(etgt - hw); absorption N -> delta(etgt + hw).
    ! rev = .true. (Auger, the time-reverse) SWAPS the Bose factors, so together
    ! with the FD occupation identity the detailed balance holds EXACTLY per
    ! sideband: f1f2(1-f3)(1-f4)*(N+1) at etgt=+hw balances (1-f1)(1-f2)f3f4*N.
    pure function shape_assist(etgt, sigma, pa, nph, hw, nbb, wrel, rev) result(shp)
        implicit none
        real(8), intent(in) :: etgt, sigma, pa
        integer, intent(in) :: nph
        real(8), intent(in) :: hw(*), nbb(*), wrel(*)
        logical, intent(in) :: rev
        real(8) :: shp, we, wa
        integer :: ip
        shp = gaussian_shape(etgt, sigma)
        if (pa <= 0d0) return
        do ip = 1, nph
            if (rev) then
                we = nbb(ip);        wa = nbb(ip) + 1d0
            else
                we = nbb(ip) + 1d0;  wa = nbb(ip)
            end if
            shp = shp + pa * wrel(ip) * (we * gaussian_shape(etgt - hw(ip), sigma) &
                                       + wa * gaussian_shape(etgt + hw(ip), sigma))
        end do
    end function shape_assist

    ! B1: precompute interk_vq over ALL SIGNED index differences
    ! d = kidx(:,i1) - kidx(:,i1p), d(i) in [-(n_i-1), n_i-1], i.e.
    ! (2n1-1)(2n2-1)(2n3-1) entries. SIGNED (not mod-n) indexing keeps the
    ! truncated 27-image umklapp sum BIT-IDENTICAL to the direct call (the
    ! image window is not translation-invariant). Lookup index (1-based):
    !   idx = 1 + (d1+n1-1) + (2n1-1)*[(d2+n2-1) + (2n2-1)*(d3+n3-1)].
    ! Rebuild once per step (lambda2 may be dynamic) -- O(8 nk) vs the former
    ! O(nk^2) evaluations per hot state per pass per kernel.
    pure subroutine build_vq_table(kn, bmat, eps_inf, qtf2, wp2, lambda2, q2reg, tab)
        implicit none
        integer, intent(in)  :: kn(3)
        real(8), intent(in)  :: bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
        real(8), intent(out) :: tab((2*kn(1)-1)*(2*kn(2)-1)*(2*kn(3)-1))
        integer :: d1, d2, d3, idx
        real(8) :: df(3)
        idx = 0
        do d3 = -(kn(3)-1), kn(3)-1
            do d2 = -(kn(2)-1), kn(2)-1
                do d1 = -(kn(1)-1), kn(1)-1
                    idx = idx + 1
                    df = (/ dble(d1)/dble(max(kn(1),1)), &
                            dble(d2)/dble(max(kn(2),1)), &
                            dble(d3)/dble(max(kn(3),1)) /)
                    tab(idx) = interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg)
                end do
            end do
        end do
    end subroutine build_vq_table

    ! A4 screening: Thomas-Fermi factor S(q) = [q/(q + q_TF)]^2 for the
    ! ACOUSTIC deformation-potential e-ph mode, tabulated over the SAME signed
    ! index differences as build_vq_table. q = the nearest-image Cartesian
    ! distance |dk + G| (min over the 27 images). q_TF from the INSTANTANEOUS
    ! free-carrier density (Debye/TF crossover, same primitives as lambda^2).
    ! q_TF <= 0 -> S = 1 (bare, the no-carrier limit). S(q=0) = 0: perfectly
    ! screened forward scattering -- exactly the small-q cut that keeps the
    ! large-E1 materials (CdS E1 = 14.5 eV [Rode 1970]) physical at
    ! n >= 1e18 cm^-3. (The maintainer-specified TF form; the full Lindhard
    ! eps_lindhard_static is available as a refinement.)
    pure subroutine build_acscreen_table(kn, bmat, qtf, tab)
        implicit none
        integer, intent(in)  :: kn(3)
        real(8), intent(in)  :: bmat(3,3), qtf
        real(8), intent(out) :: tab((2*kn(1)-1)*(2*kn(2)-1)*(2*kn(3)-1))
        integer :: d1, d2, d3, g1, g2, g3, idx
        real(8) :: df(3), qc(3), q2, qmin
        if (qtf <= 0d0) then
            tab = 1d0
            return
        end if
        idx = 0
        do d3 = -(kn(3)-1), kn(3)-1
            do d2 = -(kn(2)-1), kn(2)-1
                do d1 = -(kn(1)-1), kn(1)-1
                    idx = idx + 1
                    df = (/ dble(d1)/dble(max(kn(1),1)), &
                            dble(d2)/dble(max(kn(2),1)), &
                            dble(d3)/dble(max(kn(3),1)) /)
                    qmin = huge(1d0)
                    do g3 = -1, 1
                        do g2 = -1, 1
                            do g1 = -1, 1
                                qc(1:3) = (df(1)+g1)*bmat(1,1:3) + (df(2)+g2)*bmat(2,1:3) &
                                        + (df(3)+g3)*bmat(3,1:3)
                                q2 = dot_product(qc, qc)
                                if (q2 < qmin) qmin = q2
                            end do
                        end do
                    end do
                    qmin = sqrt(qmin)
                    tab(idx) = (qmin / (qmin + qtf))**2
                end do
            end do
        end do
    end subroutine build_acscreen_table

    ! =====================================================================
    ! GRAPHENE 2D Auger / carrier multiplication -- the [R07] branch
    ! (F. Rana, arXiv:0705.1204v2; wiki/07 sec.6). Gapless Dirac spectrum
    ! E_s(k) = s*v*|k| (a.u., hbar=1), spin x valley degeneracy g=4. The CCCV
    ! recombination rate collapses (collinear phase space, overlaps -> 1) to
    ! the 3D integral of [R07]; CVVV(n,p) = CCCV(p,n); 1/tau_r = R/min(n,p).
    ! 2D rates in a.u.^-2 * a.u.t^-1 -- there is NO single C [cm^6/s]:
    ! validation is by the cited LIFETIME targets (test_rana_2d).
    ! =====================================================================

    ! Carrier density (per area, a.u.) of one Dirac branch at quasi-Fermi mu:
    ! n(mu) = (g/2pi) int_0^inf k f(vk - mu) dk,  g = 4 (spin x K,K').
    pure function dirac_n_2d(mu, kT, v) result(n)
        implicit none
        real(8), intent(in) :: mu, kT, v
        real(8) :: n, k, dk, kmax, f
        integer :: i
        integer, parameter :: NGRID = 400
        kmax = (max(mu, 0d0) + 20d0 * kT) / v
        dk = kmax / NGRID
        n = 0d0
        do i = 1, NGRID
            k = (i - 0.5d0) * dk
            f = 1d0 / (exp(min((v * k - mu) / kT, 60d0)) + 1d0)
            n = n + k * f
        end do
        n = n * dk * 4d0 / (2d0 * pi)
    end function dirac_n_2d

    ! Kinetic-energy density (per area, a.u.) of one Dirac branch at (mu, kT),
    ! measured from the Dirac point: eps = (g/2pi) int_0^inf k (v k) f(vk - mu) dk.
    ! Intrinsic (mu = 0): eps = (2/pi)(3 zeta(3)/2)(kT)^3/v^2 (closed form, unit-tested).
    pure function dirac_e_2d(mu, kT, v) result(e)
        implicit none
        real(8), intent(in) :: mu, kT, v
        real(8) :: e, k, dk, kmax, f
        integer :: i
        integer, parameter :: NGRID = 400
        kmax = (max(mu, 0d0) + 20d0 * kT) / v
        dk = kmax / NGRID
        e = 0d0
        do i = 1, NGRID
            k = (i - 0.5d0) * dk
            f = 1d0 / (exp(min((v * k - mu) / kT, 60d0)) + 1d0)
            e = e + k * (v * k) * f
        end do
        e = e * dk * 4d0 / (2d0 * pi)
    end function dirac_e_2d

    ! =====================================================================
    ! Two-temperature model of the Coulomb sector (graphene, wiki/12 sec. 3):
    ! from the gathered electron density n, hole density p and total kinetic
    ! energy density eps (all per area, a.u.) solve for a COMMON carrier
    ! temperature kT_e and the two quasi-Fermi levels,
    !     n = n(mu_c, T_e),  p = n(mu_h, T_e),  eps = e(mu_c, T_e) + e(mu_h, T_e).
    ! At fixed (n, p) the energy is monotone increasing in T -> bisection in
    ! ln T with the density inversions nested inside. The lattice (phonon bath)
    ! stays at kT_min: T_e is READ from the distribution every ring step; the
    ! cooling toward the lattice is done by the e-ph kinetics themselves, so
    ! this is the two-temperature model without a separate T_e rate equation.
    ! Falls back to kT_min when there are no carriers or when the energy is at
    ! or below its value at the bath temperature (carriers not hotter than the
    ! lattice); clamps at kT_max. mu_h is the HOLE quasi-Fermi level measured
    ! downward (the valence-branch electron level is mu_v = -mu_h).
    ! =====================================================================
    pure subroutine dirac_fit_te(n, p, eps, v, kT_min, kT_max, kTe, mu_c, mu_h)
        implicit none
        real(8), intent(in)  :: n, p, eps, v, kT_min, kT_max
        real(8), intent(out) :: kTe, mu_c, mu_h
        real(8) :: lo, hi, mid, nn, pp
        integer :: it
        real(8), parameter :: n_eps = 1d-14
        nn = max(n, 0d0);  pp = max(p, 0d0)
        kTe = kT_min
        if (nn < n_eps .and. pp < n_eps) then
            mu_c = dirac_mu_2d(nn, kTe, v);  mu_h = dirac_mu_2d(pp, kTe, v)
            return
        end if
        if (eps <= e_of(kT_min)) then
            kTe = kT_min
        else if (eps >= e_of(kT_max)) then
            kTe = kT_max
        else
            lo = log(kT_min);  hi = log(kT_max)
            do it = 1, 60
                mid = 0.5d0 * (lo + hi)
                if (e_of(exp(mid)) < eps) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            kTe = exp(0.5d0 * (lo + hi))
        end if
        mu_c = dirac_mu_2d(nn, kTe, v)
        mu_h = dirac_mu_2d(pp, kTe, v)
    contains
        pure function e_of(kT) result(e)
            real(8), intent(in) :: kT
            real(8) :: e
            e = dirac_e_2d(dirac_mu_2d(nn, kT, v), kT, v) + dirac_e_2d(dirac_mu_2d(pp, kT, v), kT, v)
        end function e_of
    end subroutine dirac_fit_te

    ! Quasi-Fermi level of a 2D Dirac branch from its carrier density (bisection).
    pure function dirac_mu_2d(n_au, kT, v) result(mu)
        implicit none
        real(8), intent(in) :: n_au, kT, v
        real(8) :: mu, lo, hi, mid
        integer :: it
        lo = -40d0 * kT
        hi =  40d0 * kT + v * sqrt(max(pi * n_au, 0d0)) * 2d0
        do it = 1, 80
            mid = 0.5d0 * (lo + hi)
            if (dirac_n_2d(mid, kT, v) < n_au) then
                lo = mid
            else
                hi = mid
            end if
        end do
        mu = 0.5d0 * (lo + hi)
    end function dirac_mu_2d

    ! Thomas-Fermi screening vector of the Dirac gas [R07 Eq. (13)]:
    ! Q_TF = (e^2 KT)/(pi eps_inf hbar^2 v^2) log[(e^{Ef+1/KT}+1)(e^{-Ef-1/KT}+1)]
    ! with eps_inf = eps_r * eps_0. In Hartree a.u. (e = hbar = 1, eps_0 = 1/4pi)
    ! this is Q_TF = 4 KT log[...] / (eps_r v^2). mu_c / mu_v are the conduction /
    ! valence quasi-Fermi levels (Ef+1 / Ef-1 of [R07]).
    pure function rana_qtf(mu_c, mu_v, kT, v, eps_inf) result(qtf)
        implicit none
        real(8), intent(in) :: mu_c, mu_v, kT, v, eps_inf
        real(8) :: qtf, l1, l2
        ! log(e^x + 1) evaluated overflow-safely
        l1 = max( mu_c / kT, 0d0) + log(1d0 + exp(-abs( mu_c / kT)))
        l2 = max(-mu_v / kT, 0d0) + log(1d0 + exp(-abs(-mu_v / kT)))
        qtf = 4d0 * kT / (eps_inf * v * v) * (l1 + l2)
    end function rana_qtf

    ! =====================================================================
    ! 2D colmem analog (graphene, wiki/10 sec. 8.11): the memory LINE of the
    ! Coulomb (Rana) collision sector = the 2D Dirac plasmon of the
    ! instantaneous e-h plasma. Long-wavelength plasmon of the Dirac gas
    ! [Hwang & Das Sarma, PRB 75, 205418 (2007)]:  w_pl^2(q) = 2 e^2 E_F q /
    ! (kappa hbar^2)  (Gaussian). For the two-component (electron + hole)
    ! plasma the Drude weights add, and at finite T each branch's E_F is the
    ! intraband Drude weight  W(mu) = 2 k_B T ln[2 cosh(mu/2 k_B T)]
    ! [Falkovsky & Varlamov, Eur. Phys. J. B 56, 281 (2007)]  (-> |mu| when
    ! degenerate, -> 2 k_B T ln 2 for the intrinsic thermal plasma). It is
    ! evaluated at the collision's own screening momentum q = Q_TF [R07 Eq.
    ! (13)] -- the momentum-transfer scale of the screened matrix element --
    ! so 1/w_pl is the build-up time of screening (Haug & Jauho, Quantum
    ! Kinetics, the "build-up of screening" chapter): the memory time of the
    ! Coulomb collision, exactly as the phonon energy is for e-ph.
    ! Hartree a.u. (e = hbar = 1, eps_inf = eps_r/4pi absorbed as in rana_qtf):
    !     w_pl^2 = 2 (W_c + W_v) Q_TF / eps_r.
    ! No new free parameters: T, eps_r and the mu's are the Rana channel's own.
    ! =====================================================================
    pure function dirac_plasmon_2d(mu_c, mu_v, kT, v, eps_r) result(wpl)
        implicit none
        real(8), intent(in) :: mu_c, mu_v, kT, v, eps_r
        real(8) :: wpl, wc, wv, qtf
        wc  = drude_weight(mu_c)
        wv  = drude_weight(mu_v)
        qtf = rana_qtf(mu_c, mu_v, kT, v, eps_r)
        wpl = sqrt(max(2d0 * (wc + wv) * qtf / eps_r, 0d0))
    contains
        ! 2 kT ln[2 cosh(mu/2kT)] = |mu| + 2kT ln(1 + e^{-|mu|/kT}), overflow-safe
        pure function drude_weight(mu) result(w)
            real(8), intent(in) :: mu
            real(8) :: w
            w = abs(mu) + 2d0 * kT * log(1d0 + exp(-abs(mu) / kT))
        end function drude_weight
    end function dirac_plasmon_2d

    ! CCCV Auger recombination rate per area, [R07 Eq. (14)] (the paper's main
    ! result; verified against the journal text supplied by the maintainer):
    !   R = (1/(hbar^2 v)) int_0^inf dk1/2pi int_0^inf dk2/2pi int_{k2}^inf dQ/2pi
    !       |M(k1,k2,Q)|^2 * sqrt((k1+Q)(Q-k2) k1 k2)
    !       * [1-f_-1(Q-k2)] [1-f_+1(k1+Q)] f_+1(k1) f_+1(k2)
    ! (the sqrt of the collinear collapse is in the NUMERATOR). The matrix
    ! element [R07 Eqs. (10)-(12)], overlaps = 1 on the collinear line:
    !   |M|^2 = |M_d|^2 + |M_e|^2 + |M_d - M_e|^2
    !   M_d = e^2/(2 eps_inf (Q + Q_TF)),  M_e = e^2/(2 eps_inf (|Q+k1-k2| + Q_TF))
    ! which in Hartree a.u. (e^2 = 1, eps_inf = eps_r/(4 pi)) is
    !   M_d = 2 pi/(eps_r (Q + Q_TF)) etc.
    ! reverse=.true. gives the GENERATION (impact-ionization CVCC) partner
    ! [R07 Eq. (17)]: occupations [1-f_c(k1)][1-f_c(k2)] f_v(Q-k2) f_c(Q+k1);
    ! at equilibrium (single mu for both branches) generation = recombination
    ! exactly (detailed balance -- unit-tested).
    ! f_v is the VALENCE-branch electron occupation: f_v(k) = f((-v*k - mu_v)/kT).
    ! Validated against the cited [R07] benchmarks in test_rana_2d: minority
    ! lifetime ~1.1 ps at p = 1e12 cm^-2 / 300 K / eps_r = 10, tau_r > 1 ps
    ! below 1e12 (all T), > 5 ps below 1e11.
    function rana_rcccv(mu_c, mu_v, kT, v, eps_inf, reverse) result(R)
        implicit none
        real(8), intent(in) :: mu_c, mu_v, kT, v, eps_inf
        logical, intent(in) :: reverse
        real(8) :: R
        integer, parameter :: N1 = 48, NQ = 64
        real(8) :: kmax, dk, dq, k1, k2, Q, qtf
        real(8) :: fc1, fc2, fcp, fvp, occ, md, me, m2, w
        integer :: i1, i2, iq

        qtf  = rana_qtf(mu_c, mu_v, kT, v, eps_inf)
        kmax = (max(mu_c, -mu_v, 0d0) + 18d0 * kT) / v
        dk   = kmax / N1
        R    = 0d0
        do i1 = 1, N1
            k1 = (i1 - 0.5d0) * dk
            fc1 = fdocc( (v * k1 - mu_c) / kT )
            do i2 = 1, N1
                k2 = (i2 - 0.5d0) * dk
                fc2 = fdocc( (v * k2 - mu_c) / kT )
                dq = (2d0 * kmax) / NQ
                do iq = 1, NQ
                    Q = k2 + (iq - 0.5d0) * dq
                    ! occupations of the final states: CB electron at k1+Q,
                    ! VB hole left at Q-k2 (valence electron REMOVED there)
                    fcp = fdocc( (v * (k1 + Q) - mu_c) / kT )
                    fvp = fdocc( (-v * (Q - k2) - mu_v) / kT )
                    if (reverse) then
                        occ = fvp * fcp * (1d0 - fc1) * (1d0 - fc2)
                    else
                        occ = (1d0 - fvp) * (1d0 - fcp) * fc1 * fc2
                    end if
                    if (occ < 1d-30) cycle
                    md = 2d0 * pi / (eps_inf * (Q + qtf))
                    me = 2d0 * pi / (eps_inf * (abs(k1 - k2 + Q) + qtf))
                    m2 = (md - me)**2 + md*md + me*me
                    ! [R07 Eq. (14)]: the collinear-collapse sqrt is a NUMERATOR factor
                    w  = m2 * sqrt(max((k1 + Q) * (Q - k2) * k1 * k2, 0d0))
                    R  = R + w * occ * dq
                end do
            end do
        end do
        R = R * dk * dk / (2d0 * pi)**3 / (v)
    contains
        pure function fdocc(x) result(f)
            real(8), intent(in) :: x
            real(8) :: f
            f = 1d0 / (exp(min(max(x, -60d0), 60d0)) + 1d0)
        end function fdocc
    end function rana_rcccv

    ! =====================================================================
    ! 2D Rana Auger/CM as a CPTP population channel (the wiki/00 TODO-1 wiring).
    ! ---------------------------------------------------------------------
    ! Net pair relaxation of the Dirac-cone populations by the [R07] rates:
    !   R = R_CCCV + R_CVVV  (recombination, Eq. 14),
    !   G = G_CVCC + G_VCCC  (generation / carrier multiplication, Eq. 17),
    ! evaluated on the INSTANTANEOUS quasi-Fermi levels mu_c(n), mu_v(p) of the
    ! gathered conduction/valence populations (dirac_mu_2d inversion). The net
    ! rate R - G is applied as a uniform-fractional population transfer
    ! CB -> VB (R > G, net recombination) or VB -> CB (G > R, net carrier
    ! multiplication), which is exactly the "R - G = (n - n0)/tau_r" relaxation
    ! of the pair density: at equilibrium R = G identically (detailed balance,
    ! unit-tested in test_rana_2d), so dpop == 0 -- the fixed point.
    !
    ! CPTP by construction:
    !   * trace: sum(dpop) = -dN + dN = 0 exactly (removal ~ f/sum(f) balances
    !     addition ~ room/sum(room));
    !   * bounds: dN is saturated below both the available source population
    !     and the available destination phase space by the smooth cap
    !     dN = cap*(1 - exp(-dN_lin/cap)), cap = min(avail, room) -- the same
    !     1-exp form as the other channels (never overshoots, linear for small
    !     rates); per-state amounts dN*f/avail <= f and dN*room_a/room <= room_a.
    ! Energy bookkeeping: the released/absorbed energy goes to the third
    ! carrier IN the R07 integrand; at this collapsed rate-model level the
    ! populations move without an explicit hot tail (e-ph + the ring channels
    ! thermalize) -- exactly the wiring the maintainer specified in wiki/00.
    !
    ! Units: f = per-k populations (sum_k f/nk = electrons per cell);
    ! area = 2D cell area [a.u.^2]; the R07 rates are per area per a.u.t, so
    ! dN_lin = |R-G| * area * nk * tau in population units. rnet_out [a.u.^-2
    ! a.u.t^-1] is returned for the tau_r diagnostic print.
    ! =====================================================================
    subroutine rana_auger_dpop(nk, nba, eval, f, occ_max, iv, ic, area, kT, vf, &
                               eps_r, tau, dpop, rnet_out, e_src_lo, e_src_hi, n2d_in, p2d_in)
        implicit none
        integer, intent(in)  :: nk, nba, iv, ic
        real(8), intent(in)  :: eval(nba, nk), f(nba, nk), occ_max, area, kT, vf, eps_r, tau
        real(8), intent(out) :: dpop(nba, nk), rnet_out
        ! Cost-preserving full-basis mask (see ii_interk_dpop). This mean-field
        ! channel is O(nk*nba) (no O(nk^2) source loop), so the window restores
        ! no cost here -- it is threaded for INTERFACE uniformity and semantic
        ! consistency: when a frozen active window is set, the pair population is
        ! only drained from / filled into in-window states (trace-exact within
        ! the window). The bulk R07 rate itself (n2d/p2d) stays the physical
        ! full-basis density. No-op when absent (default all-active; graphene
        ! always runs all-active, so this path is dormant there).
        real(8), intent(in), optional :: e_src_lo, e_src_hi
        ! 2D colmem analog (wiki/10 sec. 8.11): memory-filtered SOURCE densities.
        ! When present, the R07 rates (quasi-Fermi levels, R and G) are evaluated
        ! on these instead of the instantaneous sums, so the fast virtual share
        ! of the Houston populations does not feed the Coulomb collision; the
        ! transfer stencils (avail/room) still use the instantaneous f, so the
        ! CPTP bounds are unchanged.
        real(8), intent(in), optional :: n2d_in, p2d_in
        real(8) :: n2d, p2d, mu_c, mu_v, rrec, rgen, rnet
        real(8) :: avail, room, dn_lin, cap, dn
        real(8) :: ec_bar, ev_bar, e_rel, we, sw, dn2, wgt
        integer :: a, ik
        logical :: do_mask
        logical, allocatable :: inwin(:, :)
        real(8), parameter :: n_eps = 1d-14

        dpop = 0d0
        rnet_out = 0d0
        if (iv < 1 .or. ic > nba .or. ic <= iv) return
        if (area <= 0d0 .or. kT <= 0d0 .or. vf <= 0d0) return

        do_mask = present(e_src_lo) .and. present(e_src_hi)
        if (do_mask) then
            allocate(inwin(nba, nk))
            inwin(:, :) = (eval(:, :) >= e_src_lo .and. eval(:, :) <= e_src_hi)
        end if

        ! Gathered CB electron / VB hole sheet densities [a.u.^-2].
        n2d = sum(f(ic:nba, :)) / (dble(nk) * area)
        p2d = (dble(iv) * occ_max - sum(f(1:iv, :)) / dble(nk)) / area
        p2d = max(p2d, 0d0)
        if (present(n2d_in)) n2d = max(n2d_in, 0d0)
        if (present(p2d_in)) p2d = max(p2d_in, 0d0)
        if (n2d < n_eps .and. p2d < n_eps) return

        ! Instantaneous quasi-Fermi levels and the [R07] rates (CCCV + CVVV
        ! recombination; their reverse generation partners).
        mu_c =  dirac_mu_2d(n2d, kT, vf)
        mu_v = -dirac_mu_2d(p2d, kT, vf)
        rrec = rana_rcccv( mu_c,  mu_v, kT, vf, eps_r, .false.) &
             + rana_rcccv(-mu_v, -mu_c, kT, vf, eps_r, .false.)
        rgen = rana_rcccv( mu_c,  mu_v, kT, vf, eps_r, .true.) &
             + rana_rcccv(-mu_v, -mu_c, kT, vf, eps_r, .true.)
        rnet = rrec - rgen
        rnet_out = rnet
        if (abs(rnet) < 1d-300) return

        dn_lin = abs(rnet) * area * dble(nk) * tau      ! population units
        if (rnet > 0d0) then
            ! net RECOMBINATION: remove electrons from the CB (~f), fill VB
            ! holes (~room). avail/room in the same population units.
            if (do_mask) then
                avail = sum(f(ic:nba, :),           mask=inwin(ic:nba, :))
                room  = sum(occ_max - f(1:iv, :),   mask=inwin(1:iv, :))
            else
                avail = sum(f(ic:nba, :))
                room  = dble(iv) * occ_max * dble(nk) - sum(f(1:iv, :))
            end if
        else
            ! net GENERATION (carrier multiplication): remove from the VB,
            ! promote into empty CB phase space.
            if (do_mask) then
                avail = sum(f(1:iv, :),             mask=inwin(1:iv, :))
                room  = sum(occ_max - f(ic:nba, :), mask=inwin(ic:nba, :))
            else
                avail = sum(f(1:iv, :))
                room  = dble(nba - ic + 1) * occ_max * dble(nk) - sum(f(ic:nba, :))
            end if
        end if
        cap = min(avail, room)
        if (cap < n_eps) return
        dn = cap * (1d0 - exp(-dn_lin / cap))           ! smooth CPTP saturation

        if (rnet > 0d0) then
            do ik = 1, nk
                do a = ic, nba
                    if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                    dpop(a, ik) = dpop(a, ik) - dn * f(a, ik) / avail
                end do
                do a = 1, iv
                    if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                    dpop(a, ik) = dpop(a, ik) + dn * (occ_max - f(a, ik)) / room
                end do
            end do
        else
            do ik = 1, nk
                do a = 1, iv
                    if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                    dpop(a, ik) = dpop(a, ik) - dn * f(a, ik) / avail
                end do
                do a = ic, nba
                    if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                    dpop(a, ik) = dpop(a, ik) + dn * (occ_max - f(a, ik)) / room
                end do
            end do
        end if

        ! ---- A6: energy bookkeeping -- the pair energy goes to a THIRD carrier.
        ! Recombination (CCCV): the mean electronic energy released per pair,
        ! E_rel = Ec_bar (f-weighted CB mean, the removed distribution) minus
        ! Ev_bar (room-weighted VB mean, the refilled states), is absorbed by a
        ! second CB electron -> shuffle dn2 of CB population UPWARD by E_rel
        ! with a thermal-width Gaussian energy match. Generation (CVCC) is the
        ! mirror: the created pair's energy is TAKEN from a hot electron ->
        ! shuffle DOWNWARD. Trace-neutral by construction; per-state caps keep
        ! the map CPTP; skipped gracefully when the target energy has no phase
        ! space inside the active window (the energy leaves the basis -- the
        ! honest statement of a finite band set).
        ! energy-match width: thermal, but never narrower than ~3 mean level
        ! spacings of the discrete CB spectrum (else the Gaussian falls between
        ! grid levels and the shuffle silently vanishes on coarse meshes).
        we = max(kT, 1d-4, &
                 3d0 * (maxval(eval(ic:nba,:)) - minval(eval(ic:nba,:))) &
                     / dble(max((nba - ic + 1) * nk, 1)))
        if (rnet > 0d0) then
            if (do_mask) then
                ec_bar = sum(eval(ic:nba,:) * f(ic:nba,:), mask=inwin(ic:nba,:)) / max(avail, n_eps)
                ev_bar = sum(eval(1:iv,:) * (occ_max - f(1:iv,:)), mask=inwin(1:iv,:)) / max(room, n_eps)
            else
                ec_bar = sum(eval(ic:nba,:) * f(ic:nba,:)) / max(avail, n_eps)
                ev_bar = sum(eval(1:iv,:) * (occ_max - f(1:iv,:))) / max(room, n_eps)
            end if
            e_rel  = ec_bar - ev_bar
            sw = 0d0
            do ik = 1, nk
                do a = ic, nba
                    if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                    sw = sw + (occ_max - f(a,ik)) &
                            * exp(-0.5d0*((eval(a,ik) - (ec_bar + e_rel))/we)**2)
                end do
            end do
            dn2 = min(dn, 0.999d0*max(avail - dn, 0d0), &
                      sw * max(1d0 - dn/max(room,n_eps), 0d0))
            if (dn2 > n_eps .and. sw > n_eps) then
                do ik = 1, nk
                    do a = ic, nba
                        if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                        wgt = (occ_max - f(a,ik)) &
                            * exp(-0.5d0*((eval(a,ik) - (ec_bar + e_rel))/we)**2)
                        dpop(a,ik) = dpop(a,ik) - dn2 * f(a,ik) / avail &
                                                + dn2 * wgt / sw
                    end do
                end do
            end if
        else
            ! generation: created CB electrons land ~room-weighted (ec_bar) and
            ! the absorbed energy e_rel comes from electrons near ec_bar+e_rel.
            if (do_mask) then
                ec_bar = sum(eval(ic:nba,:) * (occ_max - f(ic:nba,:)), mask=inwin(ic:nba,:)) / max(room, n_eps)
                ev_bar = sum(eval(1:iv,:) * f(1:iv,:), mask=inwin(1:iv,:)) / max(avail, n_eps)
            else
                ec_bar = sum(eval(ic:nba,:) * (occ_max - f(ic:nba,:))) / max(room, n_eps)
                ev_bar = sum(eval(1:iv,:) * f(1:iv,:)) / max(avail, n_eps)
            end if
            e_rel  = ec_bar - ev_bar
            sw = 0d0
            do ik = 1, nk
                do a = ic, nba
                    if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                    sw = sw + f(a,ik) &
                            * exp(-0.5d0*((eval(a,ik) - (ec_bar + e_rel))/we)**2)
                end do
            end do
            dn2 = min(dn, sw, 0.999d0*max(room - dn, 0d0))
            if (dn2 > n_eps .and. sw > n_eps) then
                do ik = 1, nk
                    do a = ic, nba
                        if (do_mask) then; if (.not. inwin(a, ik)) cycle; end if
                        wgt = f(a,ik) &
                            * exp(-0.5d0*((eval(a,ik) - (ec_bar + e_rel))/we)**2)
                        dpop(a,ik) = dpop(a,ik) - dn2 * wgt / sw &
                                                + dn2 * (occ_max - f(a,ik)) / max(room,n_eps)
                    end do
                end do
            end if
        end if
        if (do_mask) deallocate(inwin)
    end subroutine rana_auger_dpop

    ! Bulk plasmon frequency squared [a.u.]: omega_p^2 = 4 pi n/(eps_inf m*).
    ! (the high-frequency plasmon screens with eps_inf, not eps_0). [Mahan]
    pure function plasmon_freq2(n_au, eps_inf, mstar) result(wp2)
        real(8), intent(in) :: n_au, eps_inf, mstar
        real(8) :: wp2
        if (n_au <= 0d0) then
            wp2 = 0d0
        else
            wp2 = 4d0 * PI * n_au / (max(eps_inf, 1d-12) * max(mstar, 1d-12))
        end if
    end function plasmon_freq2

    ! Coupled LO-phonon-plasmon (LOPC) branches, model (c, GaAs only):
    !   omega_{L+/-}^2 = 1/2 [ (wp^2+wLO^2) +/- sqrt((wp^2+wLO^2)^2 - 4 wp^2 wTO^2) ].
    ! Vieta: wLp2+wLm2 = wp^2+wLO^2, wLp2*wLm2 = wp^2 wTO^2. The branches
    ! anticross at wp = wLO. [Varga PR 137, A1896; Mooradian-McWhorter PR 177, 1231]
    pure subroutine lopc_branches(wp2, wLO2, wTO2, wLp2, wLm2)
        real(8), intent(in)  :: wp2, wLO2, wTO2
        real(8), intent(out) :: wLp2, wLm2
        real(8) :: b, disc
        b = wp2 + wLO2
        disc = sqrt(max(b * b - 4d0 * wp2 * wTO2, 0d0))
        wLp2 = 0.5d0 * (b + disc)
        wLm2 = 0.5d0 * (b - disc)
    end subroutine lopc_branches

    ! Part C3: energy-conserving final-state weights (broadened-delta energy
    ! bins; NO tetrahedra, NO Monte-Carlo). For each candidate energy e_cand(i),
    !   w(i) = gaussian_bin(e_cand(i) - e_target, sigma)   if |dE| <= max_window
    !        = 0                                           otherwise,
    ! and wsum = sum(w). max_window realizes the energy-windowed expanding-radius
    ! search (the caller grows it with the hot-carrier excess energy). The caller
    ! normalizes the chosen transitions by wsum. Deterministic (CPTP-safe).
    ! [Stobbe PRB 49, 4494 (1994); Kunikiyo JAP 75, 297 (1994)]
    pure subroutine energy_partner_weights(e_target, e_cand, n_cand, sigma, &
                                           max_window, w, wsum)
        integer, intent(in)  :: n_cand
        real(8), intent(in)  :: e_target, e_cand(n_cand), sigma, max_window
        real(8), intent(out) :: w(n_cand), wsum
        integer :: i
        real(8) :: dE
        wsum = 0d0
        do i = 1, n_cand
            dE = e_cand(i) - e_target
            if (abs(dE) <= max_window) then
                w(i) = gaussian_bin(dE, sigma)
            else
                w(i) = 0d0
            end if
            wsum = wsum + w(i)
        end do
    end subroutine energy_partner_weights

    ! =====================================================================
    ! Part F -- carrier-carrier (e-e/e-h) thermalization map and the
    ! number-and-energy-conserving Fermi-Dirac fit it relaxes toward.
    ! =====================================================================

    pure function fermi_dirac(beta, e, mu) result(f)
        real(8), intent(in) :: beta, e, mu
        real(8) :: f, x
        x = beta * (e - mu)
        if (x > 7d2) then
            f = 0d0
        else if (x < -7d2) then
            f = 1d0
        else
            f = 1d0 / (exp(x) + 1d0)
        end if
    end function fermi_dirac

    ! For fixed beta, bisect mu so sum_a f_FD(eps_a) = ntot; return E and ftgt.
    subroutine e_at_beta(nlev, eps, ntot, beta, Eout, mu, ftgt)
        integer, intent(in)  :: nlev
        real(8), intent(in)  :: eps(nlev), ntot, beta
        real(8), intent(out) :: Eout, mu, ftgt(nlev)
        real(8) :: mlo, mhi, mmid, nsum
        integer :: it, a
        mlo = minval(eps) - 1d0; mhi = maxval(eps) + 1d0
        do it = 1, 100
            mmid = 0.5d0 * (mlo + mhi)
            nsum = 0d0
            do a = 1, nlev
                nsum = nsum + fermi_dirac(beta, eps(a), mmid)
            end do
            if (nsum < ntot) then
                mlo = mmid
            else
                mhi = mmid
            end if
        end do
        mu = 0.5d0 * (mlo + mhi)
        Eout = 0d0
        do a = 1, nlev
            ftgt(a) = fermi_dirac(beta, eps(a), mu)
            Eout = Eout + eps(a) * ftgt(a)
        end do
    end subroutine e_at_beta

    ! Fit (beta, mu) so sum_a f_FD = ntot and sum_a eps_a f_FD = etot, by nested
    ! bisection (outer beta log-scale, inner mu). ok=.false. if etot is outside
    ! [E(T=0), E(T=inf)] (not Fermi-Dirac representable, e.g. population inversion).
    subroutine fit_fermi_dirac(nlev, eps, ntot, etot, beta, mu, ftgt, ok)
        integer, intent(in)  :: nlev
        real(8), intent(in)  :: eps(nlev), ntot, etot
        real(8), intent(out) :: beta, mu, ftgt(nlev)
        logical, intent(out) :: ok
        real(8) :: blo, bhi, bmid, Elo, Ehi, Emid, mtmp
        integer :: it
        ok = .false.
        blo = 1d-4; bhi = 1d4
        call e_at_beta(nlev, eps, ntot, blo, Elo, mtmp, ftgt)
        call e_at_beta(nlev, eps, ntot, bhi, Ehi, mtmp, ftgt)
        if (etot > Elo + 1d-12 .or. etot < Ehi - 1d-12) return
        mu = mtmp; bmid = sqrt(blo * bhi)
        do it = 1, 100
            bmid = sqrt(blo * bhi)
            call e_at_beta(nlev, eps, ntot, bmid, Emid, mu, ftgt)
            if (Emid > etot) then
                blo = bmid
            else
                bhi = bmid
            end if
            if (abs(Emid - etot) < 1d-13 * max(1d0, abs(etot))) exit
        end do
        beta = bmid
        ok = .true.
    end subroutine fit_fermi_dirac

    ! Carrier-carrier CPTP relaxation map (Part F): relax the adiabatic
    ! populations toward the number- and energy-matched Fermi-Dirac and damp the
    ! coherences (EID):  rho -> (1-a) rho + a diag(occ f_FD), a = 1-exp(-nu tau),
    ! a convex combination of identity and a constant-state channel -> CPTP.
    ! Conserves Tr rho (number) and sum_a eps_a rho_aa (energy) exactly. A no-op
    ! when the population set is empty/full or not Fermi-Dirac representable.
    !
    ! TRACE-EXACT FORM: the FD target is fitted to the [0,1]-CLAMPED populations
    ! (the fit needs a representable set), so the update is written as a pure
    ! TRANSFER against that same clamped set,
    !     rho_aa += alpha * occ * (ftgt(a) - f_clamped(a)),   sum_a(...) = 0,
    ! which equals the convex mix above whenever the raw diagonal is in [0,occ].
    ! The old form mixed the RAW diagonal with the clamped-fit target; on a
    ! frozen active window whose Houston diagonal dips below 0 (PSD loss at the
    ! active<->frozen boundary) that mismatch CREATED trace at rate
    ! alpha*(occ*ntot_clamped - Tr rho) EVERY sub-step -- the monotone
    ! electrons 8.000 -> 8.012 drift seen with yn_sbe_eeh + frozen core.
    ! alpha_out reports the applied mixing weight (0 when the map was a no-op)
    ! so the caller can Kraus-extend the EID coherence damping, sqrt(1-alpha),
    ! to the active<->frozen coherence blocks it owns.
    subroutine carrier_carrier_relax(nlev, rho, eps, occ, nu, tau, alpha_out)
        integer,    intent(in)    :: nlev
        complex(8), intent(inout) :: rho(nlev, nlev)
        real(8),    intent(in)    :: eps(nlev), occ, nu, tau
        real(8),    intent(out)   :: alpha_out
        real(8) :: f(nlev), ftgt(nlev), ntot, etot, alpha, beta, mu
        integer :: a, b
        logical :: ok
        alpha_out = 0d0
        ntot = 0d0; etot = 0d0
        do a = 1, nlev
            f(a) = min(max(real(rho(a, a)) / occ, 0d0), 1d0)
            ntot = ntot + f(a)
            etot = etot + eps(a) * f(a)
        end do
        if (ntot < 1d-9 .or. (dble(nlev) - ntot) < 1d-9) return
        call fit_fermi_dirac(nlev, eps, ntot, etot, beta, mu, ftgt, ok)
        if (.not. ok) return
        alpha = 1d0 - exp(-nu * tau)
        alpha_out = alpha
        do b = 1, nlev
            do a = 1, nlev
                if (a == b) then
                    rho(a, a) = rho(a, a) + cmplx(alpha * occ * (ftgt(a) - f(a)), 0d0, 8)
                else
                    rho(a, b) = (1d0 - alpha) * rho(a, b)
                end if
            end do
        end do
    end subroutine carrier_carrier_relax

    ! =====================================================================
    ! INTER-K e-ph (the "super-mode ring" intervalley channel, Part C5/D).
    ! ---------------------------------------------------------------------
    ! The k-LOCAL apply_eph_relaxation relaxes a carrier to the energy-matched
    ! partner BAND at the SAME k -- correct only when the final valley folds onto
    ! a same-k band (folded cell) or for intra-valley polar-optical. On the
    ! PRIMITIVE cell the intervalley final states live at DIFFERENT k, so the
    ! search must run over all (k,band). This pure routine takes the GATHERED
    ! Houston spectrum and returns the net diagonal population change; the caller
    ! gathers eval/f (all-gather / ring) and applies dpop. Enabled when the ring
    ! (yn_sbe_superres) is on -- "if the ring is on, inter-k goes through it".
    !
    ! EXACTLY trace-conserving (CPTP-safe): each source (ik,a) loses
    !   out_tot = f*(1-exp(-Gamma_out*tau)) <= f   (no negativity),
    ! distributed to destinations proportional to their partial rate, so
    ! sum(dpop) = 0 by construction. Detailed balance (emission fe / absorption
    ! fa from the Bose factor) + Pauli blocking (1 - f_dest/occ_max) at the
    ! destination. Intra-k (jq==ik) is included automatically -- this SUBSUMES
    ! the intra-k channel when the ring is on. O(nk^2 nba^2 nph) all-pairs (same
    ! order as the Coulomb all-pairs sum it rides alongside).
    subroutine eph_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, evbm, &
                               nph, hw, wrel, nb_bose, kt_bath, db_realized, &
                               nu_sat, nu_eps0, nu_n, &
                               sigma, tau, dpop, gout, kidx, kn, pol_tab, pol_norm, ip_polar, &
                               ac_tab, ip_ac, ib_scale, e_src_lo, e_src_hi)
        implicit none
        integer, intent(in)  :: nk, nba, nph
        real(8), intent(in)  :: eval(nba, nk), f(nba, nk), occ_max, a2half
        real(8), intent(in)  :: ecbm, evbm, hw(nph), wrel(nph), nb_bose(nph)
        ! bath temperature [Ha]: the detailed-balance factor is evaluated at the
        ! REALIZED transfer, so it cannot be folded into nb_bose(hw_p) beforehand.
        ! nb_bose is still used for the mode banner and for the sigma -> 0 limit.
        real(8), intent(in)  :: kt_bath
        ! .true.  -- weight each transition by the detailed-balance factor of the
        !            energy it ACTUALLY transfers (eph_thermal_split_de) and use a
        !            pair-symmetric nu, so FD(kt_bath) is stationary for any sigma.
        ! .false. -- the historical behaviour: one split per MODE, from N_B(hw_p),
        !            and nu taken from the source alone. Correct only while
        !            sigma << hw_p; kept as the default so the gapped materials
        !            (Si, CdS, GaAs), whose optical modes satisfy that, are
        !            bit-identical. Set .true. for the gapless 2D materials, where
        !            the acoustic mode is 5 meV against a 0.1 eV search width.
        logical, intent(in)  :: db_realized
        real(8), intent(in)  :: nu_sat, nu_eps0, nu_n, sigma, tau   ! nu_n = saturation exponent
        real(8), intent(out) :: dpop(nba, nk)
        real(8), intent(out), optional :: gout(nba, nk)  ! total out-rate Gamma_out per source (coherence damping)
        ! A3 (optional): screened Frohlich 1/q^2 weight for the POLAR-LO mode
        ! ip_polar -- multiply its partial rates by pol_tab(q)/pol_norm (the B1
        ! signed-difference table; unit average preserves the nu_sat scale).
        integer, intent(in), optional :: kidx(3, nk), kn(3), ip_polar
        real(8), intent(in), optional :: pol_tab(*), pol_norm
        ! A4-CdS: Thomas-Fermi screening of the ACOUSTIC deformation mode --
        ! multiply mode ip_ac's partial rates by S(q) = [q/(q+q_TF)]^2 from
        ! build_acscreen_table (maintainer-critical for CdS E1 = 14.5 eV
        ! [Rode 1970]: the bare DP channel is unphysically strong at
        ! n >= 1e18 cm^-3; screening cuts the small-q part).
        real(8), intent(in), optional :: ac_tab(*)
        integer, intent(in), optional :: ip_ac
        ! Calibration knob (default 1): multiplies the partial rates of the
        ! INTERBAND pairs -- source and destination Houston branches on opposite
        ! sides of mid-gap. Those pairs carry the phonon-assisted band-to-band
        ! (BTBT / dressing-conversion) channel, whose model scale nu_sat is the
        ! cited REAL-carrier intervalley rate; the off-shell suppression of the
        ! eph matrix element between virtual states is not in the model, so the
        ! channel's ABSOLUTE rate is an upper estimate. Calibrated against the
        ! Keldysh/Hurkx brackets (samples/x12 rate_benchmark). Intraband pairs
        ! (real-carrier cooling / heating) are untouched.
        real(8), intent(in), optional :: ib_scale
        ! Cost-preserving full-basis mask: skip SOURCE dressed states whose energy is
        ! outside [e_src_lo, e_src_hi] (the frozen-window range). The dressed basis is
        ! still full (untruncated), only the scattering is restricted -- restores the
        ! O(nk^2 * n_active) cost without the truncation over-generation (wiki/06 sec.6).
        real(8), intent(in), optional :: e_src_lo, e_src_hi
        logical :: do_mask
        integer :: ik, jq, a, b, ip, ipol, ipac
        real(8) :: eps_kin, nu_a, nu_b, fe, fa, dE, shp, th, blk, gam, gamtot, out_tot
        real(8) :: delta, eps_kin_b
        real(8) :: gpart(nba, nk), emid, ibs
        logical :: src_cb
        real(8), parameter :: occ_eps = 1d-12

        ibs = 1d0
        if (present(ib_scale)) ibs = ib_scale
        do_mask = present(e_src_lo) .and. present(e_src_hi)
        emid = 0.5d0 * (ecbm + evbm)
        dpop = 0d0
        ipol = 0
        if (present(ip_polar)) then
            if (present(pol_tab) .and. present(pol_norm) .and. &
                present(kidx) .and. present(kn)) then
                if (ip_polar >= 1 .and. pol_norm > 0d0) ipol = ip_polar
            end if
        end if
        ipac = 0
        if (present(ip_ac)) then
            if (present(ac_tab) .and. present(kidx) .and. present(kn)) then
                if (ip_ac >= 1) ipac = ip_ac
            end if
        end if
        if (present(gout)) gout = 0d0
        ! Each (a, ik) source is independent: partial rates target shared sinks
        ! dpop(b, jq), hence the array reduction; gout(:, ik) is owner-written.
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp   private(ik, a, ip, jq, b, eps_kin, eps_kin_b, nu_a, nu_b, fe, fa, dE, &
        !$omp           delta, shp, th, blk, gam, gamtot, out_tot, gpart, src_cb) &
        !$omp   reduction(+:dpop)
        do ik = 1, nk
            do a = 1, nba
                if (f(a, ik) < occ_eps) cycle
                if (do_mask) then
                    if (eval(a, ik) < e_src_lo .or. eval(a, ik) > e_src_hi) cycle
                end if
                ! carrier kinetic energy from the nearest band edge (restore A^2/2,
                ! the k-independent Houston offset; it cancels in energy MATCHING).
                ! kinetic energy from the nearest band edge. a2half (the dropped
                ! Houston A^2/2) is NOT restored here: it is a GLOBAL scalar that
                ! cancels against the equally-shifted band edge (ecbm/evbm), and
                ! the field-heating is ALREADY carried by the Houston eigenvalue
                ! eval. Restoring it double-counted the field -> at strong drive
                ! a2half (tens of eV) swamped the real band energy and made every
                ! carrier "hot" (spurious e-ph/II/Auger). Fixed 2026-07-12.
                eps_kin = max(eval(a, ik) - ecbm, evbm - eval(a, ik), 0d0)
                nu_a = nu_saturation(eps_kin, nu_sat, nu_eps0, nu_n)
                if (nu_a * tau < 1d-14) cycle
                src_cb = (eval(a, ik) > emid)

                gpart = 0d0
                gamtot = 0d0
                do ip = 1, nph
                    if (.not. db_realized) call eph_thermal_split(nb_bose(ip), fe, fa)
                    do jq = 1, nk
                        do b = 1, nba
                            if (jq == ik .and. b == a) cycle
                            delta = abs(eval(b, jq) - eval(a, ik))
                            dE = abs(delta - hw(ip))
                            if (db_realized) then
                                shp = gaussian_shape(dE, sigma)
                                if (shp <= 0d0) cycle
                                ! detailed balance on the realized transfer: this is
                                ! what makes FD(kt_bath) stationary for every pair the
                                ! broadening admits, whatever sigma is.
                                call eph_thermal_split_de(delta, kt_bath, fe, fa)
                                th = merge(fe, fa, eval(b, jq) < eval(a, ik))
                                if (th <= 0d0) cycle
                                ! ... and a pair-symmetric nu, so the two directions of
                                ! one pair share a prefactor (nu(eps) alone differs by a
                                ! factor 2 across 0.4-0.8 eV at eps0 = 0.8 eV).
                                eps_kin_b = max(eval(b, jq) - ecbm, evbm - eval(b, jq), 0d0)
                                nu_b = nu_saturation(eps_kin_b, nu_sat, nu_eps0, nu_n)
                                nu_b = sqrt(nu_a * nu_b)
                            else
                                th = merge(fe, fa, eval(b, jq) < eval(a, ik))
                                if (th <= 0d0) cycle
                                shp = gaussian_shape(dE, sigma)
                                if (shp <= 0d0) cycle
                                nu_b = nu_a
                            end if
                            blk = min(max(1d0 - f(b, jq) / occ_max, 0d0), 1d0)
                            gam = nu_b * wrel(ip) * th * shp * blk
                            ! gap-straddling pair: the phonon-assisted BTBT /
                            ! dressing-conversion channel (see ib_scale above)
                            if (src_cb .neqv. (eval(b, jq) > emid)) gam = gam * ibs
                            ! A3: polar-LO Frohlich q-weight (screened 1/q^2)
                            if (ip == ipol) gam = gam * pol_tab(1 &
                                + (kidx(1,ik)-kidx(1,jq)+kn(1)-1) &
                                + (2*kn(1)-1)*((kidx(2,ik)-kidx(2,jq)+kn(2)-1) &
                                + (2*kn(2)-1)*(kidx(3,ik)-kidx(3,jq)+kn(3)-1))) / pol_norm
                            ! A4: TF-screened acoustic mode [q/(q+q_TF)]^2
                            if (ip == ipac) gam = gam * ac_tab(1 &
                                + (kidx(1,ik)-kidx(1,jq)+kn(1)-1) &
                                + (2*kn(1)-1)*((kidx(2,ik)-kidx(2,jq)+kn(2)-1) &
                                + (2*kn(2)-1)*(kidx(3,ik)-kidx(3,jq)+kn(3)-1)))
                            gpart(b, jq) = gpart(b, jq) + gam
                            gamtot = gamtot + gam
                        end do
                    end do
                end do

                if (present(gout)) gout(a, ik) = gamtot
                if (gamtot * tau < 1d-14) cycle
                out_tot = f(a, ik) * (1d0 - exp(-gamtot * tau))
                dpop(a, ik) = dpop(a, ik) - out_tot
                do jq = 1, nk
                    do b = 1, nba
                        if (gpart(b, jq) > 0d0) &
                            dpop(b, jq) = dpop(b, jq) + out_tot * gpart(b, jq) / gamtot
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine eph_interk_dpop

    ! =====================================================================
    ! INTER-K (momentum-conserving) impact ionization through the ring.
    ! ---------------------------------------------------------------------
    ! The TRUE 2-particle event: a hot conduction e- (k1, band ih) and a valence
    ! e- (k2, iv) -> two conduction e- at (k1', ic) and (k2', ic) leaving a hole,
    ! with crystal momentum k1 + k2 = k1' + k2' (mod G) [k2' from mp_partner_triple]
    ! and energy E(k1,ih)+E(k2,iv) = E(k1',ic)+E(k2',ic) [broadened Fermi golden
    ! rule]. The threshold magnitude g0 = pref*(eps_kin - E_th)^expo is the cited
    ! Stobbe-fit rate; the screened Coulomb |V(q)|^2 (q = k1-k1') and the energy
    ! delta SHAPE which momentum-conserving final config it goes to. EXACTLY
    ! trace-conserving: each event writes -amt,-amt,+amt,+amt (primary out +
    ! valence out = the two conduction gains), so sum(dpop)=0 by construction.
    ! Primary out capped at f(ih,k1)*(1-exp(-Gamma*tau)) <= f (no negativity).
    ! klut(0:nk-1): flattened MP lookup, triple (m1,m2,m3) -> ik via
    ! lidx = m1 + n1*(m2 + n2*m3). The caller precomputes it once.
    ! i1_lo/i1_hi (optional): restrict the OUTER source loop to a k subrange.
    ! The kernel is exactly additive over i1 (each i1 only accumulates its own
    ! -amt/+amt quadruples), so ranks can each run a disjoint subrange and
    ! comm_summation the dpop -- O(nk^3/P) instead of every rank redoing the
    ! full O(nk^3) sum. Omitting them keeps the full serial loop (unit tests).
    subroutine ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                              pref, expo, iv, ic, kidx, kn, klut, &
                              bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dpop, &
                              i1_lo, i1_hi, opts, e_src_lo, e_src_hi)
        implicit none
        integer, intent(in)  :: nk, nba, iv, ic, kidx(3, nk), kn(3), klut(0:nk-1)
        real(8), intent(in)  :: eval(nba, nk), f(nba, nk), occ_max, a2half
        real(8), intent(in)  :: ecbm, eth, pref, expo, sigma, tau
        real(8), intent(in)  :: bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
        real(8), intent(out) :: dpop(nba, nk)
        integer, intent(in), optional :: i1_lo, i1_hi
        type(t_ring_opts), intent(in), optional :: opts
        ! Cost-preserving full-basis mask: skip a hot SOURCE (electron ih in
        ! [ic,nba] or deep hole ih in [1,iv]) whose Houston energy eval(ih,i1)
        ! is outside [e_src_lo, e_src_hi] (the frozen active window). The dressed
        ! basis stays full; only the O(nk^2) 2-particle scattering is restricted
        ! to in-window sources -- restores cost without the truncation
        ! over-generation. No-op when the pair is absent (default all-active).
        real(8), intent(in), optional :: e_src_lo, e_src_hi
        logical :: do_mask
        integer :: i1, i1p, i2, ih, ipass, jj, m2p(3), d, i1s, i1e, nphl
        real(8) :: ekin, dd, g0, etgt, vq, df(3), shp, pauli, gpart, gamtot, out_tot, amt
        logical :: havetab
        real(8) :: vfloor, fk, pa, prefh, evbm, room
        real(8), allocatable :: hwl(:), nbbl(:), wrell(:)
        real(8), parameter :: occ_eps = 1d-12

        dpop = 0d0
        if (iv < 1 .or. ic > nba .or. ic <= iv) return
        i1s = 1;  if (present(i1_lo)) i1s = max(i1_lo, 1)
        i1e = nk; if (present(i1_hi)) i1e = min(i1_hi, nk)
        do_mask = present(e_src_lo) .and. present(e_src_hi)
        call ring_opts_unpack(opts, havetab, vfloor, fk, pa, nphl, prefh, evbm, &
                              hwl, nbbl, wrell)

        ! Independent (i1, ih) hot sources; shared sinks -> array reduction.
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp   private(i1, ih, ipass, i1p, i2, jj, m2p, d, ekin, dd, g0, etgt, &
        !$omp           vq, df, shp, pauli, gpart, gamtot, out_tot, amt) &
        !$omp   reduction(+:dpop)
        do i1 = i1s, i1e
            do ih = ic, nba
                if (f(ih, i1) < occ_eps) cycle
                if (do_mask) then
                    if (eval(ih,i1) < e_src_lo .or. eval(ih,i1) > e_src_hi) cycle
                end if
                ekin = eval(ih, i1) - ecbm   ! a2half NOT restored (cancels vs shifted CBM; see e-ph note)
                dd = ekin - eth
                ! A5: Franz-Keldysh field softening -- softplus with the
                ! electro-optic width (fk = hbar*theta): -> max(dd,0) as fk -> 0.
                if (fk > 0d0) dd = fk * log(1d0 + exp(min(dd / fk, 4d1)))
                if (dd <= 0d0) cycle
                g0 = pref * dd ** expo            ! cited Stobbe-fit total magnitude

                ! Pass 1 (ipass=1): accumulate gamtot; Pass 2: distribute out_tot.
                gamtot = 0d0
                out_tot = 0d0
                do ipass = 1, 2
                    if (ipass == 2) then
                        if (gamtot * tau < 1d-14) exit
                        out_tot = f(ih, i1) * (1d0 - exp(-gamtot * tau))
                    end if
                    do i1p = 1, nk
                        ! transferred momentum q = k1 - k1' (+G): CDRB-screened
                        ! Cartesian-metric weight with the 27-image umklapp sum
                        ! (B1: table lookup over the signed index difference)
                        if (havetab) then
                            vq = opts%vq_tab(1 + (kidx(1,i1)-kidx(1,i1p)+kn(1)-1) &
                               + (2*kn(1)-1)*((kidx(2,i1)-kidx(2,i1p)+kn(2)-1) &
                               + (2*kn(2)-1)*(kidx(3,i1)-kidx(3,i1p)+kn(3)-1)))
                        else
                            do d = 1, 3
                                df(d) = dble(kidx(d, i1) - kidx(d, i1p)) / dble(max(kn(d), 1))
                            end do
                            vq = interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg)
                        end if
                        if (vq < vfloor) cycle          ! B3 windowing (0 = off)
                        do i2 = 1, nk
                            call mp_partner_triple(kidx(:,i1), kidx(:,i2), kidx(:,i1p), kn, m2p)
                            jj = klut(m2p(1) + kn(1) * (m2p(2) + kn(2) * m2p(3))) ! O(1) k2'
                            if (jj < 1) cycle
                            etgt = eval(ih,i1) + eval(iv,i2) - eval(ic,i1p) - eval(ic,jj)
                            shp = shape_assist(etgt, sigma, pa, nphl, hwl, nbbl, wrell, .false.)
                            if (shp <= 0d0) cycle
                            pauli = (f(iv,i2) / occ_max) &
                                  * min(max(1d0 - f(ic,i1p)/occ_max, 0d0), 1d0) &
                                  * min(max(1d0 - f(ic,jj )/occ_max, 0d0), 1d0)
                            gpart = g0 * vq * shp * pauli
                            if (ipass == 1) then
                                gamtot = gamtot + gpart
                            else
                                amt = out_tot * gpart / gamtot
                                dpop(ih, i1)  = dpop(ih, i1)  - amt   ! hot primary leaves
                                dpop(ic, i1p) = dpop(ic, i1p) + amt   ! primary relaxed
                                dpop(iv, i2)  = dpop(iv, i2)  - amt   ! valence e- -> hole
                                dpop(ic, jj)  = dpop(ic, jj)  + amt   ! promoted to conduction
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! ================= A2: HOLE-INITIATED impact ionization (hhe) =========
        ! A hot HOLE deep in the valence band relaxes toward the VBM and the
        ! released energy ionizes a pair. Electron picture: an electron falls
        ! from the top valence (iv, k1') into the deep hole (ih, k1); a valence
        ! electron (iv, k2) is promoted to the conduction band (ic, k2').
        ! Momentum: k1' + k2 = k1 + k2' (mod G)  ->  partner(i1p, i2, i1).
        ! Rate scale = pref * (Cp/Cn) (registry-cited [L90]/[S14] ratio),
        ! threshold on the HOLE kinetic energy evbm - E. Trace-exact stencil.
        if (prefh > 0d0) then
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp   private(i1, ih, ipass, i1p, i2, jj, m2p, d, room, ekin, dd, g0, &
        !$omp           etgt, vq, df, shp, pauli, gpart, gamtot, out_tot, amt) &
        !$omp   reduction(+:dpop)
        do i1 = i1s, i1e
            do ih = 1, iv
                room = occ_max - f(ih, i1)          ! deep-hole capacity
                if (room < occ_eps) cycle
                if (do_mask) then
                    if (eval(ih,i1) < e_src_lo .or. eval(ih,i1) > e_src_hi) cycle
                end if
                ekin = evbm - eval(ih, i1)   ! a2half NOT restored (cancels vs shifted VBM)
                dd = ekin - eth
                if (fk > 0d0) dd = fk * log(1d0 + exp(min(dd / fk, 4d1)))
                if (dd <= 0d0) cycle
                g0 = prefh * dd ** expo
                gamtot = 0d0
                out_tot = 0d0
                do ipass = 1, 2
                    if (ipass == 2) then
                        if (gamtot * tau < 1d-14) exit
                        out_tot = room * (1d0 - exp(-gamtot * tau))
                    end if
                    do i1p = 1, nk
                        if (havetab) then
                            vq = opts%vq_tab(1 + (kidx(1,i1)-kidx(1,i1p)+kn(1)-1) &
                               + (2*kn(1)-1)*((kidx(2,i1)-kidx(2,i1p)+kn(2)-1) &
                               + (2*kn(2)-1)*(kidx(3,i1)-kidx(3,i1p)+kn(3)-1)))
                        else
                            do d = 1, 3
                                df(d) = dble(kidx(d, i1) - kidx(d, i1p)) / dble(max(kn(d), 1))
                            end do
                            vq = interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg)
                        end if
                        if (vq < vfloor) cycle
                        do i2 = 1, nk
                            call mp_partner_triple(kidx(:,i1p), kidx(:,i2), kidx(:,i1), kn, m2p)
                            jj = klut(m2p(1) + kn(1) * (m2p(2) + kn(2) * m2p(3)))
                            if (jj < 1) cycle
                            etgt = eval(iv,i1p) + eval(iv,i2) - eval(ih,i1) - eval(ic,jj)
                            shp = shape_assist(etgt, sigma, pa, nphl, hwl, nbbl, wrell, .false.)
                            if (shp <= 0d0) cycle
                            pauli = (f(iv,i1p) / occ_max) * (f(iv,i2) / occ_max) &
                                  * min(max(1d0 - f(ic,jj)/occ_max, 0d0), 1d0)
                            gpart = g0 * vq * shp * pauli
                            if (ipass == 1) then
                                gamtot = gamtot + gpart
                            else
                                amt = out_tot * gpart / gamtot
                                dpop(ih, i1)  = dpop(ih, i1)  + amt   ! deep hole filled
                                dpop(iv, i1p) = dpop(iv, i1p) - amt   ! hole surfaces at the VBM
                                dpop(iv, i2)  = dpop(iv, i2)  - amt   ! pair: hole created
                                dpop(ic, jj)  = dpop(ic, jj)  + amt   ! pair: electron created
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end if
    end subroutine ii_interk_dpop

    ! =====================================================================
    ! INTER-K (momentum-conserving) AUGER RECOMBINATION through the ring --
    ! the EXACT TIME-REVERSE of ii_interk_dpop (detailed balance).
    ! ---------------------------------------------------------------------
    ! The reverse 2-particle event: two conduction e- at (k1',ic) and (k2',ic)
    ! + a hole at (k2,iv) -> one e- recombines into the hole and the released
    ! energy promotes the other to the hot state (k1,ih). SAME quadruples,
    ! SAME weights g0*|V(q)|^2*delta_sigma as the impact-ionization kernel (no
    ! new constant: the rate scale IS the cited Stobbe/Keldysh II magnitude --
    ! Auger and II share |M|^2, only the occupation factors swap [Rana 2007;
    ! Kioupakis 2015]), REVERSED occupation product
    !   pauli_rev = (1 - f_v/occ) * (f_c1'/occ) * (f_c2'/occ)
    ! and NEGATED dpop signs: +hot, -c1', +valence, -c2'. Trace-conserving by
    ! construction (sum(dpop)=0). The hot-state gain is capped at
    ! (occ - f(ih,k1))*(1-exp(-Gamma*tau)) (no over-filling), mirroring the
    ! II primary-out cap. DETAILED BALANCE: for Fermi-Dirac occupations and an
    ! energy-conserving quadruple, f1 f2 (1-f3)(1-f4) = (1-f1)(1-f2) f3 f4, so
    ! in the linear (Gamma*tau -> 0) regime the net II + Auger dpop vanishes
    ! identically -- the equilibrium-fixed-point unit test.
    ! =====================================================================
    ! i1_lo/i1_hi: same optional outer-loop subrange as ii_interk_dpop (the
    ! kernel is additive over i1 -- MPI ranks sum disjoint subranges).
    subroutine auger_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                                 pref, expo, iv, ic, kidx, kn, klut, &
                                 bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dpop, &
                                 i1_lo, i1_hi, opts, e_src_lo, e_src_hi)
        implicit none
        integer, intent(in)  :: nk, nba, iv, ic, kidx(3, nk), kn(3), klut(0:nk-1)
        real(8), intent(in)  :: eval(nba, nk), f(nba, nk), occ_max, a2half
        real(8), intent(in)  :: ecbm, eth, pref, expo, sigma, tau
        real(8), intent(in)  :: bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
        real(8), intent(out) :: dpop(nba, nk)
        integer, intent(in), optional :: i1_lo, i1_hi
        type(t_ring_opts), intent(in), optional :: opts
        ! Cost-preserving full-basis mask (see ii_interk_dpop): skip the hot
        ! recombination SOURCE -- the promoted-to state ih in [ic,nba] (forward
        ! Auger) or the deep valence electron ih in [1,iv] (reverse hole-Auger)
        ! -- when eval(ih,i1) is outside the frozen active window. No-op when
        ! absent (default all-active).
        real(8), intent(in), optional :: e_src_lo, e_src_hi
        logical :: do_mask
        integer :: i1, i1p, i2, ih, ipass, jj, m2p(3), d, i1s, i1e, nphl
        real(8) :: ekin, dd, g0, etgt, vq, df(3), shp, pauli, gpart, gamtot, in_tot, amt, room
        logical :: havetab
        real(8) :: vfloor, fk, pa, prefh, evbm, out_tot
        real(8), allocatable :: hwl(:), nbbl(:), wrell(:)
        real(8), parameter :: occ_eps = 1d-12

        dpop = 0d0
        if (iv < 1 .or. ic > nba .or. ic <= iv) return
        ! empty conduction band -> Auger needs two occupied CB sources; no-op.
        if (maxval(f(ic:nba, :)) < occ_eps) return
        i1s = 1;  if (present(i1_lo)) i1s = max(i1_lo, 1)
        i1e = nk; if (present(i1_hi)) i1e = min(i1_hi, nk)
        do_mask = present(e_src_lo) .and. present(e_src_hi)
        call ring_opts_unpack(opts, havetab, vfloor, fk, pa, nphl, prefh, evbm, &
                              hwl, nbbl, wrell)

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp   private(i1, ih, ipass, i1p, i2, jj, m2p, d, room, ekin, dd, g0, &
        !$omp           etgt, vq, df, shp, pauli, gpart, gamtot, in_tot, amt) &
        !$omp   reduction(+:dpop)
        do i1 = i1s, i1e
            do ih = ic, nba
                room = occ_max - f(ih, i1)          ! empty hot-state phase space
                if (room < occ_eps) cycle
                if (do_mask) then
                    if (eval(ih,i1) < e_src_lo .or. eval(ih,i1) > e_src_hi) cycle
                end if
                ekin = eval(ih, i1) - ecbm   ! a2half NOT restored (cancels vs shifted CBM; see e-ph note)
                dd = ekin - eth
                if (fk > 0d0) dd = fk * log(1d0 + exp(min(dd / fk, 4d1)))
                if (dd <= 0d0) cycle
                g0 = pref * dd ** expo              ! same II magnitude (shared |M|^2)

                gamtot = 0d0
                in_tot = 0d0
                do ipass = 1, 2
                    if (ipass == 2) then
                        if (gamtot * tau < 1d-14) exit
                        in_tot = room * (1d0 - exp(-gamtot * tau))
                    end if
                    do i1p = 1, nk
                        ! SAME CDRB-screened umklapp weight as the II kernel
                        ! (shared |M|^2 -> detailed balance preserved exactly)
                        if (havetab) then
                            vq = opts%vq_tab(1 + (kidx(1,i1)-kidx(1,i1p)+kn(1)-1) &
                               + (2*kn(1)-1)*((kidx(2,i1)-kidx(2,i1p)+kn(2)-1) &
                               + (2*kn(2)-1)*(kidx(3,i1)-kidx(3,i1p)+kn(3)-1)))
                        else
                            do d = 1, 3
                                df(d) = dble(kidx(d, i1) - kidx(d, i1p)) / dble(max(kn(d), 1))
                            end do
                            vq = interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg)
                        end if
                        if (vq < vfloor) cycle          ! B3 (0 = off)
                        do i2 = 1, nk
                            call mp_partner_triple(kidx(:,i1), kidx(:,i2), kidx(:,i1p), kn, m2p)
                            jj = klut(m2p(1) + kn(1) * (m2p(2) + kn(2) * m2p(3)))
                            if (jj < 1) cycle
                            etgt = eval(ih,i1) + eval(iv,i2) - eval(ic,i1p) - eval(ic,jj)
                            shp = shape_assist(etgt, sigma, pa, nphl, hwl, nbbl, wrell, .true.)
                            if (shp <= 0d0) cycle
                            ! REVERSED occupations: hole present at (iv,i2), both
                            ! conduction sources occupied at (ic,i1p) and (ic,jj).
                            pauli = min(max(1d0 - f(iv,i2)/occ_max, 0d0), 1d0) &
                                  * (f(ic,i1p) / occ_max) &
                                  * (f(ic,jj ) / occ_max)
                            gpart = g0 * vq * shp * pauli
                            if (ipass == 1) then
                                gamtot = gamtot + gpart
                            else
                                amt = in_tot * gpart / gamtot
                                dpop(ih, i1)  = dpop(ih, i1)  + amt   ! promoted to hot
                                dpop(ic, i1p) = dpop(ic, i1p) - amt   ! promoted-source leaves
                                dpop(iv, i2)  = dpop(iv, i2)  + amt   ! hole filled (recombination)
                                dpop(ic, jj)  = dpop(ic, jj)  - amt   ! recombining e- leaves
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! ============ A2 reverse: hole-Auger (the hhe time-reverse) ===========
        ! An electron IN the deep valence state (ih, k1) is excited to the top
        ! valence (iv, k1') by the recombination of e(ic, k2') with h(iv, k2) --
        ! reversed occupations, negated stencil, same quadruple/weights as the
        ! forward hole channel (detailed balance per sideband via rev=.true.).
        if (prefh > 0d0) then
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp   private(i1, ih, ipass, i1p, i2, jj, m2p, d, ekin, dd, g0, etgt, &
        !$omp           vq, df, shp, pauli, gpart, gamtot, out_tot, amt) &
        !$omp   reduction(+:dpop)
        do i1 = i1s, i1e
            do ih = 1, iv
                if (f(ih, i1) < occ_eps) cycle      ! deep electron present
                if (do_mask) then
                    if (eval(ih,i1) < e_src_lo .or. eval(ih,i1) > e_src_hi) cycle
                end if
                ekin = evbm - eval(ih, i1)   ! a2half NOT restored (cancels vs shifted VBM)
                dd = ekin - eth
                if (fk > 0d0) dd = fk * log(1d0 + exp(min(dd / fk, 4d1)))
                if (dd <= 0d0) cycle
                g0 = prefh * dd ** expo
                gamtot = 0d0
                out_tot = 0d0
                do ipass = 1, 2
                    if (ipass == 2) then
                        if (gamtot * tau < 1d-14) exit
                        out_tot = f(ih, i1) * (1d0 - exp(-gamtot * tau))
                    end if
                    do i1p = 1, nk
                        if (havetab) then
                            vq = opts%vq_tab(1 + (kidx(1,i1)-kidx(1,i1p)+kn(1)-1) &
                               + (2*kn(1)-1)*((kidx(2,i1)-kidx(2,i1p)+kn(2)-1) &
                               + (2*kn(2)-1)*(kidx(3,i1)-kidx(3,i1p)+kn(3)-1)))
                        else
                            do d = 1, 3
                                df(d) = dble(kidx(d, i1) - kidx(d, i1p)) / dble(max(kn(d), 1))
                            end do
                            vq = interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg)
                        end if
                        if (vq < vfloor) cycle
                        do i2 = 1, nk
                            call mp_partner_triple(kidx(:,i1p), kidx(:,i2), kidx(:,i1), kn, m2p)
                            jj = klut(m2p(1) + kn(1) * (m2p(2) + kn(2) * m2p(3)))
                            if (jj < 1) cycle
                            etgt = eval(iv,i1p) + eval(iv,i2) - eval(ih,i1) - eval(ic,jj)
                            shp = shape_assist(etgt, sigma, pa, nphl, hwl, nbbl, wrell, .true.)
                            if (shp <= 0d0) cycle
                            pauli = min(max(1d0 - f(iv,i1p)/occ_max, 0d0), 1d0) &
                                  * min(max(1d0 - f(iv,i2 )/occ_max, 0d0), 1d0) &
                                  * (f(ic,jj) / occ_max)
                            gpart = g0 * vq * shp * pauli
                            if (ipass == 1) then
                                gamtot = gamtot + gpart
                            else
                                amt = out_tot * gpart / gamtot
                                dpop(ih, i1)  = dpop(ih, i1)  - amt   ! deep e- excited away
                                dpop(iv, i1p) = dpop(iv, i1p) + amt   ! lands at the VBM
                                dpop(iv, i2)  = dpop(iv, i2)  + amt   ! hole filled (recomb.)
                                dpop(ic, jj)  = dpop(ic, jj)  - amt   ! recombining e- leaves
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        end if
    end subroutine auger_interk_dpop

    ! =====================================================================
    ! Monkhorst-Pack momentum-conservation index map (for the inter-k / nonlocal
    ! collision channels: impact ionization, e-e). A regular MP mesh has points
    ! k_m = (2m - n + 1)/(2n) per dimension (m = 0..n-1), so crystal-momentum
    ! conservation k1 + k2 - k1' (mod G) is EXACT integer arithmetic on the index
    ! triples: m2' = mod(m1 + m2 - m1', n). mp_grid_triple recovers the triple
    ! from the reduced coordinate and reports the residual (0 for a true MP point,
    ! so the caller can gate the channel off on a non-MP / symmetry-reduced grid).
    ! =====================================================================
    pure subroutine mp_grid_triple(kred, n, m, resid)
        real(8), intent(in)  :: kred(3)
        integer, intent(in)  :: n(3)
        integer, intent(out) :: m(3)
        real(8), intent(out) :: resid     ! max wrapped |kred - reconstructed MP point|
        integer :: d
        real(8) :: km
        resid = 0d0
        do d = 1, 3
            if (n(d) <= 0) then
                m(d) = 0
                cycle
            end if
            m(d) = modulo(nint(kred(d) * n(d) + (n(d) - 1) * 0.5d0), n(d))
            km = (2d0 * m(d) - n(d) + 1d0) / (2d0 * n(d))
            resid = max(resid, abs(modulo(kred(d) - km + 0.5d0, 1d0) - 0.5d0))
        end do
    end subroutine mp_grid_triple

    ! Momentum-conserving partner index triple: m2' = mod(m1 + m2 - m1', n).
    pure subroutine mp_partner_triple(m1, m2, m1p, n, m2p)
        integer, intent(in)  :: m1(3), m2(3), m1p(3), n(3)
        integer, intent(out) :: m2p(3)
        integer :: d
        do d = 1, 3
            m2p(d) = modulo(m1(d) + m2(d) - m1p(d), max(n(d), 1))
        end do
    end subroutine mp_partner_triple

    ! =====================================================================
    ! Velocity-gauge basis sufficiency / N_b convergence primitives.
    ! These quantify whether the band count N_b carried into the dynamics is
    ! large enough -- a correctness axis SEPARATE from the plane-wave cutoff,
    ! and NOT fixed by rotating into the Houston basis (which only diagonalizes
    ! the already-truncated H_VG^(N_b) = P_Nb H_VG P_Nb). See the wiki page
    ! "VG Basis Sufficiency & N_b Convergence".
    ! [Hylleraas-Undheim Z.Phys.65,759(1930); MacDonald PR 43,830(1933);
    !  Wismer-Yakovlev PRB 97,144302(2018)]
    ! =====================================================================

    ! Criterion (c): dimensionless adiabatic admixture of band c into band a
    !   eta_ac = A_max |pi_ac| / |eps_a - eps_c|.
    ! The basis is safe for level a when eta to the FIRST DISCARDED band << 1;
    ! eta >~ 1 means the admixture is nonperturbative and c cannot be dropped.
    pure function vg_eta_admixture(A_max, pi_ac, gap_ac) result(eta)
        real(8), intent(in) :: A_max, pi_ac, gap_ac
        real(8) :: eta
        if (abs(gap_ac) <= 0d0) then
            eta = huge(1d0)
        else
            eta = A_max * abs(pi_ac) / abs(gap_ac)
        end if
    end function vg_eta_admixture

    ! Second-order (Rayleigh-Schrodinger) truncation shift of a retained level a
    ! from the n DISCARDED bands c > N_b:
    !   delta eps_a = sum_c A^2 |pi_ac|^2 / (eps_a - eps_c).
    ! This is exactly the error the Houston basis inherits (it diagonalizes the
    ! projected operator, so this shift is NOT recovered). Valid for eta << 1.
    pure function vg_trunc_shift2(A, pi_ac, eps_a, eps_disc, n) result(dshift)
        integer, intent(in) :: n
        real(8), intent(in) :: A, pi_ac(n), eps_a, eps_disc(n)
        real(8) :: dshift
        integer :: c
        real(8) :: denom
        dshift = 0d0
        do c = 1, n
            denom = eps_a - eps_disc(c)
            if (abs(denom) > 0d0) dshift = dshift + A*A*pi_ac(c)*pi_ac(c) / denom
        end do
    end function vg_trunc_shift2

    ! Criterion (b): relative L2 convergence error between an observable computed
    ! with N_b and with N_b+Delta bands (current, HHG spectrum, carrier number):
    !   eps_conv = ||O_large - O_small|| / ||O_large||.
    pure function vg_conv_error(o_small, o_large, n) result(eps_conv)
        integer, intent(in) :: n
        real(8), intent(in) :: o_small(n), o_large(n)
        real(8) :: eps_conv, num, den
        integer :: i
        num = 0d0; den = 0d0
        do i = 1, n
            num = num + (o_large(i) - o_small(i))**2
            den = den + o_large(i)**2
        end do
        if (den <= 0d0) then
            eps_conv = 0d0
        else
            eps_conv = sqrt(num / den)
        end if
    end function vg_conv_error

    ! Criterion (a): does the top retained adiabatic band's peak occupation
    ! P_top = max_k rho~_{Nb,Nb}(k) exceed the tolerance (default 1e-3)? A .true.
    ! means population reached the basis edge -> enlarge N_b.
    pure function vg_ptop_exceeds(ptop, thr) result(over)
        real(8), intent(in) :: ptop, thr
        logical :: over
        over = (ptop > thr)
    end function vg_ptop_exceeds

    ! =====================================================================
    ! Collisional-memory (non-Markovian) dephasing lines.
    ! [MT99] Meier & Tannor, JCP 111, 3365 (1999); [B25] RPP 88, 070501
    ! (2025); maintainer decision wiki/10 sec. 8.6: the production
    ! decoherence is COLLISIONAL -- the Markovian exp(-gout*tau/2) factor
    ! of the e-ph Lindblad is the zero-memory limit of the channel's own
    ! bath correlation kernel. The kernel is built VERBATIM from the cited
    ! phonon table (no new free parameters):
    !
    !   phi(tau) = (1/A) sum_p wrel_p [ (N_p+1) e^{-i w_p tau}
    !                                 +  N_p    e^{+i w_p tau} ] e^{-tau/tau_c}
    !
    ! = Lorentzian lines [MT99 Eq. (15)] at the cited mode energies w_p,
    ! emission/absorption weights (N_p+1)/N_p (detailed balance), common
    ! width 1/tau_c (default sigma_E -- the time-domain completion of the
    ! sigma_E-broadened golden rule already in use). The anchor
    ! A = Re sum_j c_j/mu_j normalizes the zero-frequency response to 1:
    ! a slow (adiabatically following) Houston coherence is damped at
    ! EXACTLY the channel's Markovian rate; only the response to
    ! sub-correlation-time modulation (the field-driven dressing) changes
    ! -- the bath cannot follow it [B25 Fig 5(b)].
    ! Lines: j=(p,+): c = wrel_p (N_p+1)/A, mu = 1/tau_c + i w_p
    !        j=(p,-): c = wrel_p  N_p   /A, mu = 1/tau_c - i w_p
    ! Frequency response (tested): R(w) = Re sum_j c_j/(mu_j + i w);
    ! R(0) = 1 exactly, R -> 0 for |w| >> w_p + 1/tau_c.
    ! =====================================================================

    subroutine colmem_lines(nph, hw, wrel, nbose, tauc, nl, cl, mul)
        integer, intent(in) :: nph
        real(8), intent(in) :: hw(nph), wrel(nph), nbose(nph), tauc
        integer, intent(out) :: nl
        complex(8), intent(out) :: cl(2*nph), mul(2*nph)
        real(8) :: anchor
        integer :: p, j

        nl = 0
        do p = 1, nph
            if (wrel(p) <= 0d0) cycle
            nl = nl + 1
            cl(nl)  = cmplx(wrel(p) * (nbose(p) + 1d0), 0d0, 8)
            mul(nl) = cmplx(1d0 / tauc, hw(p), 8)
            if (nbose(p) > 0d0) then
                nl = nl + 1
                cl(nl)  = cmplx(wrel(p) * nbose(p), 0d0, 8)
                mul(nl) = cmplx(1d0 / tauc, -hw(p), 8)
            end if
        end do
        if (nl == 0) return
        anchor = 0d0
        do j = 1, nl
            anchor = anchor + real(cl(j) / mul(j))
        end do
        if (anchor <= 0d0) error stop 'colmem_lines: non-positive Markov anchor'
        do j = 1, nl
            cl(j) = cl(j) / anchor
        end do
    end subroutine colmem_lines

    ! Population-sector memory filter (wiki/10 sec. 8.8): one discrete step of
    ! the line convolution on a POPULATION f, returning the memory-filtered
    ! ftil that feeds the ring collision kernels. Auxiliary fields z_j evolve
    ! as z <- z e^{-mu tau} + f tau; ftil = Re sum_j c_j z_j / gnorm with the
    ! DISCRETE normalization gnorm = Re sum_j c_j tau/(1-e^{-mu_j tau}), so a
    ! CONSTANT f is a machine-exact fixed point (ftil = f) at ANY tau -- the
    ! calibrated Markovian rates are untouched for slow populations. Fast
    ! (sub-correlation-time) modulation -- the A^2(t) dressing breathing --
    ! filters out of the collision SOURCE: the time-domain ICFE statement
    ! that the phonon bath cannot follow sub-cycle population modulation.
    ! Initialize z with colmem_pop_init so ftil(0) = f(0) exactly.
    pure subroutine colmem_pop_filter(nl, cl, mul, tau, f, z, ftil)
        integer, intent(in) :: nl
        complex(8), intent(in) :: cl(nl), mul(nl)
        real(8), intent(in) :: tau, f
        complex(8), intent(inout) :: z(nl)
        real(8), intent(out) :: ftil
        real(8) :: gnorm
        integer :: j
        ftil = 0d0
        gnorm = 0d0
        do j = 1, nl
            z(j) = z(j) * exp(-mul(j) * tau) + f * tau
            ftil = ftil + real(cl(j) * z(j))
            gnorm = gnorm + real(cl(j) * tau / (1d0 - exp(-mul(j) * tau)))
        end do
        ftil = ftil / max(gnorm, 1d-30)
    end subroutine colmem_pop_filter

    ! Fixed-point initialization of the filter fields: z_j = f tau/(1-e^{-mu_j tau})
    ! makes the very first filtered value equal f exactly (and a constant f
    ! stays a fixed point of colmem_pop_filter thereafter).
    pure subroutine colmem_pop_init(nl, mul, tau, f, z)
        integer, intent(in) :: nl
        complex(8), intent(in) :: mul(nl)
        real(8), intent(in) :: tau, f
        complex(8), intent(out) :: z(nl)
        integer :: j
        do j = 1, nl
            z(j) = f * tau / (1d0 - exp(-mul(j) * tau))
        end do
    end subroutine colmem_pop_init

    ! =====================================================================
    ! Option A (wiki/10 sec. 3A): DRESSED-REFERENCE carrier measure.
    ! The ring dissipators already act in the instantaneous dressed basis
    ! (H_VG + A.p (+ Sigma^HF) is diagonalized every step -- the ring-basis
    ! audit, wiki/10 sec. 8.10), but they MEASURED carriers against the
    ! static {occ, 0} reference. The field-rotated GROUND STATE then looks
    ! like carriers: its Houston-frame populations f0_a = occ sum_{v<=nv}
    ! |W(v,a)|^2 breathe with A^2 -- and the wiki/00 diagnostic shows the
    ! actual sub-cycle state tracks exactly this rotation background
    ! (corr 0.99, 99.9% return). The dressed reference subtracts it:
    !
    !   delta0_a = occ sum_{v<=nv} |W(v,a)|^2 - occ [a <= nv],
    !   f_eff    = f - delta0
    !
    ! Properties (tested): sum_a delta0_a = 0 exactly (unitarity -> the
    ! carrier measure is trace-neutral); A -> 0 => W -> 1 => delta0 = 0
    ! (post-pulse carriers counted exactly); a FROZEN (non-reacting) state
    ! f = f0 gives f_eff = {occ,0} -> the channels see equilibrium, zero
    ! fabrication; an ADIABATIC state f = {occ,0} gives negative upper
    ! excess -> clamped, also protected. Only genuine (non-adiabatic,
    ! collision-created or multiphoton) population survives in f_eff.
    ! Composes with the collisional-memory filters: subtract BEFORE the
    ! gather/gate/filter, so the linear pop filter smooths f - delta0.
    ! =====================================================================
    ! f0_in (optional): the actual INITIAL diagonal occupation of the active
    ! bands at this k. Present for a doped / finite-temperature sheet, where the
    ! reference is no longer "nv full bands, the rest empty":
    !     delta_a = sum_b f0_b |W(b,a)|^2 - f0_a
    ! i.e. the adiabatically rotated initial state minus the initial state --
    ! still trace-neutral (sum_a delta_a = 0 by unitarity) and still zero at
    ! A -> 0. Without it the integer-filling form is used (identical result when
    ! f0 = {occ, ..., occ, 0, ..., 0}).
    pure subroutine dressed_ref_delta(nba, nv, occ, W, delta, f0_in)
        integer, intent(in) :: nba, nv
        real(8), intent(in) :: occ
        complex(8), intent(in) :: W(nba, nba)
        real(8), intent(out) :: delta(nba)
        real(8), intent(in), optional :: f0_in(nba)
        integer :: a, v, b
        real(8) :: f0
        if (present(f0_in)) then
            do a = 1, nba
                f0 = 0d0
                do b = 1, nba
                    f0 = f0 + f0_in(b) * abs(W(b, a))**2
                end do
                delta(a) = f0 - f0_in(a)
            end do
            return
        end if
        do a = 1, nba
            f0 = 0d0
            do v = 1, min(nv, nba)
                f0 = f0 + abs(W(v, a))**2
            end do
            delta(a) = occ * f0
            if (a <= nv) delta(a) = delta(a) - occ
        end do
    end subroutine dressed_ref_delta

    ! =====================================================================
    ! Velocity-gauge f-sum-rule completeness of a truncated band basis
    ! (wiki/12 sec. 6). J = Tr[(p + A) rho]/V carries the DIAMAGNETIC term
    ! A N_e/V for every electron; for a filled band it is cancelled by the
    ! paramagnetic (interband) response only with a COMPLETE basis:
    !   S_n^a(k) = sum_{m/=n} 2|p_nm^a|^2/(eps_m - eps_n) = 1 - d2eps_n/dk_a2,
    ! whose full-band average is 1. With nb bands S < 1 and the fraction
    !   eta_a = 1 - <S_n^a>_{k, n occupied}
    ! of A N_e/V survives as a spurious REACTIVE current ~ A = E/omega
    ! (negligible in the near-IR, a plasma mirror at THz: graphene nb = 2 gives
    ! eta = 0.29, R = 0.85 at 100 kV/cm, 3.4 THz). The missing bands lie far
    ! above (Delta >> hbar w, A.p) and respond adiabatically, so subtracting
    ! eta_a (N_e/V) A_a(t) from the current is exact to linear order in A and
    ! removes the artifact at any nb (calc_current_bloch, yn_sbe_vg_sumrule).
    ! Inputs: eigen(nb,nk) [Ha], p(nb,nb,3,nk) velocity matrix elements [a.u.],
    ! occ(nb,nk) ground-state occupations, w(nk) k-weights. Degenerate pairs
    ! (|de| < de_min, e.g. the exactly gapless K point) are skipped.
    ! =====================================================================
    pure subroutine vg_sumrule_eta(nk, nb, eigen, p, occ, w, de_min, eta, scap)
        implicit none
        integer, intent(in) :: nk, nb
        real(8), intent(in) :: eigen(nb, nk), occ(nb, nk), w(nk), de_min
        complex(8), intent(in) :: p(nb, nb, 3, nk)
        real(8), intent(out) :: eta(3), scap(3)
        real(8) :: norm, de
        integer :: ik, n, m, a
        scap = 0d0;  norm = 0d0
        do ik = 1, nk
            do n = 1, nb
                if (occ(n, ik) <= 0d0) cycle
                norm = norm + w(ik) * occ(n, ik)
                do m = 1, nb
                    if (m == n) cycle
                    de = eigen(m, ik) - eigen(n, ik)
                    if (abs(de) < de_min) cycle
                    do a = 1, 3
                        scap(a) = scap(a) + w(ik) * occ(n, ik) * 2d0 * abs(p(n, m, a, ik))**2 / de
                    end do
                end do
            end do
        end do
        if (norm > 0d0) scap = scap / norm
        eta = 1d0 - scap
    end subroutine vg_sumrule_eta

    ! Adiabatic ground-state current of the truncated velocity-gauge Hamiltonian
    ! H_k(A) = eps + A.p at ONE k (wiki/12 sec. 6a). Input: its eigenvectors W
    ! (columns, ascending eigenvalues) and the ground-state occupations occ(n)
    ! of the n-th lowest level (adiabatic continuation from A = 0);
    !     jgs_a = sum_n occ_n <phi_n(A)| p_a + A_a |phi_n(A)>
    !           = d/dA_a sum_n occ_n [E_n(A) + A^2/2]          (Hellmann-Feynman).
    ! A uniform static A is a pure gauge: in a COMPLETE basis this vanishes after
    ! the BZ sum for every A; the truncated basis leaves the diamagnetic
    ! remainder eta N_e A/V at linear order (vg_sumrule_eta) plus a nonlinear
    ! part once A ~ the k-distance to the band-touching points. calc_current_bloch
    ! subtracts the BZ sum of this from the SBE current: the k-trace of rho is
    ! conserved, so the artifact lives only in this ground-state sum and the
    ! subtraction is exact for ANY population, with no adjustable parameter.
    pure subroutine vg_gs_current_k(nb, W, p, A, occ, jgs)
        implicit none
        integer, intent(in) :: nb
        complex(8), intent(in) :: W(nb, nb), p(nb, nb, 3)
        real(8), intent(in) :: A(3), occ(nb)
        real(8), intent(out) :: jgs(3)
        integer :: n, ia, i, j
        complex(8) :: acc
        jgs = 0d0
        do n = 1, nb
            if (occ(n) <= 0d0) cycle
            do ia = 1, 3
                acc = (0d0, 0d0)
                do j = 1, nb
                    do i = 1, nb
                        acc = acc + conjg(W(i, n)) * p(i, j, ia) * W(j, n)
                    end do
                end do
                jgs(ia) = jgs(ia) + occ(n) * (real(acc) + A(ia))
            end do
        end do
    end subroutine vg_gs_current_k

    ! Steady-state damping response of the line set at coherence modulation
    ! frequency w (a.u.): R(0) = 1 (the Markov anchor) by construction.
    pure function colmem_response(nl, cl, mul, w) result(r)
        integer, intent(in) :: nl
        complex(8), intent(in) :: cl(nl), mul(nl)
        real(8), intent(in) :: w
        real(8) :: r
        integer :: j
        r = 0d0
        do j = 1, nl
            r = r + real(cl(j) / (mul(j) + cmplx(0d0, w, 8)))
        end do
    end function colmem_response

    ! =====================================================================
    ! SFSB non-Markovian heat bath (strong-field spin-boson model).
    ! [B25] Boroumand, Thorpe, Bart, Parks, Toutounji, Vampa, Brabec, Wang,
    ! Rep. Prog. Phys. 88, 070501 (2025) -- transcribed from the journal PDF
    ! (wiki/10 sec. 6). The bath enters the driven two-band dynamics ONLY
    ! through the correlation function [B25 Eq. (5)]
    !
    !   C(t) = int_{-inf}^{inf} dw W(w) [ i sin(wt) - (1-cos(wt)) coth(hw/2kT) ]
    !
    ! W(w) is the odd continuum spectral weight of the shift operator
    ! D = exp{-sum_q g_q(b_q^+ - b_q)/(h w_q)} [B25 Eq. (2)]: the (g_q/h w_q)^2
    ! mapping gives W(w) = jo*g(|w|)/w with a dimensionless coupling jo and a
    ! cutoff profile g (g(0)=1). NORMALIZATION ANCHOR: for ANY such W the high-T
    ! limit is C -> -t/T2 with T2 = hbar/(2 pi kB T jo) -- exactly the printed
    ! [B25 sec. 2] Debye/relaxation-time-approximation anchor (verified in
    ! test_bath_corr). Profiles implemented (the letter's supplement with the
    ! full model list is not available -- STRICT provenance, wiki/10 sec. 6):
    !   'ohmic':  g(w) = exp(-w/wc)          (used for ALL numerics in [B25])
    !   'debye':  g(w) = wc^2/(w^2+wc^2)     (Drude-Lorentz; RTA anchor model)
    !   'rta':    C(t) = -t/T2 exactly       (the Markovian reference [B25 Fig 2])
    ! Closed-form anchors (tested): Im C_ohmic = 2 jo atan(wc t) exact;
    ! Re C_ohmic(T=0) = -jo ln(1+wc^2 t^2); Im C_debye = pi jo (1-exp(-wc t)).
    ! All quantities in Hartree a.u. (hbar = kB = 1; kT in Ha).
    ! =====================================================================

    ! High-temperature dephasing time of the jo-coupled bath:
    ! T2 = hbar/(2 pi kB T jo)  [B25 sec. 2, printed formula], a.u.
    pure function bath_t2_high_t(kT, jo) result(t2)
        real(8), intent(in) :: kT, jo
        real(8) :: t2
        t2 = 1d0 / (2d0 * PI * max(kT, 1d-300) * max(jo, 1d-300))
    end function bath_t2_high_t

    ! coth(x) with the series/asymptotic switches that keep it accurate and
    ! finite over the full bath integration range (x > 0).
    pure function bath_coth(x) result(c)
        real(8), intent(in) :: x
        real(8) :: c
        if (x > 19d0) then
            c = 1d0
        else if (x < 1d-4) then
            c = 1d0 / x + x / 3d0
        else
            c = 1d0 / tanh(x)
        end if
    end function bath_coth

    ! Tabulate C(tau) on tau = m*dtau, m = 0..nt [B25 Eq. (5)].
    ! kT <= 0 means T = 0 (coth -> 1). jo <= 0 gives C = 0 (bath off).
    ! model = 'ohmic' | 'debye' | 'rta'; for 'rta' the optional t2_rta (> 0)
    ! overrides the high-T Debye T2 derived from (kT, jo).
    ! The w-integral is composite Simpson on [0, wmax] with the w -> 0 limit
    ! at the first node; wmax/resolution chosen so the closed-form anchors
    ! reproduce to <= 1e-4 relative (test_bath_corr).
    subroutine bath_corr_table(nt, dtau, kT, jo, wc, model, ctab, t2_rta)
        integer, intent(in) :: nt
        real(8), intent(in) :: dtau, kT, jo, wc
        character(*), intent(in) :: model
        complex(8), intent(out) :: ctab(0:nt)
        real(8), intent(in), optional :: t2_rta
        real(8) :: t2, taumax, dom, wmax, w, tau, sn, sh
        real(8), allocatable :: wgt(:), cth(:), omg(:)
        real(8) :: re, im
        integer :: nw, j, m

        ctab(:) = (0d0, 0d0)
        if (jo <= 0d0 .and. trim(model) /= 'rta') return

        select case (trim(model))
        case ('rta')
            t2 = bath_t2_high_t(kT, jo)
            if (present(t2_rta)) then
                if (t2_rta > 0d0) t2 = t2_rta
            end if
            do m = 0, nt
                ctab(m) = cmplx(-m * dtau / t2, 0d0, 8)
            end do
            return
        case ('ohmic')
            wmax = 45d0 * wc
        case ('debye')
            ! the Drude-Lorentz 1/w^2 tail decays slowly; extend past the
            ! thermal crossover so the truncated tail is O((wc/wmax)^2)
            wmax = max(300d0 * wc, 40d0 * max(kT, 0d0))
        case default
            error stop 'bath_corr_table: unknown bath model (ohmic|debye|rta)'
        end select

        taumax = max(nt * dtau, dtau)
        dom = min(2d0 * PI / (20d0 * taumax), wc / 24d0)
        nw = int(wmax / dom) + 1
        nw = min(max(nw + mod(nw, 2), 8), 4000000)   ! even panel count, capped
        dom = wmax / nw

        allocate(wgt(0:nw), cth(0:nw), omg(0:nw))
        do j = 0, nw
            omg(j) = j * dom
            ! Simpson: 1,4,2,...,2,4,1 times dom/3; fold in 2*jo*g(w)/w
            if (j == 0 .or. j == nw) then
                wgt(j) = dom / 3d0
            else if (mod(j, 2) == 1) then
                wgt(j) = 4d0 * dom / 3d0
            else
                wgt(j) = 2d0 * dom / 3d0
            end if
            if (j > 0) then
                w = omg(j)
                if (trim(model) == 'ohmic') then
                    wgt(j) = wgt(j) * 2d0 * jo * exp(-w / wc) / w
                else
                    wgt(j) = wgt(j) * 2d0 * jo * wc**2 / ((w**2 + wc**2) * w)
                end if
                if (kT > 0d0) then
                    cth(j) = bath_coth(w / (2d0 * kT))
                else
                    cth(j) = 1d0
                end if
            else
                wgt(0) = wgt(0) * 2d0 * jo    ! limit node: g(0)/w folded below
                cth(0) = 1d0
            end if
        end do

        !$omp parallel do private(m, tau, re, im, j, sn, sh) schedule(static)
        do m = 0, nt
            tau = m * dtau
            ! w -> 0 limit of g/w * [i sin - (1-cos) coth] = i tau - kT tau^2
            re = -wgt(0) * merge(max(kT, 0d0) * tau**2, 0d0, kT > 0d0)
            im =  wgt(0) * tau
            do j = 1, nw
                sn = sin(omg(j) * tau)
                sh = sin(0.5d0 * omg(j) * tau)
                im = im + wgt(j) * sn
                re = re - wgt(j) * 2d0 * sh * sh * cth(j)
            end do
            ctab(m) = cmplx(re, im, 8)
        end do
        !$omp end parallel do

        deallocate(wgt, cth, omg)
    end subroutine bath_corr_table

    ! Second-order (Dyson) conduction population of a driven two-level pair
    ! with the bath memory kernel [B25 Eq. (3)]:
    !
    !   nc(t) = (1/2) Re int^t dt1 int^t1 dt2 Om*(t1) Om(t2)
    !                                  exp[ i S(t1,t2) + C(t1-t2) ]
    !
    ! S(t1,t2) = int_{t2}^{t1} Es dtau, Es = sqrt(dE^2 + |Om|^2) the
    ! Stark-shifted gap [B25 after Eq. (4)]. Om(t) = 2 d(K_t) . E(t) is the
    ! generalized Rabi frequency (a.u., e = hbar = 1), sampled by the caller
    ! along the K_t = K + A(t) trajectory. The phase factorizes across steps,
    ! exp[iS] = e^{i th(t1)} e^{-i th(t2)}; the bath kernel exp[C(t1-t2)] does
    ! NOT -- that non-factorizable factor IS the memory (non-Markovian) part,
    ! so the inner t2 integral is a true history sum (Volterra), truncatable
    ! at nwin steps once |exp(Re C)| has decayed. C = 0 recovers textbook
    ! second-order perturbation theory; C = -tau/T2 recovers the relaxation
    ! time approximation (both asserted in test_sfsb_kernel).
    ! nwin <= 0 or >= nt means the full history.
    pure subroutine sfsb_nc_series(nt, dt, om, es, ctab, nwin, nc)
        integer, intent(in) :: nt, nwin
        real(8), intent(in) :: dt
        complex(8), intent(in) :: om(0:nt)
        real(8), intent(in) :: es(0:nt)
        complex(8), intent(in) :: ctab(0:nt)
        real(8), intent(out) :: nc(0:nt)
        real(8), allocatable :: th(:)
        complex(8), allocatable :: phi(:), kexp(:)
        complex(8) :: inner
        real(8) :: g_prev, g_cur
        integer :: i, m, lo, win

        allocate(th(0:nt), phi(0:nt), kexp(0:nt))
        win = nwin
        if (win <= 0 .or. win > nt) win = nt

        th(0) = 0d0
        do i = 1, nt
            th(i) = th(i - 1) + 0.5d0 * (es(i - 1) + es(i)) * dt
        end do
        do i = 0, nt
            phi(i)  = om(i) * exp(cmplx(0d0, -th(i), 8))
            kexp(i) = exp(ctab(i))
        end do

        nc(0) = 0d0
        g_prev = 0d0
        do i = 1, nt
            lo = max(0, i - win)
            ! trapezoid over t2 in [t_lo, t_i]; kexp(0) = exp(C(0)) = 1
            inner = 0.5d0 * (phi(lo) * kexp(i - lo) + phi(i))
            do m = lo + 1, i - 1
                inner = inner + phi(m) * kexp(i - m)
            end do
            inner = inner * dt
            g_cur = 0.5d0 * real(conjg(om(i)) * exp(cmplx(0d0, th(i), 8)) * inner)
            nc(i) = nc(i - 1) + 0.5d0 * (g_prev + g_cur) * dt
            g_prev = g_cur
        end do

        deallocate(th, phi, kexp)
    end subroutine sfsb_nc_series

end module sbe_superres_ssbe
