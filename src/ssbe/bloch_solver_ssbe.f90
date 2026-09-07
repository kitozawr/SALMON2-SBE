module bloch_solver_ssbe
    use math_constants, only: pi, zi
    use phys_constants, only: au_ev
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast
    use gs_info_ssbe
    use util_ssbe, only: split_range
    implicit none

    private
    public :: s_sbe_bloch_solver, init_sbe_bloch_solver, calc_current_bloch, &
              dt_evolve_bloch_cf4, calc_trace, calc_energy, calc_bloch_population_k, &
              calc_unfolded_population_k, calc_diabatic_population_k, &
              calc_diabatic_unfolded_population_k, calc_intraband_current_houston, &
              calc_nex_nonad

    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
        logical :: flag_vnl_correction

        ! Frozen core handling
        logical, allocatable :: is_active(:) ! .true. if band is active, .false. if frozen
        integer :: n_active_bands = 0
        integer, allocatable :: active_idx(:)  ! Mapping: 1..n_active -> global band index

        ! Houston-basis branch (wave-packet) positions X_a(k,t) used by the
        ! Kuhn-Zurek/Caldeira-Leggett dephasing kernel. By the origin-shift
        ! invariance of the dephasing map (only differences X_a-X_b enter),
        ! the choice X_a(0)=0 is physically irrelevant; it merely fixes a
        ! reproducible reference for restarts.
        real(8), allocatable :: X_branch(:, :)  ! (1:nb, ik_min:ik_max)

        ! Kuhn-Zurek/Caldeira-Leggett decoherence: lambda = kB*T / tau_m
        real(8) :: lambda_decoh = 0d0
        logical :: flag_decoh   = .false.

        ! k-local impact-ionization Lindblad channel (Stobbe rate fit):
        !   gamma(e_kin) = P (e_kin - E_th)^4 Theta(e_kin - E_th),
        ! e_kin measured from the global field-free CBM. Each event is split
        ! (Hartree-Fock / two-particle closure) into two frozen-rate
        ! amplitude-damping channels: primary relaxation h -> h' and cold-pair
        ! creation v1 -> c1, applied in the same Houston basis as the
        ! Kuhn-Zurek dephasing (no extra ZHEEV).
        logical :: flag_impact  = .false.
        real(8) :: ii_pref_au   = 0d0   ! P in 1/(Ha^a a.u.time)
        real(8) :: ii_exponent  = 4d0   ! fit exponent a (4 GaAs Stobbe, 2 Si Keldysh)
        real(8) :: ii_eth_au    = 0d0   ! threshold E_th [Ha]
        real(8) :: ii_ramp_au   = 0d0   ! linear Theta-smoothing width [Ha]
        real(8) :: ii_ecbm_au   = 0d0   ! global CBM of the field-free bands [Ha]
        real(8) :: ii_eg_au     = 0d0   ! band gap E_g [Ha] (primary loses E_g)
        integer :: nv_act       = 0     ! valence branches inside the active subspace
        real(8) :: occ_max      = 2d0   ! 2 (scalar bands) / 1 (spinor bands)
        real(8), allocatable :: f0_act(:, :)   ! (nba, nk) INITIAL occupation of the active bands (doped/thermal sheet): dressed-reference baseline

        ! Bandgap-renormalization-coupled impact-ionization threshold (Part C7).
        ! When the excited carrier density exceeds bgr_n_gate, the II threshold
        ! shrinks: E_th(t) = E_th0 - |K n(t)^(1/3)|. [Vashishta-Kalia PRB 25, 6492]
        logical :: flag_bgr     = .false.
        real(8) :: ii_eth0_au   = 0d0   ! fixed input threshold E_th0 [Ha]
        real(8) :: bgr_n_gate   = 5d18  ! gate density [cm^-3]
        real(8) :: bgr_coeff    = 1.9d-8 ! K [eV cm]
        integer :: homo_idx     = 0     ! HOMO band index (valence/conduction split)
        real(8) :: au_dens_cm3  = 0d0   ! a.u.^-3 -> cm^-3 number-density conversion
        ! Cost-preserving mask for the FULL-basis dressed ring (yn_sbe_full_dressed):
        ! the dressed projection is truncation-free (all bands), but the ring kernels
        ! scatter only the dressed states whose energy is in the frozen-window range
        ! [ring_e_lo, ring_e_hi] (a.u.) -- restores the O(nk^2 * n_active) cost while
        ! keeping the correct (untruncated) dressed basis (wiki/06 sec.6).
        logical :: ring_mask_on = .false.
        real(8) :: ring_e_lo    = -1d30
        real(8) :: ring_e_hi    =  1d30
        ! Dissipator sub-cycling (Part C8): split the dissipative half-step into
        ! diss_subcycle CPTP sub-steps when the collision rate is fast vs dt.
        real(8) :: eph_numax_au = 0d0   ! estimated peak e-ph rate [1/a.u.time]

        ! Carrier-carrier (e-e/e-h) thermalization (Part F): CPTP relaxation of
        ! the adiabatic populations toward a Fermi-Dirac with the SAME number
        ! and energy, plus coherence damping (EID). Off unless yn_sbe_eeh='y'.
        logical :: flag_eeh   = .false.
        real(8) :: eeh_nu_au  = 0d0     ! carrier-carrier rate [1/a.u.time]

        ! Auger recombination (Sec 13): density-gated, number-conserving CPTP
        ! channel. Per-carrier rate gamma = auger_c_au * n^2 (n in cm^-3), so the
        ! total recombination rate is R = C n^3. Inert below auger_n_gate_cm3.
        logical :: flag_auger      = .false.
        real(8) :: auger_c_au      = 0d0   ! C[cm^6/s]*(au_fs*1e-15): gamma=auger_c_au*n_cm3^2
        real(8) :: auger_n_gate_cm3 = 0d0  ! activation density [cm^-3]
        real(8) :: auger_eg_au     = 0d0   ! band gap E_g [Ha] (hot-carrier target offset)
        real(8) :: n_exc_cm3       = 0d0   ! running excited-carrier density [cm^-3]

        ! 2D Rana Auger / carrier multiplication (graphene [R07]; wiki/07 sec.6):
        ! net CPTP pair relaxation R - G of the Dirac-cone populations on the
        ! instantaneous quasi-Fermi levels (rana_auger_dpop). Ring-gated (needs
        ! the global population gather, like graphene e-ph). Replaces BOTH the
        ! k-local C n^3 Auger and the gap-threshold ring Auger for graphene.
        logical :: flag_rana2d  = .false.
        real(8) :: rana_vf_au   = 0d0   ! Dirac velocity [a.u.] (1e8 cm/s [R07])
        real(8) :: rana_eps_r   = 0d0   ! background eps_r (R07 Fig.4; substrate)
        real(8) :: rana_kt_au   = 0d0   ! carrier/bath temperature [Ha]
        real(8) :: rana_area_au = 0d0   ! 2D cell area [a.u.^2]

        ! Approved-improvements bundle (wiki/00 2026-07-04, PR-approved):
        real(8) :: ring_vq_floor = 0d0   ! B3: relative vq floor (0 = off)
        logical :: flag_ii_fk    = .false. ! A5: FK-softened II threshold
        real(8) :: ii_fk_mu      = 0d0   !     reduced mass [m_e]
        real(8) :: ii_phassist   = 0d0   ! A1: phonon-assisted sideband strength
        logical :: flag_ii_holes = .false. ! A2: hole-initiated II/Auger
        real(8) :: ii_cpcn       = 0d0   !     cited Cp/Cn ratio
        ! C1 ledger: CUMULATIVE per-channel conduction-population and energy
        ! change (1 = e-ph, 2 = II, 3 = ring Auger, 4 = 2D Rana). Identical on
        ! every rank (kernels are deterministic on the gathered state).
        real(8) :: led_dn(4) = 0d0
        real(8) :: led_de(4) = 0d0
        integer :: eph_ip_ac = 0        ! A4: index of the appended acoustic mode (0 = none)

        ! Nonlocal impact ionization (Part C4): the hot electron ionizes a
        ! valence partner drawn from the WHOLE BZ (momentum exchange), so the
        ! partner-population / Pauli factors use the global BZ-averaged active-
        ! band occupations gathered once per step (rides the ring/all-gather).
        ! The full momentum-resolved final-state version is the refinement.
        logical :: flag_nl_ii = .false.
        real(8), allocatable :: glob_occ(:)   ! (nba) BZ-averaged active-band occupation

        ! Sublattice (band-unfolding) awareness of the impact-ionization
        ! channel. When the unfolding weights gs%unfold_w are available
        ! (non-zero), the k-local two-particle event is resolved per FCC
        ! sublattice so that the secondary pair is created in the SAME
        ! primitive-cell sector (same primitive crystal momentum) as the
        ! primary -- restoring primitive momentum conservation that the
        ! cubic folding hides, and removing the spurious cross-sublattice
        ! ("false generation") events. Reduces identically to the folded
        ! channel when the weights are absent.
        logical :: flag_unfold_ii = .false.

        ! Coulomb (time-dependent Hartree-Fock / exchange) renormalization,
        ! Golde-Kira-Meier-Koch SBE (Phys. Status Solidi B 248, 863 (2011),
        ! Eqs. 4-5). The exchange (Fock) self-energy
        !   Sigma_nm(k) = - sum_{q/=k} V(k-q) rho_nm(q)
        ! is a NON-k-local mean field (couples all k). Added to H_VG it
        ! reproduces both the renormalized single-particle energies
        ! (diagonal: eps~ = eps - sum_q V_{k-q} f_q) and the renormalized
        ! Rabi frequency (off-diagonal: Omega = d.E + sum_q V_{k-q} p_q),
        ! with the (1-f_e-f_h) Pauli factor coming for free from the
        ! von Neumann commutator. Frozen over a dt step (mean-field predictor).
        logical :: flag_coulomb = .false.
        logical :: flag_hf_subproj = .false. ! project Sigma^HF onto FCC sublattice blocks
        logical :: flag_coset_proj = .false. ! project the momentum coupling p block-diagonal over cosets

        ! Population-relaxing electron-phonon Lindblad (Part C5, super-mode).
        ! k-local skeleton: each adiabatic level relaxes toward an energy-matched
        ! (+-hw) partner at the saturated rate nu(eps_kin), emission ~ (N_B+1),
        ! absorption ~ N_B, Pauli-blocked. CPTP (amplitude damping). Off unless
        ! yn_sbe_eph='y'. Drives THz bleaching (collision-rate saturation).
        logical :: flag_eph     = .false.
        real(8) :: eph_nusat_au = 0d0   ! saturation rate nu_sat [1/a.u.time]
        real(8) :: eph_eps0_au  = 0d0   ! nu(eps) onset eps_0 [Ha]
        real(8) :: eph_kt_au    = 0d0   ! phonon-bath kT [Ha]: detailed balance on the
                                        ! REALIZED transfer needs kT, not only N_B(hw_p)
        ! Enforce detailed balance on the realized transfer (and a pair-symmetric nu).
        ! Set from the material: ON for the gapless 2D Dirac materials, where the
        ! appended acoustic mode is ~5 meV against a ~0.1 eV search width and the
        ! historical per-mode split runs the carriers to thousands of K with no field.
        ! OFF elsewhere, so the gapped materials stay bit-identical to their published
        ! validations -- the same violation exists there, but sigma/hw is O(1) rather
        ! than O(20), so it is a small correction rather than a runaway.
        logical :: eph_db_realized = .false.
        real(8) :: eph_n        = 2d0   ! nu(eps) shape exponent
        real(8) :: eph_ib_scale = 1d0   ! gap-straddling (BTBT) rate calibration factor
        real(8) :: eph_sigma_au = 0d0   ! energy-bin width sigma_E [Ha]
        real(8) :: eph_ecbm_au  = 0d0   ! field-free CBM [Ha]
        real(8) :: eph_evbm_au  = 0d0   ! field-free VBM [Ha]
        ! Ring virtual-transient gate (real-vs-virtual separation): the ring
        ! POPULATION kernels see the persistent Houston floor f_gate -- drops
        ! instantly with f, rises toward f with time constant ring_gate_tau_au
        ! ~ 2*pi/Egap (the energy-time virtuality scale). Sub-lifetime LZ /
        ! dressing transients never enter the transfer rates; they still lose
        ! coherence at the full nu (the gout/T2 role is untouched).
        real(8) :: ring_gate_tau_au = 0d0            ! 0 = gate off
        real(8), allocatable :: f_ring_gate(:, :)    ! (nba, nk) persistent floor
        integer :: eph_nph      = 0     ! number of phonon modes in the table
        real(8), allocatable :: eph_hw(:)   ! phonon energies hw_p [Ha]
        real(8), allocatable :: eph_nb(:)   ! Bose factors N_B(hw_p, T_ph)
        real(8), allocatable :: eph_wrel(:) ! relative golden-rule weights (sum=1)

        ! Collisional-memory (non-Markovian) dephasing of the e-ph ring gout
        ! (wiki/10 sec. 8.6, maintainer-approved): the instantaneous
        ! exp(-gout*tau/2) is replaced by a memory convolution whose kernel is
        ! built VERBATIM from the cited phonon table (Lorentzian lines at the
        ! mode energies, (N+1)/N thermal weights, width 1/tau_c = sigma_E by
        ! default -- colmem_lines in sbe_superres). Auxiliary coherence fields
        ! zmem live elementwise in the Houston frame (upper triangle a<b),
        ! attached to sorted-branch indices like every other ring quantity.
        logical :: flag_colmem = .false.
        integer :: colmem_nl   = 0
        complex(8), allocatable :: colmem_c(:), colmem_mu(:)
        complex(8), allocatable :: zmem(:, :, :, :)  ! (nba,nba,nl,ik_min:ik_max)
        ! population-sector memory filter (wiki/10 sec. 8.8): the ring
        ! collision kernels read the memory-filtered f instead of the
        ! instantaneous Houston populations. Rank-identical (filters the
        ! GATHERED f_all, like the ring gate).
        logical :: flag_colmem_pop = .false.
        complex(8), allocatable :: zpop(:, :, :)     ! (nba, nl, nk)
        ! Option A (wiki/10 sec. 3A/8.10): the ring channels measure carriers
        ! against the FIELD-ROTATED ground state (dressed reference) instead
        ! of the static {occ,0} -- the rotation background delta0 is
        ! subtracted from the Houston populations before the gather.
        logical :: flag_dressed_ref = .false.
        ! graphene 2D collisional-memory analog (wiki/10 sec. 8.11): the Rana
        ! (Coulomb) SOURCE densities n2d/p2d are memory-filtered with the 2D
        ! Dirac-plasmon line w_pl(n,p; Q_TF) of the instantaneous e-h plasma
        ! (two lines +/- w_pl x {n, p}); the e-ph sectors keep the phonon lines.
        logical :: flag_colmem_2d = .false.
        real(8) :: colmem_tauc_au = 0d0
        complex(8) :: zpl(2, 2) = (0d0, 0d0)        ! (line, {n, p}) auxiliary fields
        logical :: zpl_init = .false., zpl_printed = .false.
        ! two-temperature Coulomb sector (wiki/12 sec. 3): carrier temperature and
        ! quasi-Fermi levels fitted from the moments of the gathered Dirac populations
        ! (diagnostics for *_sbe_te.data; kT_e drives the Rana R-G and the plasmon line)
        logical :: flag_rana_te = .false.
        ! velocity-gauge f-sum-rule correction (wiki/12 sec. 6): the fraction eta of
        ! the diamagnetic current A N_e/V the truncated basis cannot cancel
        logical :: flag_vg_sumrule = .false.
        real(8) :: vg_eta(3) = 0d0, vg_scap(3) = 0d0, vg_ne_cell = 0d0
        real(8) :: te_au = 0d0, te_mu_c = 0d0, te_mu_h = 0d0, te_n2d = 0d0, te_p2d = 0d0

        real(8) :: coul_pref     = 0d0  ! strength * 4 pi / (eps * Omega_cell * Nk)
        real(8) :: coul_screen2  = 0d0  ! kappa^2 [1/Bohr^2] (Yukawa regulariser)
        ! A7: 2D-sheet exchange kernel (graphene): V_2D = 2 pi/(eps A Nk (q+kappa))
        logical :: coul_2d       = .false.
        real(8) :: coul_pref2d   = 0d0  ! strength * 2 pi / (eps * A_cell * Nk)
        real(8) :: coul_screen1  = 0d0  ! kappa [1/Bohr] (2D linear screening)
        integer :: icomm         = 0    ! MPI communicator (for the non-local exchange sum)
        complex(8), allocatable :: sigma_hf(:, :, :) ! (nba, nba, ik_min:ik_max) exchange Sigma

        ! Ring/pipeline MPI (Part D): systolic-ring all-pairs pass replacing the
        ! all-gather for the non-local sums in super-mode. Memory O(N_k/P + one
        ! transit block) instead of O(N_k). [Plimpton JCP 117, 1 (1995)]
        logical :: flag_ring = .false.
        integer :: irank = 0, nproc = 1
        integer, allocatable :: itbl_min(:), itbl_max(:)  ! k-partition per rank (0:nproc-1)

        ! Monkhorst-Pack momentum-conservation map (built lazily on first use by
        ! the momentum-resolved nonlocal impact ionization). kmap_ok=.false. on a
        ! non-MP / symmetry-reduced grid -> the momentum-resolved II is gated off.
        logical :: kmap_built = .false.
        logical :: kmap_ok    = .false.
        integer :: kmap_n(3)  = 0
        integer, allocatable :: kmap_idx(:, :)   ! (3, nk) 0-based MP triple
        integer, allocatable :: kmap_lut(:)      ! (0:nk-1) flattened triple -> ik
    end type

    !=========================================================================
    ! CF4 (commutator-free Magnus 4) + Suzuki-Yoshida composition constants
    !=========================================================================
    ! Two-point Gauss-Legendre nodes on [0,1]: c = 1/2 -+ sqrt(3)/6
    real(8), parameter :: cf4_c1 = 0.21132486540518713d0
    real(8), parameter :: cf4_c2 = 0.78867513459481287d0
    ! CF4 combination weights: alpha = 1/4 -+ sqrt(3)/6
    real(8), parameter :: cf4_alpha1 =  0.53867513459481287d0
    real(8), parameter :: cf4_alpha2 = -0.03867513459481287d0
    ! Suzuki-Yoshida triple-jump constants (4th-order composition of a
    ! 2nd-order base scheme): p1 + p2 + p1 = 1
    real(8), parameter :: yoshida_p1 =  1.35120719196d0
    real(8), parameter :: yoshida_p2 = -1.70241438392d0

contains

! Build the per-material phonon table for the e-ph channel: energies hw_p [Ha],
! Bose factors N_B(hw_p, T), and normalized relative golden-rule weights
! w_p ~ D_p^2/hw_p (the common pi/rho factors cancel in the per-material
! normalization). Si: 6 intervalley (g/f) phonons. GaAs: the Frohlich LO plus
! 5 intervalley; the polar LO is given the dominant weight (sum of the
! intervalley weights) -- a documented skeleton choice (the full Frohlich asinh
! energy-dependence is a later refinement). [tables + sources in 02_constants]
! Abort with a guiding message when a material-dependent channel is requested
! for a material that is not in the registry. Lists the supported names so the
! fix (add a `case` in get_material_params + MAT_SUPPORTED) is obvious.
subroutine stop_unknown_material(name, channel)
    use sbe_superres_ssbe, only: MAT_SUPPORTED
    implicit none
    character(*), intent(in) :: name, channel
    write(*, '(a)') '# ERROR: material "'//trim(name)//'" is not in the SBE material registry,'
    write(*, '(a)') '#        but '//trim(channel)//' needs its constants.'
    write(*, '(a)') '#        Supported: '//MAT_SUPPORTED
    write(*, '(a)') '#        Add a case to get_material_params() in sbe_superres_ssbe.f90.'
    error stop 'unknown epm_material for an SBE channel (see message above)'
end subroutine stop_unknown_material


! Abort when a channel is enabled for a material that has NO cited constants for
! it (provenance gate). Constants are never transferred from another material:
! no source => forbidden. The fix is to add a CITED constant + flip the *_ok
! flag in get_material_params, or to disable the channel for this material.
subroutine stop_forbidden_channel(name, channel)
    implicit none
    character(*), intent(in) :: name, channel
    write(*, '(a)') '# ERROR: '//trim(channel)//' is FORBIDDEN for material "'//trim(name)//'":'
    write(*, '(a)') '#        no cited CdS-/material-specific constants exist for it, and none'
    write(*, '(a)') '#        may be transferred from another material (no source = not valid).'
    write(*, '(a)') '#        Disable the channel for this material, or add a CITED constant and'
    write(*, '(a)') '#        set the matching *_ok flag in get_material_params().'
    error stop 'forbidden SBE channel for this material (see message above)'
end subroutine stop_forbidden_channel


subroutine init_eph_phonon_table(sbe, mp, kT_au, ac_qtyp_au, ac_xi_ev)
    use sbe_superres_ssbe, only: bose_factor, mev_to_ha, s_material_params
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_material_params),  intent(in)    :: mp
    real(8),                  intent(in)    :: kT_au
    ! A4: when > 0, append a GRID-RESOLVED quasi-elastic acoustic deformation
    ! mode: q = ac_qtyp_au (one k-grid spacing -- the smallest hop the mesh
    ! supports), hw_ac = c_s q, effective coupling D_ac(q) = Xi_d q (the
    ! long-wavelength deformation-potential matrix element). Constants cited
    ! per material in the registry (Si: Jacoboni-Reggiani; GaAs: Fischetti-
    ! Laux; graphene: Hwang-Das Sarma). Removes the sub-optical-phonon
    ! cooling freeze-out.
    real(8), intent(in), optional :: ac_qtyp_au
    ! substrate-dependent DP override (sbe_eph_ac_xi_ev; graphene use-case)
    real(8), intent(in), optional :: ac_xi_ev
    integer :: np, p, nadd
    real(8) :: wsum, hw_ac_mev, wraw_ac, q_cm, cs_au, xi_ev
    real(8), parameter :: BOHR_CM = 0.52917721067d-8, V_AU_CMPS = 2.18769126364d8

    np = mp%eph_nph
    nadd = 0
    hw_ac_mev = 0d0; wraw_ac = 0d0
    if (present(ac_qtyp_au)) then
        xi_ev = mp%eph_ac_xi_ev
        if (present(ac_xi_ev)) then
            if (ac_xi_ev > 0d0) xi_ev = ac_xi_ev     ! substrate override
        end if
        if (ac_qtyp_au > 0d0 .and. xi_ev > 0d0) then
            nadd = 1
            cs_au = mp%eph_ac_cs_cmps / V_AU_CMPS
            hw_ac_mev = cs_au * ac_qtyp_au * 27.211386245988d3   ! Ha -> meV
            q_cm = ac_qtyp_au / BOHR_CM
            ! raw weight in the table convention D[1e8 eV/cm]^2 / hw[meV]
            wraw_ac = (xi_ev * q_cm / 1d8)**2 / max(hw_ac_mev, 1d-12)
        end if
    end if
    sbe%eph_nph = np + nadd
    sbe%eph_kt_au = kT_au
    ! Detailed balance on the REALIZED transfer, for every material. The split has to
    ! be taken at the energy a pair actually moves, not at hw_p, whenever the matching
    ! width is not small against the phonon energies -- and it never is: sigma/hw runs
    ! 3.2..20 for silicon (0.2 eV against a 13.25 meV appended acoustic mode carrying
    ! 23.4 % of the weight) just as it runs ~20 for graphene. Restricting this to the
    ! 2D Dirac registry was a staging decision, not a physical one; wiki/12 has the
    ! measurement that closed it, including why narrowing sigma is not an alternative
    ! (it removes the leak and the useful relaxation in equal measure).
    sbe%eph_db_realized = .true.
    allocate(sbe%eph_hw(np+nadd), sbe%eph_nb(np+nadd), sbe%eph_wrel(np+nadd))

    wsum = sum(mp%eph_wraw(1:np)) + wraw_ac
    do p = 1, np
        sbe%eph_hw(p)   = mev_to_ha(mp%eph_hw_mev(p))
        sbe%eph_wrel(p) = mp%eph_wraw(p) / max(wsum, 1d-300)
        sbe%eph_nb(p)   = bose_factor(sbe%eph_hw(p), kT_au)
    end do
    if (nadd == 1) then
        sbe%eph_hw(np+1)   = mev_to_ha(hw_ac_mev)
        sbe%eph_wrel(np+1) = wraw_ac / max(wsum, 1d-300)
        sbe%eph_nb(np+1)   = min(bose_factor(sbe%eph_hw(np+1), kT_au), 1d12)
    end if
end subroutine init_eph_phonon_table


subroutine init_sbe_bloch_solver(sbe, gs, nb_sbe, icomm, verbose)
    use util_ssbe
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast, &
                             comm_get_max, comm_sync_all
    use salmon_global, only: frozen_core_threshold_ev, frozen_free_threshold_ev, &
                             yn_sbe_full_dressed, &
                             sbe_decoh_temperature_k, sbe_decoh_tau_m_fs, yn_sbe_spinor, &
                             yn_sbe_impact_ionization, sbe_ii_prefactor, &
                             sbe_ii_threshold_ev, sbe_ii_ramp_ev, &
                             sbe_ii_form, sbe_ii_exponent, &
                             yn_sbe_coulomb, sbe_coulomb_epsilon, &
                             sbe_coulomb_strength, sbe_coulomb_screen_au, &
                             yn_sbe_hf_sublattice_proj, yn_sbe_coset_proj, &
                             yn_sbe_eph, sbe_eph_temperature_k, sbe_eph_nu_sat, &
                             sbe_eph_eps0_ev, sbe_eph_n, sbe_search_sigma_e_ev, &
                             sbe_ring_gate_fs, sbe_eph_interband_scale, &
                             yn_sbe_bgr_threshold, sbe_bgr_n_gate, sbe_bgr_coeff, &
                             yn_sbe_superres, yn_sbe_eeh, sbe_eeh_nu_sat, epm_material, &
                             yn_sbe_auger, sbe_auger_c_cm6s, sbe_auger_n_gate_cm3, &
                             sbe_ring_vq_floor, yn_sbe_ii_fk_soften, sbe_ii_fk_mu, &
                             sbe_ii_phassist, yn_sbe_ii_holes, yn_sbe_eph_acoustic, &
                             sbe_eph_ac_xi_ev, &
                             yn_sbe_colmem, sbe_colmem_tau_fs, yn_sbe_colmem_pop, &
                             yn_sbe_dressed_ref, yn_sbe_rana_te, yn_sbe_vg_sumrule, &
                             num_kgrid
    use sbe_superres_ssbe, only: bose_factor, s_material_params, colmem_lines, colmem_response, &
                                 get_material_params, MAT_SUPPORTED, vg_sumrule_eta, vg_gs_current_k
    use math_constants, only: pi
    use phys_constants, only: au_fs, kB_au, au_ev
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    ! Print the channel banners (default .true.). The multiscale driver sets it
    ! .false. for all but one macropoint so the diagnostics appear once, not once
    ! per macropoint group.
    logical, intent(in), optional :: verbose
    integer :: ik, ib, nk_proc, irank, nproc, ierr, count_active
    integer, allocatable :: itbl_min(:), itbl_max(:)
    real(8) :: eigen_ev, fermi_energy_ev
    integer :: homo_idx, lumo_idx
    integer :: ik_gamma
    real(8) :: dmin, dk
    real(8), allocatable :: eig_gamma(:)
    logical :: lprint
    type(s_material_params) :: mp
    character(20) :: ii_form_eff
    real(8) :: ii_exp_eff, ii_pref_eff, ii_thr_eff, coul_eps_eff
    real(8) :: colmem_tauc

    call comm_get_groupinfo(icomm, irank, nproc)

    ! Banner-print gate: only the group root, and only when verbose (or absent).
    lprint = (irank == 0)
    if (present(verbose)) lprint = lprint .and. verbose

    ! Per-material constants come from the single registry in sbe_superres_ssbe.
    ! Each channel auto-selects through `mp`; a value the user set explicitly in
    ! the namelist (non-sentinel) always overrides the material default below.
    mp = get_material_params(epm_material)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)
    ! store the partition + rank info for the ring-pipeline (Part D)
    sbe%irank = irank
    sbe%nproc = nproc
    ! Communicator the k-points are distributed over (the per-macropoint group
    ! icomm_macro in multiscale, or the world comm in single-cell realtime). Set
    ! unconditionally here: every collective in the new channels reduces over it
    ! (Coulomb all-gather/ring AND the nonlocal-II / BGR global reductions), so
    ! it must be valid even when Coulomb is off. (Previously it was set only in
    ! the Coulomb block, leaving icomm=MPI_COMM_NULL for a BGR/nl-II-only run.)
    sbe%icomm = icomm
    allocate(sbe%itbl_min(0:nproc-1), sbe%itbl_max(0:nproc-1))
    sbe%itbl_min = itbl_min
    sbe%itbl_max = itbl_max

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
    
    ! =========================================================================
    ! ИНИЦИАЛИЗАЦИЯ rho: ИСПОЛЬЗУЕМ gs%occup (КАК В ОРИГИНАЛЕ!)
    ! =========================================================================
    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)  ! ← КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ
        end do
    end do
    
    sbe%flag_vnl_correction = .false.

    ! =========================================================================
    ! Houston-basis branch positions X_a(k,t): explicit zero-init.
    ! By the invariance theorem of the dephasing kernel under a global shift
    ! of the X_a origin (only X_a-X_b enters exp[-lambda(X_a-X_b)^2 tau]),
    ! the choice X_a(t=0) = 0 carries no physical content; it only fixes a
    ! reproducible convention for restarts.
    ! =========================================================================
    allocate(sbe%X_branch(1:sbe%nb, sbe%ik_min:sbe%ik_max))
    sbe%X_branch = 0d0

    ! =========================================================================
    ! Kuhn-Zurek/Caldeira-Leggett decoherence strength: lambda = kB*T / tau_m
    ! Both temperature and relaxation time must be positive to enable it;
    ! otherwise the dephasing map reduces identically to the identity (D=0),
    ! which is trivially CPTP.
    ! =========================================================================
    if (sbe_decoh_temperature_k > 0d0 .and. sbe_decoh_tau_m_fs > 0d0) then
        sbe%lambda_decoh = kB_au * sbe_decoh_temperature_k / (sbe_decoh_tau_m_fs / au_fs)
        sbe%flag_decoh   = .true.
    else
        sbe%lambda_decoh = 0d0
        sbe%flag_decoh   = .false.
    end if

    ! Kuhn-Zurek (single-particle wave-packet) dephasing is UNPHYSICAL for gapless
    ! Dirac carriers: in graphene coherence loss is an intrinsically many-body
    ! effect, not a wave-packet-separation kernel. Forbid the combination rather
    ! than silently producing a meaningless decoherence rate.
    if (sbe%flag_decoh .and. trim(epm_material) == 'graphene') then
        write(*, '(a)') '# ERROR: Kuhn-Zurek decoherence (sbe_decoh_*) is not valid for graphene.'
        write(*, '(a)') '#        Gapless Dirac coherence loss is a many-body effect, not a'
        write(*, '(a)') '#        single-particle wave-packet dephasing. Disable sbe_decoh_temperature_k'
        write(*, '(a)') '#        / sbe_decoh_tau_m_fs for graphene.'
        error stop 'Kuhn-Zurek decoherence forbidden for graphene (many-body coherence)'
    end if

    ! =========================================================================
    ! ЛОГИКА is_active (Frozen Core)
    ! =========================================================================
    
    ! 1. Calculate Fermi Energy (Assumes closed-shell / even number of electrons)

    if (yn_sbe_spinor == 'y') then
        ! Spinor bands: one electron per band, ne occupied bands
        homo_idx = gs%ne
    else
        homo_idx = gs%ne / 2
        if (mod(gs%ne, 2) /= 0 .and. irank == 0) then
            write(*, '(a)') 'WARNING: Odd number of electrons. Fermi energy assumes closed-shell.'
        end if
    end if
    lumo_idx = homo_idx + 1

    ! 1b. Anchor the frozen-core window at the GAMMA point, NOT at k-index 1.
    !     The band structure is complex and the k-ordering is grid-dependent: on
    !     a 9^3 MP grid k=1 is the CORNER (-4/9,-4/9,-4/9), not Gamma, so keying
    !     the window off "the first k-point" is fragile. Every rank finds Gamma
    !     the same way -- the reduced-coordinate point closest to the origin
    !     (exactly (0,0,0) on odd grids) -- from the REPLICATED gs%kpoint, so
    !     ik_gamma is identical on all ranks.
    ik_gamma = 1
    dmin = huge(1d0)
    do ik = 1, sbe%nk
        dk = gs%kpoint(1, ik)**2 + gs%kpoint(2, ik)**2 + gs%kpoint(3, ik)**2
        if (dk < dmin) then
            dmin = dk
            ik_gamma = ik
        end if
    end do

    ! 1c. Distribute the Gamma-point band energies to EVERY rank ("gather all /
    !     разошлёт Г точку на всех"): rank 0 copies gs%eigen(:, ik_gamma) and
    !     broadcasts the nb-number reference vector. The ENTIRE window is then
    !     derived from this single, bit-identical vector, so the active-band set
    !     ("линии") is the same on every rank by construction -- no per-rank or
    !     per-k divergence, whatever the band complexity.
    allocate(eig_gamma(1:sbe%nb))
    if (irank == 0) eig_gamma(1:sbe%nb) = gs%eigen(1:sbe%nb, ik_gamma)
    call comm_bcast(eig_gamma, icomm, 0)

    ! gs%eigen is stored in atomic units (Hartree). Convert to eV here so the
    ! frozen-core window thresholds (frozen_core/free_threshold_ev, genuine eV
    ! inputs) are compared in eV as named. The dynamics always use gs%eigen in
    ! a.u. directly; only the active-band selection and the printout use eV.
    fermi_energy_ev = ((eig_gamma(homo_idx) + eig_gamma(lumo_idx)) * 0.5d0) * au_ev

    ! 2. Initialize active bands array
    allocate(sbe%is_active(1:sbe%nb))
    sbe%is_active = .false.
    sbe%n_active_bands = 0

    ! 3. Build the active-band mask from the broadcast Gamma reference on every
    !    rank. Identical input vector -> identical mask; the guard in 5a verifies.
    do ib = 1, sbe%nb
        eigen_ev = eig_gamma(ib) * au_ev
        ! Note: Ensure frozen_core_threshold_ev is negative if it represents a window below E_F
        if (eigen_ev > fermi_energy_ev + frozen_core_threshold_ev .and. &
            eigen_ev < fermi_energy_ev + frozen_free_threshold_ev) then
            sbe%is_active(ib) = .true.
        end if
    end do
    deallocate(eig_gamma)

    ! 3b. FULL-BASIS dressed projection (default, wiki/06 sec.6): a NARROWED frozen
    !     window truncates the ring's H_VG diagonalisation to the active subspace,
    !     dropping the A.p coupling to the frozen bands -- the field-dressing that
    !     belongs in the frozen bands then piles into the active conduction states
    !     and the dissipators over-scatter it into real carriers (measured x10^4-10^6
    !     over-generation at sub-gap THz; Si 4^3 all-active 2.2e18 vs frozen 1.2e22).
    !     The reversible VG unitary already runs on the full basis; the DRESSED
    !     PROJECTION the dissipators use must be truncation-free too. Default 'y'
    !     forces the ring to the full band basis (all bands active for dissipation);
    !     the narrowed frozen window is honoured only under the explicit opt-out
    !     yn_sbe_full_dressed='n' (fast, but over-generates at sub-gap fields).
    if (yn_sbe_full_dressed == 'y') then
        ! Cost-preserving: keep the correct FULL dressed basis (is_active = all), but
        ! if the user narrowed the frozen window, scatter only the dressed states in
        ! that ENERGY range -- the ring kernels skip out-of-window sources, restoring
        ! the O(nk^2 * n_active) cost without truncating the dressed projection.
        if (.not. all(sbe%is_active)) then
            sbe%ring_mask_on = .true.
            sbe%ring_e_lo = (fermi_energy_ev + frozen_core_threshold_ev) / au_ev
            sbe%ring_e_hi = (fermi_energy_ev + frozen_free_threshold_ev) / au_ev
            if (irank == 0) &
                write(*, '(a,f8.2,a,f8.2,a)') '# yn_sbe_full_dressed=y: FULL-basis dressed '// &
                    'projection; ring dissipates only the [', &
                    fermi_energy_ev + frozen_core_threshold_ev, ',', &
                    fermi_energy_ev + frozen_free_threshold_ev, '] eV window (cost-preserving).'
        end if
        sbe%is_active(:) = .true.
    end if

    ! 4. Derive the active-band count from the mask (single source of truth).
    sbe%n_active_bands = count(sbe%is_active)

    ! 5a. Cross-rank consistency guard (compiler-agnostic: only MPI_Allreduce +
    !     Fortran error stop). Every rank built the mask from its OWN replicated
    !     gs%eigen; if any rank's copy differs (a failed gs broadcast on a
    !     non-synchronized start), the counts diverge -- catch it here with a
    !     clean, collective abort instead of a lone segfault deep in the time
    !     loop. min == max  <=>  all ranks agree on the count.
    block
        integer :: n_hi, n_lo
        n_hi =  sbe%n_active_bands
        n_lo = -sbe%n_active_bands
        call comm_get_max(n_hi, icomm)
        call comm_get_max(n_lo, icomm)
        n_lo = -n_lo
        if (n_hi /= n_lo) then
            if (irank == 0) then
                write(*, '(a)')          '# ERROR: frozen-core active-band set differs across MPI ranks.'
                write(*, '(a,i0,a,i0,a)') '#        n_active_bands range over ranks = [', n_lo, ', ', n_hi, '].'
                write(*, '(a)')          '#        Cause: non-synchronized distributed start or a failed'
                write(*, '(a)')          '#        broadcast of gs%eigen. Re-run; if it recurs, check the'
                write(*, '(a)')          '#        launcher pins all ranks before the solver init.'
            end if
            call comm_sync_all(icomm)
            error stop 'frozen-core active-band set inconsistent across MPI ranks'
        end if
    end block

    ! 5b. A frozen window that selects fewer than two bands (no valence OR no
    !     conduction level) cannot support any dynamics: the mean active level
    !     spacing is zero (division), the ring has nothing to scatter, and the
    !     run dies downstream with an opaque "process killed". Abort here with an
    !     actionable message and the numbers that produced the empty window.
    if (sbe%n_active_bands < 2) then
        if (irank == 0) then
            write(*, '(a,i0,a)') '# ERROR: frozen-core window selected ', sbe%n_active_bands, &
                ' active band(s) -- need >= 2 (>= 1 valence + >= 1 conduction).'
            write(*, '(a,f10.4,a)') '#        Fermi (from k-point 1) = ', fermi_energy_ev, ' eV.'
            write(*, '(a,f10.4,a,f10.4,a)') '#        active window = [', &
                fermi_energy_ev + frozen_core_threshold_ev, ', ', &
                fermi_energy_ev + frozen_free_threshold_ev, '] eV.'
            write(*, '(a)') '#        Widen frozen_core_threshold_ev / frozen_free_threshold_ev, or check'
            write(*, '(a)') '#        that the GS *_eigen.data is the intended file (units / k-ordering).'
        end if
        call comm_sync_all(icomm)
        error stop 'frozen-core window selected < 2 active bands (empty/degenerate window)'
    end if

    ! 6. Build active_idx (1..n_active -> global band index) from the mask.
    if (sbe%n_active_bands > 0) then
        allocate(sbe%active_idx(sbe%n_active_bands))
        count_active = 0
        do ib = 1, sbe%nb
            if (sbe%is_active(ib)) then
                count_active = count_active + 1
                sbe%active_idx(count_active) = ib
            end if
        end do
        ! Belt-and-suspenders: the mask walk MUST consume exactly n_active_bands
        ! slots. If it does not, the mask was mutated between step 5 and here --
        ! abort rather than propagate an out-of-range active_idx.
        if (count_active /= sbe%n_active_bands) then
            write(*, '(a,i0,a,i0,a)') '# ERROR: active_idx fill (', count_active, &
                ') /= n_active_bands (', sbe%n_active_bands, ') -- corrupt frozen-core mask.'
            error stop 'frozen-core active_idx inconsistent with mask'
        end if
    else
        ! Modern Fortran handles zero-sized arrays natively.
        allocate(sbe%active_idx(0))
    end if

    ! =========================================================================
    ! k-local impact-ionization channel (optional, yn_sbe_impact_ionization)
    ! =========================================================================
    sbe%occ_max = merge(1d0, 2d0, yn_sbe_spinor == 'y')

    ! Active-subspace valence-branch count (gap edge: v1 = nv_act, c1 = nv_act+1
    ! in energy-ordered active indexing). Needed by BOTH the impact-ionization
    ! and the Auger channels, so compute it unconditionally here (not inside the
    ! impact block -- an Auger-only run must still have it set).
    sbe%nv_act = 0
    do ib = 1, sbe%n_active_bands
        if (sbe%active_idx(ib) <= homo_idx) sbe%nv_act = sbe%nv_act + 1
    end do

    ! Initial active-band occupations. Only stored when they are NOT the plain
    ! integer filling (a doped / finite-temperature sheet, wiki/12 sec. 4a): the
    ! dressed reference then measures carriers against this baseline instead of
    ! "nv full, the rest empty". Absent => byte-identical to the previous code.
    if (sbe%n_active_bands > 0) then
        block
            integer :: ika, iba
            real(8) :: dmax
            dmax = 0d0
            do ika = 1, gs%nk
                do iba = 1, sbe%n_active_bands
                    dmax = max(dmax, abs(gs%occup(sbe%active_idx(iba), ika) &
                                       - gs%occup_ref(sbe%active_idx(iba), ika)))
                end do
            end do
            if (dmax > 1d-12) then
                allocate(sbe%f0_act(sbe%n_active_bands, gs%nk))
                do ika = 1, gs%nk
                    do iba = 1, sbe%n_active_bands
                        sbe%f0_act(iba, ika) = gs%occup(sbe%active_idx(iba), ika)
                    end do
                end do
                if (lprint) write(*, '(a,es10.3,a)') &
                    '# dressed reference: doped/thermal baseline stored (max |f0 - integer filling| = ', &
                    dmax, ')'
            end if
        end block
    end if

    ! Super-mode ring flag -- set UNCONDITIONALLY and BEFORE every channel block:
    ! it gates the inter-k (ring) variants of e-ph, impact ionization and Auger,
    ! and the ring-vs-allgather choice inside Coulomb HF. (BUG FIXED 2026-07-02:
    ! this used to be set inside the Coulomb block only, so with yn_sbe_coulomb='n'
    ! the inter-k channels were SILENTLY OFF even with yn_sbe_superres='y' --
    ! runs fell back to the k-local / BZ-averaged paths without warning.)
    sbe%flag_ring = (yn_sbe_superres == 'y')

    sbe%flag_impact = (yn_sbe_impact_ionization == 'y')
    if (sbe%flag_impact .and. .not. mp%ii_ok) &
        call stop_forbidden_channel(epm_material, 'impact ionization (yn_sbe_impact_ionization)')
    if (sbe%flag_impact) then
        ! Fit form/exponent/prefactor/threshold default to the material registry
        ! ('auto' / negative sentinels); an explicit namelist value overrides it.
        ! GaAs Stobbe quartic (a=4) [Stobbe-Redmer-Schattke, PRB 49, 4494 (1994)];
        ! Si Keldysh quadratic (a=2) [Keldysh, JETP 21, 1135 (1965); Cartier APL
        ! 62, 3339 (1993)]; Si full-band a=4.6 [Kamakura JAP 75, 3500 (1994)].
        if (trim(sbe_ii_form) == 'auto') then
            if (.not. mp%found) call stop_unknown_material(epm_material, 'impact ionization (sbe_ii_form=auto)')
            ii_form_eff = mp%ii_form
        else
            ii_form_eff = sbe_ii_form
        end if
        ii_exp_eff  = merge(sbe_ii_exponent,     mp%ii_exponent,     sbe_ii_exponent     > 0d0)
        ii_pref_eff = merge(sbe_ii_prefactor,    mp%ii_prefactor,    sbe_ii_prefactor    > 0d0)
        ii_thr_eff  = merge(sbe_ii_threshold_ev, mp%ii_threshold_ev, sbe_ii_threshold_ev >= 0d0)
        ! A material whose registry prefactor is a sentinel (no cited value, e.g.
        ! CdS: the II prefactor is a fit parameter) requires an explicit
        ! sbe_ii_prefactor -- never invent or borrow one.
        if (ii_pref_eff <= 0d0) then
            write(*, '(a)') '# ERROR: impact ionization for material "'//trim(epm_material)// &
                            '" has no cited prefactor (it is a fit parameter).'
            write(*, '(a)') '#        Set sbe_ii_prefactor explicitly in &sbe.'
            error stop 'impact-ionization prefactor not cited for this material; set sbe_ii_prefactor'
        end if
        if ((sbe_ii_exponent <= 0d0 .or. sbe_ii_prefactor <= 0d0 .or. &
             sbe_ii_threshold_ev < 0d0) .and. .not. mp%found) &
            call stop_unknown_material(epm_material, 'impact ionization (sentinel default)')
        sbe%ii_exponent = ii_exp_eff
        ! Prefactor P [s^-1 eV^-a] -> [1/(Ha^a a.u.time)]:
        ! rate_au = P * t_au[s] * (dE[Ha] * au_ev)^a
        sbe%ii_pref_au = ii_pref_eff * (au_fs * 1d-15) * au_ev**ii_exp_eff
        sbe%ii_eth_au  = ii_thr_eff / au_ev
        sbe%ii_ramp_au = sbe_ii_ramp_ev / au_ev
        ! Global CBM of the field-free band structure (kinetic-energy zero of
        ! the Stobbe fit) and the gap lost by the primary electron per event.
        if (homo_idx + 1 > gs%nb) stop "impact ionization: no conduction bands"
        sbe%ii_ecbm_au = minval(gs%eigen(homo_idx + 1, :))
        sbe%ii_eg_au   = gs%eg_au
        ! (sbe%nv_act -- the gap-edge valence-branch count -- is set
        !  unconditionally above, before the channel blocks.)
        ! Sublattice resolution is enabled iff the unfolding weights were
        ! loaded (gs%unfold_w not all zero). When absent, the channel falls
        ! back to the original folded (single-pool) treatment.
        sbe%flag_unfold_ii = .false.
        if (allocated(gs%unfold_w)) then
            if (maxval(abs(gs%unfold_w)) > 1d-12) sbe%flag_unfold_ii = .true.
        end if
        if (lprint) then
            write(*, '(a)') '# impact ionization (k-local Lindblad, Stobbe fit) enabled:'
            write(*, '(a,ES12.5,a)') '#   P     = ', sbe%ii_pref_au, ' 1/(Ha^4 a.u.t)'
            write(*, '(a,ES12.5,a,ES12.5,a)') '#   E_th  = ', sbe%ii_eth_au, ' Ha, ramp = ', sbe%ii_ramp_au, ' Ha'
            write(*, '(a,ES12.5,a,ES12.5,a)') '#   E_CBM = ', sbe%ii_ecbm_au, ' Ha, E_g  = ', sbe%ii_eg_au, ' Ha'
            write(*, '(a,i4,a,i4)') '#   active valence branches = ', sbe%nv_act, ' / ', sbe%n_active_bands
            if (sbe%flag_unfold_ii) then
                write(*, '(a)') '#   sublattice-resolved (unfolding weights loaded): primitive'
                write(*, '(a)') '#   momentum conservation enforced over 4 FCC sublattices'
            else
                write(*, '(a)') '#   folded (single-pool) treatment: no unfolding weights found'
            end if
            if (sbe%nv_act < 1 .or. sbe%nv_act >= sbe%n_active_bands) &
                write(*, '(a)') '#   WARNING: no valence or no conduction branches active -- channel is inert'
        end if
        ! Nonlocal impact ionization (Part C4): partner sourced from the whole BZ
        sbe%flag_nl_ii = (yn_sbe_superres == 'y')
        if (sbe%n_active_bands > 0) allocate(sbe%glob_occ(sbe%n_active_bands))
        if (sbe%flag_nl_ii) sbe%glob_occ = 0d0
        if (sbe%flag_nl_ii .and. lprint) &
            write(*, '(a)') '#   nonlocal mode: valence partner drawn from the whole BZ (momentum exchange)'
    end if

    ! =========================================================================
    ! Coulomb (time-dependent Hartree-Fock / exchange) renormalization
    ! [Golde-Kira-Meier-Koch, Phys. Status Solidi B 248, 863 (2011)].
    ! Exchange self-energy Sigma_nm(k) = - sum_{q/=k} V(k-q) rho_nm(q) with the
    ! statically-screened Coulomb V(q) = strength * 4 pi / (eps Omega Nk (q^2+kappa^2)).
    ! Stored in the velocity-gauge Bloch basis and added to H_VG; the Houston
    ! basis (eigenbasis of H_VG + Sigma) the dissipators use becomes the
    ! Coulomb-renormalized one. Off in the no-Coulomb limit (yn_sbe_coulomb='n').
    ! =========================================================================
    sbe%flag_coulomb = (yn_sbe_coulomb == 'y')
    if (sbe%flag_coulomb .and. sbe_coulomb_epsilon <= 0d0 .and. .not. mp%coulomb_ok) &
        call stop_forbidden_channel(epm_material, 'Coulomb HF (yn_sbe_coulomb, no cited dielectric)')
    if (sbe%flag_coulomb) then
        ! Background dielectric defaults to the material registry (sentinel <=0);
        ! an explicit sbe_coulomb_epsilon overrides it.
        if (sbe_coulomb_epsilon > 0d0) then
            coul_eps_eff = sbe_coulomb_epsilon
        else if (mp%found) then
            coul_eps_eff = mp%eps0
        else
            call stop_unknown_material(epm_material, 'Coulomb HF (sbe_coulomb_epsilon default)')
        end if
        ! Discrete exchange-sum prefactor: continuum (1/(2pi)^3) int d^3q  ->
        ! (1/(Omega_cell Nk)) sum_grid; V(q) = 4 pi / (eps q^2) in a.u.
        sbe%coul_pref = sbe_coulomb_strength * 4d0 * pi &
                      / (max(coul_eps_eff, 1d-12) * gs%volume * dble(gs%nk))
        sbe%coul_screen2 = sbe_coulomb_screen_au**2
        ! A7: 2D-sheet kernel for slab materials (graphene): V_2D(q) =
        ! 2 pi/(eps A Nk (q+kappa)), q = the IN-PLANE |k-q| (b3 = the vacuum
        ! axis). Unlocks the Dirac-velocity / excitonic HF renormalization.
        sbe%coul_2d = mp%found .and. mp%coulomb_2d
        if (sbe%coul_2d) then
            block
                real(8) :: b3len, area2d
                b3len = sqrt(dot_product(gs%b_matrix(3,1:3), gs%b_matrix(3,1:3)))
                if (b3len < 1d-12 .or. &
                    abs(gs%b_matrix(3,1)) + abs(gs%b_matrix(3,2)) > 1d-8*b3len) then
                    write(*,'(a)') '# ERROR: 2D Sigma^HF needs a slab cell (b3 || z'
                    write(*,'(a)') '#        out-of-plane) and the b-matrix header.'
                    error stop 'A7: 2D Coulomb needs a 2D slab geometry'
                end if
                area2d = gs%volume * b3len / (2d0 * pi)
                sbe%coul_pref2d = sbe_coulomb_strength * 2d0 * pi &
                                / (max(coul_eps_eff, 1d-12) * area2d * dble(gs%nk))
                sbe%coul_screen1 = max(sbe_coulomb_screen_au, 0d0)
            end block
        end if
        ! sbe%icomm already set above (used by all collectives, not just Coulomb)
        if (sbe%n_active_bands > 0) &
            allocate(sbe%sigma_hf(1:sbe%n_active_bands, 1:sbe%n_active_bands, &
                                  sbe%ik_min:sbe%ik_max))
        ! Folding fix: project Sigma^HF block-diagonally onto the 4 FCC
        ! sublattices (zero the spurious inter-sublattice exchange that the
        ! cubic-cell folding creates). Requires the unfold weights and a
        ! folded cubic cell; otherwise inert. [Popescu-Zunger, PRB 85, 085201]
        sbe%flag_hf_subproj = (yn_sbe_hf_sublattice_proj == 'y')
        if (sbe%flag_hf_subproj) then
            if (.not. allocated(gs%unfold_w)) then
                sbe%flag_hf_subproj = .false.
            else if (maxval(abs(gs%unfold_w)) <= 1d-12) then
                sbe%flag_hf_subproj = .false.
            end if
        end if
        if (lprint) then
            write(*, '(a)') '# Coulomb HF (Golde-Kira-Meier-Koch SBE) enabled:'
            if (sbe%coul_2d) then
                write(*, '(a,ES12.5)') '#   A7 2D-SHEET kernel: 2pi*str/(eps*A*Nk) = ', sbe%coul_pref2d
            else
                write(*, '(a,ES12.5)') '#   exchange prefactor 4pi*str/(eps*Omega*Nk) = ', sbe%coul_pref
            end if
            write(*, '(a,f8.3,a,ES12.5,a)') '#   eps = ', coul_eps_eff, &
                ', screening kappa^2 = ', sbe%coul_screen2, ' 1/Bohr^2'
            write(*, '(a)') '#   NOTE: non-k-local mean field, O(Nk^2) per step (frozen over dt)'
            if (sbe%flag_hf_subproj) then
                write(*, '(a)') '#   sublattice projection ON: inter-sublattice exchange zeroed'
            else
                write(*, '(a)') '#   sublattice projection OFF (no unfold weights / disabled)'
            end if
        end if
    end if

    ! =========================================================================
    ! Coset block-diagonal projection of the FIELD coupling (momentum matrix p).
    ! A translationally invariant perturbation conserves primitive crystal
    ! momentum, so <coset s|p|coset s'> = 0 for s/=s'; the EPM eigenvector
    ! mixing at folded-valley degeneracies makes these spurious (here ~0.7x the
    ! intra-coset coupling), artificially hybridizing the valleys. Projecting p
    ! block-diagonal over the cosets (off-diagonal elements x sum_s w_s(i)w_s(j),
    ! same soft projector as the HF fix) keeps rho block-diagonal and restores
    ! the per-valley (primitive) Zener physics. Needs the unfold weights; inert
    ! otherwise. Applied to p_active in the propagator (hence H_VG, the Houston
    ! basis and the branch velocity) -- the core of the dynamics.
    sbe%flag_coset_proj = (yn_sbe_coset_proj == 'y')
    if (sbe%flag_coset_proj) then
        if (.not. allocated(gs%unfold_w)) then
            sbe%flag_coset_proj = .false.
        else if (maxval(abs(gs%unfold_w)) <= 1d-12) then
            sbe%flag_coset_proj = .false.
        end if
    end if
    if (lprint .and. sbe%flag_coset_proj) &
        write(*, '(a)') '# Coset projection ON: inter-coset momentum coupling p block-diagonalized (folding fix)'

    ! =========================================================================
    ! Population-relaxing electron-phonon Lindblad (Part C5, super-mode).
    ! k-local skeleton with a single effective optical phonon. Off by default.
    ! [Jacoboni-Reggiani RMP 55, 645 (1983); nu saturation: Meng et al.,
    !  PRB 91, 075201 (2015); Fischetti-Laux PRB 38, 9721 (1988)]
    ! =========================================================================
    sbe%flag_eph = (yn_sbe_eph == 'y')
    if (sbe%flag_eph .and. .not. mp%eph_ok) &
        call stop_forbidden_channel(epm_material, 'electron-phonon (yn_sbe_eph, no cited nu_sat)')
    if (sbe%flag_eph) then
        ! The phonon table and the nu_sat default both come from the material
        ! registry -- the channel cannot run without one.
        if (.not. mp%found) call stop_unknown_material(epm_material, 'electron-phonon (yn_sbe_eph)')
        if (yn_sbe_eph_acoustic == 'y') then
            ! A4: cited-constants gate + the grid-resolved acoustic q
            if (.not. mp%found .or. mp%eph_ac_xi_ev <= 0d0) &
                call stop_forbidden_channel(epm_material, &
                    'acoustic e-ph (yn_sbe_eph_acoustic, no cited Xi_d/c_s; CdS: piezo pending)')
            block
                real(8) :: qt, bl
                integer :: iid
                qt = huge(1d0)
                do iid = 1, 3
                    bl = sqrt(dot_product(gs%b_matrix(iid,1:3), gs%b_matrix(iid,1:3)))
                    if (num_kgrid(iid) > 1) qt = min(qt, bl / dble(num_kgrid(iid)))
                end do
                if (qt >= huge(1d0)) qt = 0d0
                call init_eph_phonon_table(sbe, mp, kB_au * sbe_eph_temperature_k, &
                                            ac_qtyp_au=qt, ac_xi_ev=sbe_eph_ac_xi_ev)
                sbe%eph_ip_ac = sbe%eph_nph          ! the appended acoustic mode
                if (lprint) then
                    write(*,'(a,f7.3,a)') &
                        '#   A4 acoustic mode appended: hw_ac = ', &
                        sbe%eph_hw(sbe%eph_nph)*au_ev*1d3, ' meV (grid-resolved q)'
                    if (sbe_eph_ac_xi_ev > 0d0) write(*,'(a,f6.2,a)') &
                        '#   A4 DP OVERRIDE from input: Xi_d = ', sbe_eph_ac_xi_ev, &
                        ' eV (substrate-dependent)'
                    write(*,'(a)') '#   A4 screening: TF factor [q/(q+q_TF)]^2 from the'
                    write(*,'(a)') '#   instantaneous carrier density (mandatory small-q cut).'
                end if
            end block
        else
            call init_eph_phonon_table(sbe, mp, kB_au * sbe_eph_temperature_k)
        end if
        ! Saturation rate (overall magnitude cap): material default if not set.
        if (sbe_eph_nu_sat > 0d0) then
            sbe%eph_nusat_au = sbe_eph_nu_sat * (au_fs * 1d-15)
        else
            sbe%eph_nusat_au = mp%eph_nu_sat_si * (au_fs * 1d-15)
        end if
        sbe%eph_eps0_au = sbe_eph_eps0_ev / au_ev
        sbe%eph_n       = sbe_eph_n
        sbe%eph_ib_scale = sbe_eph_interband_scale
        if (abs(sbe%eph_ib_scale - 1d0) > 1d-12 .and. lprint) &
            write(*, '(a,es10.3,a)') '#   interband (BTBT) rate scale = ', &
                sbe%eph_ib_scale, ' (calibration, sbe_eph_interband_scale)'
        if (sbe_search_sigma_e_ev > 0d0) then
            sbe%eph_sigma_au = sbe_search_sigma_e_ev / au_ev
        else
            sbe%eph_sigma_au = 0.2d0 / au_ev     ! grid-matched default (Stobbe 0.2 eV)
        end if
        if (homo_idx + 1 <= gs%nb) then
            sbe%eph_ecbm_au = minval(gs%eigen(homo_idx + 1, :))
        end if
        sbe%eph_evbm_au = maxval(gs%eigen(homo_idx, :))
        sbe%eph_numax_au = 2d0 * sbe%eph_nusat_au   ! peak-rate estimate for sub-cycling
        if (lprint) then
            write(*, '(a)') '# electron-phonon population relaxation (Part C5) enabled:'
            write(*, '(a,i2,a,ES12.5,a)') '#   ', sbe%eph_nph, &
                ' phonon modes; nu_sat = ', sbe%eph_nusat_au, ' 1/a.u.t'
            do ib = 1, sbe%eph_nph
                write(*, '(a,i2,a,f7.2,a,f6.3,a,f6.3)') '#   mode ', ib, ': hw = ', &
                    sbe%eph_hw(ib)*au_ev*1d3, ' meV, N_B = ', sbe%eph_nb(ib), &
                    ', weight = ', sbe%eph_wrel(ib)
            end do
            write(*, '(a)') '#   k-local skeleton; CPTP amplitude damping; toggle Kuhn-Zurek off'
        end if
    end if

    ! =========================================================================
    ! Collisional-memory (non-Markovian) dephasing of the e-ph ring gout
    ! (wiki/10 sec. 8.6, maintainer-approved 2026-07-20). The kernel lines are
    ! read from the cited phonon table just built above -- no new constants.
    ! =========================================================================
    ! =========================================================================
    ! Velocity-gauge f-sum rule (wiki/12 sec. 6a): measure on the GS data how much
    ! of the diamagnetic current A N_e/V the nb-band basis cancels (diagnostic
    ! banner). With yn_sbe_vg_sumrule = 'y' calc_current_bloch subtracts the
    ! adiabatic ground-state current of the same truncated H_k(A(t)) -- the
    ! pure-gauge identity, parameter-free and exact at every A (its linear
    ! limit is this eta N_e A/V).
    ! =========================================================================
    sbe%flag_vg_sumrule = (yn_sbe_vg_sumrule == 'y')
    block
        real(8) :: eta3(3), scap3(3), wsum
        call vg_sumrule_eta(gs%nk, sbe%nb, gs%eigen(1:sbe%nb, :), gs%p_tm_matrix(1:sbe%nb, 1:sbe%nb, :, :), &
                            gs%occup_ref(1:sbe%nb, :), gs%kweight, 1d-3, eta3, scap3)  ! undoped reference; skip |de| < 27 meV (the exact Dirac point)
        wsum = sum(gs%kweight)
        sbe%vg_eta = eta3;  sbe%vg_scap = scap3
        sbe%vg_ne_cell = sum(gs%kweight * sum(gs%occup_ref(1:sbe%nb, :), dim=1)) / max(wsum, 1d-300)
        if (lprint) then
            write(*, '(a,3f8.4,a,3f8.4,a)') '# VG f-sum rule (occupied bands, this nb): captured oscillator'// &
                ' strength S = (', scap3, ' ), uncancelled diamagnetic fraction eta = (', eta3, ' )'
            if (sbe%flag_vg_sumrule) then
                write(*, '(a)') '#   pure-gauge restoration ON (yn_sbe_vg_sumrule): J(t) -= adiabatic ground-state'// &
                    ' current of the same truncated H_k(A(t)) [= eta N_e A/V at linear order, exact beyond]'
            else if (maxval(abs(eta3(1:2))) > 0.02d0) then
                write(*, '(a)') '#   NOTE: |eta| > 2 % in-plane -> a spurious reactive current ~ A = E/omega'// &
                    ' (decisive at THz); raise nstate or set yn_sbe_vg_sumrule = ''y'''
            end if
        end if
    end block

    sbe%flag_colmem     = (yn_sbe_colmem == 'y')
    sbe%flag_colmem_pop = (yn_sbe_colmem_pop == 'y')
    sbe%flag_dressed_ref = (yn_sbe_dressed_ref == 'y')
    if (sbe%flag_dressed_ref .and. .not. sbe%flag_ring) then
        if (lprint) write(*, '(a)') '# ERROR: yn_sbe_dressed_ref applies to the ring channels' // &
            ' -- enable yn_sbe_superres.'
        error stop 'dressed-reference carrier measure requires the ring'
    end if
    if (sbe%flag_dressed_ref .and. lprint) &
        write(*, '(a)') '# Option A dressed-reference carrier measure enabled:' // &
            ' ring channels read f - delta0[U(A)] (wiki/10 sec. 3A/8.10)'
    if (sbe%flag_colmem .or. sbe%flag_colmem_pop) then
        if (.not. (sbe%flag_eph .and. sbe%flag_ring)) then
            if (lprint) write(*, '(a)') '# ERROR: yn_sbe_colmem/_pop ride the e-ph ring' // &
                ' -- enable yn_sbe_eph and yn_sbe_superres.'
            error stop 'collisional-memory dephasing requires eph + ring'
        end if
        if (trim(epm_material) == 'graphene') then
            ! 2D colmem analog (maintainer 2026-09-04, wiki/10 sec. 8.11): the
            ! 2026-07-20 exclusion is lifted for the COLLISIONAL memory (it was a
            ! guard, not physics -- Kuhn-Zurek single-particle dephasing stays
            ! forbidden). The e-ph sectors use the graphene phonon table (E2g/
            ! A1' + acoustic) exactly like every material; the Coulomb (Rana)
            ! SOURCE gets the 2D Dirac-plasmon memory line (flag_colmem_2d,
            ! set in the Rana block below once that channel is configured).
            if (lprint) write(*, '(a)') '# graphene: 2D collisional-memory analog (wiki/10 sec. 8.11)'
        end if
        if (sbe_colmem_tau_fs > 0d0) then
            colmem_tauc = sbe_colmem_tau_fs / au_fs
        else
            colmem_tauc = 1d0 / sbe%eph_sigma_au      ! hbar/sigma_E default
        end if
        sbe%colmem_tauc_au = colmem_tauc
        allocate(sbe%colmem_c(2 * sbe%eph_nph), sbe%colmem_mu(2 * sbe%eph_nph))
        call colmem_lines(sbe%eph_nph, sbe%eph_hw, sbe%eph_wrel, sbe%eph_nb, &
                          colmem_tauc, sbe%colmem_nl, sbe%colmem_c, sbe%colmem_mu)
        if (sbe%colmem_nl < 1) error stop 'collisional-memory: no kernel lines'
        if (sbe%flag_colmem) then
            allocate(sbe%zmem(sbe%n_active_bands, sbe%n_active_bands, sbe%colmem_nl, &
                              sbe%ik_min:sbe%ik_max))
            sbe%zmem = (0d0, 0d0)
        end if
        if (lprint) then
            if (sbe%flag_colmem) &
                write(*, '(a)') '# collisional-memory dephasing (non-Markovian e-ph gout) enabled:'
            if (sbe%flag_colmem_pop) &
                write(*, '(a)') '# collisional-memory POPULATION filter (ring kernels read the' // &
                    ' memory-filtered f) enabled:'
            write(*, '(a,i3,a,f8.3,a)') '#   ', sbe%colmem_nl, &
                ' kernel lines from the cited phonon table; tau_c = ', &
                colmem_tauc * au_fs, ' fs'
            write(*, '(a,f10.6,a,es10.3)') '#   Markov anchor R(0) = ', &
                colmem_response(sbe%colmem_nl, sbe%colmem_c, sbe%colmem_mu, 0d0), &
                ';  R(10 w_max) = ', colmem_response(sbe%colmem_nl, sbe%colmem_c, &
                sbe%colmem_mu, 10d0 * maxval(sbe%eph_hw))
        end if
    end if

    ! =========================================================================
    ! Carrier-carrier (e-e/e-h) thermalization channel (Part F). Rate scale only;
    ! the per-step (mu,T) Fermi-Dirac fit is done inside the channel.
    ! [rate scale: Goodnick-Lugli PRB 37, 2578; Fischetti-Laux PRB 38, 9721]
    ! =========================================================================
    sbe%flag_eeh = (yn_sbe_eeh == 'y')
    if (sbe%flag_eeh .and. sbe_eeh_nu_sat <= 0d0 .and. .not. mp%eeh_ok) &
        call stop_forbidden_channel(epm_material, 'carrier-carrier (yn_sbe_eeh, no cited rate)')
    if (sbe%flag_eeh) then
        if (sbe_eeh_nu_sat > 0d0) then
            sbe%eeh_nu_au = sbe_eeh_nu_sat * (au_fs * 1d-15)
        else
            sbe%eeh_nu_au = 1.0d14 * (au_fs * 1d-15)   ! 1e13-1e14 s^-1 scale
        end if
        if (lprint) then
            write(*, '(a)') '# carrier-carrier (e-e/e-h) thermalization (Part F) enabled:'
            write(*, '(a,ES12.5,a)') '#   nu_cc = ', sbe%eeh_nu_au, &
                ' 1/a.u.t; CPTP relax to Fermi-Dirac (conserves number+energy)'
        end if
    end if

    ! Ring virtual-transient gate (real-vs-virtual separation, see the struct
    ! comment). Auto time constant = 2*pi/Egap -- a Houston population must
    ! outlive the interband virtuality time before the population kernels may
    ! scatter it; sub-cycle LZ/dressing transients only dephase (gout/T2).
    sbe%ring_gate_tau_au = 0d0
    if (sbe_ring_gate_fs > 0d0) then
        sbe%ring_gate_tau_au = sbe_ring_gate_fs / au_fs
    else if (sbe_ring_gate_fs < 0d0 .and. homo_idx + 1 <= gs%nb) then
        block
            real(8) :: egap_au
            egap_au = minval(gs%eigen(homo_idx + 1, :)) - maxval(gs%eigen(homo_idx, :))
            if (egap_au > 1d-6) sbe%ring_gate_tau_au = 2d0 * pi / egap_au
        end block
    end if
    if (sbe%ring_gate_tau_au > 0d0 .and. lprint) &
        write(*, '(a,f8.3,a)') '# ring virtual-transient gate: tau = ', &
            sbe%ring_gate_tau_au * au_fs, ' fs (population kernels read the'// &
            ' persistent Houston floor; coherence/T2 rates unchanged;'// &
            ' sbe_ring_gate_fs=0 disables)'

    ! =========================================================================
    ! Auger recombination (Sec 13): density-gated, number-conserving CPTP
    ! channel. A conduction electron recombines with a valence hole and the
    ! released gap energy promotes a second conduction electron to a hot state
    ! (gap-edge mean-field closure). Per-carrier rate gamma = C n^2 (R = C n^3).
    ! Provenance-gated: the material must supply a cited C, or the user must set
    ! sbe_auger_c_cm6s explicitly. NO material currently ships a verified default
    ! C (the former CdS "Haury 1998" coefficient was fabricated and removed), so a
    ! plain yn_sbe_auger='y' run aborts unless sbe_auger_c_cm6s is given. Cited
    ! per-material coefficients (GaAs/Si/graphene) are the subject of the
    ! nonlocal-Auger task -- see wiki/07_nonlocal_auger.md.
    ! =========================================================================
    sbe%flag_auger = (yn_sbe_auger == 'y')
    sbe%flag_rana2d = sbe%flag_auger .and. mp%found .and. mp%auger_2d_rana
    if (sbe%flag_rana2d) then
        ! ---- 2D RANA Auger/CM (graphene [R07], wiki/07 sec.6): the cited
        ! gapless collinear CCCV/CVVV recombination + CVCC generation, applied
        ! as the net CPTP relaxation R - G of the Dirac-cone pair density on
        ! the instantaneous quasi-Fermi levels. NO C [cm^6/s] exists or is
        ! used; the rate comes from the R07 integrals (lifetime benchmarks
        ! unit-tested in test_rana_2d). Ring-gated like graphene e-ph: the
        ! quasi-Fermi inversion needs the GLOBAL gathered populations.
        if (.not. sbe%flag_ring) then
            write(*, '(a)') '# ERROR: the graphene 2D Rana Auger/CM channel needs the global'
            write(*, '(a)') '#        population gather: enable yn_sbe_superres=''y'' together'
            write(*, '(a)') '#        with yn_sbe_auger=''y'' (same ring gate as graphene e-ph).'
            error stop '2D Rana Auger requires the ring (yn_sbe_superres)'
        end if
        sbe%rana_vf_au = mp%rana_vf_au
        sbe%rana_eps_r = mp%rana_eps_r
        if (sbe_coulomb_epsilon > 0d0) sbe%rana_eps_r = sbe_coulomb_epsilon
        sbe%rana_kt_au = kB_au * sbe_eph_temperature_k
        ! 2D colmem analog (wiki/10 sec. 8.11): with the population memory
        ! filter on, the Rana SOURCE densities are filtered with the 2D
        ! Dirac-plasmon line instead of the phonon lines (apply_ring_channels).
        sbe%flag_colmem_2d = sbe%flag_colmem_pop
        ! two-temperature Coulomb sector (wiki/12 sec. 3): R - G at the carrier T_e
        sbe%flag_rana_te = (yn_sbe_rana_te == 'y')
        sbe%te_au = sbe%rana_kt_au
        ! 2D cell area A = V*|b3|/(2 pi): exact when b3 is the out-of-plane
        ! (vacuum) axis, i.e. b3 || z -- true for the graphene datasets.
        block
            real(8) :: b3len
            b3len = sqrt(dot_product(gs%b_matrix(3, 1:3), gs%b_matrix(3, 1:3)))
            if (b3len < 1d-12) then
                write(*, '(a)') '# ERROR: 2D Rana Auger needs the reciprocal vectors in the'
                write(*, '(a)') '#        GS k.data header (# b1/b2/b3) to compute the cell area.'
                error stop '2D Rana Auger: b_matrix missing from the GS dataset'
            end if
            if (abs(gs%b_matrix(3, 1)) + abs(gs%b_matrix(3, 2)) > 1d-8 * b3len) then
                write(*, '(a)') '# ERROR: 2D Rana Auger needs an out-of-plane third axis'
                write(*, '(a)') '#        (b3 || z, in-plane a1/a2 + vacuum a3).'
                error stop '2D Rana Auger: cell is not a 2D slab geometry'
            end if
            sbe%rana_area_au = gs%volume * b3len / (2d0 * pi)
        end block
        if (lprint) then
            write(*, '(a)') '# Auger/CM: graphene 2D RANA channel [R07] enabled (ring-gated):'
            write(*, '(a)') '#   net CPTP pair relaxation R - G on instantaneous quasi-Fermi'
            write(*, '(a)') '#   levels; R = G at equilibrium (detailed balance, unit-tested).'
            write(*, '(a,f8.4,a,f6.1,a,f8.2,a)') '#   v_F = ', sbe%rana_vf_au, &
                ' a.u., eps_r = ', sbe%rana_eps_r, ', T = ', &
                sbe%rana_kt_au / kB_au, ' K (e-ph bath)'
            write(*, '(a,f10.3,a)') '#   2D cell area = ', sbe%rana_area_au, ' a.u.^2'
            write(*, '(a)') '#   k-local C n^3 and gap-threshold ring Auger: OFF (2D branch).'
            if (sbe%flag_colmem_2d) write(*, '(a)') '#   2D colmem analog: R - G evaluated on n, p'// &
                ' memory-filtered with the 2D Dirac-plasmon line w_pl(n,p; Q_TF) (wiki/10 sec. 8.11)'
            if (sbe%flag_rana_te) write(*, '(a)') '#   two-temperature model: R - G, Q_TF and the plasmon'// &
                ' line at the carrier T_e fitted from (n, p, energy) of the cone populations;'// &
                ' lattice = e-ph bath T (cooling through the phonon channel; wiki/12 sec. 3)'
        end if
    else if (sbe%flag_auger .and. sbe%flag_ring) then
        ! ---- NONLOCAL (ring) Auger: the exact time-reverse of the nonlocal
        ! impact ionization -- same |M|^2 kernel (momentum map, screened |V(q)|^2,
        ! broadened energy delta, cited Stobbe/Keldysh magnitude), swapped Fermi
        ! factors [detailed balance; Rana 2007; Kioupakis 2015]. NO separate
        ! coefficient C is needed or used: the rate scale IS the II kernel's.
        ! It therefore REQUIRES the impact-ionization channel configured (the
        ! shared constants live in sbe%ii_*): enable yn_sbe_impact_ionization.
        if (.not. sbe%flag_impact) then
            write(*, '(a)') '# ERROR: the nonlocal (ring) Auger is the time-reverse of the'
            write(*, '(a)') '#        nonlocal impact ionization and shares its kernel constants.'
            write(*, '(a)') '#        Enable yn_sbe_impact_ionization=''y'' together with'
            write(*, '(a)') '#        yn_sbe_auger=''y'' when yn_sbe_superres=''y''.'
            error stop 'ring Auger requires the impact-ionization channel (shared kernel)'
        end if
        sbe%auger_eg_au = gs%eg_au
        if (lprint) then
            write(*, '(a)') '# Auger recombination: NONLOCAL inter-k (ring) mode enabled:'
            write(*, '(a)') '#   time-reverse of the nonlocal impact ionization (same |M|^2,'
            write(*, '(a)') '#   swapped Fermi factors, detailed balance) -- no separate C.'
            write(*, '(a)') '#   k-local C*n^3 channel gated OFF (no double count).'
        end if
    else if (sbe%flag_auger) then
        ! ---- k-LOCAL (C n^3) Auger: needs a verified coefficient. ----
        if (.not. mp%found) call stop_unknown_material(epm_material, 'Auger (yn_sbe_auger)')
        if (.not. mp%auger_ok .and. sbe_auger_c_cm6s <= 0d0) &
            call stop_forbidden_channel(epm_material, 'Auger recombination (yn_sbe_auger, no cited C)')
        ! C [cm^6/s]: explicit override, else the cited material default.
        if (sbe_auger_c_cm6s > 0d0) then
            sbe%auger_c_au = sbe_auger_c_cm6s * (au_fs * 1d-15)
        else
            sbe%auger_c_au = mp%auger_c_cm6s * (au_fs * 1d-15)
        end if
        if (sbe_auger_n_gate_cm3 > 0d0) then
            sbe%auger_n_gate_cm3 = sbe_auger_n_gate_cm3
        else
            sbe%auger_n_gate_cm3 = mp%auger_n_gate_cm3
        end if
        sbe%auger_eg_au = gs%eg_au
        if (homo_idx + 1 > gs%nb) stop "Auger: no conduction bands"
        if (lprint) then
            write(*, '(a)') '# Auger recombination (Sec 13, density-gated CPTP) enabled:'
            write(*, '(a,ES12.5,a,ES12.5,a)') '#   C = ', &
                sbe%auger_c_au / (au_fs * 1d-15), ' cm^6/s, n_gate = ', &
                sbe%auger_n_gate_cm3, ' cm^-3'
            write(*, '(a,ES12.5,a)') '#   E_g = ', sbe%auger_eg_au, &
                ' Ha; number-conserving (recombination + hot-carrier promotion)'
        end if
    end if

    ! =========================================================================
    ! Bandgap-renormalization-coupled impact-ionization threshold (Part C7).
    ! Bookkeeping (homo_idx, a.u.->cm^-3 density conversion) is always stored;
    ! the threshold only moves when yn_sbe_bgr_threshold='y' AND impact
    ! ionization is on AND the carrier density exceeds the gate. [Vashishta-Kalia]
    ! =========================================================================
    sbe%homo_idx    = homo_idx
    sbe%au_dens_cm3 = 1d24 / (0.52917721067d0)**3     ! a.u.^-3 -> cm^-3 (Bohr in Angstrom)
    sbe%ii_eth0_au  = sbe%ii_eth_au                   ! fixed reference threshold
    ! ---- approved-improvements parameter parsing (wiki/00 2026-07-04) ----
    sbe%ring_vq_floor = max(sbe_ring_vq_floor, 0d0)
    sbe%flag_ii_fk = (yn_sbe_ii_fk_soften == 'y')
    if (sbe%flag_ii_fk) then
        if (sbe_ii_fk_mu <= 0d0) then
            write(*,'(a)') '# ERROR: yn_sbe_ii_fk_soften needs an explicit reduced mass'
            write(*,'(a)') '#        sbe_ii_fk_mu > 0 (hbar*theta = (F^2/2mu)^(1/3)).'
            error stop 'A5: sbe_ii_fk_mu required'
        end if
        sbe%ii_fk_mu = sbe_ii_fk_mu
    end if
    sbe%ii_phassist = max(sbe_ii_phassist, 0d0)
    if (sbe%ii_phassist > 0d0 .and. sbe%eph_nph <= 0) then
        ! A1 sidebands need the cited phonon table even when e-ph is off
        if (.not. mp%found .or. mp%eph_nph <= 0) &
            call stop_forbidden_channel(epm_material, &
                'phonon-assisted II/Auger (sbe_ii_phassist, no cited phonon table)')
        call init_eph_phonon_table(sbe, mp, kB_au * sbe_eph_temperature_k)
    end if
    sbe%flag_ii_holes = (yn_sbe_ii_holes == 'y')
    if (sbe%flag_ii_holes) then
        if (.not. mp%found .or. mp%ii_cpcn <= 0d0) &
            call stop_forbidden_channel(epm_material, &
                'hole-initiated II/Auger (yn_sbe_ii_holes, no cited Cp/Cn)')
        if (.not. sbe%flag_ring .or. .not. (sbe%flag_impact .or. sbe%flag_auger)) then
            write(*,'(a)') '# ERROR: yn_sbe_ii_holes rides the nonlocal ring II/Auger:'
            write(*,'(a)') '#        needs yn_sbe_superres + impact_ionization (and/or auger).'
            error stop 'A2: hole channel needs the ring II/Auger'
        end if
        sbe%ii_cpcn = mp%ii_cpcn
        if (lprint) write(*,'(a,f7.3,a)') '# A2 hole-initiated II/Auger enabled: Cp/Cn = ', &
            sbe%ii_cpcn, ' (cited, registry)'
    end if
    ! C3: energy-conservation broadening vs the actual grid level spacing
    if ((sbe%flag_eph .or. sbe%flag_impact .or. sbe%flag_auger) .and. sbe%flag_ring) then
        block
            real(8) :: spac
            integer :: ib2
            spac = 0d0
            do ib2 = 1, sbe%n_active_bands - 1
                spac = spac + abs(gs%eigen(sbe%active_idx(ib2+1), 1) - gs%eigen(sbe%active_idx(ib2), 1))
            end do
            spac = spac / dble(max(sbe%n_active_bands - 1, 1))
            if (lprint .and. sbe%eph_sigma_au > 0d0) then
                if (sbe%eph_sigma_au < 0.1d0 * spac .or. sbe%eph_sigma_au > 5d0 * spac) &
                    write(*,'(a,es10.2,a,es10.2,a)') &
                        '# NOTE (C3): delta_sigma = ', sbe%eph_sigma_au, &
                        ' Ha vs mean active level spacing ', spac, &
                        ' Ha -- consider retuning sbe_search_sigma_e_ev to the grid.'
            end if
        end block
    end if

    sbe%flag_bgr    = (yn_sbe_bgr_threshold == 'y') .and. sbe%flag_impact
    ! GUARD -- BGR and Sigma^HF are mutually exclusive (wiki/07 Sec 0.2b): both
    ! renormalise the gap with carrier density. With Coulomb HF on, the Houston
    ! eigenvalues already carry the dynamic gap shrinkage (diagonal
    ! eps~ = eps - sum_q V_{k-q} f_q), so ALSO lowering the II threshold by the
    ! Vashishta-Kalia n^(1/3) law would count the same physics twice. BGR is the
    ! cheap STAND-IN for the HF shift when Coulomb is off -- never both.
    if (sbe%flag_bgr .and. sbe%flag_coulomb) then
        write(*, '(a)') '# ERROR: yn_sbe_bgr_threshold=''y'' together with yn_sbe_coulomb=''y'''
        write(*, '(a)') '#        double-counts the density-driven gap renormalisation:'
        write(*, '(a)') '#        Sigma^HF already shifts the Houston eigenvalues the impact-'
        write(*, '(a)') '#        ionization threshold is measured against. Use BGR only as the'
        write(*, '(a)') '#        stand-in when Coulomb HF is off (see wiki/07 Sec 0.2b).'
        error stop 'BGR threshold + Coulomb HF are mutually exclusive (gap double-count)'
    end if
    if (sbe%flag_bgr) then
        sbe%bgr_n_gate = sbe_bgr_n_gate
        sbe%bgr_coeff  = sbe_bgr_coeff
        if (lprint) then
            write(*, '(a)') '# BGR-coupled impact-ionization threshold (Part C7) enabled:'
            write(*, '(a,ES12.5,a,ES12.5,a)') '#   gate n = ', sbe%bgr_n_gate, &
                ' cm^-3, K = ', sbe%bgr_coeff, ' eV cm'
        end if
    end if

    ! 7. Diagnostic Print
    if (lprint) then
        write(*, '(a)') '=========================================='
        write(*, '(a)') ' SBE real-time (velocity gauge) -- run configuration'
        write(*, '(a, i6, a, i4, a, i4)') '   k-points = ', sbe%nk, &
            ',   bands = ', sbe%nb, ',   active = ', sbe%n_active_bands
        if (yn_sbe_spinor == 'y') then
            write(*, '(a, f4.1)') '   basis    = spinor (spin-orbit split),  occ/band = ', sbe%occ_max
        else
            write(*, '(a, f4.1)') '   basis    = scalar (no spin-orbit),      occ/band = ', sbe%occ_max
        end if
        if (gs%have_unfold) then
            write(*, '(a, i1, a)') '   cell     = FOLDED supercell (', gs%n_coset, &
                '-coset unfold map present)'
        else
            write(*, '(a)') '   cell     = PRIMITIVE (no folding, no unfold map)'
        end if
        write(*, '(a)') '------------------------------------------'
        write(*, '(a)') 'DIAGNOSTIC: Frozen Core Check'
        write(*, '(a, f8.2, a)') '  frozen_core_threshold_ev = ', frozen_core_threshold_ev, ' eV'
        write(*, '(a, f8.2, a)') '  frozen_free_threshold_ev = ', frozen_free_threshold_ev, ' eV'
        write(*, '(a, f12.4, a)') '  Fermi energy (eV)      = ', fermi_energy_ev, ' eV'
        write(*, '(a, i6, a, 3f8.4, a)') '  Gamma reference: k-point ', ik_gamma, &
            ' at (', gs%kpoint(1, ik_gamma), gs%kpoint(2, ik_gamma), gs%kpoint(3, ik_gamma), &
            ') [reduced]'
        write(*, '(a, i4, a, i4)') '  n_active_bands         = ', sbe%n_active_bands, ' / ', sbe%nb
        ! What 'active' means depends on yn_sbe_full_dressed. With the default
        ! 'y' the dressed PROJECTION is truncation-free -> is_active is forced
        ! all-true (n_active_bands = nb) and a narrowed frozen window survives
        ! ONLY as the ring SOURCE mask [ring_e_lo, ring_e_hi]. Printing just
        ! is_active then shows "active = T" for every band and hides the window
        ! entirely -- so report the source mask separately (it is the thing that
        ! actually skips work).
        if (sbe%ring_mask_on) then
            write(*, '(a)') '  yn_sbe_full_dressed = y  -> is_active is ALL bands by design'
            write(*, '(a)') '     (the dressed projection must be truncation-free; wiki/06 sec.6).'
            write(*, '(a,f8.2,a,f8.2,a)') '  ring SOURCE mask (cost only) = [', &
                sbe%ring_e_lo * au_ev, ',', sbe%ring_e_hi * au_ev, '] eV:'
            write(*, '(a,i4,a,i4,a)') '     dissipator sources are skipped outside it -- ', &
                count(gs%eigen(1:sbe%nb, ik_gamma) >= sbe%ring_e_lo .and. &
                      gs%eigen(1:sbe%nb, ik_gamma) <= sbe%ring_e_hi), ' / ', sbe%nb, &
                ' bands in-window at Gamma.'
            write(*, '(a)') '     NOTE: the runtime mask tests the INSTANTANEOUS Houston'
            write(*, '(a)') '     eigenvalue per k, so the Gamma count below is indicative only.'
        end if
        write(*, '(a)') '----------------------------------------'
        if (sbe%ring_mask_on) then
            write(*, '(a)') '  Band energies (at Gamma) relative to Fermi level'
            write(*, '(a)') '  [active = in the dressed projection; src = usable as a dissipator source]:'
        else
            write(*, '(a)') '  Band energies (at Gamma) relative to Fermi level:'
        end if

        do ib = 1, min(sbe%nb, 100)  ! Print first 100 bands
            eigen_ev = gs%eigen(ib, ik_gamma) * au_ev
            if (sbe%ring_mask_on) then
                write(*, '(a, i3, a, f10.4, a, f8.2, a, l1, a, l1)') &
                    '    Band ', ib, ': E = ', eigen_ev, ' eV, E-E_F = ', &
                    (eigen_ev - fermi_energy_ev), ' eV, active = ', sbe%is_active(ib), &
                    ', src = ', (gs%eigen(ib, ik_gamma) >= sbe%ring_e_lo .and. &
                                 gs%eigen(ib, ik_gamma) <= sbe%ring_e_hi)
            else
                write(*, '(a, i3, a, f10.4, a, f8.2, a, l1)') &
                    '    Band ', ib, ': E = ', eigen_ev, ' eV, E-E_F = ', &
                    (eigen_ev - fermi_energy_ev), ' eV, active = ', sbe%is_active(ib)
            end if
        end do
        
        if (sbe%nb > 100) write(*, '(a)') '    ... (more bands)'
        write(*, '(a)') '=========================================='
    end if

end subroutine

subroutine calc_current_bloch(sbe, gs, Ac, jmat, icomm)
    use eigen_lapack, only: eigen_zheev
    use sbe_superres_ssbe, only: vg_gs_current_k
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm
    integer :: ik, idir, ib, jb, nb
    complex(8) :: tmp1(1:3), tmp(1:3), v_mat(1:sbe%nb, 1:sbe%nb)
    complex(8) :: trace_val

    nb = sbe%nb
    tmp1 = 0d0

    !$omp parallel do default(shared) private(ik, idir, ib, jb, v_mat, trace_val) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            v_mat = gs%p_tm_matrix(:, :, idir, ik)
            do ib = 1, nb
                v_mat(ib, ib) = v_mat(ib, ib) + Ac(idir)
            end do
            if (sbe%flag_vnl_correction) then
                v_mat = v_mat + gs%rvnl_tm_matrix(:, :, idir, ik)
            endif

            trace_val = 0d0
            do ib = 1, nb
                do jb = 1, nb
                    trace_val = trace_val + v_mat(ib, jb) * sbe%rho(jb, ib, ik)
                end do
            end do
            tmp1(idir) = tmp1(idir) + gs%kweight(ik) * trace_val
        end do
    end do
    !$omp end parallel do

    call comm_summation(tmp1, tmp, 3, icomm)
    jmat(:) = real(tmp(:)) / (sum(gs%kweight) * gs%volume)

    ! Velocity-gauge pure-gauge (f-sum-rule) restoration, wiki/12 sec. 6a: subtract
    ! the adiabatic ground-state current of the SAME truncated H_k(A(t)) = eps + A.p,
    ! evaluated on the UNDOPED T -> 0 filling gs%occup_ref -- NOT on the doped/thermal
    ! occupation: the artifact belongs to the filled sea, while the partially filled
    ! band's adiabatic current IS the physical Drude response (wiki/12 sec. 4a).
    ! A uniform A is a pure gauge, so a complete basis gives zero here for every A;
    ! the nb-band basis leaves eta N_e A/V at linear order (banner at start-up)
    ! and a nonlinear remainder once A ~ the k-distance to the Dirac points --
    ! the linear static form over-corrects there and, with the self-consistent
    ! sheet field, an over-corrected (anti-inductive) sheet runs away. Exact for
    ! any population (Tr_k rho is conserved), one ZHEEV per k per call, no parameter.
    if (sbe%flag_vg_sumrule) then
        block
            complex(8) :: Hk(nb, nb), Wk(nb, nb), pk(nb, nb, 3)
            real(8) :: evk(nb), jgs(3)
            tmp1 = 0d0
            !$omp parallel do default(shared) private(ik, Hk, Wk, pk, evk, jgs) reduction(+:tmp1)
            do ik = sbe%ik_min, sbe%ik_max
                call vg_coupling_k(sbe, gs, ik, nb, pk)
                call build_HVG(nb, gs%eigen(1:nb, ik), pk, Ac, Hk)
                call eigen_zheev(Hk, evk, Wk)
                call vg_gs_current_k(nb, Wk, pk, Ac, gs%occup_ref(1:nb, ik), jgs)
                tmp1(:) = tmp1(:) + gs%kweight(ik) * jgs(:)
            end do
            !$omp end parallel do
            call comm_summation(tmp1, tmp, 3, icomm)
            jmat(:) = jmat(:) - real(tmp(:)) / (sum(gs%kweight) * gs%volume)
        end block
    end if
end subroutine calc_current_bloch


! The field-coupling momentum block the propagator uses at k (p_tm + optional
! nonlocal correction, optional coset projection) -- the pure-gauge reference in
! calc_current_bloch must be built from the SAME coupling as the dynamics.
subroutine vg_coupling_k(sbe, gs, ik, nb, pk)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: ik, nb
    complex(8), intent(out) :: pk(nb, nb, 3)
    integer :: i, j, idir, s
    real(8) :: pcoset
    pk(:, :, :) = gs%p_tm_matrix(1:nb, 1:nb, :, ik)
    if (sbe%flag_vnl_correction) pk(:, :, :) = pk(:, :, :) + gs%rvnl_tm_matrix(1:nb, 1:nb, :, ik)
    if (sbe%flag_coset_proj) then
        do idir = 1, 3
            do j = 1, nb
                do i = 1, nb
                    if (i == j) cycle
                    pcoset = 0d0
                    do s = 1, 4
                        pcoset = pcoset + gs%unfold_w(s, i, ik) * gs%unfold_w(s, j, ik)
                    end do
                    pk(i, j, idir) = pk(i, j, idir) * pcoset
                end do
            end do
        end do
    end if
end subroutine vg_coupling_k


function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication, only: comm_get_groupinfo, comm_summation
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    real(8) :: tr

    integer :: ik, ib
    real(8) :: tmp, tmp1

    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, nb_max
            tmp1 = tmp1 + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    tr = tmp / sum(gs%kweight)

    return
end function calc_trace


function calc_energy(sbe, gs, Ac, icomm) result(energy)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    real(8), intent(in) :: Ac(1:3)
    integer :: ik, ib, jb, idir
    real(8) :: tmp1, tmp, energy
    
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    tmp1 = tmp1 &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_mod_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            tmp1 = tmp1 &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
                & + 0.5 * dot_product(Ac, Ac) &
                & ) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    energy = tmp / sum(gs%kweight)

    return
end function calc_energy


!=============================================================================
! CF4 (commutator-free Magnus, 4th order) propagator on Gauss-Legendre nodes,
! composed via the Suzuki-Yoshida triple-jump for the unitary part, combined
! with a strictly CPTP Kuhn-Zurek/Caldeira-Leggett dephasing map applied
! through Strang splitting with an exact Hadamard/Gaussian kernel in the
! instantaneous (Houston) eigenbasis:
!
!   rho(t+h) = D(h/2) o [ S2(p1 h) o S2(p2 h) o S2(p1 h) ] o D(h/2) [rho(t)]
!
! IMPORTANT: the Suzuki-Yoshida composition wraps ONLY the unitary CF4
! sub-steps S2(.), never the dephasing map D(.). A unitary step run with a
! negative sub-step (p2 h < 0) is simply a unitary rotation run backwards in
! time -- always valid. A dephasing step run with tau < 0 would replace the
! Hadamard/Gaussian kernel exp[-lambda (X_a-X_b)^2 tau] by its reciprocal,
! exp[+lambda (X_a-X_b)^2 |tau|], which is not positive semi-definite (it
! fails the Schoenberg/Bochner criterion for an RBF kernel and the Schur
! product theorem for the Hadamard map), and would break completely positive
! trace preservation. Hence D(h/2) is applied twice, with tau=+h/2>0 each
! time, by Strang splitting around the (always-safe) unitary composition.
!=============================================================================

subroutine dt_evolve_bloch_cf4(sbe, gs, t_start, dt, Ac_begin, Ac_end)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: t_start  ! time at the beginning of the step, rho(t_start) -> rho(t_start+dt)
    real(8), intent(in) :: dt
    real(8), intent(in) :: Ac_begin(1:3)  ! external vector potential A(t_start)
    real(8), intent(in) :: Ac_end(1:3)    ! external vector potential A(t_start+dt)

    real(8) :: tau_sub(3), t_sub(3)
    real(8) :: t_node(2, 3), s_node
    real(8) :: Ac_node(1:3, 2, 3)
    integer :: isub
    real(8) :: pcoset, cross_damp

    integer :: ik, nb, nba, i, j, idir, in, im

    complex(8), allocatable :: p_active(:, :, :)
    complex(8), allocatable :: rho_a(:, :)
    complex(8), allocatable :: H1f(:, :), H2f(:, :), HVG(:, :)   ! H1f/H2f: FULL-basis VG unitary
    real(8),    allocatable :: eigen_active(:)
    real(8),    allocatable :: V_begin(:), V_end(:), X_a(:)
    real(8),    allocatable :: w_act_sub(:, :)   ! (4, nba) field-free sublattice weights of active bands
    complex(8), allocatable :: p_k_full(:, :, :)
    complex(8), allocatable :: rho_n_full(:, :)
    integer :: s

    nb  = sbe%nb
    nba = sbe%n_active_bands

    !-------------------------------------------------------------------------
    ! The external field is known only at the step endpoints (the analytic
    ! pulse in realtime_ssbe, or the macroscopic Maxwell field in
    ! multiscale_ssbe -- both callers supply Ac(t_start) and Ac(t_start+dt)).
    ! CF4(Gauss-Legendre)+Yoshida needs A at several intermediate sub-nodes;
    ! we obtain them by linear interpolation in time,
    !   A(t_start + s*dt) = (1-s) Ac_begin + s Ac_end,   s in [0,1],
    ! which is consistent with the existing multiscale convention (compare
    ! the "linear interpolation for A(t+dt/2)" used by the previous ETDRK4
    ! step) and strictly more accurate than the old approach of treating A as
    ! constant over the whole step. These nodes are identical for every
    ! k-point, so they are evaluated once before the OpenMP/k-point loop.
    !-------------------------------------------------------------------------
    tau_sub(1) = yoshida_p1 * dt
    tau_sub(2) = yoshida_p2 * dt
    tau_sub(3) = yoshida_p1 * dt

    t_sub(1) = t_start
    t_sub(2) = t_sub(1) + tau_sub(1)
    t_sub(3) = t_sub(2) + tau_sub(2)
    ! t_sub(3) + tau_sub(3) = t_start + dt, since p1 + p2 + p1 = 1

    do isub = 1, 3
        t_node(1, isub) = t_sub(isub) + cf4_c1 * tau_sub(isub)
        t_node(2, isub) = t_sub(isub) + cf4_c2 * tau_sub(isub)

        s_node = (t_node(1, isub) - t_start) / dt
        Ac_node(:, 1, isub) = (1d0 - s_node) * Ac_begin + s_node * Ac_end

        s_node = (t_node(2, isub) - t_start) / dt
        Ac_node(:, 2, isub) = (1d0 - s_node) * Ac_begin + s_node * Ac_end
    end do

    ! Coulomb HF mean field, frozen over this dt step: form the non-k-local
    ! exchange self-energy Sigma(k) from rho(t_start) BEFORE the (OpenMP) k-loop
    ! (it is an MPI all-gather + convolution, must run collectively once).
    if (sbe%flag_coulomb .and. nba > 0) call compute_coulomb_selfenergy(sbe, gs)

    ! BGR-coupled impact-ionization threshold (Part C7): shift E_th with the
    ! running excited-carrier density (global reduction; once per step).
    if (sbe%flag_bgr) call update_bgr_threshold(sbe, gs)

    ! Auger recombination (Sec 13): the per-carrier rate gamma = C n^2 needs the
    ! running excited-carrier density n(t); compute it once per step (same global
    ! reduction as BGR) and store it for apply_auger_recombination.
    ! (the 2D Rana branch computes its own sheet densities from the gather;
    !  the 3D volume density n_exc_cm3 feeds only the k-local C n^3 rate)
    if (sbe%flag_auger .and. .not. sbe%flag_rana2d) call update_excited_density(sbe, gs)

    ! Nonlocal impact ionization (Part C4): gather the BZ-averaged active-band
    ! occupations once per step (the valence partner is sourced from anywhere).
    if (sbe%flag_nl_ii .and. nba > 0) call gather_global_occupation(sbe, gs)

    !$omp parallel default(shared) &
    !$omp    private(ik, i, j, idir, in, im, isub, s, pcoset, cross_damp) &
    !$omp    private(p_active, rho_a, H1f, H2f, HVG, eigen_active, V_begin, V_end, X_a, w_act_sub) &
    !$omp    private(p_k_full, rho_n_full)

    if (nba > 0) then
        allocate(p_active(nba, nba, 3), rho_a(nba, nba), HVG(nba, nba))
        allocate(eigen_active(nba), V_begin(nba), V_end(nba), X_a(nba))
        allocate(w_act_sub(4, nba))
    end if
    allocate(p_k_full(nb, nb, 1:3), rho_n_full(nb, nb))
    allocate(H1f(nb, nb), H2f(nb, nb))   ! FULL-basis velocity-gauge unitary

    !$omp do
    do ik = sbe%ik_min, sbe%ik_max

        ! ---- FULL-basis field coupling + density matrix (the VG unitary basis)
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) &
            p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)
        rho_n_full(:, :) = sbe%rho(:, :, ik)

        ! Coset block-diagonal projection of the field coupling (folding fix):
        ! suppress the spurious inter-coset momentum matrix elements. Off-diagonal
        ! only; intra-band (i=j) velocity untouched. On the FULL basis.
        if (sbe%flag_coset_proj) then
            do idir = 1, 3
                do j = 1, nb
                    do i = 1, nb
                        if (i == j) cycle
                        pcoset = 0d0
                        do s = 1, 4
                            pcoset = pcoset + gs%unfold_w(s, i, ik) * gs%unfold_w(s, j, ik)
                        end do
                        p_k_full(i, j, idir) = p_k_full(i, j, idir) * pcoset
                    end do
                end do
            end do
        end if

        ! ---- active Houston window: quantities the dissipator needs ----------
        ! FROZEN CORE is a Houston/dissipator-cost reduction (fewer ZHEEV + a
        ! smaller ring O(nk^3 x nba)), NOT a velocity-gauge cutoff: the reversible
        ! VG unitary below runs on the FULL band basis so a large A(t) can push
        ! population up/down through the frozen bands and bring it back (VG basis
        ! sufficiency). The dissipators act only in this active window.
        if (nba > 0) then
            do idir = 1, 3
                do j = 1, nba
                    im = sbe%active_idx(j)
                    do i = 1, nba
                        in = sbe%active_idx(i)
                        p_active(i, j, idir) = p_k_full(in, im, idir)
                    end do
                end do
            end do
            do i = 1, nba
                eigen_active(i) = gs%eigen(sbe%active_idx(i), ik)
                X_a(i) = sbe%X_branch(sbe%active_idx(i), ik)
            end do
            if (sbe%flag_impact .and. sbe%flag_unfold_ii) then
                do i = 1, nba
                    do s = 1, 4
                        w_act_sub(s, i) = gs%unfold_w(s, sbe%active_idx(i), ik)
                    end do
                end do
            else if (sbe%flag_impact) then
                w_act_sub = 0d0
            end if
        end if
        V_begin = 0d0; V_end = 0d0

        ! ================= Strang split:  D(h/2) U(h) D(h/2) =================

        ! Step 1: D(h/2) -- dissipative half-step, ACTIVE Houston window only.
        ! Truncate rho to the active block, transform+dissipate inside
        ! houston_dissipate, glue the evolved block back into the full matrix.
        if (nba > 0 .and. (sbe%flag_decoh .or. sbe%flag_impact .or. &
                           sbe%flag_eph .or. sbe%flag_eeh .or. sbe%flag_auger)) then
            do j = 1, nba; im = sbe%active_idx(j)
                do i = 1, nba; in = sbe%active_idx(i)
                    rho_a(i, j) = rho_n_full(in, im)
                end do; end do
            call build_HVG(nba, eigen_active, p_active, Ac_begin, HVG)
            if (sbe%flag_coulomb) HVG = HVG + sbe%sigma_hf(:, :, ik)
            call houston_dissipate(sbe, nba, rho_a, HVG, p_active, Ac_begin, X_a, &
                                   0.5d0 * dt, V_begin, w_act_sub, cross_damp)
            do j = 1, nba; im = sbe%active_idx(j)
                do i = 1, nba; in = sbe%active_idx(i)
                    rho_n_full(in, im) = rho_a(i, j)
                end do; end do
            ! CP extension of the block dissipators to the frozen sector: the
            ! active<->frozen coherences must carry the same loss-Kraus factor
            ! as the block coherences, or rho loses PSD (|rho_af|^2 > f_a f_f)
            ! and the Houston diagonal goes negative on later steps.
            if (cross_damp < 1d0 .and. nba < nb) then
                do im = 1, nb
                    if (sbe%is_active(im)) cycle
                    do i = 1, nba; in = sbe%active_idx(i)
                        rho_n_full(in, im) = rho_n_full(in, im) * cross_damp
                        rho_n_full(im, in) = rho_n_full(im, in) * cross_damp
                    end do
                end do
            end if
        end if

        ! Step 2: S4 unitary = S2(p1 h) o S2(p2 h) o S2(p1 h) on the FULL basis.
        ! Each S2(tau) is a CF4 (two-exponential) commutator-free Magnus step on
        ! the two Gauss-Legendre nodes; the middle negative-tau Yoshida jump is an
        ! exact backward-time rotation. Sigma^HF (active block) is embedded into
        ! the full H_VG at the active positions.
        do isub = 1, 3
            call build_HVG(nb, gs%eigen(:, ik), p_k_full, Ac_node(:, 1, isub), H1f)
            call build_HVG(nb, gs%eigen(:, ik), p_k_full, Ac_node(:, 2, isub), H2f)
            if (sbe%flag_coulomb .and. nba > 0) then
                do j = 1, nba; im = sbe%active_idx(j)
                    do i = 1, nba; in = sbe%active_idx(i)
                        H1f(in, im) = H1f(in, im) + sbe%sigma_hf(i, j, ik)
                        H2f(in, im) = H2f(in, im) + sbe%sigma_hf(i, j, ik)
                    end do; end do
            end if
            call cf4_unitary_step(nb, rho_n_full, H1f, H2f, tau_sub(isub))
        end do

        ! Step 3: D(h/2) -- dissipative half-step (post-field), same window.
        if (nba > 0 .and. (sbe%flag_decoh .or. sbe%flag_impact .or. &
                           sbe%flag_eph .or. sbe%flag_eeh .or. sbe%flag_auger)) then
            do j = 1, nba; im = sbe%active_idx(j)
                do i = 1, nba; in = sbe%active_idx(i)
                    rho_a(i, j) = rho_n_full(in, im)
                end do; end do
            call build_HVG(nba, eigen_active, p_active, Ac_end, HVG)
            if (sbe%flag_coulomb) HVG = HVG + sbe%sigma_hf(:, :, ik)
            call houston_dissipate(sbe, nba, rho_a, HVG, p_active, Ac_end, X_a, &
                                   0.5d0 * dt, V_end, w_act_sub, cross_damp)
            do j = 1, nba; im = sbe%active_idx(j)
                do i = 1, nba; in = sbe%active_idx(i)
                    rho_n_full(in, im) = rho_a(i, j)
                end do; end do
            ! CP extension to the frozen sector (see the first D(h/2) above).
            if (cross_damp < 1d0 .and. nba < nb) then
                do im = 1, nb
                    if (sbe%is_active(im)) cycle
                    do i = 1, nba; in = sbe%active_idx(i)
                        rho_n_full(in, im) = rho_n_full(in, im) * cross_damp
                        rho_n_full(im, in) = rho_n_full(im, in) * cross_damp
                    end do
                end do
            end if
        end if

        ! Branch-position update via the midpoint (average endpoint) velocity --
        ! 4th-order-consistent with CF4 (both D steps read the pre-step X_a).
        if (nba > 0) then
            do i = 1, nba
                sbe%X_branch(sbe%active_idx(i), ik) = &
                    X_a(i) + 0.5d0 * (V_begin(i) + V_end(i)) * dt
            end do
        end if

        sbe%rho(:, :, ik) = rho_n_full(:, :)

        ! Hermiticity (numerical safeguard). NO freeze-reset of the inactive
        ! bands: they now evolve reversibly under the full-basis VG unitary (they
        ! hold the field-dressed virtual population and give it back); only the
        ! dissipators skip them. Population output still reads the gap-edge bands.
        do j = 1, nb; do i = 1, nb
            sbe%rho(i, j, ik) = 0.5d0 * (sbe%rho(i, j, ik) + conjg(sbe%rho(j, i, ik)))
        end do; end do

    end do
    !$omp end do

    if (nba > 0) then
        deallocate(p_active, rho_a, HVG, eigen_active, V_begin, V_end, X_a, w_act_sub)
    end if
    deallocate(p_k_full, rho_n_full, H1f, H2f)
    !$omp end parallel

    ! NOTE: a frozen band legitimately holds population -- the full-basis unitary
    ! lets carriers TUNNEL / field-couple from the active window up into the
    ! frozen bands and back (that reversible field-dressed population is exactly
    ! the basis-sufficiency the frozen scheme preserves; the current captures it
    ! and it correctly does NOT enter the active-window Sigma^HF / dissipators).
    ! So a nonzero deviation of a frozen band from its ground occupation is
    ! EXPECTED and is not an error -- no diagnostic is raised for it.

    ! Inter-k e-ph through the super-mode ring: once per step on the post-step
    ! B2: ALL nonlocal ring channels (inter-k e-ph; nonlocal II + its Auger
    ! time-reverse; graphene 2D Rana) through ONE shared Houston pass + gather
    ! -- one ZHEEV per k per step instead of up to three; every channel sees
    ! the same pre-step populations (Strang-consistent, first order in dt like
    ! before). MPI-collective -> outside the OpenMP region, every rank calls.
    ! |E(t)| = |Ac_end - Ac_begin|/dt feeds the A5 Franz-Keldysh threshold.
    if (sbe%flag_ring .and. (sbe%flag_eph .or. sbe%flag_impact .or. sbe%flag_auger)) &
        call apply_ring_channels(sbe, gs, Ac_end, &
             sqrt(dot_product(Ac_end - Ac_begin, Ac_end - Ac_begin)) / dt, dt)

end subroutine dt_evolve_bloch_cf4


! Part C4: BZ-averaged occupation of each active band (global reduction), used
! by the nonlocal impact ionization so the valence partner / Pauli factors are
! sourced from the whole BZ rather than the local k-point (momentum exchange).
subroutine gather_global_occupation(sbe, gs)
    use communication, only: comm_summation
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: nba, ik, a
    real(8), allocatable :: loc(:), glob(:)
    nba = sbe%n_active_bands
    allocate(loc(nba), glob(nba))
    loc = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do a = 1, nba
            loc(a) = loc(a) + real(sbe%rho(sbe%active_idx(a), sbe%active_idx(a), ik)) * gs%kweight(ik)
        end do
    end do
    call comm_summation(loc, glob, nba, sbe%icomm)
    sbe%glob_occ(1:nba) = glob(1:nba) / sum(gs%kweight)
    deallocate(loc, glob)
end subroutine gather_global_occupation


! Part C7: update the impact-ionization threshold from the running excited
! carrier density n(t). n = (excited electrons per cell, BZ-averaged) / V_cell,
! converted to cm^-3 (the /N_k is already in the kweight average -- see the
! density note in the wiki). Above the gate, E_th(t) = E_th0 - |K n^(1/3)|.
subroutine update_bgr_threshold(sbe, gs)
    use communication, only: comm_summation
    use sbe_superres_ssbe, only: bgr_gap_shift_ev
    use phys_constants, only: au_ev
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: ik, ib
    real(8) :: loc, glob, n_cm3, dEbgr_ev

    loc = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = sbe%homo_idx + 1, sbe%nb
            loc = loc + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    call comm_summation(loc, glob, sbe%icomm)
    n_cm3 = (glob / sum(gs%kweight)) / gs%volume * sbe%au_dens_cm3

    if (n_cm3 > sbe%bgr_n_gate) then
        dEbgr_ev = bgr_gap_shift_ev(n_cm3, sbe%bgr_coeff)        ! negative [eV]
        sbe%ii_eth_au = sbe%ii_eth0_au - abs(dEbgr_ev) / au_ev   ! threshold shrinks
    else
        sbe%ii_eth_au = sbe%ii_eth0_au
    end if
end subroutine update_bgr_threshold


! Running excited-carrier density n(t) [cm^-3] for the Auger rate (same global
! reduction / normalization as update_bgr_threshold). Stored in sbe%n_exc_cm3.
subroutine update_excited_density(sbe, gs)
    use communication, only: comm_summation
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: ik, ib
    real(8) :: loc, glob

    loc = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = sbe%homo_idx + 1, sbe%nb
            loc = loc + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    call comm_summation(loc, glob, sbe%icomm)
    sbe%n_exc_cm3 = (glob / sum(gs%kweight)) / gs%volume * sbe%au_dens_cm3
end subroutine update_excited_density


! Build the instantaneous velocity-gauge Hamiltonian in the active subspace:
!   H_VG(t) = diag(eigen) + A(t) . pi
subroutine build_HVG(nba, eigen_active, p_active, Ac, H)
    implicit none
    integer,    intent(in)  :: nba
    real(8),    intent(in)  :: eigen_active(nba)
    complex(8), intent(in)  :: p_active(nba, nba, 3)
    real(8),    intent(in)  :: Ac(3)
    complex(8), intent(out) :: H(nba, nba)
    integer :: i, idir

    H = Ac(1) * p_active(:, :, 1) + Ac(2) * p_active(:, :, 2) + Ac(3) * p_active(:, :, 3)
    do i = 1, nba
        H(i, i) = H(i, i) + eigen_active(i)
    end do
end subroutine build_HVG


!=============================================================================
! Time-dependent Hartree-Fock (exchange/Fock) Coulomb self-energy, in the
! velocity-gauge stationary-Bloch active basis [Golde-Kira-Meier-Koch, Phys.
! Status Solidi B 248, 863 (2011), Eqs. 4-5]:
!
!   Sigma_nm(k) = - sum_{q/=k} V(k-q) rho_nm(q),
!   V(p) = strength * 4 pi / ( eps Omega_cell Nk (|p|^2 + kappa^2) )   [a.u.]
!
! This is a NON-k-local mean field: it convolves the active-band density
! matrix with the screened Coulomb kernel over the whole BZ. The convolution
! is gauge-covariant under the uniform velocity-gauge shift k -> k - A(t)
! (A cancels in k - q), so it is evaluated directly on the grid-k rho. Adding
! Sigma to H_VG renormalizes the diagonal energies (eps~ = eps - sum_q V f_q)
! and the off-diagonal Rabi term (Omega = d.E + sum_q V p_q); the (1-f_e-f_h)
! Pauli factor then follows from the von Neumann commutator [H+Sigma, rho].
! Sigma is Hermitian (rho Hermitian, V real), so the propagation stays unitary.
!
! Equilibrium subtraction: the EPM/DFT bands carry NO explicit Coulomb
! exchange, so Sigma is built from the DEVIATION from the ground state,
! drho = rho - rho_0 (rho_0 = diag(gs%occup)). All of drho's pieces vanish at
! t=0 (no excited electrons/holes, no polarization), so Sigma(t=0)=0 and the
! equilibrium gap stays the EPM gap; only the dynamical (carrier-induced)
! renormalization is added -- exactly the f^e, f^h, p of Eqs. 4-5.
!
! MPI: rho is k-partitioned (ik_min:ik_max). We zero-pad the local block and
! comm_summation to obtain the full-BZ drho on every rank (an all-gather), then
! each rank forms Sigma for its local k. Cost O(Nk^2 nba^2) per call; evaluated
! once per dt step (the mean field is frozen over the Strang/CF4 sub-steps).
!=============================================================================
! Dispatcher: form Sigma^HF either via the all-gather (default) or the systolic
! ring (super-mode, Part D), then apply the sublattice-block projection (Part E).
subroutine compute_coulomb_selfenergy(sbe, gs)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    if (sbe%n_active_bands < 1) return
    if (sbe%flag_ring) then
        call compute_coulomb_selfenergy_ring(sbe, gs)
    else
        call compute_coulomb_selfenergy_allgather(sbe, gs)
    end if
    if (sbe%flag_hf_subproj) call apply_hf_sublattice_projection(sbe, gs)
end subroutine compute_coulomb_selfenergy


! Screened-Coulomb kernel V(k-q) [a.u.] with minimum-image |k-q| and the q=k
! self-term excluded (returns 0). Shared by the all-gather and the ring.
pure function coulomb_kernel(sbe, gs, ik, iq) result(vkq)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info),      intent(in) :: gs
    integer,                  intent(in) :: ik, iq
    real(8) :: vkq, dr(3), dkx, dky, dkz, q2
    dr(1) = gs%kpoint(1, ik) - gs%kpoint(1, iq)
    dr(2) = gs%kpoint(2, ik) - gs%kpoint(2, iq)
    dr(3) = gs%kpoint(3, ik) - gs%kpoint(3, iq)
    dr(1) = dr(1) - dnint(dr(1))
    dr(2) = dr(2) - dnint(dr(2))
    dr(3) = dr(3) - dnint(dr(3))
    dkx = dr(1)*gs%b_matrix(1,1) + dr(2)*gs%b_matrix(2,1) + dr(3)*gs%b_matrix(3,1)
    dky = dr(1)*gs%b_matrix(1,2) + dr(2)*gs%b_matrix(2,2) + dr(3)*gs%b_matrix(3,2)
    dkz = dr(1)*gs%b_matrix(1,3) + dr(2)*gs%b_matrix(2,3) + dr(3)*gs%b_matrix(3,3)
    if (sbe%coul_2d) then
        ! A7: 2D sheet -- IN-PLANE |q|, V ~ 1/(q + kappa) (linear 2D screening)
        q2 = dkx*dkx + dky*dky
        if (q2 < 1d-12) then
            vkq = 0d0
        else
            vkq = sbe%coul_pref2d / (sqrt(q2) + sbe%coul_screen1)
        end if
        return
    end if
    q2 = dkx*dkx + dky*dky + dkz*dkz
    if (q2 < 1d-12) then
        vkq = 0d0                          ! exclude q = k (the V(0) self term)
    else
        vkq = sbe%coul_pref / (q2 + sbe%coul_screen2)
    end if
end function coulomb_kernel


! All-gather implementation: zero-pad the local drho and comm_summation to the
! full-BZ drho on every rank, then form Sigma for the local k. Memory O(Nk).
subroutine compute_coulomb_selfenergy_allgather(sbe, gs)
    use communication, only: comm_summation
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: nba, ik, iq, i, j, in, im
    real(8) :: vkq
    complex(8), allocatable :: rho_loc(:, :, :), rho_all(:, :, :)

    nba = sbe%n_active_bands
    allocate(rho_loc(nba, nba, sbe%nk), rho_all(nba, nba, sbe%nk))
    rho_loc = (0d0, 0d0)
    do ik = sbe%ik_min, sbe%ik_max
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                rho_loc(i, j, ik) = sbe%rho(in, im, ik)
            end do
            rho_loc(j, j, ik) = rho_loc(j, j, ik) - cmplx(gs%occup(im, ik), 0d0, 8)
        end do
    end do
    call comm_summation(rho_loc, rho_all, nba * nba * sbe%nk, sbe%icomm)

    !$omp parallel do default(shared) schedule(dynamic) private(ik, iq, i, j, vkq)
    do ik = sbe%ik_min, sbe%ik_max
        sbe%sigma_hf(:, :, ik) = (0d0, 0d0)
        do iq = 1, sbe%nk
            vkq = coulomb_kernel(sbe, gs, ik, iq)
            if (vkq == 0d0) cycle
            do j = 1, nba
                do i = 1, nba
                    sbe%sigma_hf(i, j, ik) = sbe%sigma_hf(i, j, ik) - vkq * rho_all(i, j, iq)
                end do
            end do
        end do
    end do
    !$omp end parallel do
    deallocate(rho_loc, rho_all)
end subroutine compute_coulomb_selfenergy_allgather


! Ring/pipeline implementation (Part D): each rank holds its local drho block
! plus ONE transit buffer. Blocks circulate around the ring (comm_exchange =
! MPI_Sendrecv) in nproc steps; at each hop the rank contracts its local k
! against the transit q-block into Sigma. After nproc-1 hops every block has
! visited, completing the full O(Nk^2) sum. Memory O(Nk/P + one block) -- does
! NOT grow with P. Bit-identical to the all-gather (same kernel, same order of
! the q-blocks). One fused pass: extra nonlocal accumulators (II, e-ph, e-e)
! can be added here without new communication. [Plimpton JCP 117, 1 (1995)]
subroutine compute_coulomb_selfenergy_ring(sbe, gs)
    use communication, only: comm_exchange
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: nba, nproc, maxn, hop, src, q0, ntr, ik, jq, iq, i, j, in, im
    integer :: idest, isrc
    real(8) :: vkq
    complex(8), allocatable :: transit(:, :, :), recvbuf(:, :, :)

    nba = sbe%n_active_bands
    nproc = sbe%nproc
    maxn = 0
    do i = 0, nproc - 1
        maxn = max(maxn, sbe%itbl_max(i) - sbe%itbl_min(i) + 1)
    end do
    allocate(transit(nba, nba, maxn), recvbuf(nba, nba, maxn))

    ! local drho block packed into transit[1:nloc]
    transit = (0d0, 0d0)
    do ik = sbe%ik_min, sbe%ik_max
        jq = ik - sbe%ik_min + 1
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                transit(i, j, jq) = sbe%rho(in, im, ik)
            end do
            transit(j, j, jq) = transit(j, j, jq) - cmplx(gs%occup(im, ik), 0d0, 8)
        end do
    end do
    do ik = sbe%ik_min, sbe%ik_max
        sbe%sigma_hf(:, :, ik) = (0d0, 0d0)
    end do

    idest = mod(sbe%irank + 1, nproc)
    isrc  = mod(sbe%irank - 1 + nproc, nproc)
    src   = sbe%irank
    do hop = 0, nproc - 1
        q0  = sbe%itbl_min(src)
        ntr = sbe%itbl_max(src) - sbe%itbl_min(src) + 1
        !$omp parallel do default(shared) schedule(dynamic) private(ik, jq, iq, i, j, vkq)
        do ik = sbe%ik_min, sbe%ik_max
            do jq = 1, ntr
                iq = q0 + jq - 1
                vkq = coulomb_kernel(sbe, gs, ik, iq)
                if (vkq == 0d0) cycle
                do j = 1, nba
                    do i = 1, nba
                        sbe%sigma_hf(i, j, ik) = sbe%sigma_hf(i, j, ik) - vkq * transit(i, j, jq)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        if (hop < nproc - 1) then
            call comm_exchange(transit, idest, recvbuf, isrc, 1, sbe%icomm)
            transit = recvbuf
            src = mod(src - 1 + nproc, nproc)
        end if
    end do
    deallocate(transit, recvbuf)
end subroutine compute_coulomb_selfenergy_ring


! Part E sublattice-block projection of Sigma^HF (shared by all-gather & ring).
subroutine apply_hf_sublattice_projection(sbe, gs)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: nba, ik, i, j, in, im, s
    real(8) :: proj
    nba = sbe%n_active_bands
    !$omp parallel do default(shared) private(ik, i, j, in, im, s, proj)
    do ik = sbe%ik_min, sbe%ik_max
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                if (i == j) cycle
                in = sbe%active_idx(i)
                proj = 0d0
                do s = 1, 4
                    proj = proj + gs%unfold_w(s, in, ik) * gs%unfold_w(s, im, ik)
                end do
                sbe%sigma_hf(i, j, ik) = sbe%sigma_hf(i, j, ik) * proj
            end do
        end do
    end do
    !$omp end parallel do
end subroutine apply_hf_sublattice_projection


! Single CF4 (commutator-free Magnus, 4th order) sub-step of length tau,
! evaluated on the two Gauss-Legendre Hamiltonians H1=H(t+c1*tau), H2=H(t+c2*tau):
!   Omega1 = tau (alpha1 H1 + alpha2 H2),  Omega2 = tau (alpha2 H1 + alpha1 H2)
!   rho <- exp(-i Omega2) exp(-i Omega1) rho exp(+i Omega1) exp(+i Omega2)
! Implemented as two successive exact unitary rotations (each built from an
! eigendecomposition of the Hermitian generator, so no Pade/Krylov truncation
! error is introduced -- the propagator is exactly unitary to machine precision).
subroutine cf4_unitary_step(nba, rho, H1, H2, tau)
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: H1(nba, nba), H2(nba, nba)
    real(8),    intent(in)    :: tau
    complex(8) :: Omega(nba, nba)

    Omega = tau * (cf4_alpha1 * H1 + cf4_alpha2 * H2)
    call apply_unitary_rotation(nba, rho, Omega)

    Omega = tau * (cf4_alpha2 * H1 + cf4_alpha1 * H2)
    call apply_unitary_rotation(nba, rho, Omega)
end subroutine cf4_unitary_step


! Apply rho -> U rho U^dagger with U = exp(-i*Omega) for Hermitian Omega,
! computed exactly via eigendecomposition Omega = W diag(lambda) W^dagger:
!   U rho U^dagger = W [ exp(-i lambda_i) (W^dagger rho W)_ij exp(+i lambda_j) ] W^dagger
subroutine apply_unitary_rotation(nba, rho, Omega)
    use eigen_lapack, only: eigen_zheev
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: Omega(nba, nba)

    real(8)    :: evals(nba)
    complex(8) :: W(nba, nba), t1(nba, nba), t2(nba, nba)
    integer :: i, j

    call eigen_zheev(Omega, evals, W)

    call ZGEMM('C', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), W,  nba, rho, nba, cmplx(0d0, 0d0, 8), t1, nba)
    call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), t1, nba, W,   nba, cmplx(0d0, 0d0, 8), t2, nba)

    do j = 1, nba
        do i = 1, nba
            t2(i, j) = t2(i, j) * exp(cmplx(0d0, -(evals(i) - evals(j)), 8))
        end do
    end do

    call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), W,  nba, t2, nba, cmplx(0d0, 0d0, 8), t1, nba)
    call ZGEMM('N', 'C', nba, nba, nba, cmplx(1d0, 0d0, 8), t1, nba, W,  nba, cmplx(0d0, 0d0, 8), rho, nba)
end subroutine apply_unitary_rotation


! Strang dissipative half-step in the instantaneous Houston (adiabatic) basis.
! One ZHEEV serves both dissipative channels:
!   1) diagonalize the instantaneous H_VG(t) -> Houston basis U, {eps_a}
!   2) rotate rho~ = U^dagger rho U
!   3) [flag_decoh ] rho~_ab <- exp[-lambda (X_a - X_b)^2 * tau] * rho~_ab
!      (Hadamard product with a PSD Gram/RBF matrix => exactly CPTP, by
!      Schoenberg/Bochner positive-definiteness + the Schur product theorem)
!   4) [flag_impact] k-local impact-ionization jump channels (see
!      apply_impact_ionization below) -- each an exact amplitude-damping map,
!      CPTP for any tau >= 0
!   5) rotate back rho = U rho~ U^dagger
! Also returns the instantaneous branch (group) velocities in the field
! polarization direction, V_a = [(U^dagger pi U)_aa . e_hat] + (A . e_hat)
! (computed only when the Kuhn-Zurek dephasing needs the branch positions).
subroutine houston_dissipate(sbe, nba, rho, H, p_active, Ac, X, tau, V, w_act_sub, &
                             cross_damp)
    use eigen_lapack, only: eigen_zheev
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: H(nba, nba)
    complex(8), intent(in)    :: p_active(nba, nba, 3)
    real(8),    intent(in)    :: Ac(3)
    real(8),    intent(in)    :: X(nba)
    real(8),    intent(in)    :: tau
    real(8),    intent(out)   :: V(nba)
    real(8),    intent(in)    :: w_act_sub(4, nba)  ! field-free sublattice weights of the active bands
    ! Kraus factor the CALLER must apply to the active<->frozen coherence
    ! blocks (which this block routine cannot see). For the carrier-carrier
    ! convex-mix channel rho -> (1-a) rho + a diag(...) the CP extension to
    ! the frozen sector damps every active<->frozen coherence by the SAME
    ! (1-a) the within-block coherences get -- accumulated over sub-steps.
    ! Without it the frozen coherences outlive the block population they
    ! belong to, |rho_af|^2 > f_a*f_f, rho loses PSD and the Houston diagonal
    ! goes negative on later steps (the nelec/nhole < 0 pathology).
    real(8),    intent(out)   :: cross_damp

    real(8)    :: evals(nba), ehat(3), Ac_norm
    complex(8) :: W(nba, nba), t1(nba, nba), t2(nba, nba)
    real(8)    :: wsub_branch(4, nba), wcoef, tau_sub_d, alpha_cc
    integer :: i, j, idir, a, s, m_sub, isub_d

    cross_damp = 1d0
    call eigen_zheev(H, evals, W)

    ! rho~ = U^dagger rho U
    call ZGEMM('C', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), W,  nba, rho, nba, cmplx(0d0, 0d0, 8), t1, nba)
    call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), t1, nba, W,   nba, cmplx(0d0, 0d0, 8), t2, nba)

    ! Exact Hadamard/Gaussian dephasing kernel (PSD for tau >= 0)
    if (sbe%flag_decoh) then
        do j = 1, nba
            do i = 1, nba
                t2(i, j) = t2(i, j) * exp(-sbe%lambda_decoh * (X(i) - X(j))**2 * tau)
            end do
        end do
    end if

    ! k-local impact-ionization (threshold-gated) + population-relaxing e-ph.
    ! These amplitude-damping channels are applied via Strang sub-cycling
    ! (Part C8): split the half-step into m_sub CPTP sub-steps so that, when the
    ! collision rate is fast vs the step (nu_max*tau >~ 0.2), the operator-split
    ! error stays small and the per-sub-step re-read of the populations protects
    ! the Pauli factors (a built-in predictor-corrector). Each sub-step is CPTP,
    ! so positivity is never threatened. m_sub = 1 unless e-ph is active, so
    ! impact-ionization-only runs are byte-for-byte unchanged.
    if (sbe%flag_impact .or. sbe%flag_eph .or. sbe%flag_eeh .or. sbe%flag_auger) then
        ! field-aware Houston-branch sublattice weights for the unfolding-aware
        ! impact ionization (computed once; field-frozen over the half-step)
        if (sbe%flag_impact .and. sbe%flag_unfold_ii) then
            wsub_branch = 0d0
            do a = 1, nba
                do i = 1, nba
                    wcoef = real(W(i, a))**2 + aimag(W(i, a))**2
                    do s = 1, 4
                        wsub_branch(s, a) = wsub_branch(s, a) + wcoef * w_act_sub(s, i)
                    end do
                end do
            end do
        else
            wsub_branch = 0d0
        end if

        m_sub = 1
        ! e-ph sub-cycling only for the k-LOCAL channel; when the ring is on the
        ! inter-k e-ph runs once per step in apply_eph_interk_ring (outside).
        if (sbe%flag_eph .and. .not. sbe%flag_ring) &
            m_sub = max(m_sub, min(20, max(1, ceiling(10d0 * sbe%eph_numax_au * tau))))
        if (sbe%flag_eeh) &
            m_sub = max(m_sub, min(20, max(1, ceiling(10d0 * sbe%eeh_nu_au * tau))))
        if (sbe%flag_auger .and. .not. sbe%flag_ring) &
            m_sub = max(m_sub, min(20, max(1, ceiling(10d0 * &
                    sbe%auger_c_au * sbe%n_exc_cm3**2 * tau))))
        tau_sub_d = tau / dble(m_sub)
        do isub_d = 1, m_sub
            ! k-LOCAL impact ionization only when the ring is OFF; with the ring
            ! the momentum-conserving inter-k II runs once per step (outside).
            if (sbe%flag_impact .and. .not. sbe%flag_ring) &
                call apply_impact_ionization(sbe, nba, t2, evals, Ac, tau_sub_d, &
                                             wsub_branch, sbe%flag_unfold_ii)
            if (sbe%flag_eph .and. .not. sbe%flag_ring) &
                call apply_eph_relaxation(sbe, nba, t2, evals, Ac, tau_sub_d)
            if (sbe%flag_eeh) then
                call apply_carrier_carrier(sbe, nba, t2, evals, tau_sub_d, alpha_cc)
                ! (1-a): the same EID factor the within-block coherences get --
                ! the convex mix of identity and the (replace-block, keep-frozen)
                ! channel, both CPTP, so the composite is CPTP by construction.
                cross_damp = cross_damp * max(1d0 - alpha_cc, 0d0)
            end if
            ! k-LOCAL (C n^3) Auger only when the ring is OFF; with the ring the
            ! momentum-conserving inter-k Auger (the II time-reverse) runs once
            ! per step in apply_ii_interk_ring (outside).
            if (sbe%flag_auger .and. .not. sbe%flag_ring) &
                call apply_auger_recombination(sbe, nba, t2, evals, tau_sub_d)
        end do
    end if

    ! rho = U rho~ U^dagger
    call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), W,  nba, t2, nba, cmplx(0d0, 0d0, 8), t1, nba)
    call ZGEMM('N', 'C', nba, nba, nba, cmplx(1d0, 0d0, 8), t1, nba, W,  nba, cmplx(0d0, 0d0, 8), rho, nba)

    ! Branch velocities, projected on the polarization direction of A(t);
    ! only the Kuhn-Zurek branch positions consume them.
    V = 0d0
    if (sbe%flag_decoh) then
        Ac_norm = sqrt(dot_product(Ac, Ac))
        if (Ac_norm > 1.0d-12) then
            ehat = Ac / Ac_norm
        else
            ehat = (/ 1d0, 0d0, 0d0 /)
        end if

        do idir = 1, 3
            call ZGEMM('C', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), W,  nba, p_active(:, :, idir), nba, &
                       cmplx(0d0, 0d0, 8), t1, nba)
            call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), t1, nba, W, nba, cmplx(0d0, 0d0, 8), t2, nba)
            do i = 1, nba
                V(i) = V(i) + ehat(idir) * (real(t2(i, i)) + Ac(idir))
            end do
        end do
    end if
end subroutine houston_dissipate


!=============================================================================
! k-local impact ionization (Stobbe-Redmer-Schattke rate fit, PRB 49, 4494)
! in the Houston basis. The quartic two-particle event
!   A_h = sqrt(gamma_St) c+_h' c+_c1 c_v1 c_h
! is closed k-locally (no momentum transfer) and Hartree-Fock-factorized
! (two-particle closure, Rosati-Iotti-Dolcini-Rossi PRB 90, 125140) into two
! effective single-particle amplitude-damping channels with FROZEN scalar
! rates (partner population and Pauli blockers enter as factors clamped to
! [0,1], guaranteeing Gamma >= 0 => each map is exactly CPTP):
!   rel : branch h  -> h'  with Gamma_rel  = gamma_St * f_v1 (1-f_c1)(1-f_h')
!   pair: branch v1 -> c1  with Gamma_pair = gamma_St * f_h  (1-f_c1)(1-f_h')
! "Cold pair": secondaries are born at the band-edge branches v1/c1 with no
! kinetic energy; the primary drops to the conduction branch closest to
! eps_h - E_g. The hot set is gated on the kinetic energy from the global CBM,
!   eps_kin = eps_h + |A|^2/2 - E_CBM
! (the A^2/2 scalar dropped in H_VG is restored here, where the comparison
! against the field-free constant E_CBM requires it; by the Houston identity
! this equals E_h(k+A) - E_CBM -- the scale on which the Stobbe fit is
! defined). For most k-points and times the gate is empty and the channel
! costs O(N_C) comparisons -- the "rare impact events" mechanism.
!
! SUBLATTICE (BAND-UNFOLDING) RESOLUTION (use_unfold = .true.)
! ------------------------------------------------------------
! Under cubic folding, one grid k hides 4 distinct primitive crystal
! momenta -- the 4 FCC sublattices. Genuine impact ionization conserves
! the PRIMITIVE crystal momentum, so the strict k-local event must close
! WITHIN a single sublattice: the primary (h -> h') and the secondary pair
! (v1 -> c1) must all sit in the same primitive-cell sector. The folded
! treatment above ignores this and lets a primary on sublattice s create
! its secondary pair at the GLOBAL band edges v1/c1 (which may belong to a
! different sublattice) -- a momentum-non-conserving "false generation"
! event that inflates the carrier multiplication rate.
!
! With the spectral weights w_s(a) (sum_s w_s(a) = 1, projected onto the
! Houston branches) we resolve the event per sublattice s = 1..4:
!   * the secondary pair is created at the s-RESOLVED band edges
!     v1(s), c1(s) (the valence/conduction branches whose weight is
!     dominant on s, taken at the gap edge),
!   * the primary's receiving branch h'(s) is the conduction branch with
!     weight on s closest to eps_h - E_g,
!   * both channel rates are scaled by w_s(h) -- the primary's fraction on
!     s -- so that summing over s reproduces the primary's total
!     relaxation rate. In the no-folding limit (one sublattice, w = 1)
!     v1(s)=v1, c1(s)=c1, h'(s)=h' and the scheme reduces EXACTLY to the
!     folded channel above.
! Every sub-channel is still an exact CPTP amplitude-damping map (rates
! >= 0 by construction), so the composition remains CPTP.
!=============================================================================
subroutine apply_impact_ionization(sbe, nba, rho_ad, evals, Ac, tau, wsub, use_unfold)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho_ad(nba, nba)   ! adiabatic-basis rho~
    real(8),    intent(in)    :: evals(nba)         ! Houston branch energies
    real(8),    intent(in)    :: Ac(3)
    real(8),    intent(in)    :: tau
    real(8),    intent(in)    :: wsub(4, nba)       ! per-branch sublattice weights (sum_s = 1)
    logical,    intent(in)    :: use_unfold

    integer :: iv1, ic1, ih, ihp, a, s
    integer :: iv1_sub(4), ic1_sub(4), ihp_s
    real(8) :: a2half, ekin, d, gam, etgt, f, ws, wbest
    real(8) :: pv1, phh, bc1, bhp, g_rel, g_pair
    real(8), parameter :: occ_eps = 1d-12
    real(8), parameter :: w_min   = 0.05d0   ! ignore sublattices the primary barely touches
    real(8), parameter :: w_tie   = 1d-12

    if (sbe%nv_act < 1 .or. sbe%nv_act >= nba) return
    iv1 = sbe%nv_act        ! topmost valence branch (energy-ordered)
    ic1 = sbe%nv_act + 1    ! lowest conduction branch
    f = sbe%occ_max
    a2half = 0.5d0 * dot_product(Ac, Ac)

    !-------------------------------------------------------------------------
    ! Folded (single-pool) treatment: no unfolding weights available
    !-------------------------------------------------------------------------
    if (.not. use_unfold) then
        do ih = ic1, nba
            ekin = evals(ih) - sbe%ii_ecbm_au   ! a2half NOT restored (cancels vs shifted CBM)
            d = ekin - sbe%ii_eth_au
            if (d <= 0d0) cycle
            if (real(rho_ad(ih, ih)) < occ_eps) cycle
            gam = sbe%ii_pref_au * d**sbe%ii_exponent
            if (sbe%ii_ramp_au > 0d0 .and. d < sbe%ii_ramp_au) gam = gam * d / sbe%ii_ramp_au
            etgt = evals(ih) - sbe%ii_eg_au
            ihp = ic1
            do a = ic1, nba
                if (abs(evals(a) - etgt) < abs(evals(ihp) - etgt)) ihp = a
            end do
            if (ihp == ih) cycle
            ! valence partner + conduction blocker: global BZ average in the
            ! nonlocal mode (Part C4, momentum exchange), else local-k.
            if (sbe%flag_nl_ii) then
                pv1 = min(max(sbe%glob_occ(iv1) / f, 0d0), 1d0)
                bc1 = min(max(1d0 - sbe%glob_occ(ic1) / f, 0d0), 1d0)
            else
                pv1 = min(max(real(rho_ad(iv1, iv1)) / f, 0d0), 1d0)
                bc1 = min(max(1d0 - real(rho_ad(ic1, ic1)) / f, 0d0), 1d0)
            end if
            phh = min(max(real(rho_ad(ih,  ih )) / f, 0d0), 1d0)
            bhp = min(max(1d0 - real(rho_ad(ihp, ihp)) / f, 0d0), 1d0)
            g_rel  = gam * pv1 * bc1 * bhp
            g_pair = gam * phh * bc1 * bhp
            call apply_damping_channel(nba, rho_ad, ih,  ihp, g_rel,  tau)
            call apply_damping_channel(nba, rho_ad, iv1, ic1, g_pair, tau)
        end do
        return
    end if

    !-------------------------------------------------------------------------
    ! Sublattice-resolved (primitive-momentum-conserving) treatment.
    ! Precompute the s-resolved band edges once per k:
    !   v1(s) = valence branch dominant on s, taken at the TOP of the valence
    !   c1(s) = conduction branch dominant on s, taken at the BOTTOM of the
    !           conduction. The energy tie-break recovers the global edges
    !           (nv_act, nv_act+1) when all weights are equal (no folding).
    !-------------------------------------------------------------------------
    do s = 1, 4
        iv1_sub(s) = iv1
        wbest = -1d0
        do a = 1, sbe%nv_act
            if (wsub(s, a) > wbest + w_tie .or. &
                (abs(wsub(s, a) - wbest) <= w_tie .and. evals(a) > evals(iv1_sub(s)))) then
                wbest = wsub(s, a); iv1_sub(s) = a
            end if
        end do
        ic1_sub(s) = ic1
        wbest = -1d0
        do a = ic1, nba
            if (wsub(s, a) > wbest + w_tie .or. &
                (abs(wsub(s, a) - wbest) <= w_tie .and. evals(a) < evals(ic1_sub(s)))) then
                wbest = wsub(s, a); ic1_sub(s) = a
            end if
        end do
    end do

    do ih = ic1, nba
        ! Threshold gate on the kinetic energy from the field-free CBM
        ekin = evals(ih) - sbe%ii_ecbm_au   ! a2half NOT restored (cancels vs shifted CBM)
        d = ekin - sbe%ii_eth_au
        if (d <= 0d0) cycle
        if (real(rho_ad(ih, ih)) < occ_eps) cycle

        ! Stobbe rate, with optional linear ramp over the fit resolution
        gam = sbe%ii_pref_au * d**sbe%ii_exponent
        if (sbe%ii_ramp_au > 0d0 .and. d < sbe%ii_ramp_au) gam = gam * d / sbe%ii_ramp_au

        etgt = evals(ih) - sbe%ii_eg_au

        ! Distribute the event across the sublattices the primary occupies
        do s = 1, 4
            ws = wsub(s, ih)
            if (ws < w_min) cycle

            ! Receiving branch of the primary on sublattice s: conduction
            ! branch with weight on s closest to eps_h - E_g
            ihp_s = ic1_sub(s)
            do a = ic1, nba
                if (wsub(s, a) < w_min) cycle
                if (abs(evals(a) - etgt) < abs(evals(ihp_s) - etgt)) ihp_s = a
            end do
            if (ihp_s == ih) cycle   ! no receiving state below: event suppressed

            ! Frozen scalar factors at the s-resolved branches, clamped to
            ! [0,1]; the leading w_s(h) splits the primary across sublattices
            ! (sum_s w_s(h) = 1 => total primary relaxation preserved)
            if (sbe%flag_nl_ii) then
                pv1 = min(max(sbe%glob_occ(iv1_sub(s)) / f, 0d0), 1d0)
                bc1 = min(max(1d0 - sbe%glob_occ(ic1_sub(s)) / f, 0d0), 1d0)
            else
                pv1 = min(max(real(rho_ad(iv1_sub(s), iv1_sub(s))) / f, 0d0), 1d0)
                bc1 = min(max(1d0 - real(rho_ad(ic1_sub(s), ic1_sub(s))) / f, 0d0), 1d0)
            end if
            phh = min(max(real(rho_ad(ih,         ih        )) / f, 0d0), 1d0)
            bhp = min(max(1d0 - real(rho_ad(ihp_s,      ihp_s     )) / f, 0d0), 1d0)
            g_rel  = gam * ws * pv1 * bc1 * bhp
            g_pair = gam * ws * phh * bc1 * bhp

            ! Sequential exact amplitude-damping maps (each CPTP for tau>=0)
            call apply_damping_channel(nba, rho_ad, ih,         ihp_s,      g_rel,  tau)
            call apply_damping_channel(nba, rho_ad, iv1_sub(s), ic1_sub(s), g_pair, tau)
        end do
    end do
end subroutine apply_impact_ionization


!=============================================================================
! Population-relaxing electron-phonon Lindblad (Part C5), k-local skeleton, in
! the Houston/adiabatic basis. Unlike Kuhn-Zurek pure dephasing this channel
! RELAXES populations (Gamma_aa != 0) -- it cools hot carriers toward the band
! edges and is what reproduces THz bleaching (collision-rate saturation).
!
! Model (single effective optical phonon hw): each adiabatic level a with a
! carrier present relaxes to the energy-matched partner level b ~ eps_a -/+ hw:
!   emission   a -> b (b below, eps_a - eps_b ~ hw): rate nu(eps_a)(N_B+1)
!   absorption a -> b (b above, eps_b - eps_a ~ hw): rate nu(eps_a) N_B
! weighted by a Gaussian energy-conservation shape exp(-dE^2/2 sigma^2) and the
! target Pauli blocker (1 - rho_bb/f) clamped to [0,1]. nu(eps) is the smooth
! saturating collision rate; eps = carrier kinetic energy from the nearest band
! edge (electron above CBM or hole below VBM). Each transfer is the exact CPTP
! amplitude-damping map amp_damp_channel -> the whole channel is CPTP.
! [Jacoboni-Reggiani RMP 55, 645 (1983); nu sat: Meng et al. PRB 91, 075201]
!=============================================================================
subroutine apply_eph_relaxation(sbe, nba, rho_ad, evals, Ac, tau)
    use sbe_superres_ssbe, only: nu_saturation, gaussian_shape, amp_damp_channel, &
                                 eph_thermal_split, eph_thermal_split_de
    implicit none
    type(s_sbe_bloch_solver), intent(in)    :: sbe
    integer,                  intent(in)    :: nba
    complex(8),               intent(inout) :: rho_ad(nba, nba)
    real(8),                  intent(in)    :: evals(nba)
    real(8),                  intent(in)    :: Ac(3)
    real(8),                  intent(in)    :: tau

    integer :: ia, ib, b_em, b_ab, ip
    real(8) :: f, a2half, eps_kin, nu, nu_p, hw, sig, ekin, fe, fa
    real(8) :: best_em, best_ab, dE, shp, gam, blk
    real(8), parameter :: occ_eps = 1d-12

    f = sbe%occ_max
    sig = sbe%eph_sigma_au
    a2half = 0.5d0 * dot_product(Ac, Ac)

    do ia = 1, nba
        if (real(rho_ad(ia, ia)) < occ_eps) cycle
        ! carrier kinetic energy from the nearest band edge (electron or hole);
        ! restore the dropped A^2/2 (Houston identity), as in impact ionization.
        ekin = evals(ia)   ! a2half NOT restored (cancels vs shifted band edge)
        eps_kin = max(ekin - sbe%eph_ecbm_au, sbe%eph_evbm_au - ekin, 0d0)
        ! saturating collision rate = total magnitude cap for this level
        nu = nu_saturation(eps_kin, sbe%eph_nusat_au, sbe%eph_eps0_au, sbe%eph_n)
        if (nu * tau < 1d-14) cycle

        ! sum over phonon modes: each mode p relaxes the carrier to the partner
        ! level energy-matched to eps_a -/+ hw_p, split into emission/absorption
        ! (detailed balance), weighted by the mode weight w_p (sum_p w_p = 1) so
        ! the total channel rate stays ~ nu(eps).
        do ip = 1, sbe%eph_nph
            hw = sbe%eph_hw(ip)
            ! best energy-matched emission (below) and absorption (above)
            b_em = 0; best_em = huge(1d0)
            b_ab = 0; best_ab = huge(1d0)
            do ib = 1, nba
                if (ib == ia) cycle
                if (evals(ib) < evals(ia)) then
                    dE = abs((evals(ia) - evals(ib)) - hw)
                    if (dE < best_em) then; best_em = dE; b_em = ib; end if
                else
                    dE = abs((evals(ib) - evals(ia)) - hw)
                    if (dE < best_ab) then; best_ab = dE; b_ab = ib; end if
                end if
            end do

            ! Detailed balance on the REALIZED transfer and a pair-symmetric nu --
            ! the same two corrections as the inter-k twin (eph_thermal_split_de).
            if (.not. sbe%eph_db_realized) call eph_thermal_split(sbe%eph_nb(ip), fe, fa)
            if (b_em > 0) then
                shp = gaussian_shape(best_em, sig)
                blk = min(max(1d0 - real(rho_ad(b_em, b_em)) / f, 0d0), 1d0)
                nu_p = nu
                if (sbe%eph_db_realized) then
                    call eph_thermal_split_de(evals(ia) - evals(b_em), sbe%eph_kt_au, fe, fa)
                    nu_p = pair_nu(sbe, nu, evals(b_em))
                end if
                gam = nu_p * sbe%eph_wrel(ip) * fe * shp * blk
                call amp_damp_channel(nba, rho_ad, ia, b_em, gam, tau)
            end if
            if (b_ab > 0) then
                nu_p = nu
                if (sbe%eph_db_realized) then
                    call eph_thermal_split_de(evals(b_ab) - evals(ia), sbe%eph_kt_au, fe, fa)
                    nu_p = pair_nu(sbe, nu, evals(b_ab))
                end if
                if (fa > 0d0) then
                    shp = gaussian_shape(best_ab, sig)
                    blk = min(max(1d0 - real(rho_ad(b_ab, b_ab)) / f, 0d0), 1d0)
                    gam = nu_p * sbe%eph_wrel(ip) * fa * shp * blk
                    call amp_damp_channel(nba, rho_ad, ia, b_ab, gam, tau)
                end if
            end if
        end do
    end do
end subroutine apply_eph_relaxation


! Pair-symmetric collision prefactor: nu depends on the carrier's distance from the
! nearest band edge, so nu(eps_a) alone gives a pair two different rates for its two
! directions and breaks detailed balance independently of the Bose factors. The
! geometric mean is the symmetric interpolation that leaves nu unchanged when the two
! energies coincide.
pure function pair_nu(sbe, nu_a, eval_b) result(nu_pair)
    use sbe_superres_ssbe, only: nu_saturation
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    real(8), intent(in) :: nu_a, eval_b
    real(8) :: nu_pair, eps_b
    eps_b = max(eval_b - sbe%eph_ecbm_au, sbe%eph_evbm_au - eval_b, 0d0)
    nu_pair = sqrt(max(nu_a, 0d0) * nu_saturation(eps_b, sbe%eph_nusat_au, &
                                                  sbe%eph_eps0_au, sbe%eph_n))
end function pair_nu


!=============================================================================
! B2: UNIFIED nonlocal ring channels -- ONE shared Houston pass + gather for
! inter-k e-ph, nonlocal II + its Auger time-reverse, and the graphene 2D Rana
! Auger/CM. One ZHEEV per k per step (was up to three); every channel computes
! its dpop from the SAME gathered pre-step populations and is applied
! sequentially through the SHARED basis (each apply re-extracts the CURRENT
! rho, so CPTP composition is exact). Also hosts, per the approved wiki/00
! plan: B1 (precomputed vq table, bit-identical), B3 (vq windowing floor),
! A1 (phonon-assisted sidebands), A2 (hole-initiated channel), A3 (screened
! Frohlich 1/q^2 weight for the polar-LO e-ph mode), A5 (Franz-Keldysh-
! softened II threshold from |E(t)|), A8 (carrier-temperature-aware
! lambda^2(n, T_c)), and C1 (the per-channel ledger).
!=============================================================================
subroutine apply_ring_channels(sbe, gs, Ac, efield_au, tau)
    use sbe_superres_ssbe, only: eph_interk_dpop, ii_interk_dpop, auger_interk_dpop, &
                                 rana_auger_dpop, mp_grid_triple, get_material_params, &
                                 s_material_params, build_vq_table, build_acscreen_table, &
                                 t_ring_opts, debye_kappa2, tf_kappa2_degenerate, fit_fermi_dirac, &
                                 colmem_pop_filter, colmem_pop_init, dressed_ref_delta, &
                                 colmem_lines, bose_factor, dirac_mu_2d, dirac_plasmon_2d, &
                                 dirac_fit_te, dirac_n_2d
    use eigen_lapack, only: eigen_zheev
    use communication, only: comm_summation
    use salmon_global, only: num_kgrid, epm_material, sbe_eph_temperature_k
    use math_constants, only: pi
    use phys_constants, only: kB_au
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    real(8),                  intent(in)    :: Ac(3), efield_au, tau

    type(s_material_params) :: mp
    type(t_ring_opts) :: opts
    integer :: nba, nk, ik, i, j, in, im, idir, a, s, m(3), iv, ic, lidx, ntab
    real(8) :: a2half, pcoset, resid, maxresid, sig, rnet
    real(8) :: eps_inf, qtf2, wp2, lambda2, q2reg, n_val, kf, blen2, kt_au
    real(8),    allocatable :: eval_loc(:,:), f_loc(:,:), eval_all(:,:), f_all(:,:)
    real(8),    allocatable :: dpop(:,:), dpop2(:,:), dpop_loc(:,:), gout(:,:), actab(:)
    integer :: ipol_use, ipac_use
    real(8) :: kappa2_c, pnorm
    ! 2D colmem analog (graphene): raw / memory-filtered Rana source densities
    real(8)    :: n_raw, p_raw, n_fil, p_fil, wpl_au, mu_c2, mu_v2
    real(8)    :: e_raw, kt_use, ed_au, mu_h2
    complex(8) :: cl2(2), mul2(2)
    integer    :: nl2
    complex(8), allocatable :: U_loc(:,:,:)
    real(8)    :: eigen_active(sbe%n_active_bands), evals(sbe%n_active_bands)
    complex(8) :: p_active(sbe%n_active_bands, sbe%n_active_bands, 3)
    complex(8) :: HVG(sbe%n_active_bands, sbe%n_active_bands)
    complex(8) :: W(sbe%n_active_bands, sbe%n_active_bands)
    complex(8) :: t1(sbe%n_active_bands, sbe%n_active_bands)
    complex(8) :: rad(sbe%n_active_bands, sbe%n_active_bands)

    nba = sbe%n_active_bands
    if (nba <= 0) return
    nk = sbe%nk
    iv = sbe%nv_act
    ic = sbe%nv_act + 1
    a2half = 0.5d0 * dot_product(Ac, Ac)
    sig    = max(sbe%eph_sigma_au, 2d-3)
    kappa2_c = 0d0

    allocate(eval_loc(nba,nk), f_loc(nba,nk), eval_all(nba,nk), f_all(nba,nk))
    allocate(dpop(nba,nk), U_loc(nba,nba, sbe%ik_min:sbe%ik_max))
    eval_loc = 0d0; f_loc = 0d0

    ! ---- shared pass 1: Houston spectrum + populations, gathered over all k --
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            do j = 1, nba
                im = sbe%active_idx(j)
                do i = 1, nba
                    in = sbe%active_idx(i)
                    p_active(i, j, idir) = gs%p_tm_matrix(in, im, idir, ik)
                    if (sbe%flag_vnl_correction) &
                        p_active(i, j, idir) = p_active(i, j, idir) + gs%rvnl_tm_matrix(in, im, idir, ik)
                end do
            end do
        end do
        if (sbe%flag_coset_proj) then
            do j = 1, nba
                im = sbe%active_idx(j)
                do i = 1, nba
                    if (i == j) cycle
                    in = sbe%active_idx(i)
                    pcoset = 0d0
                    do s = 1, 4
                        pcoset = pcoset + gs%unfold_w(s, in, ik) * gs%unfold_w(s, im, ik)
                    end do
                    p_active(i, j, 1:3) = p_active(i, j, 1:3) * pcoset
                end do
            end do
        end if
        do i = 1, nba
            eigen_active(i) = gs%eigen(sbe%active_idx(i), ik)
        end do
        call build_HVG(nba, eigen_active, p_active, Ac, HVG)
        if (sbe%flag_coulomb) HVG = HVG + sbe%sigma_hf(:, :, ik)
        call eigen_zheev(HVG, evals, W)
        U_loc(:,:,ik) = W
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                rad(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        call ZGEMM('C','N', nba,nba,nba, cmplx(1d0,0d0, 8), W,nba, rad,nba, cmplx(0d0,0d0, 8), t1,nba)
        call ZGEMM('N','N', nba,nba,nba, cmplx(1d0,0d0, 8), t1,nba, W,nba, cmplx(0d0,0d0, 8), rad,nba)
        do a = 1, nba
            eval_loc(a,ik) = evals(a)
            f_loc(a,ik)    = real(rad(a,a))
        end do
        ! ---- Option A: dressed-reference carrier measure (wiki/10 sec. 3A) --
        ! Subtract the field-rotated-GS background delta0 (trace-neutral,
        ! vanishes at A = 0) so the channels see only genuine excess; the
        ! sub-cycle state tracks this rotation to corr 0.99 (wiki/00), so
        ! the dominant fabrication seed is removed at the SOURCE. Composes
        ! with the ring gate and the collisional-memory pop filter (both
        ! act downstream on the gathered f). Clamped to [0, occ] after the
        ! gather (an adiabatic state gives negative upper excess -> 0).
        if (sbe%flag_dressed_ref) then
            block
                real(8) :: dref(nba)
                if (allocated(sbe%f0_act)) then
                    call dressed_ref_delta(nba, sbe%nv_act, sbe%occ_max, U_loc(:,:,ik), dref, sbe%f0_act(:, ik))
                else
                    call dressed_ref_delta(nba, sbe%nv_act, sbe%occ_max, U_loc(:,:,ik), dref)
                end if
                f_loc(:, ik) = f_loc(:, ik) - dref(:)
            end block
        end if
    end do
    call comm_summation(eval_loc, eval_all, nba*nk, sbe%icomm)
    call comm_summation(f_loc,    f_all,    nba*nk, sbe%icomm)
    if (sbe%flag_dressed_ref) f_all = min(max(f_all, 0d0), sbe%occ_max)

    ! ---- virtual-transient gate: population kernels see REAL carriers only --
    ! f_gate is the persistent Houston floor: it drops with f instantly (what
    ! left cannot be scattered) but rises toward f only with the virtuality
    ! time constant 2*pi/Egap -- so sub-lifetime LZ / field-dressing transients
    ! (which the reversible unitary takes back) never feed the population
    ! transfer, the free-carrier screen, or the FD fit. They still LOSE
    ! COHERENCE at the full nu through the gout/T2 exponential (whose rate is
    ! per-carrier, not population-weighted) -- e-ph keeps its T2 role for the
    ! virtual electrons; only the real ones are transported. Deterministic and
    ! rank-identical (f_all is). Not checkpointed: after a restart the floor
    ! re-seeds from the instantaneous f (one-gate-time transient).
    if (sbe%ring_gate_tau_au > 0d0) then
        if (.not. allocated(sbe%f_ring_gate)) then
            allocate(sbe%f_ring_gate(nba, nk))
            sbe%f_ring_gate = f_all
        end if
        block
            real(8) :: gfac
            gfac = 1d0 - exp(-tau / sbe%ring_gate_tau_au)
            sbe%f_ring_gate = min(f_all, sbe%f_ring_gate + &
                                  (f_all - sbe%f_ring_gate) * gfac)
        end block
        f_all = sbe%f_ring_gate
    end if

    ! ---- collisional-memory POPULATION filter (wiki/10 sec. 8.8) -----------
    ! The Stark-dressed diagonal carries a virtual share ~ |eE.d/E_g|^2 ~ I
    ! that the Markovian golden-rule kernels convert at a linear rate (the
    ! surviving ~ I tail after the coherence-sector colmem fix). Remedy: the
    ! kernels read the MEMORY-FILTERED populations -- auxiliary fields z_j
    ! per (state, k) convolve f with the SAME cited phonon lines; a constant
    ! f is a machine-exact fixed point (discrete anchor, calibrated rates
    ! untouched), while the A^2(t) dressing breathing (2*w_laser >> phonon
    ! lines) filters out of the collision SOURCE -- the time-domain ICFE.
    ! Filters the gathered f_all like the ring gate => rank-identical; all
    ! downstream consumers (kernels, free-carrier screen, FD fit) see the
    ! filtered f; application to rho keeps the REAL populations + limiter.
    ! Not checkpointed: after a restart z re-seeds from the instantaneous f
    ! (one-memory-time transient), same policy as the ring gate.
    !
    ! 2D colmem analog (graphene, wiki/10 sec. 8.11): capture the Rana SOURCE
    ! densities BEFORE the phonon-line population filter -- after the dressed
    ! reference / ring gate, so the basis-level subtraction composes with the
    ! time filter -- they get their own Coulomb-sector memory line (channel 4).
    n_raw = 0d0;  p_raw = 0d0;  e_raw = 0d0
    if ((sbe%flag_colmem_2d .or. sbe%flag_rana_te) .and. iv >= 1 .and. ic <= nba) then
        n_raw = sum(f_all(ic:nba, :)) / (dble(nk) * sbe%rana_area_au)
        p_raw = max((dble(iv) * sbe%occ_max - sum(f_all(1:iv, :)) / dble(nk)) &
                    / sbe%rana_area_au, 0d0)
        ! kinetic energy density of electrons + holes from the Dirac point E_D
        ! (midpoint of the instantaneous band edges: exact by e-h symmetry)
        ed_au = 0.5d0 * (minval(eval_all(ic, :)) + maxval(eval_all(iv, :)))
        e_raw = (sum(f_all(ic:nba, :) * (eval_all(ic:nba, :) - ed_au)) &
               + sum((sbe%occ_max - f_all(1:iv, :)) * (ed_au - eval_all(1:iv, :)))) &
                / (dble(nk) * sbe%rana_area_au)
        e_raw = max(e_raw, 0d0)
    end if
    if (sbe%flag_colmem_pop) then
        if (.not. allocated(sbe%zpop)) then
            allocate(sbe%zpop(nba, sbe%colmem_nl, nk))
            do ik = 1, nk
                do a = 1, nba
                    call colmem_pop_init(sbe%colmem_nl, sbe%colmem_mu, tau, &
                                         f_all(a, ik), sbe%zpop(a, :, ik))
                end do
            end do
        end if
        block
            real(8) :: ftil
            do ik = 1, nk
                do a = 1, nba
                    call colmem_pop_filter(sbe%colmem_nl, sbe%colmem_c, sbe%colmem_mu, &
                                           tau, f_all(a, ik), sbe%zpop(a, :, ik), ftil)
                    f_all(a, ik) = min(max(ftil, 0d0), sbe%occ_max)
                end do
            end do
        end block
    end if

    ! ---- II/Auger screening context + the B1 vq table (also reused by A3) ---
    mp = get_material_params(epm_material)
    if ((sbe%flag_impact .or. sbe%flag_auger) .and. .not. sbe%flag_rana2d &
        .or. (sbe%flag_eph .and. mp%found .and. mp%eph_polar) &
        .or. (sbe%flag_eph .and. sbe%eph_ip_ac > 0)) then
        eps_inf = merge(mp%eps_inf, 1d0, mp%found)
        n_val   = dble(gs%ne) / max(gs%volume, 1d-30)
        kf      = (3d0 * pi * pi * n_val) ** (1d0 / 3d0)
        qtf2    = 4d0 * kf / pi
        wp2     = 4d0 * pi * n_val
        lambda2 = 0d0
        q2reg   = huge(1d0)
        do idir = 1, 3
            blen2 = dot_product(gs%b_matrix(idir, 1:3), gs%b_matrix(idir, 1:3))
            q2reg = min(q2reg, blen2 / dble(max(num_kgrid(idir), 1))**2)
        end do
        q2reg = 0.25d0 * q2reg
        ! Shared instantaneous free-carrier screen kappa^2 = min(Debye, TF) on
        ! the gathered CB density with the FD-fit carrier temperature (A8).
        ! Consumers: the GaAs lambda^2(n,T_c) in |V(q)|^2 (registry-gated) and
        ! the MANDATORY A4 acoustic TF factor [q/(q+q_TF)]^2 (all materials).
        if ((mp%found .and. mp%dyn_lambda_ok .or. sbe%eph_ip_ac > 0) .and. ic <= nba) then
            block
                real(8) :: n_exc_au, ntot, etot, beta, muf
                real(8), allocatable :: ecb(:), fcb(:), ftgt(:)
                logical :: okf
                logical, save :: lam_printed = .false.
                integer :: ncb, aa, kk, ii2
                n_exc_au = sum(f_all(ic:nba, :)) / (dble(nk) * max(gs%volume, 1d-30))
                if (n_exc_au > 0d0) then
                    kt_au = kB_au * max(sbe_eph_temperature_k, 1d0)
                    ncb = (nba - ic + 1) * nk
                    allocate(ecb(ncb), fcb(ncb), ftgt(ncb))
                    ii2 = 0
                    ntot = 0d0; etot = 0d0
                    do kk = 1, nk
                        do aa = ic, nba
                            ii2 = ii2 + 1
                            ecb(ii2) = eval_all(aa, kk)
                            fcb(ii2) = min(max(f_all(aa, kk) / sbe%occ_max, 0d0), 1d0)
                            ntot = ntot + fcb(ii2)
                            etot = etot + ecb(ii2) * fcb(ii2)
                        end do
                    end do
                    if (ntot > 1d-9 .and. dble(ncb) - ntot > 1d-9) then
                        call fit_fermi_dirac(ncb, ecb, ntot, etot, beta, muf, ftgt, okf)
                        if (okf .and. beta > 1d-12) kt_au = max(kt_au, 1d0 / beta)
                    end if
                    kappa2_c = min(debye_kappa2(n_exc_au, mp%eps0, kt_au), &
                                   tf_kappa2_degenerate(n_exc_au, mp%eps0))
                    if (mp%dyn_lambda_ok) lambda2 = kappa2_c
                    if (.not. lam_printed .and. kappa2_c > 1d-4 .and. sbe%irank == 0) then
                        write(*, '(a,es12.4,a,es12.4,a,f8.1,a)') &
                            '# free-carrier screen active: kappa^2 = ', &
                            kappa2_c, ' a.u. at n_exc = ', n_exc_au * sbe%au_dens_cm3, &
                            ' cm^-3, T_c = ', kt_au / kB_au, ' K (Debye/TF; -> lambda^2'// &
                            ' and/or the acoustic [q/(q+q_TF)]^2 cut)'
                        lam_printed = .true.
                    end if
                end if
            end block
        end if
        ! B1: the signed-difference vq table (bit-identical to the direct call)
        ntab = (2*num_kgrid(1)-1) * (2*num_kgrid(2)-1) * (2*num_kgrid(3)-1)
        allocate(opts%vq_tab(ntab))
        call build_vq_table(num_kgrid, gs%b_matrix, eps_inf, qtf2, wp2, lambda2, &
                            q2reg, opts%vq_tab)
        opts%use_tab  = .true.
        opts%vq_floor = sbe%ring_vq_floor * maxval(opts%vq_tab)     ! B3 (0 = off)
    end if

    ! ---- channel 1: inter-k e-ph (+ A3 screened Frohlich weight, polar LO) --
    if (sbe%flag_eph .and. sbe%eph_nph > 0) then
        allocate(gout(nba, nk))
        ! q-resolved weights: A3 polar Frohlich (unit-average-normalized -> the
        ! cited nu_sat total rate is preserved, only the q-distribution becomes
        ! small-q-peaked) and the A4 acoustic TF screen [q/(q+q_TF)]^2 from the
        ! instantaneous carrier density (MANDATORY small-q cut -- CdS E1=14.5 eV
        ! [Rode 1970] would otherwise blow up at n >= 1e18 cm^-3). Neutral
        ! placeholders (index 0 / size 1) when a weight is inactive.
        ipol_use = 0;  pnorm = 0d0;  ipac_use = 0
        if ((mp%found .and. mp%eph_polar .and. allocated(opts%vq_tab)) &
            .or. sbe%eph_ip_ac > 0) then
            if (.not. sbe%kmap_built) call build_kmap(sbe, gs)
        end if
        if (sbe%kmap_ok) then
            if (mp%found .and. mp%eph_polar .and. allocated(opts%vq_tab)) then
                ipol_use = 1
                pnorm = sum(opts%vq_tab) / dble(size(opts%vq_tab))
            end if
            if (sbe%eph_ip_ac > 0) ipac_use = sbe%eph_ip_ac
        end if
        if (ipac_use > 0) then
            allocate(actab((2*num_kgrid(1)-1)*(2*num_kgrid(2)-1)*(2*num_kgrid(3)-1)))
            call build_acscreen_table(num_kgrid, gs%b_matrix, &
                                      sqrt(max(kappa2_c, 0d0)), actab)
        else
            allocate(actab(1));  actab = 1d0
        end if
        if (ipol_use > 0 .or. ipac_use > 0) then
            if (ipol_use > 0) then
                call eph_interk_dpop(nk, nba, eval_all, f_all, sbe%occ_max, a2half, &
                         sbe%eph_ecbm_au, sbe%eph_evbm_au, sbe%eph_nph, &
                         sbe%eph_hw(1:sbe%eph_nph), sbe%eph_wrel(1:sbe%eph_nph), &
                         sbe%eph_nb(1:sbe%eph_nph), sbe%eph_kt_au, sbe%eph_db_realized, &
                         sbe%eph_nusat_au, sbe%eph_eps0_au, &
                         sbe%eph_n, sbe%eph_sigma_au, tau, dpop, gout, &
                         kidx=sbe%kmap_idx, kn=sbe%kmap_n, pol_tab=opts%vq_tab, &
                         pol_norm=pnorm, ip_polar=1, ac_tab=actab, ip_ac=ipac_use, &
                         ib_scale=sbe%eph_ib_scale, e_src_lo=sbe%ring_e_lo, e_src_hi=sbe%ring_e_hi)
            else
                call eph_interk_dpop(nk, nba, eval_all, f_all, sbe%occ_max, a2half, &
                         sbe%eph_ecbm_au, sbe%eph_evbm_au, sbe%eph_nph, &
                         sbe%eph_hw(1:sbe%eph_nph), sbe%eph_wrel(1:sbe%eph_nph), &
                         sbe%eph_nb(1:sbe%eph_nph), sbe%eph_kt_au, sbe%eph_db_realized, &
                         sbe%eph_nusat_au, sbe%eph_eps0_au, &
                         sbe%eph_n, sbe%eph_sigma_au, tau, dpop, gout, &
                         kidx=sbe%kmap_idx, kn=sbe%kmap_n, &
                         ac_tab=actab, ip_ac=ipac_use, ib_scale=sbe%eph_ib_scale, &
                         e_src_lo=sbe%ring_e_lo, e_src_hi=sbe%ring_e_hi)
            end if
        else
            call eph_interk_dpop(nk, nba, eval_all, f_all, sbe%occ_max, a2half, &
                     sbe%eph_ecbm_au, sbe%eph_evbm_au, sbe%eph_nph, &
                     sbe%eph_hw(1:sbe%eph_nph), sbe%eph_wrel(1:sbe%eph_nph), &
                     sbe%eph_nb(1:sbe%eph_nph), sbe%eph_kt_au, sbe%eph_db_realized, &
                     sbe%eph_nusat_au, sbe%eph_eps0_au, &
                     sbe%eph_n, sbe%eph_sigma_au, tau, dpop, gout, &
                     ib_scale=sbe%eph_ib_scale, e_src_lo=sbe%ring_e_lo, e_src_hi=sbe%ring_e_hi)
        end if
        deallocate(actab)
        call ring_ledger(sbe, 1, nba, nk, ic, eval_all, dpop)
        call ring_apply_dpop(sbe, gs, U_loc, dpop, tau, gout)
        deallocate(gout)
    end if

    ! ---- channels 2+3: nonlocal II + ring Auger (gap materials) -------------
    if ((sbe%flag_impact .or. sbe%flag_auger) .and. .not. sbe%flag_rana2d) then
        if (.not. sbe%kmap_built) call build_kmap(sbe, gs)
        if (sbe%kmap_ok) then
            ! A5: Franz-Keldysh electro-optic width from the instantaneous field
            opts%fk_theta = 0d0
            if (sbe%flag_ii_fk .and. efield_au > 0d0) &
                opts%fk_theta = (efield_au**2 / (2d0 * sbe%ii_fk_mu)) ** (1d0/3d0)
            ! A1: phonon-assisted sidebands from the cited table
            if (sbe%ii_phassist > 0d0 .and. sbe%eph_nph > 0) then
                opts%phassist = sbe%ii_phassist
                opts%nph  = sbe%eph_nph
                opts%hw   = sbe%eph_hw(1:sbe%eph_nph)
                opts%nbb  = sbe%eph_nb(1:sbe%eph_nph)
                opts%wrel = sbe%eph_wrel(1:sbe%eph_nph)
            end if
            ! A2: hole-initiated channel, cited Cp/Cn scale; Houston VBM
            if (sbe%flag_ii_holes) then
                opts%pref_h = sbe%ii_cpcn * sbe%ii_pref_au
                opts%evbm   = maxval(eval_all(iv, :))
            end if
            allocate(dpop_loc(nba, nk), dpop2(nba, nk))
            dpop_loc = 0d0
            if (sbe%flag_impact) then
                call ii_interk_dpop(nk, nba, eval_all, f_all, sbe%occ_max, a2half, &
                        sbe%ii_ecbm_au, sbe%ii_eth_au, sbe%ii_pref_au, sbe%ii_exponent, &
                        iv, ic, sbe%kmap_idx, sbe%kmap_n, sbe%kmap_lut, &
                        gs%b_matrix, eps_inf, qtf2, wp2, lambda2, q2reg, sig, tau, dpop_loc, &
                        i1_lo=sbe%ik_min, i1_hi=sbe%ik_max, opts=opts, &
                        e_src_lo=sbe%ring_e_lo, e_src_hi=sbe%ring_e_hi)
                call comm_summation(dpop_loc, dpop, nba*nk, sbe%icomm)
                call ring_ledger(sbe, 2, nba, nk, ic, eval_all, dpop)
            else
                dpop = 0d0
            end if
            if (sbe%flag_auger) then
                dpop_loc = 0d0
                call auger_interk_dpop(nk, nba, eval_all, f_all, sbe%occ_max, a2half, &
                        sbe%ii_ecbm_au, sbe%ii_eth_au, sbe%ii_pref_au, sbe%ii_exponent, &
                        iv, ic, sbe%kmap_idx, sbe%kmap_n, sbe%kmap_lut, &
                        gs%b_matrix, eps_inf, qtf2, wp2, lambda2, q2reg, sig, tau, dpop_loc, &
                        i1_lo=sbe%ik_min, i1_hi=sbe%ik_max, opts=opts, &
                        e_src_lo=sbe%ring_e_lo, e_src_hi=sbe%ring_e_hi)
                call comm_summation(dpop_loc, dpop2, nba*nk, sbe%icomm)
                call ring_ledger(sbe, 3, nba, nk, ic, eval_all, dpop2)
                dpop = dpop + dpop2
            end if
            call ring_apply_dpop(sbe, gs, U_loc, dpop, tau)
            deallocate(dpop_loc, dpop2)
        end if
    end if

    ! ---- channel 4: graphene 2D Rana Auger/CM -------------------------------
    if (sbe%flag_rana2d .and. iv >= 1 .and. ic <= nba) then
        kt_use = sbe%rana_kt_au
        if (sbe%flag_rana_te) then
            ! two-temperature model (wiki/12 sec. 3): a COMMON carrier temperature and
            ! the quasi-Fermi levels from the moments (n, p, eps) of the cone
            ! populations; the lattice stays at the phonon-bath T and cools the
            ! carriers through the e-ph channel. T_e is READ from the distribution
            ! each ring step (no separate T_e rate equation).
            ! No carriers -> no temperature: below 1e-3 of the lattice balance density
            ! the three-moment fit is ill-conditioned (a handful of virtual remnants
            ! read as 1e4 K) and the rates do not depend on T_e anyway; hold the bath.
            if (n_raw + p_raw < 1d-3 * dirac_n_2d(0d0, sbe%rana_kt_au, sbe%rana_vf_au)) then
                kt_use = sbe%rana_kt_au;  mu_c2 = 0d0;  mu_h2 = 0d0
            else
                call dirac_fit_te(n_raw, p_raw, e_raw, sbe%rana_vf_au, sbe%rana_kt_au, 2d0, &
                                  kt_use, mu_c2, mu_h2)
            end if
            sbe%te_au = kt_use;  sbe%te_mu_c = mu_c2;  sbe%te_mu_h = mu_h2
            sbe%te_n2d = n_raw;  sbe%te_p2d = p_raw
        end if
        if (sbe%flag_colmem_2d) then
            ! 2D colmem analog: memory-filter the SOURCE densities with the 2D
            ! Dirac-plasmon line of the instantaneous e-h plasma, w_pl(n,p) at
            ! the collision's own screening momentum Q_TF [R07 Eq. 13] -- the
            ! build-up time of screening is the memory time of the Coulomb
            ! collision. Lines +/- w_pl with the (N+1)/N thermal split and the
            ! common width 1/tau_c; DISCRETE anchor -> a constant density is a
            ! machine-exact fixed point (the calibrated R07 rates are untouched
            ! for slow populations) while the 2*w_laser breathing of the
            ! virtual share filters out of the collision source. Rank-identical
            ! (built from the gathered f_all); z re-seeds after a restart.
            mu_c2 =  dirac_mu_2d(max(n_raw, 1d-30), kt_use, sbe%rana_vf_au)
            mu_v2 = -dirac_mu_2d(max(p_raw, 1d-30), kt_use, sbe%rana_vf_au)
            wpl_au = max(dirac_plasmon_2d(mu_c2, mu_v2, kt_use, sbe%rana_vf_au, &
                                          sbe%rana_eps_r), 1d-6)
            call colmem_lines(1, (/ wpl_au /), (/ 1d0 /), &
                              (/ min(bose_factor(wpl_au, kt_use), 1d12) /), &
                              sbe%colmem_tauc_au, nl2, cl2, mul2)
            if (.not. sbe%zpl_init) then
                call colmem_pop_init(nl2, mul2(1:nl2), tau, n_raw, sbe%zpl(1:nl2, 1))
                call colmem_pop_init(nl2, mul2(1:nl2), tau, p_raw, sbe%zpl(1:nl2, 2))
                sbe%zpl_init = .true.
            end if
            call colmem_pop_filter(nl2, cl2(1:nl2), mul2(1:nl2), tau, n_raw, sbe%zpl(1:nl2, 1), n_fil)
            call colmem_pop_filter(nl2, cl2(1:nl2), mul2(1:nl2), tau, p_raw, sbe%zpl(1:nl2, 2), p_fil)
            ! one-time diagnostic once a meaningful density exists (1e9 cm^-2)
            if (sbe%irank == 0 .and. .not. sbe%zpl_printed .and. n_raw > 2.8d-8) then
                write(*, '(a,f8.2,a,es10.3,a,es10.3,a)') '# 2D colmem analog: w_pl = ', &
                    wpl_au * au_ev * 1d3, ' meV at n_2d = ', n_raw / (0.52917721067d-8)**2, &
                    ' cm^-2 (raw) -> ', max(n_fil, 0d0) / (0.52917721067d-8)**2, &
                    ' cm^-2 (memory-filtered source)'
                sbe%zpl_printed = .true.
            end if
            call rana_auger_dpop(nk, nba, eval_all, f_all, sbe%occ_max, iv, ic, &
                                 sbe%rana_area_au, kt_use, sbe%rana_vf_au, &
                                 sbe%rana_eps_r, tau, dpop, rnet, &
                                 e_src_lo=sbe%ring_e_lo, e_src_hi=sbe%ring_e_hi, &
                                 n2d_in=n_fil, p2d_in=p_fil)
        else
            call rana_auger_dpop(nk, nba, eval_all, f_all, sbe%occ_max, iv, ic, &
                                 sbe%rana_area_au, kt_use, sbe%rana_vf_au, &
                                 sbe%rana_eps_r, tau, dpop, rnet, &
                                 e_src_lo=sbe%ring_e_lo, e_src_hi=sbe%ring_e_hi)
        end if
        call ring_ledger(sbe, 4, nba, nk, ic, eval_all, dpop)
        call ring_apply_dpop(sbe, gs, U_loc, dpop, tau)
    end if

    deallocate(eval_loc, f_loc, eval_all, f_all, dpop, U_loc)
end subroutine apply_ring_channels


! C1: accumulate the per-channel ledger -- the conduction-population change
! (pair creation > 0 / recombination < 0) and the population-weighted energy
! change of channel ich in {1 e-ph, 2 II, 3 ring Auger, 4 Rana}. dpop is the
! full-BZ change (identical on every rank), so the ledger needs no reduction.
subroutine ring_ledger(sbe, ich, nba, nk, ic, eval_all, dpop)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    integer, intent(in) :: ich, nba, nk, ic
    real(8), intent(in) :: eval_all(nba, nk), dpop(nba, nk)
    if (ic <= nba) sbe%led_dn(ich) = sbe%led_dn(ich) + sum(dpop(ic:nba, :)) / dble(nk)
    sbe%led_de(ich) = sbe%led_de(ich) + sum(eval_all * dpop) / dble(nk)
end subroutine ring_ledger


! Shared pass 2: apply a Houston-diagonal population change dpop to the LOCAL
! rho(k) through the SHARED basis U_loc, re-extracting the CURRENT rho (so
! sequential channel application composes exactly). Coherence damping:
! sqrt(f_new/f_old) for population-losing levels (the exact amplitude-damping
! Kraus factor) and, when gout is present (e-ph), the out-rate exponential
! exp(-(g_a+g_b) tau / 2) -- the two conventions of the pre-B2 routines.
subroutine ring_apply_dpop(sbe, gs, U_loc, dpop, tau, gout)
    use communication, only: comm_get_min
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    complex(8), intent(in) :: U_loc(sbe%n_active_bands, sbe%n_active_bands, sbe%ik_min:sbe%ik_max)
    real(8),    intent(in) :: dpop(sbe%n_active_bands, sbe%nk), tau
    real(8),    intent(in), optional :: gout(sbe%n_active_bands, sbe%nk)
    integer :: ik, i, j, in, im, a, b, nba, jl
    real(8) :: fold, fnew, scal, grate
    complex(8) :: mem_c, cnew
    real(8)    :: damp(sbe%n_active_bands)
    real(8)    :: dampf(sbe%n_active_bands)
    complex(8) :: W(sbe%n_active_bands, sbe%n_active_bands)
    complex(8) :: t1(sbe%n_active_bands, sbe%n_active_bands)
    complex(8) :: rad(sbe%n_active_bands, sbe%n_active_bands)
    complex(8) :: Dk(sbe%n_active_bands, sbe%n_active_bands)
    complex(8) :: cvec(sbe%n_active_bands), ctmp
    complex(8) :: rho0(sbe%n_active_bands, sbe%n_active_bands)  ! dressed background (Houston basis), Option A coherence sector
    real(8), allocatable :: fold_loc(:,:)
    logical, save :: limiter_printed = .false.
    real(8), parameter :: dtol = 1d-14
    ! Absolute one-step occupancy slack for the limiter ratios: a state whose
    ! REMAINING capacity is pure roundoff (a full valence band at t~0, or a
    ! frozen-window diagonal a hair outside [0,occ]) must not stall the WHOLE
    ! BZ's dissipation through the global min over capacity/flux -- with the
    ! slack the per-state overdraw is bounded by captol (trace stays exact,
    ! positivity is violated by at most captol) while any MATERIAL overdraw is
    ! still scaled away.
    real(8), parameter :: captol = 1d-12

    nba = sbe%n_active_bands

    ! ---- CPTP LIMITER (pass A): the per-state dpop of a nonlocal channel is a
    ! SUM over many independent quadruples computed against the step-START
    ! populations; within one (large) step several sources can pile onto the
    ! same sink and overdraw it (f + dpop < 0 or > occ_max). Truncating that
    ! overdraw per-state (the old max(...,0)) silently DESTROYS the trace: the
    ! clipped loss keeps its paired gains => net particle creation (the
    ! runaway "electrons 8 -> 16" seen with all channels + large dt). The CPTP
    ! repair: scale the WHOLE dpop field by the largest s in [0,1] that keeps
    ! every state inside [0, occ_max]. s*dpop still sums to 0 exactly (trace
    ! preserved), all bounds hold, and s -> 1 for a properly resolved dt.
    allocate(fold_loc(nba, sbe%ik_min:sbe%ik_max))
    scal = 1d0
    do ik = sbe%ik_min, sbe%ik_max
        W = U_loc(:,:,ik)
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                rad(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        call ZGEMM('C','N', nba,nba,nba, cmplx(1d0,0d0, 8), W,nba, rad,nba, cmplx(0d0,0d0, 8), t1,nba)
        call ZGEMM('N','N', nba,nba,nba, cmplx(1d0,0d0, 8), t1,nba, W,nba, cmplx(0d0,0d0, 8), rad,nba)
        do a = 1, nba
            fold = min(max(real(rad(a,a)), 0d0), sbe%occ_max)
            fold_loc(a, ik) = real(rad(a,a))
            if (dpop(a,ik) < -dtol) scal = min(scal, (fold + captol) / (-dpop(a,ik)))
            if (dpop(a,ik) >  dtol) scal = min(scal, (sbe%occ_max - fold + captol) / dpop(a,ik))
        end do
    end do
    scal = max(scal, 0d0)
    if (sbe%nproc > 1) call comm_get_min(scal, sbe%icomm)
    if (scal < 0.999d0 .and. .not. limiter_printed .and. sbe%irank == 0) then
        write(*,'(a,es10.3,a)') '# ring CPTP limiter engaged: dpop scaled by s = ', &
            scal, ' (population flux exceeds one-step capacity -- consider a smaller dt)'
        limiter_printed = .true.
    end if

    ! ---- pass B: apply the (scaled) transfer in the same Houston basis ------
    do ik = sbe%ik_min, sbe%ik_max
        W = U_loc(:,:,ik)
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                rad(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        call ZGEMM('C','N', nba,nba,nba, cmplx(1d0,0d0, 8), W,nba, rad,nba, cmplx(0d0,0d0, 8), t1,nba)
        call ZGEMM('N','N', nba,nba,nba, cmplx(1d0,0d0, 8), t1,nba, W,nba, cmplx(0d0,0d0, 8), rad,nba)
        do a = 1, nba
            fold = fold_loc(a, ik)
            ! TRACE-EXACT diagonal update: sum_a (scal*dpop) = scal*0 = 0, so the
            ! written diagonal preserves Tr(rad) EXACTLY. The old max(.,0) clip used
            ! the RAW fold while the limiter `scal` was built from the [0,occ]-clamped
            ! fold; when a FROZEN-window active subblock (a principal submatrix with
            ! the active<->frozen coherences dropped) has a Houston population dip a
            ! hair below 0, that clip silently CREATED population every step (a ~0.5%
            ! electron leak at strong field). In the all-active case the full rho is
            ! PSD so raw fold in [0,occ] and this is a no-op. Positivity is preserved
            ! to numerical tolerance; trace is now preserved exactly.
            fnew = fold + scal * dpop(a,ik)
            if (dpop(a,ik) < 0d0 .and. fold > 1d-12) then
                damp(a) = sqrt(max(fnew, 0d0) / fold)   ! clamp the coherence ratio only (avoid sqrt<0)
            else
                damp(a) = 1d0
            end if
            ! gout is the coherence out-rate paired with the population loss;
            ! scale it consistently with the limited population transfer.
            ! With collisional MEMORY on (wiki/10 sec. 8.6) the within-block
            ! gout damping is applied through the memory convolution below
            ! instead of this instantaneous exponential; the active<->frozen
            ! extension keeps the Markovian factor (dampf) -- those coherences
            ! are far off-shell and get the conservative CP-consistent damping.
            dampf(a) = damp(a)
            if (present(gout)) then
                if (sbe%flag_colmem) then
                    dampf(a) = dampf(a) * exp(-0.5d0 * scal * gout(a, ik) * tau)
                else
                    damp(a) = damp(a) * exp(-0.5d0 * scal * gout(a, ik) * tau)
                    dampf(a) = damp(a)
                end if
            end if
            rad(a,a) = cmplx(fnew, 0d0, 8)
        end do
        ! Option A, COHERENCE sector (снятие одёжки): damp toward the reversible
        ! dressed background rho0(a,b) = occ sum_{v<=nv} conj(W(v,a)) W(v,b) -- the
        ! field-free GS carried into the instantaneous dressed (Houston) basis --
        ! instead of toward 0. The static dressed_ref only subtracted the DIAGONAL
        ! of rho0 (f0_a) from the population measure; the reversible dressing lives
        ! equally in the OFF-DIAGONAL coherences, and damping those realifies it
        ! (the dominant, dt-divergent fabrication -- wiki/06 sec.6-iv). Damping the
        ! EXCESS coherence rad-rho0 preserves the reversible dressing exactly:
        !   rho_ab -> rho0_ab + (rho_ab - rho0_ab) damp_a damp_b.
        ! rho0 is Hermitian PSD with |rho0_ab|^2 <= f0_a f0_b, A->0 => W->1 =>
        ! rho0 -> diag(occ,0) (off-diagonals vanish) => byte-identical to the old
        ! damp-to-0 path. Trace untouched (diagonal not modified here).
        if (sbe%flag_dressed_ref) then
            call ZGEMM('C','N', nba,nba, min(sbe%nv_act,nba), &
                       cmplx(sbe%occ_max,0d0,8), W,nba, W,nba, cmplx(0d0,0d0,8), rho0,nba)
            do b = 1, nba
                do a = 1, nba
                    if (a /= b) rad(a,b) = rho0(a,b) + (rad(a,b) - rho0(a,b)) * damp(a) * damp(b)
                end do
            end do
        else
            do b = 1, nba
                do a = 1, nba
                    if (a /= b) rad(a,b) = rad(a,b) * damp(a) * damp(b)
                end do
            end do
        end if
        ! ---- collisional-memory dephasing (wiki/10 sec. 8.6) -----------------
        ! Replaces the instantaneous exp(-(g_a+g_b) tau/2): the Houston-frame
        ! coherence drives auxiliary fields z_j (one per kernel line, decay
        ! e^{-mu_j tau} + source), and is damped by the convolution
        ! rho_ab -= (g_a+g_b)/2 * tau * sum_j c_j z_j. A slow (adiabatic)
        ! coherence is damped at exactly the Markovian rate (anchor R(0)=1);
        ! sub-correlation-time (field-driven) modulation is protected -- the
        ! phonon bath cannot follow it. z is attached to the sorted Houston
        ! branch indices, first-order-consistent like every ring quantity.
        ! Trace untouched (diagonal not modified); Hermiticity by mirroring.
        if (sbe%flag_colmem .and. present(gout)) then
            do b = 2, nba
                do a = 1, b - 1
                    mem_c = (0d0, 0d0)
                    do jl = 1, sbe%colmem_nl
                        mem_c = mem_c + sbe%colmem_c(jl) * sbe%zmem(a, b, jl, ik)
                    end do
                    grate = 0.5d0 * scal * (gout(a, ik) + gout(b, ik))
                    cnew = rad(a, b) - grate * tau * mem_c
                    do jl = 1, sbe%colmem_nl
                        sbe%zmem(a, b, jl, ik) = sbe%zmem(a, b, jl, ik) &
                            * exp(-sbe%colmem_mu(jl) * tau) + rad(a, b) * tau
                    end do
                    rad(a, b) = cnew
                    rad(b, a) = conjg(cnew)
                end do
            end do
        end if
        call ZGEMM('N','N', nba,nba,nba, cmplx(1d0,0d0, 8), W,nba, rad,nba, cmplx(0d0,0d0, 8), t1,nba)
        call ZGEMM('N','C', nba,nba,nba, cmplx(1d0,0d0, 8), t1,nba, W,nba, cmplx(0d0,0d0, 8), rad,nba)
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                sbe%rho(in, im, ik) = rad(i, j)
            end do
        end do
        ! ---- CP extension to the FROZEN sector -------------------------------
        ! The loss-Kraus K0 = W diag(damp) W^dagger that damped the WITHIN-block
        ! coherences must also damp the active<->frozen coherence blocks (the
        ! full-space Kraus is K0 (+) 1_frozen; the gain operators have no frozen
        ! part). Skipping it leaves |rho_af|^2 > f_a f_f after a population
        ! loss, the full rho stops being PSD, and the next step's Houston
        ! diagonal dips below 0 -- the accumulating negative-population
        ! pathology (nelec/nhole < 0) seen on frozen windows at large dt, which
        ! the eeh channel then converted into a trace leak.
        if (nba < sbe%nb .and. any(dampf < 1d0)) then
            do j = 1, nba
                t1(:, j) = W(:, j) * dampf(j)
            end do
            call ZGEMM('N','C', nba,nba,nba, cmplx(1d0,0d0, 8), t1,nba, W,nba, &
                       cmplx(0d0,0d0, 8), Dk,nba)
            do im = 1, sbe%nb
                if (sbe%is_active(im)) cycle
                do i = 1, nba
                    cvec(i) = sbe%rho(sbe%active_idx(i), im, ik)
                end do
                do i = 1, nba
                    ctmp = cmplx(0d0, 0d0, 8)
                    do j = 1, nba
                        ctmp = ctmp + Dk(i, j) * cvec(j)
                    end do
                    sbe%rho(sbe%active_idx(i), im, ik) = ctmp
                    sbe%rho(im, sbe%active_idx(i), ik) = conjg(ctmp)
                end do
            end do
        end if
    end do
    deallocate(fold_loc)
end subroutine ring_apply_dpop


! MP momentum-conservation index map for the ring II/Auger (built once).
subroutine build_kmap(sbe, gs)
    use sbe_superres_ssbe, only: mp_grid_triple
    use salmon_global, only: num_kgrid
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: ik, m(3), lidx, nk
    real(8) :: resid, maxresid
    nk = sbe%nk
    allocate(sbe%kmap_idx(3, nk), sbe%kmap_lut(0:nk-1))
    sbe%kmap_n = num_kgrid
    sbe%kmap_lut = 0
    maxresid = 0d0
    do ik = 1, nk
        call mp_grid_triple(gs%kpoint(:, ik), num_kgrid, m, resid)
        sbe%kmap_idx(:, ik) = m
        maxresid = max(maxresid, resid)
        lidx = m(1) + num_kgrid(1) * (m(2) + num_kgrid(2) * m(3))
        if (lidx >= 0 .and. lidx <= nk - 1) sbe%kmap_lut(lidx) = ik
    end do
    sbe%kmap_ok = (maxresid < 1d-6) .and. &
                  (num_kgrid(1) * num_kgrid(2) * num_kgrid(3) == nk) .and. &
                  (minval(sbe%kmap_lut) >= 1)
    sbe%kmap_built = .true.
    if (sbe%irank == 0 .and. .not. sbe%kmap_ok) &
        write(*,'(a)') '# NOTE: momentum-resolved nonlocal II disabled '// &
                       '(k-grid is not a regular MP mesh)'
end subroutine build_kmap


!=============================================================================
! Part F: carrier-carrier (e-e/e-h) thermalization as a CPTP channel in the
! Houston/adiabatic basis. Carrier-carrier scattering conserves total carrier
! NUMBER and ENERGY (it thermalizes the distribution to a hot Fermi-Dirac and
! produces excitation-induced dephasing, but does NOT relax energy to the
! lattice). We realize this by relaxing the adiabatic populations toward the
! Fermi-Dirac with the SAME number and energy, with coherence damping (EID):
!
!   rho~ -> (1 - alpha) rho~ + alpha * diag(occ * f_FD),  alpha = 1 - exp(-nu_cc tau),
!
! a convex combination of the identity and the constant-state channel
! diag(occ f_FD) -> EXACTLY CPTP. Because f_FD is fitted to the current (N,E),
! both Tr rho~ (number) and sum_a eps_a rho~_aa (energy) are conserved exactly.
! This is the intra-k thermalization model (the standard relaxation-time /
! BGK carrier-carrier closure); the full inter-k momentum-resolved collision
! integral rides the ring (Part D) and is the documented refinement.
! [Taj-Rossi PRA 78, 052113; Rosati et al. PRB 90, 125140; Goodnick-Lugli
!  PRB 37, 2578; conserves N and E -- validation invariants]
!=============================================================================
subroutine apply_carrier_carrier(sbe, nba, rho_ad, evals, tau, alpha_out)
    use sbe_superres_ssbe, only: carrier_carrier_relax
    implicit none
    type(s_sbe_bloch_solver), intent(in)    :: sbe
    integer,                  intent(in)    :: nba
    complex(8),               intent(inout) :: rho_ad(nba, nba)
    real(8),                  intent(in)    :: evals(nba), tau
    real(8),                  intent(out)   :: alpha_out
    call carrier_carrier_relax(nba, rho_ad, evals, sbe%occ_max, sbe%eeh_nu_au, tau, &
                               alpha_out)
end subroutine apply_carrier_carrier


!=============================================================================
! Auger recombination (Sec 13): density-gated, number-conserving CPTP channel.
! Gap-edge mean-field closure: a conduction electron (lowest CB branch ic1)
! recombines with a valence hole (top VB branch iv1) -- destroying an e-h pair
! -- and the released gap energy E_g promotes a SECOND ic1 electron to the
! conduction state ic_hot energy-matched to E(ic1)+E_g. Both transfers run at
! the per-carrier Auger rate gamma = C n^2 (so the recombination rate is C n^3),
! gated to switch on only above n_gate. Realized as two amplitude-damping maps
! (amp_damp_channel), each trace-preserving -> TOTAL carrier number conserved;
! the Pauli factors (CB electron present, VB hole present, hot target empty) are
! normalized by occ_max and clamped to [0,1], so the map is exactly CPTP and
! recombination stops as the holes fill (hv1 -> 0). Energy is conserved to the
! mean-field (HF-factorization) order, like the impact-ionization channel.
!
! NOTE: this acts on the HOUSTON/adiabatic (real-carrier) populations, not the
! virtual driving polarization, and the rate is C n^2 with a typically tiny C
! -- so Auger is a RARE event that only becomes visible at very high real
! carrier density / strong fields. Its job here is to be present and exactly
! CPTP, not to dominate the dynamics. The coefficient C must come from a cited
! per-material value or an explicit sbe_auger_c_cm6s -- no material ships a
! verified default today (the former CdS "Haury 1998" C was fabricated/removed).
! [GKLS: Taj-Rossi PRA 78, 052113 (2008); cited C: see wiki/07_nonlocal_auger.md]
!=============================================================================
subroutine apply_auger_recombination(sbe, nba, rho_ad, evals, tau)
    use sbe_superres_ssbe, only: amp_damp_channel
    implicit none
    type(s_sbe_bloch_solver), intent(in)    :: sbe
    integer,                  intent(in)    :: nba
    complex(8),               intent(inout) :: rho_ad(nba, nba)
    real(8),                  intent(in)    :: evals(nba)
    real(8),                  intent(in)    :: tau

    integer :: iv1, ic1, ic_hot, a
    real(8) :: gamma0, f, etgt, fc1, hv1, bhot, g_rec, g_prom

    if (sbe%nv_act < 1 .or. sbe%nv_act >= nba) return
    if (sbe%n_exc_cm3 < sbe%auger_n_gate_cm3) return        ! density gate
    iv1 = sbe%nv_act          ! top valence branch (energy-ordered)
    ic1 = sbe%nv_act + 1      ! lowest conduction branch
    f = sbe%occ_max

    ! per-carrier Auger rate gamma = C n^2 [1/a.u.t]
    gamma0 = sbe%auger_c_au * sbe%n_exc_cm3**2
    if (gamma0 * tau < 1d-14) return

    ! hot-carrier target: conduction state energy-matched to E(ic1) + E_g
    etgt = evals(ic1) + sbe%auger_eg_au
    ic_hot = ic1
    do a = ic1, nba
        if (abs(evals(a) - etgt) < abs(evals(ic_hot) - etgt)) ic_hot = a
    end do

    ! Pauli factors (normalized by occ_max, clamped to [0,1])
    fc1  = min(max(real(rho_ad(ic1,   ic1  )) / f, 0d0), 1d0)   ! CB electron present
    hv1  = min(max(1d0 - real(rho_ad(iv1, iv1)) / f, 0d0), 1d0) ! VB hole present
    bhot = min(max(1d0 - real(rho_ad(ic_hot, ic_hot)) / f, 0d0), 1d0) ! hot target empty

    g_rec  = gamma0 * fc1 * hv1            ! recombination: needs CB e- + VB hole
    g_prom = gamma0 * fc1 * hv1 * bhot     ! promotion: + hot target empty (Pauli)
    if (ic_hot /= ic1) &
        call amp_damp_channel(nba, rho_ad, ic1, ic_hot, g_prom, tau)
    call amp_damp_channel(nba, rho_ad, ic1, iv1, g_rec, tau)
end subroutine apply_auger_recombination


! Exact finite-time map of a single Lindblad jump channel L = sqrt(Gamma)|f><i|
! (amplitude damping i -> f) in the basis where i, f are basis states; O(N):
!   rho_ii -> e^{-G tau} rho_ii
!   rho_ff -> rho_ff + (1 - e^{-G tau}) rho_ii
!   rho_ib -> e^{-G tau/2} rho_ib  (and rho_bi), for all b /= i
! Trace-preserving and completely positive for any Gamma, tau >= 0. The jump
! L rho L^dagger = Gamma rho_ii |f><f| feeds only the f-diagonal: no new
! coherences are created, while the populations and coherences of the source
! branch are damped -- ionization is itself a decoherence channel.
subroutine apply_damping_channel(nba, rho, i_src, i_dst, gamma_ch, tau)
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    integer,    intent(in)    :: i_src, i_dst
    real(8),    intent(in)    :: gamma_ch, tau

    real(8) :: g, gh, transfer
    integer :: b

    if (gamma_ch * tau < 1d-14) return
    g  = exp(-gamma_ch * tau)
    gh = sqrt(g)
    transfer = (1d0 - g) * real(rho(i_src, i_src))

    do b = 1, nba
        if (b == i_src) cycle
        rho(i_src, b) = gh * rho(i_src, b)
        rho(b, i_src) = gh * rho(b, i_src)
    end do
    rho(i_src, i_src) = g * rho(i_src, i_src)
    rho(i_dst, i_dst) = rho(i_dst, i_dst) + transfer
end subroutine apply_damping_channel


! Population of band `ib_target` resolved per k-point, in the instantaneous
! Houston (adiabatic) eigenbasis -- the SAME basis the propagator and the
! dissipation half-step use, so an unexcited adiabatically-following state
! reports ZERO conduction population (no spurious gauge offset).
!
! In the Velocity Gauge (VG) the SBE propagates rho(k,t) at the fixed grid
! canonical crystal momentum k; the physical (kinetic) momentum of the carrier
! is k + A(t) (A in a.u., e/hbar = 1). The state it occupies is therefore the
! eigenstate of the instantaneous VG Hamiltonian
!   H_VG(k,t) = H_0(k) + A(t)·p(k) ≈ H_0(k + A(t))            [build_HVG]
! -- the field-free Hamiltonian at the SHIFTED kinetic momentum k + A. The
! population at k must project onto THIS basis (the +A·p Houston basis), not
! the opposite-shifted H_0(k - A) = H_0 - A·p: projecting onto the wrong-sign
! basis leaves an unexcited valence state with a spurious, reversible CB weight
! ~ (2 A·p / E_g)^2 (an O(A^2) offset that grows with the field envelope and is
! largest where the interband coupling/gap is large -- e.g. the folded L-valley),
! masking the genuine non-adiabatic (Zener/multiphoton) excitation.
! We diagonalise H_VG to get U (the Houston rotation), apply a greedy bipartite
! match on |U_ij| to correct the energy-sort ambiguity of ZHEEV at near-
! degeneracies, then form rho_houston = U_sorted^dagger rho_VG U_sorted and
! return Re(rho_houston[ia_target, ia_target]). The Coulomb HF self-energy is
! added when active so the basis matches houston_dissipate exactly.
subroutine calc_bloch_population_k(sbe, gs, Ac, ib_target, pop_k, icomm)
    use eigen_lapack, only: eigen_zheev
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(in)  :: Ac(3)
    integer,                  intent(in)  :: ib_target
    real(8),                  intent(out) :: pop_k(1:sbe%nk)
    integer,                  intent(in)  :: icomm

    integer :: nba, ia_target, ik, i, j, in, im, n_done, best_i, best_j
    real(8)  :: curr_max
    integer,  allocatable :: zone_map(:)
    logical,  allocatable :: row_used(:), col_used(:)
    real(8),    allocatable :: pop_local(:), evals(:), p_k_full(:,:,:), eigen_a(:)
    complex(8), allocatable :: H(:,:), W(:,:), W_sorted(:,:), t1(:,:), t2(:,:), rho_a(:,:)

    nba = sbe%n_active_bands

    ! Locate ib_target in the active subspace; return all-zero if it is frozen.
    ia_target = 0
    do i = 1, nba
        if (sbe%active_idx(i) == ib_target) then
            ia_target = i
            exit
        end if
    end do
    if (ia_target == 0) then
        pop_k = 0d0
        return
    end if

    allocate(pop_local(1:sbe%nk))
    allocate(evals(nba), H(nba,nba), W(nba,nba), W_sorted(nba,nba), t1(nba,nba), t2(nba,nba))
    allocate(rho_a(nba,nba), p_k_full(sbe%nb, sbe%nb, 3), eigen_a(nba))
    allocate(zone_map(nba), row_used(nba), col_used(nba))
    pop_local = 0d0

    do ik = sbe%ik_min, sbe%ik_max
        ! Build the full p matrix exactly as the propagator does.
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) &
            p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)

        ! Restrict H_VG and rho to the active subspace (mirrors dt_evolve_bloch_cf4).
        do i = 1, nba
            eigen_a(i) = gs%eigen(sbe%active_idx(i), ik)
        end do
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                ! H_VG(k,t) = H_0(k) + A·p  (instantaneous Houston basis, k' = k + A)
                H(i, j) = Ac(1)*p_k_full(in,im,1) &
                          + Ac(2)*p_k_full(in,im,2) &
                          + Ac(3)*p_k_full(in,im,3)
                rho_a(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        do i = 1, nba
            H(i, i) = H(i, i) + eigen_a(i)
        end do
        ! Match houston_dissipate's basis exactly when the HF mean field is on.
        if (sbe%flag_coulomb) H(:, :) = H(:, :) + sbe%sigma_hf(:, :, ik)

        ! Diagonalize H_VG = H_0(k+A): H = W Lambda W^dagger  (LAPACK ZHEEV)
        call eigen_zheev(H, evals, W)

        ! Overlap-tracking permutation (greedy bipartite match on |W_ij|).
        ! H_0(k) = diag(eigen) => reference U_0 = I => overlap S = W directly.
        row_used = .false.
        col_used = .false.
        zone_map  = 0
        do n_done = 1, nba
            curr_max = -1d0
            best_i = 1; best_j = 1
            do j = 1, nba
                if (col_used(j)) cycle
                do i = 1, nba
                    if (row_used(i)) cycle
                    if (abs(W(i, j)) > curr_max) then
                        curr_max = abs(W(i, j))
                        best_i   = i
                        best_j   = j
                    end if
                end do
            end do
            zone_map(best_j) = best_i
            row_used(best_i) = .true.
            col_used(best_j) = .true.
        end do

        ! Reorder columns: W_sorted(:, zone_map(j)) = W(:, j)
        do j = 1, nba
            W_sorted(:, zone_map(j)) = W(:, j)
        end do

        ! rho_crystal = U_sorted^dagger rho_VG U_sorted
        call ZGEMM('C', 'N', nba, nba, nba, cmplx(1d0,0d0, 8), W_sorted, nba, rho_a,    nba, cmplx(0d0,0d0, 8), t1, nba)
        call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0,0d0, 8), t1,       nba, W_sorted, nba, cmplx(0d0,0d0, 8), t2, nba)
        pop_local(ik) = real(t2(ia_target, ia_target))
    end do

    call comm_summation(pop_local, pop_k, sbe%nk, icomm)
    deallocate(pop_local, evals, H, W, W_sorted, t1, t2, rho_a, p_k_full, eigen_a)
    deallocate(zone_map, row_used, col_used)
end subroutine calc_bloch_population_k


! Real-carrier (diabatic) population of band ib_target, per k-point, in the
! FIXED field-free Bloch basis -- the k-resolved analogue of the standard
! excited-electron count n_ex (_sbe_nex.data, calc_trace over the conduction
! bands). rho is stored in this basis, so NO field-dependent projection is
! applied: pop = Re(rho_{ib,ib}(k)). This drops the reversible virtual
! polarization (the A^2(t) "breathing" that the instantaneous-Houston
! projection carries) and reports only the real promoted carriers. The
! velocity-gauge diagonal population is the diabatic transition probability:
! it accumulates monotonically and freezes when the field passes, matching
! n_ex exactly when summed over the conduction bands.
subroutine calc_diabatic_population_k(sbe, ib_target, pop_k, icomm)
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    integer,                  intent(in)  :: ib_target
    real(8),                  intent(out) :: pop_k(1:sbe%nk)
    integer,                  intent(in)  :: icomm

    integer :: ik
    real(8), allocatable :: pop_local(:)

    allocate(pop_local(1:sbe%nk))
    pop_local = 0d0
    if (ib_target >= 1 .and. ib_target <= sbe%nb) then
        do ik = sbe%ik_min, sbe%ik_max
            pop_local(ik) = real(sbe%rho(ib_target, ib_target, ik))
        end do
    end if
    call comm_summation(pop_local, pop_k, sbe%nk, icomm)
    deallocate(pop_local)
end subroutine calc_diabatic_population_k


! NON-ADIABATIC (real) excited density: the excitation measured against the
! INSTANTANEOUS DRESSED ground state, i.e. the population that leaks into the
! dressed CONDUCTION states. Per k we diagonalise the active-window velocity-
! gauge Hamiltonian H_VG(k, A(t)) = H_0 + A.p (+ Sigma^HF) -> dressed eigenstates
! (ZHEEV, ascending energy: the lowest nv_act are the dressed valence), transform
! rho into that basis, and sum the diagonal over the dressed conduction states.
!
! This is exactly zero when the system follows the field adiabatically (rho = the
! dressed ground state), so it DROPS the reversible virtual "dressing" (the A^2(t)
! breathing that the fixed-basis conduction sum _sbe_nex.data carries, ~5x at the
! pulse peak) and reports only the genuinely non-adiabatic promoted carriers.
!
! Returns TWO BZ-averaged counts per cell (the caller divides by the cell volume
! to get a density, matching _sbe_nex.data):
!   nex_proj = sum over the dressed CONDUCTION diagonal of W^dag rho W. This still
!              carries the field-rotated-GS leakage of the field-free ground state
!              into the dressed conduction manifold.
!   nex_dref = the SAME dressed-conduction sum with that leakage removed via the
!              Option-A dressed-reference delta0 (dressed_ref_delta), clamped to
!              [0, occ] per (band,k) -- EXACTLY the carrier measure the ring
!              dissipators see (wiki/10 sec.3A). This is the density that drives
!              the density-dependent rates (e-ph nu, screening, FD fit), so it is
!              the physically relevant "real carrier" tracer.
subroutine calc_nex_nonad(sbe, gs, Ac, icomm, nex_proj, nex_dref)
    use eigen_lapack, only: eigen_zheev
    use communication, only: comm_summation
    use sbe_superres_ssbe, only: dressed_ref_delta
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(in)  :: Ac(3)
    integer,                  intent(in)  :: icomm
    real(8),                  intent(out) :: nex_proj, nex_dref

    integer :: nba, nv, ik, i, j, in, im, c
    real(8) :: acc_p, acc_d, fd, dsum(2), gsum(2)
    real(8),    allocatable :: evals(:), p_k_full(:, :, :), eigen_a(:), dref(:)
    complex(8), allocatable :: H(:, :), W(:, :), t1(:, :), t2(:, :), rho_a(:, :)

    ! FULL-basis dressed projection (wiki/06 sec.6): diagonalise H_VG on ALL nb
    ! bands, NOT the truncated frozen active window. Diagonalising on the active
    ! subspace drops the A.p coupling to the frozen bands, so the field-dressing
    ! that belongs in the frozen bands piles into the active conduction states and
    ! inflates nex by ~x300 on a frozen window (measured: Si 4^3, frozen 8-band vs
    ! all-active 16-band, same coherent rho -> 9.7e19 vs 3.3e17). The measure must
    ! be truncation-free even when the DISSIPATORS act only on the active window.
    ! (For an all-active run active_idx(i)=i and nv_act=homo_idx, so this is
    !  byte-identical to the old active-window form.)
    nba = sbe%nb
    nv  = sbe%homo_idx
    nex_proj = 0d0
    nex_dref = 0d0
    ! Need at least one dressed valence and one dressed conduction level.
    if (nba < 1 .or. nv < 1 .or. nv >= nba) return

    dsum(:) = 0d0
    !$omp parallel default(shared) &
    !$omp    private(ik, i, j, in, im, c, acc_p, acc_d, fd, evals, H, W, t1, t2, rho_a, p_k_full, eigen_a, dref)
    allocate(evals(nba), H(nba, nba), W(nba, nba), t1(nba, nba), t2(nba, nba), rho_a(nba, nba))
    allocate(p_k_full(sbe%nb, sbe%nb, 3), eigen_a(nba), dref(nba))
    !$omp do reduction(+: dsum)
    do ik = sbe%ik_min, sbe%ik_max
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) &
            p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)
        do i = 1, nba
            eigen_a(i) = gs%eigen(i, ik)
        end do
        do j = 1, nba
            im = j
            do i = 1, nba
                in = i
                H(i, j) = Ac(1) * p_k_full(in, im, 1) &
                        + Ac(2) * p_k_full(in, im, 2) &
                        + Ac(3) * p_k_full(in, im, 3)
                rho_a(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        do i = 1, nba
            H(i, i) = H(i, i) + eigen_a(i)
        end do
        if (sbe%flag_coulomb) H(:, :) = H(:, :) + sbe%sigma_hf(:, :, ik)
        ! H = W diag(evals) W^dagger, evals ascending -> dressed valence = 1..nv
        call eigen_zheev(H, evals, W)
        ! rho_dressed = W^dagger rho_a W ; its diagonal is the dressed occupation
        call ZGEMM('C', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), W,  nba, rho_a, nba, cmplx(0d0, 0d0, 8), t1, nba)
        call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0, 0d0, 8), t1, nba, W,     nba, cmplx(0d0, 0d0, 8), t2, nba)
        ! Option-A dressed-reference delta0 (same call the ring uses).
        if (allocated(sbe%f0_act)) then
            call dressed_ref_delta(nba, nv, sbe%occ_max, W, dref, sbe%f0_act(:, ik))
        else
            call dressed_ref_delta(nba, nv, sbe%occ_max, W, dref)
        end if
        acc_p = 0d0
        acc_d = 0d0
        do c = nv + 1, nba
            fd = real(t2(c, c))
            acc_p = acc_p + fd
            acc_d = acc_d + min(max(fd - dref(c), 0d0), sbe%occ_max)
        end do
        dsum(1) = dsum(1) + acc_p * gs%kweight(ik)
        dsum(2) = dsum(2) + acc_d * gs%kweight(ik)
    end do
    !$omp end do
    deallocate(evals, H, W, t1, t2, rho_a, p_k_full, eigen_a, dref)
    !$omp end parallel
    call comm_summation(dsum, gsum, 2, icomm)
    nex_proj = gsum(1) / sum(gs%kweight)
    nex_dref = gsum(2) / sum(gs%kweight)
end subroutine calc_nex_nonad


! Population of the PHYSICAL lowest conduction band (CB1) of each folded
! primitive BZ point, resolved per cubic k-point.
!
! In the cubic-supercell EPM dataset every band at k_sc belongs to exactly
! one FCC sublattice s (= primitive point k_prim = k_sc + G0(s); the folding
! is exact, see the unfold map written by `epm_gaas_reference.py unfoldmap`).
! Energy-ordered supercell branch indices therefore mix DIFFERENT physical
! primitive bands from k to k. This routine projects rho into the crystal
! gauge exactly like calc_bloch_population_k and then accumulates, for each
! sublattice s, the population of its physical CB1 level:
!   spinor input: primitive ranks nv_prim+1 and nv_prim+2 (Kramers pair,
!                 spins summed); scalar input: rank nv_prim+1.
! pop_lev(L, s, ik) is the crystal-gauge population of the L-th physical primitive
! level (spins summed) at the folded primitive point k_prim = k_sc(ik) + G0(s):
!   L = 1 -> VB-1 (second valence from top), 2 -> VB (top valence),
!   L = 3 -> CB1 (lowest conduction),         4 -> CB2 (second conduction).
! Band assignment comes entirely from the EPM `unfoldmap` (gs%unfold_sub/prim).
! Energies for the spectral plot are taken from *_bandpath.data, so none are
! emitted here.
subroutine calc_unfolded_population_k(sbe, gs, Ac, pop_lev, icomm)
    use eigen_lapack, only: eigen_zheev
    use salmon_global, only: yn_sbe_spinor
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(in)  :: Ac(3)
    real(8),                  intent(out) :: pop_lev(1:4, 1:4, 1:sbe%nk)
    integer,                  intent(in)  :: icomm

    integer :: nba, ik, i, j, in, im, n_done, best_i, best_j
    integer :: isub, irank_prim, n_spin, nv_phys, pphys, off, islot, s
    real(8)  :: curr_max, popi
    integer,  allocatable :: zone_map(:)
    logical,  allocatable :: row_used(:), col_used(:)
    real(8),    allocatable :: pop_local(:, :, :), evals(:), p_k_full(:,:,:), eigen_a(:)
    complex(8), allocatable :: H(:,:), W(:,:), W_sorted(:,:), t1(:,:), t2(:,:), rho_a(:,:)

    nba = sbe%n_active_bands
    pop_lev = 0d0
    if (.not. gs%have_unfold .or. nba == 0) return

    ! Spin states per physical level in the primitive rank ordering, and the
    ! number of physical valence bands (Kramers doublets for spinor input).
    n_spin  = merge(2, 1, yn_sbe_spinor == 'y')
    nv_phys = gs%nv_prim / n_spin

    allocate(pop_local(1:4, 1:4, 1:sbe%nk))
    allocate(evals(nba), H(nba,nba), W(nba,nba), W_sorted(nba,nba), t1(nba,nba), t2(nba,nba))
    allocate(rho_a(nba,nba), p_k_full(sbe%nb, sbe%nb, 3), eigen_a(nba))
    allocate(zone_map(nba), row_used(nba), col_used(nba))
    pop_local = 0d0

    do ik = sbe%ik_min, sbe%ik_max
        ! Instantaneous Houston projection: identical to calc_bloch_population_k
        ! (H_VG = H_0(k) + A·p, the +A·p basis the propagator populates).
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) &
            p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)

        do i = 1, nba
            eigen_a(i) = gs%eigen(sbe%active_idx(i), ik)
        end do
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                H(i, j) = Ac(1)*p_k_full(in,im,1) &
                          + Ac(2)*p_k_full(in,im,2) &
                          + Ac(3)*p_k_full(in,im,3)
                rho_a(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        do i = 1, nba
            H(i, i) = H(i, i) + eigen_a(i)
        end do
        if (sbe%flag_coulomb) H(:, :) = H(:, :) + sbe%sigma_hf(:, :, ik)

        call eigen_zheev(H, evals, W)

        row_used = .false.
        col_used = .false.
        zone_map  = 0
        do n_done = 1, nba
            curr_max = -1d0
            best_i = 1; best_j = 1
            do j = 1, nba
                if (col_used(j)) cycle
                do i = 1, nba
                    if (row_used(i)) cycle
                    if (abs(W(i, j)) > curr_max) then
                        curr_max = abs(W(i, j))
                        best_i   = i
                        best_j   = j
                    end if
                end do
            end do
            zone_map(best_j) = best_i
            row_used(best_i) = .true.
            col_used(best_j) = .true.
        end do
        do j = 1, nba
            W_sorted(:, zone_map(j)) = W(:, j)
        end do

        call ZGEMM('C', 'N', nba, nba, nba, cmplx(1d0,0d0, 8), W_sorted, nba, rho_a,    nba, cmplx(0d0,0d0, 8), t1, nba)
        call ZGEMM('N', 'N', nba, nba, nba, cmplx(1d0,0d0, 8), t1,       nba, W_sorted, nba, cmplx(0d0,0d0, 8), t2, nba)

        ! Accumulate the spin-summed population of the top two valence and the
        ! bottom two conduction physical primitive bands of every sublattice.
        ! The cubic band's primitive rank (from the unfold map) fixes the physical
        ! slot {VB-1, VB, CB1, CB2}; its population is then DISTRIBUTED over the 4
        ! sublattices by the spectral weights gs%unfold_w (sum_s = 1). At a
        ! symmetry degeneracy this splits the population equally among the
        ! equivalent primitive points, instead of an argmax dumping it on one.
        do i = 1, nba
            in = sbe%active_idx(i)
            isub = gs%unfold_sub(in, ik)
            if (isub < 1 .or. isub > gs%n_coset) cycle
            irank_prim = gs%unfold_prim(in, ik)
            if (irank_prim < 1) cycle
            if (irank_prim <= gs%nv_prim) then
                ! Valence: physical index counted from the band bottom; offset
                ! from the top valence band (off = 0 top, -1 next).
                pphys = (irank_prim + n_spin - 1) / n_spin
                off   = pphys - nv_phys
                if (off == 0) then
                    islot = 2          ! VB (top valence)
                else if (off == -1) then
                    islot = 1          ! VB-1
                else
                    cycle
                end if
            else
                ! Conduction: physical index counted from the band bottom.
                pphys = (irank_prim - gs%nv_prim + n_spin - 1) / n_spin
                if (pphys == 1) then
                    islot = 3          ! CB1
                else if (pphys == 2) then
                    islot = 4          ! CB2
                else
                    cycle
                end if
            end if
            popi = real(t2(i, i))
            do s = 1, gs%n_coset
                pop_local(islot, s, ik) = pop_local(islot, s, ik) &
                    & + gs%unfold_w(s, in, ik) * popi
            end do
        end do
    end do

    call comm_summation(pop_local, pop_lev, 16 * sbe%nk, icomm)
    deallocate(pop_local, evals, H, W, W_sorted, t1, t2, rho_a, p_k_full, eigen_a)
    deallocate(zone_map, row_used, col_used)
end subroutine calc_unfolded_population_k


! Real-carrier (diabatic) twin of calc_unfolded_population_k: same band ->
! {VB-1,VB,CB1,CB2} x coset distribution from the unfold map, but the per-band
! weight is the FIXED-basis occupation Re(rho_{in,in}(k)) instead of the
! instantaneous-Houston projection -- no field-dependent diagonalization, no
! A^2(t) breathing. Reports the real promoted carriers per primitive point.
subroutine calc_diabatic_unfolded_population_k(sbe, gs, pop_lev, icomm)
    use salmon_global, only: yn_sbe_spinor
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(out) :: pop_lev(1:4, 1:4, 1:sbe%nk)
    integer,                  intent(in)  :: icomm

    integer :: nba, ik, i, in, isub, irank_prim, n_spin, nv_phys, pphys, off, islot, s
    real(8) :: popi
    real(8), allocatable :: pop_local(:, :, :)

    nba = sbe%n_active_bands
    pop_lev = 0d0
    if (.not. gs%have_unfold .or. nba == 0) return

    n_spin  = merge(2, 1, yn_sbe_spinor == 'y')
    nv_phys = gs%nv_prim / n_spin

    allocate(pop_local(1:4, 1:4, 1:sbe%nk))
    pop_local = 0d0

    do ik = sbe%ik_min, sbe%ik_max
        do i = 1, nba
            in = sbe%active_idx(i)
            isub = gs%unfold_sub(in, ik)
            if (isub < 1 .or. isub > gs%n_coset) cycle
            irank_prim = gs%unfold_prim(in, ik)
            if (irank_prim < 1) cycle
            if (irank_prim <= gs%nv_prim) then
                pphys = (irank_prim + n_spin - 1) / n_spin
                off   = pphys - nv_phys
                if (off == 0) then
                    islot = 2
                else if (off == -1) then
                    islot = 1
                else
                    cycle
                end if
            else
                pphys = (irank_prim - gs%nv_prim + n_spin - 1) / n_spin
                if (pphys == 1) then
                    islot = 3
                else if (pphys == 2) then
                    islot = 4
                else
                    cycle
                end if
            end if
            popi = real(sbe%rho(in, in, ik))      ! diabatic (fixed-basis) occupation
            do s = 1, gs%n_coset
                pop_local(islot, s, ik) = pop_local(islot, s, ik) &
                    & + gs%unfold_w(s, in, ik) * popi
            end do
        end do
    end do

    call comm_summation(pop_local, pop_lev, 16 * sbe%nk, icomm)
    deallocate(pop_local)
end subroutine calc_diabatic_unfolded_population_k


! Intra-band (group-velocity) current in the instantaneous Houston basis.
! In the velocity gauge only the TOTAL current J = Tr[(p + A + v_nl) rho] is
! gauge invariant; its intra/inter split is basis dependent and is physical in
! the Houston (adiabatic) basis, where the diagonal carries the Boltzmann drift
! of each populated band and the off-diagonal the interband polarization:
!   J_intra = sum_k w_k sum_a f^H_a v^H_aa,   f^H_a = (U^dagger rho U)_aa,
!             v^H_aa = (U^dagger (p + A + v_nl) U)_aa,
! with U diagonalizing H_VG = H_0(k) + A.p (+ Sigma_HF) -- the same Houston
! basis the propagator and dissipation use. Summing J_intra + J_inter recovers
! the gauge-invariant total. [intra/inter decomposition: T. Otobe, PRB 94,
! 235152 (2016)]
subroutine calc_intraband_current_houston(sbe, gs, Ac, jmat_intra, icomm)
    use eigen_lapack, only: eigen_zheev
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(in)  :: Ac(3)
    real(8),                  intent(out) :: jmat_intra(3)
    integer,                  intent(in)  :: icomm

    integer :: nba, ik, i, j, idir, in, im, a
    real(8) :: tmp1(3), tmp(3)
    real(8),    allocatable :: evals(:), eigen_a(:)
    complex(8), allocatable :: p_k_full(:,:,:), p_active(:,:,:), H(:,:), W(:,:)
    complex(8), allocatable :: rho_a(:,:), t1(:,:), rhoH(:,:), vH(:,:)

    nba = sbe%n_active_bands
    tmp1 = 0d0
    if (nba == 0) then
        call comm_summation(tmp1, tmp, 3, icomm)
        jmat_intra = tmp / (sum(gs%kweight) * gs%volume)
        return
    end if

    allocate(evals(nba), eigen_a(nba), p_k_full(sbe%nb, sbe%nb, 3))
    allocate(p_active(nba, nba, 3), H(nba, nba), W(nba, nba))
    allocate(rho_a(nba, nba), t1(nba, nba), rhoH(nba, nba), vH(nba, nba))

    do ik = sbe%ik_min, sbe%ik_max
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) &
            p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)
        do idir = 1, 3
            do j = 1, nba
                im = sbe%active_idx(j)
                do i = 1, nba
                    in = sbe%active_idx(i)
                    p_active(i, j, idir) = p_k_full(in, im, idir)
                end do
            end do
        end do
        do i = 1, nba
            eigen_a(i) = gs%eigen(sbe%active_idx(i), ik)
        end do
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                rho_a(i, j) = sbe%rho(in, im, ik)
            end do
        end do

        ! Houston basis: diagonalize H_VG = H_0(k) + A.p (+ Sigma_HF)
        H(:, :) = Ac(1)*p_active(:,:,1) + Ac(2)*p_active(:,:,2) + Ac(3)*p_active(:,:,3)
        do i = 1, nba
            H(i, i) = H(i, i) + eigen_a(i)
        end do
        if (sbe%flag_coulomb) H(:, :) = H(:, :) + sbe%sigma_hf(:, :, ik)
        call eigen_zheev(H, evals, W)

        ! f^H = U^dagger rho U  (diagonal = Houston populations)
        call ZGEMM('C','N', nba,nba,nba, cmplx(1d0,0d0, 8), W, nba, rho_a, nba, cmplx(0d0,0d0, 8), t1, nba)
        call ZGEMM('N','N', nba,nba,nba, cmplx(1d0,0d0, 8), t1, nba, W, nba, cmplx(0d0,0d0, 8), rhoH, nba)

        do idir = 1, 3
            ! velocity operator v = p + A + v_nl (A on the diagonal), rotate to Houston
            H(:, :) = p_active(:, :, idir)
            do i = 1, nba
                H(i, i) = H(i, i) + Ac(idir)
            end do
            call ZGEMM('C','N', nba,nba,nba, cmplx(1d0,0d0, 8), W, nba, H, nba, cmplx(0d0,0d0, 8), t1, nba)
            call ZGEMM('N','N', nba,nba,nba, cmplx(1d0,0d0, 8), t1, nba, W, nba, cmplx(0d0,0d0, 8), vH, nba)
            do a = 1, nba
                tmp1(idir) = tmp1(idir) + gs%kweight(ik) * real(rhoH(a, a)) * real(vH(a, a))
            end do
        end do
    end do

    call comm_summation(tmp1, tmp, 3, icomm)
    jmat_intra(:) = tmp(:) / (sum(gs%kweight) * gs%volume)
    deallocate(evals, eigen_a, p_k_full, p_active, H, W, rho_a, t1, rhoH, vH)
end subroutine calc_intraband_current_houston


end module
