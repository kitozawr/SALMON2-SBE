!
!  test_eph_interk_cptp.f90  -  inter-k e-ph net-population map (Part C5/D ring).
!
!  Tests eph_interk_dpop: the inter-k (intervalley) e-ph relaxation that runs over
!  the GATHERED Houston spectrum (all k, all bands), enabled when the super-mode
!  ring is on. The net diagonal change dpop must:
!    1) conserve the total trace EXACTLY (sum(dpop) = 0) -- CPTP,
!    2) keep populations in [0, occ_max] after f += dpop,
!    3) move population from a carrier to the ENERGY-MATCHED partner at a
!       DIFFERENT k (true inter-k transfer), in the downhill (emission) direction,
!    4) be a no-op (dpop = 0) when no phonon energy-matches any pair (gamma = 0).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_eph_interk_cptp
    use sbe_superres_ssbe, only: eph_interk_dpop
    implicit none
    integer, parameter :: nk = 3, nba = 2, nph = 1
    real(8), parameter :: nu_n = 2d0          ! saturation exponent (real)
    integer :: nfail, a, ik
    real(8) :: eval(nba, nk), f(nba, nk), dpop(nba, nk)
    real(8) :: occ_max, a2half, ecbm, evbm, sigma, tau, total
    real(8) :: hw(nph), wrel(nph), nb_bose(nph), nu_sat, nu_eps0, kt_bath
    real(8), parameter :: TOL = 1d-12

    nfail = 0
    occ_max = 2d0; a2half = 0d0; sigma = 0.01d0; tau = 0.1d0
    hw(1) = 0.05d0; wrel(1) = 1d0; nb_bose(1) = 0.3d0
    ! kt_bath consistent with nb_bose = 0.3 at hw = 0.05 Ha; the checks below run the
    ! HISTORICAL split (db_realized = .false.), which ignores it -- the realized-transfer
    ! branch has its own test, test_eph_detailed_balance.f90.
    kt_bath = hw(1) / log(1d0 + 1d0 / nb_bose(1))
    nu_sat = 1d0; nu_eps0 = 0.1d0

    ! Band layout: a=1 valence (full), a=2 conduction. Conduction at k1 sits one
    ! phonon (0.05) ABOVE k2 -> emission k1->k2 is energy-matched (inter-k).
    eval(1, :) = (/ -0.30d0, -0.30d0, -0.30d0 /)     ! valence (full, Pauli-blocks down)
    eval(2, 1) = 0.30d0                              ! conduction @ k1 (carrier here)
    eval(2, 2) = 0.25d0                              ! conduction @ k2 = 0.30 - hw
    eval(2, 3) = 0.50d0                              ! conduction @ k3 (far, off-resonant)
    ecbm = 0.25d0; evbm = -0.30d0
    f = 0d0
    f(1, :) = 2d0                                    ! valence full
    f(2, 1) = 0.5d0                                  ! one excited carrier at (cond, k1)

    call eph_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, evbm, &
                         nph, hw, wrel, nb_bose, kt_bath, .false., nu_sat, nu_eps0, nu_n, &
                         sigma, tau, dpop)

    ! (1) trace conserved
    total = sum(dpop)
    call chk("trace conserved (sum dpop = 0)", total, 0d0, TOL)

    ! (2) populations stay in [0, occ_max]
    do ik = 1, nk
        do a = 1, nba
            if (f(a,ik)+dpop(a,ik) < -TOL .or. f(a,ik)+dpop(a,ik) > occ_max+1d-9) &
                call bad("population left [0, occ_max]")
        end do
    end do

    ! (3) inter-k transfer: source (cond,k1) loses, energy-matched (cond,k2) gains
    if (dpop(2,1) >= 0d0) call bad("source (cond,k1) did not lose population")
    if (dpop(2,2) <= 0d0) call bad("energy-matched (cond,k2) did not gain (no inter-k transfer)")
    ! the dominant transfer is k1->k2 (the resonant one), not k1->k3 (off-resonant)
    if (dpop(2,2) <= dpop(2,3)) call bad("resonant k2 should gain more than off-resonant k3")
    ! valence is full -> Pauli-blocked, should not gain
    call chk("full valence Pauli-blocked (no gain)", dpop(1,1), 0d0, 1d-9)

    ! (4) gamma = 0: a phonon that matches nothing -> exact no-op
    block
        real(8) :: hw_big(nph), dz(nba, nk)
        hw_big(1) = 10d0                             ! no pair is within ~sigma of 10 Ha
        call eph_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, evbm, &
                             nph, hw_big, wrel, nb_bose, kt_bath, .false., &
                             nu_sat, nu_eps0, nu_n, &
                             sigma, tau, dz)
        call chk("gamma=0 -> no-op (max|dpop|)", maxval(abs(dz)), 0d0, 1d-12)
    end block

    ! (5) empty system -> no-op
    block
        real(8) :: f0(nba, nk), dz(nba, nk)
        f0 = 0d0
        call eph_interk_dpop(nk, nba, eval, f0, occ_max, a2half, ecbm, evbm, &
                             nph, hw, wrel, nb_bose, kt_bath, .false., nu_sat, nu_eps0, nu_n, &
                             sigma, tau, dz)
        call chk("empty -> no-op (max|dpop|)", maxval(abs(dz)), 0d0, 1d-12)
    end block

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (inter-k e-ph net-dpop: trace-conserving, bounded, '// &
                       'inter-k resonant transfer downhill, gamma=0 no-op)'
        call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'; call exit(1)
    end if

contains
    subroutine chk(name, got, want, tol)
        character(*), intent(in) :: name
        real(8), intent(in) :: got, want, tol
        if (abs(got-want) > tol*max(1d0,abs(want))) then
            write(*,'(a,a,a,es16.8,a,es16.8)') '  FAIL: ',name,' got=',got,' want=',want
            nfail = nfail + 1
        end if
    end subroutine chk
    subroutine bad(name)
        character(*), intent(in) :: name
        write(*,'(a,a)') '  FAIL: ', name
        nfail = nfail + 1
    end subroutine bad
end program test_eph_interk_cptp
