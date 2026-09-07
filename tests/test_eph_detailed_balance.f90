!
!  test_eph_detailed_balance.f90  -  the e-ph ring must leave a thermal sheet alone.
!
!  A relaxation channel coupled to a bath at kT has exactly one fixed point: the
!  Fermi-Dirac distribution at that kT. Start there and the net population change must
!  vanish. Nothing in the CPTP tests catches a violation of this -- trace is still
!  conserved, populations still stay in [0, occ_max], transfers still go to
!  energy-matched partners -- yet a channel that fails it pumps energy into the carriers
!  with no field at all.
!
!  It did. The Gaussian energy matching has width sigma (finite mesh, finite lifetime),
!  so a source is connected to partners at |dE| anywhere within a few sigma of hw_p,
!  while the emission/absorption split was taken once per MODE from N_B(hw_p). A pair
!  that actually transfers delta was then weighted by exp(hw_p/kT) instead of
!  exp(delta/kT), and for delta >> hw_p the upward rate was too large by
!  exp((delta - hw_p)/kT). On graphene, whose appended acoustic mode is 5.4 meV against
!  a 0.1 eV search width, that ran a doped sheet from 300 K to ~2300 K in 350 fs with
!  the drive switched off (x14 README SS7.16).
!
!  Checks, on a Dirac-like spectrum held at an exact FD(300 K):
!    1) db_realized = .true.: the net energy change is < 1e-3 of what the historical
!       path produces, and shrinks as sigma does,
!    2) it is stationary in the STRONG sense -- |sum E dpop| below a tight absolute
!       bound at the production sigma = 0.1 eV,
!    3) trace is still conserved exactly (the fix must not cost CPTP),
!    4) db_realized = .false. still heats, i.e. the gate really does select the old
!       behaviour for the gapped materials that were validated with it.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_eph_detailed_balance
    use sbe_superres_ssbe, only: eph_interk_dpop, bose_factor
    implicit none
    integer, parameter :: nk = 60, nba = 2, nph = 1
    real(8), parameter :: HA_EV = 27.211386245988d0
    real(8), parameter :: KB_AU = 3.166811563d-6          ! Ha/K
    real(8) :: eval(nba, nk), f(nba, nk), dpop(nba, nk)
    real(8) :: hw(nph), wrel(nph), nb(nph)
    real(8) :: mu, kT, occ_max, e, sig, tau, nu_sat, eps0, nun
    real(8) :: dE_fix_100, dE_fix_020, dE_old_100, tr
    integer :: ik, a, nfail

    nfail = 0
    occ_max = 2d0
    mu = 0.6d0 / HA_EV                       ! E_F = 0.6 eV, as in x14
    kT = KB_AU * 300d0
    hw(1) = 5.39d-3 / HA_EV                  ! the appended acoustic mode
    wrel(1) = 1d0
    nb(1) = bose_factor(hw(1), kT)
    nu_sat = 1d0; eps0 = 0.8d0 / HA_EV; nun = 2d0
    tau = 1d-3

    do ik = 1, nk                            ! +-e_k, e_k over 0 .. 1.2 eV
        e = (dble(ik) - 0.5d0) / dble(nk) * (1.2d0 / HA_EV)
        eval(1, ik) = -e
        eval(2, ik) =  e
    end do
    do ik = 1, nk                            ! EXACT Fermi-Dirac at the bath temperature
        do a = 1, nba
            f(a, ik) = occ_max / (1d0 + exp(min(max((eval(a, ik) - mu) / kT, -60d0), 60d0)))
        end do
    end do

    sig = 0.100d0 / HA_EV
    call eph_interk_dpop(nk, nba, eval, f, occ_max, 0d0, 0d0, 0d0, &
                         nph, hw, wrel, nb, kT, .true., nu_sat, eps0, nun, sig, tau, dpop)
    dE_fix_100 = sum(eval * dpop) * HA_EV * 1d3
    tr = sum(dpop)
    call chk("trace conserved, realized-transfer balance", abs(tr), 0d0, 1d-14)

    call eph_interk_dpop(nk, nba, eval, f, occ_max, 0d0, 0d0, 0d0, &
                         nph, hw, wrel, nb, kT, .false., nu_sat, eps0, nun, sig, tau, dpop)
    dE_old_100 = sum(eval * dpop) * HA_EV * 1d3
    call chk("trace conserved, historical split", abs(sum(dpop)), 0d0, 1d-14)

    sig = 0.020d0 / HA_EV
    call eph_interk_dpop(nk, nba, eval, f, occ_max, 0d0, 0d0, 0d0, &
                         nph, hw, wrel, nb, kT, .true., nu_sat, eps0, nun, sig, tau, dpop)
    dE_fix_020 = sum(eval * dpop) * HA_EV * 1d3

    write(*, '(a)') '  a thermal sheet, left alone by the bath it is in equilibrium with:'
    write(*, '(a,es12.3,a)') '    sum E dpop, historical split, sigma = 0.1 eV : ', dE_old_100, ' meV/step'
    write(*, '(a,es12.3,a)') '    sum E dpop, realized transfer, sigma = 0.1 eV: ', dE_fix_100, ' meV/step'
    write(*, '(a,es12.3,a)') '    sum E dpop, realized transfer, sigma = 0.02  : ', dE_fix_020, ' meV/step'

    ! (1) the fix is at least three orders below the historical leak
    if (.not. (abs(dE_fix_100) < 1d-3 * abs(dE_old_100))) then
        write(*, '(a)') '  FAIL: realized-transfer balance is not << the historical leak'
        nfail = nfail + 1
    end if
    ! (2) stationary in absolute terms at the production width
    if (.not. (abs(dE_fix_100) < 1d-3)) then
        write(*, '(a)') '  FAIL: |sum E dpop| >= 1e-3 meV/step at sigma = 0.1 eV'
        nfail = nfail + 1
    end if
    ! (3) narrowing sigma cannot make it worse
    if (.not. (abs(dE_fix_020) <= abs(dE_fix_100))) then
        write(*, '(a)') '  FAIL: residual grows as sigma narrows'
        nfail = nfail + 1
    end if
    ! (4) the gate really selects the old behaviour
    if (.not. (dE_old_100 > 1d-2)) then
        write(*, '(a)') '  FAIL: db_realized=.false. no longer reproduces the historical heating'
        nfail = nfail + 1
    end if

    if (nfail == 0) then
        write(*, '(a)') 'PASS  (e-ph detailed balance: FD(T_bath) is stationary at the '// &
                        'production search width; the gate preserves the historical path)'
    else
        write(*, '(a,i0,a)') 'FAIL  (', nfail, ' check(s))'
        stop 1
    end if

contains
    subroutine chk(label, got, want, tol)
        character(*), intent(in) :: label
        real(8), intent(in) :: got, want, tol
        if (abs(got - want) > tol) then
            write(*, '(a,a,a,es12.4,a,es12.4)') '  FAIL: ', label, '  got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine chk
end program test_eph_detailed_balance
