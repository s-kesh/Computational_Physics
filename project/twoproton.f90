! All lengths are in Angstrom and energies in eV
module twoproton
        implicit none

        character, parameter :: TAB = achar(9)
        character, parameter :: N_LINE = achar(10)

        ! physical constanst
        real, parameter :: a0 = 0.529           ! unit A⁰
        real, parameter :: e = 27.2             ! unit eV hatree energy

        ! Inputs
        integer, parameter :: ncoord = 3


        contains
                ! Calculate variation a
                subroutine cala(S, a)
                        implicit none
                        real, intent(in) :: S
                        real, intent(out) :: a
                        real :: aold
                        a = 0.5
                        aold = 0
                        do
                                if (abs(a - aold) .lt. 1.0E-6)  then
                                        exit
                                endif
                                aold = a
                                a = 1. / (1. + exp(-S / aold))
                        enddo
                end subroutine cala

                ! Generate random coordinates of system
                subroutine genconfig(r1, r2)
                        implicit none
                        real, dimension(:), intent(out) :: r1, r2
                        call random_number(r1)
                        call random_number(r2)
                        r1 = r1 - 0.5
                        r2 = r2 - 0.5
!                        print *, "r1", r1
!                        print *, "r2", r2
                end subroutine genconfig

                subroutine cal_para(r1, r2, S, a, beta, alpha  &
                                , r1l, r1r, r2l, r2r, r12 &
                                , r1l_norm, r1r_norm, r2l_norm, r2r_norm, r12_norm &
                                , phi1L, phi1R, phi1, phi2L, phi2R, phi2, phi12)
                        implicit none
                        real, dimension(:), intent(in) :: r1, r2
                        real, intent(in) :: S, a, beta, alpha
                        real, dimension(:), intent(out) :: r1l, r1r, r2l, r2r, r12
                        real, intent(out) :: r1l_norm, r1r_norm, r2l_norm, r2r_norm, r12_norm
                        real, intent(out) :: phi1L, phi1R, phi2L, phi2R
                        real, intent(out) :: phi1, phi2, phi12

                        real, dimension(size(r1)) :: L, R

                        L = [-S / 2., 0., 0.]
                        R = [S / 2., 0., 0.]

                        r1l = r1 - L
                        r1r = r1 - R
                        r2l = r2 - L
                        r2r = r2 - R
                        r12 = r1 - r2


                        r1l_norm = norm2(r1l)
                        r1r_norm = norm2(r1r)
                        r2l_norm = norm2(r2l)
                        r2r_norm = norm2(r2r)
                        r12_norm = norm2(r12)

                        phi1L = exp(-r1l_norm / a)
                        phi1R = exp(-r1r_norm / a)
                        phi2L = exp(-r2l_norm / a)
                        phi2R = exp(-r2r_norm / a)

                        phi1 = phi1L + phi1R
                        phi2 = phi2L + phi2R
                        phi12 = exp( r12_norm / (alpha * (1. + beta * r12_norm)))

                end subroutine

                subroutine energyloc(r1, r2, S, a, beta, alpha, E_local, delbeta, delE)
                        implicit none
                        real, dimension(:), intent(in) :: r1, r2
                        real, intent(in) :: S, a, beta, alpha
                        real, intent(out) :: E_local, delbeta, delE

                        real, dimension(size(r1)) :: r1l, r1r, r2l, r2r, r12
                        real :: r1l_norm, r1r_norm, r2l_norm, r2r_norm, r12_norm
                        real :: phi1L, phi1R, phi2L, phi2R
                        real :: phi1, phi2, phi12


                        real, dimension(size(r1)) :: r1l_unit, r1r_unit, r2l_unit, r2r_unit, r12_unit
                        real :: kinetic , potential
                        real :: e1, e2, cor, cross
                        real :: cross1, cross2

                        call cal_para(r1, r2, S, a, beta, alpha  &
                                , r1l, r1r, r2l, r2r, r12 &
                                , r1l_norm, r1r_norm, r2l_norm, r2r_norm, r12_norm &
                                , phi1L, phi1R, phi1, phi2L, phi2R, phi2, phi12)

                        r1l_unit = r1l / r1l_norm
                        r1r_unit = r1r / r1r_norm
                        r2l_unit = r2l / r2l_norm
                        r2r_unit = r2r / r2r_norm
                        r12_unit = r12 / r12_norm

                        potential = -(1. / r1l_norm + 1. / r1r_norm &
                                + 1. / r2l_norm + 1. / r2r_norm) &
                                + 1. / r12_norm

                        e1 = (phi1L / r1l_norm + phi1R / r1r_norm) &
                                / (a * phi1)
                        e2 = (phi12 / r2l_norm + phi2R / r2r_norm) &
                                / (a * phi2)

                        cross1 = (phi1L * sum(r1l_unit * r12_unit) &
                                + phi1R * sum(r1r_unit * r12_unit)) &
                                / (phi1 * 2. * a * (1 + beta * r12_norm)**2)
                        cross2 = (phi2L * sum(r2l_unit * r12_unit) &
                                + phi2R * sum(r2r_unit * r12_unit)) &
                                / (phi2 * 2. * a * (1 + beta * r12_norm)**2)
                        cross = cross1 - cross2

                        cor = (-1. / a**2) &
                                - ((1. + 4. * beta) * r12_norm + 4.) &
                                / (4. * r12_norm * (beta * r12_norm + 1.)**4)

                        kinetic = e1 + e2 + cross + cor

                        E_local = kinetic + potential

                        delbeta = -r12_norm**2 / (alpha * (1. + beta * r12_norm)**2)
                        delE = E_local * delbeta

!                        print *, "kinetic", kinetic
!                        print *, "potential", potential
!                        print *, "Eloal", E_local
                end subroutine

                real function prob(r1, r2, S, a, beta, alpha)
                        implicit none
                        real, dimension(:), intent(in) :: r1, r2
                        real, intent(in) :: S, a, beta, alpha

                        real, dimension(size(r1)) :: r1l, r1r, r2l, r2r, r12
                        real :: r1l_norm, r1r_norm, r2l_norm, r2r_norm, r12_norm
                        real :: phi1L, phi1R, phi2L, phi2R
                        real :: phi1, phi2, phi12

                        call cal_para(r1, r2, S, a, beta, alpha  &
                                , r1l, r1r, r2l, r2r, r12 &
                                , r1l_norm, r1r_norm, r2l_norm, r2r_norm, r12_norm &
                                , phi1L, phi1R, phi1, phi2L, phi2R, phi2, phi12)

!                        print *, "r1l", r1l
!                        print *, "r1r", r1r
!                        print *, "r2l", r2l
!                        print *, "r2r", r2r
!                        print *, "r12", r12

                        prob = (phi1 * phi2 * phi12)**2
!                        print *, "phi1", phi1
!                        print *, "phi2", phi2
!                        print *, "phi12", phi12
!                        print *, "prob", prob
                end function

                subroutine metrosinglestep(r1, r2, w, S, a, beta, alpha, accept, stepsize)
                        implicit none
                        real, intent(in) :: S, a, beta, alpha, stepsize
                        real, dimension(:), intent(inout) :: r1, r2
                        real, intent(inout) :: w
                        integer, intent(inout) :: accept

                        real, dimension(ncoord) :: r1old, r2old
                        real :: wold
                        real, dimension(ncoord) :: ran_no
                        real :: ranno

                        r1old = r1
                        r2old = r2
                        wold = w

                        call random_number(ran_no)
                        r1 = r1old + stepsize * (ran_no - 0.5)
                        call random_number(ran_no)
                        r2 = r2old + stepsize * (ran_no - 0.5)

                        call random_number(ranno)
                        w = prob(r1, r2, S, a, beta, alpha)
                        if (w / wold .gt. ranno) then
                                accept = accept + 1
                        else
                                r1 = r1old
                                r2 = r2old
                                w = wold
                        endif
                end subroutine metrosinglestep

                subroutine metropolis(r1, r2, S, a, beta, alpha &
                                , accept, stepsize, Nthermal, Ntotal    &
                                , Energy, delb, dele)
                        implicit none
                        integer, intent(in) :: Nthermal, Ntotal
                        real, intent(in) :: S, a, beta, alpha, stepsize
                        real, dimension(:), intent(inout) :: r1, r2
                        integer, intent(out) :: accept
                        real, dimension(:), intent(out) :: Energy, delb, dele

                        real :: w
                        integer :: i

                        accept = 0
                        w = prob(r1, r2, S, a, beta, alpha)
                        print *, "prob", w
                        do i = 1, Nthermal
                                call metrosinglestep(r1, r2, w, S, a, beta, alpha, accept, stepsize)
                        enddo

                        accept = 0
                        call energyloc(r1, r2, S, a, beta, alpha, Energy(1), delb(1), dele(1))
                        do i = 2, Ntotal
                                call metrosinglestep(r1, r2, w, S, a, beta, alpha, accept, stepsize)
                                call energyloc(r1, r2, S, a, beta, alpha, Energy(i), delb(i), dele(i))
                        enddo
                end subroutine metropolis

                real function min_beta(beta, Energy, delb, delE)
                        implicit none
                        real, dimension(:), intent(in) :: Energy, delb, delE
                        real, intent(inout) :: beta

                        min_beta = beta - 2 * (avg(delE) - avg(Energy) * avg(delb))
                end function

                real function avg(x)
                        implicit none
                        real, dimension(:), intent(in) :: x

                        avg = sum(x) / size(x)
                end function

                real function sd(x)
                        implicit none
                        real, dimension(:), intent(in) :: x

                        sd = sqrt(avg(x * x) - avg(x)**2)
                end function

end module
