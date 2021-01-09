module twoproton
        use global
        implicit none

        contains
                real function avg(x)
                        implicit none
                        real, dimension(:), intent(in) :: x
                        avg = sum(x) / size(x)
                end function

                real function var(x)
                        implicit none
                        real, dimension(:), intent(in) :: x
                        var = sum(x * x) / size(x) - avg(x)**2
                end function

                ! Calculate variation a
                subroutine cala(S, a)
                        implicit none
                        real, intent(in) :: S
                        real, intent(out) :: a
                        real :: aold
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
                subroutine genconfig(config)
                        implicit none
                        real, dimension(:, :), intent(out) :: config
                        call random_number(config)
                        config = config - 0.5
                end subroutine genconfig

                real function min_beta(beta, Energy, delb, delE)
                        implicit none
                        real, intent(in) :: beta, Energy
                        real, dimension(:), intent(in) :: delb, delE

                        min_beta = beta - gamm * 2 *  ( avg(delE(Nthermal: Ntotal)) &
                                - Energy * avg(delb(Nthermal: Ntotal)) )
                end function

                subroutine calpara(config, r1, r2, S, a, beta  &
                                , r1l, r1r, r2l, r2r, r12 &
                                , r1lnorm, r1rnorm, r2lnorm, r2rnorm, r12norm &
                                , phi1L, phi1R, phi1, phi2L, phi2R, phi2, phi12)
                        implicit none
                        real, dimension(:, :), intent(in) :: config
                        real, intent(in) :: S, a, beta
                        real, dimension(:, :), intent(out) :: r1, r2, r1l, r1r, r2l, r2r, r12
                        real, dimension(:), intent(out) &
                                :: r1lnorm, r1rnorm, r2lnorm, r2rnorm, r12norm &
                                , phi1L, phi1R, phi2L, phi2R, phi1, phi2, phi12

                        real, dimension(ncoord, TotalWalkers) :: proton, xunit

                        xunit = 0.
                        xunit(1, :) = 1.
                        proton = 0.5 * S * xunit

                        r1 = config(1:3, :)
                        r2 = config(4:6, :)
                        r1l = r1 + proton
                        r1r = r1 - proton
                        r2l = r2 + proton
                        r2r = r2 - proton
                        r12 = r1 - r2

                        r1lnorm = norm2(r1l, dim=1)
                        r1rnorm = norm2(r1r, dim=1)
                        r2lnorm = norm2(r2l, dim=1)
                        r2rnorm = norm2(r2r, dim=1)
                        r12norm = norm2(r12, dim=1)

                        phi1L = exp(-r1lnorm / a)
                        phi1R = exp(-r1rnorm / a)
                        phi2L = exp(-r2lnorm / a)
                        phi2R = exp(-r2rnorm / a)

                        phi1 = phi1L + phi1R
                        phi2 = phi2L + phi2R
                        phi12 = exp( r12norm / (alpha * (1. + beta * r12norm)))
                end subroutine

                subroutine energyloc(config, S, a, beta, Elocal, delphi, Edelphi)
                        implicit none
                        real, dimension(:, :), intent(in) :: config
                        real, intent(in) :: S, a, beta
                        real, intent(out) :: Elocal, delphi, Edelphi

                        real, dimension(ncoord, TotalWalkers) &
                                :: r1, r1L, r1R, r2, r2L, r2R, r12 &
                                , r1Lunit, r1Runit, r2Lunit, r2Runit, r12unit
                        real, dimension(TotalWalkers) &
                                :: phi1, phi1L, phi1R, phi2, phi2L, phi2R, phi12 &
                                , r1Lnorm, r1Rnorm, r2Lnorm, r2Rnorm, r12norm
                        real, dimension(TotalWalkers) :: kinetic, potential &
                                , e1, e2, cor, cross, cross1, cross2, cross3 &
                                , Etot, delphiarray

                        call calpara(config, r1, r2, s, a, beta &
                                , r1L, r1R, r2L, r2R, r12 &
                                , r1Lnorm, r1Rnorm, r2Lnorm, r2Rnorm, r12norm &
                                , phi1L, phi1R, phi1, phi2L, phi2R, phi2, phi12)

                        r1Lunit = r1L / spread(r1Lnorm, 1, ncoord)
                        r1Runit = r1R / spread(r1Rnorm, 1, ncoord)
                        r2Lunit = r2L / spread(r2Lnorm, 1, ncoord)
                        r2Runit = r2R / spread(r2Rnorm, 1, ncoord)
                        r12unit = r12 / spread(r12norm, 1, ncoord)

                        e1 = (phi1L / r1Lnorm + phi1R / r1Rnorm) / (a * phi1)
                        e2 = (phi2L / r2Lnorm + phi2R / r2Rnorm) / (a * phi2)

                        cross1 = (phi1L * sum(r1Lunit * r12unit, dim=1) &
                                + phi1R * sum(r1Runit * r12unit , dim=1)) / phi1
                        cross2 = (phi2L * sum(r2Lunit * r12unit, dim=1) &
                                + phi2R * sum(r2Runit * r12unit , dim=1)) / phi2
                        cross3 = 1. / (2. * a * (1.+ beta * r12norm)**2)

                        cross = (cross1 - cross2) * cross3
                        cor = - ((1. + 4. * beta) * r12norm + 4.) / (4. * r12norm * (1 + beta * r12norm)**4.)

                        potential = 1. / r12norm - 1. / r1Lnorm - 1. / r1Rnorm - 1. / r2Lnorm - 1. / r2Rnorm
                        kinetic = e1 + e2 + cross + cor
                        Etot = kinetic + potential + 1. / s - 1. / a**2

                        delphiarray = -r12norm**2 / (alpha *(1 + beta * r12norm)**2)

                        Elocal = avg(Etot)
                        delphi = avg(delphiarray)
                        Edelphi = avg(Etot * delphiarray)
                end subroutine

                function prob(config, S, a, beta)
                        implicit none
                        real, dimension(:, :), intent(in) :: config
                        real, intent(in) :: S, a, beta

                        real, dimension(size(config,1) / 2, TotalWalkers) &
                                :: r1, r2, r1l, r1r, r2l, r2r, r12
                        real, dimension(TotalWalkers) :: r1lnorm, r1rnorm, r2lnorm, r2rnorm, r12norm &
                                , phi1L, phi1R, phi2L, phi2R, phi1, phi2, phi12
                        real, dimension(TotalWalkers) :: prob

                        call calpara(config, r1, r2, S, a, beta &
                                , r1l, r1r, r2l, r2r, r12 &
                                , r1lnorm, r1rnorm, r2lnorm, r2rnorm, r12norm &
                                , phi1L, phi1R, phi1, phi2L, phi2R, phi2, phi12)

                        prob = phi1 * phi2 * phi12
                end function

                subroutine metropolis(config, S, a, beta, Energy, delphi, Edelphi)
                        implicit none
                        real, intent(in) :: S, a, beta
                        real, dimension(:, :), intent(inout) :: config
                        real, dimension(:), intent(out) :: Energy, delphi, Edelphi

                        real, dimension(ncoord * 2, TotalWalkers) :: confignew, ranno
                        real, dimension(TotalWalkers) :: ratio, accept
                        real :: ran, stepsize
                        integer :: step, i

                        stepsize = 1.
                        accept = 0
                        do step = 1, Ntotal
                                call random_number(ranno)
                                confignew = config + stepsize * (ranno - 0.5)
                                ratio = (prob(confignew, S, a, beta) / prob(config, S, a, beta))**2
                                do i = 1, TotalWalkers
                                        call random_number(ran)
                                        if ( ratio(i) .gt. ran) then
                                                accept(i) = accept(i) + 1
                                                config(:, i) = confignew(:, i)
                                        end if
                                end do
                                if (step .ge. Nthermal) then
                                        call energyloc(config, S, a, beta &
                                                , Energy(step), delphi(step), Edelphi(step))
                                end if

                                if (mod(step, 100) .eq. 0) then
                                        stepsize = stepsize * (sum(accept) / (50.0 * TotalWalkers))
                                        accept = 0
                                end if
                        end do
                end subroutine metropolis
end module
