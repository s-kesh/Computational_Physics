program main
        use globaldata
        use twoproton
        implicit none

        real :: a
        real, dimension(ncoord * 2, TotalWalkers) :: configuration
        real, dimension(NTotal) :: Elocal, delphibeta, Edelphibeta
        real, dimension(Nos, minimize_steps) :: betaarray, earray, evararray
        real, dimension(Nos) :: S
        integer :: i, j

        do i = 1, Nos
                S(i) = 0.5 + (i - 1) * 0.5
        enddo

        !$OMP PARALLEL DO NUM_THREADS(14)
        do i = 1, Nos
                a = 0.5
                call cala(S(i), a)

                betaarray = 0.4
                do j = 1, minimize_steps
                        call random_seed()
                        call genconfig(configuration)
                        call metropolis(configuration, S(i), a, betaarray(i, j), Elocal, delphibeta, Edelphibeta)
                        earray(i, j) = avg(Elocal(Nthermal: NTotal))
                        evararray(i, j) = var(Elocal(Nthermal: NTotal))
                        if (j .ne. minimize_steps) &
                                betaarray(i, j + 1) = min_beta(betaarray(i, j), earray(i, j), delphibeta, Edelphibeta)
                end do
        end do
        !$OMP END PARALLEL DO

        call finish_and_print(S, betaarray, earray, evararray)
end program main
