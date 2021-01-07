program main
        use twoproton
        implicit none

        ! variational parameters
        ! Two are fixed beforehand
        ! Reducing no of variational parameters to be 2
        real, parameter :: alpha = 2., S = 1.411
        real :: a, beta

        integer, parameter :: N = 1000000
        integer, parameter :: Nth = 10000
        integer, parameter :: minimize_steps = 40
        real, dimension(ncoord) :: r1, r2
        real, dimension(:), allocatable :: Elocal, delphibeta, Edelphibeta
        real :: Energy, Energysd
        integer :: accept
        real :: step

        integer :: i


        ! Calculate variational parameter a using coulumbs cups condition
        call cala(S, a)

!        write (*, '(A, A, A, A)', advance='no') "stepsize", TAB, "Accept", N_LINE
        step = 1.55
        beta = 0.4
        do i = 1, minimize_steps
                ! Generate position of electrons randomly
                ! And calculating probabiliy and local energy of that configuration
                call random_seed()
                call genconfig(r1, r2)
                allocate(Elocal(N))
                allocate(delphibeta(N))
                allocate(Edelphibeta(N))
                call metropolis(r1, r2, S, a, beta, alpha, accept, step, Nth, N, Elocal, delphibeta, Edelphibeta)
                beta = min_beta(beta, Elocal, delphibeta, Edelphibeta)
                Energy = avg(Elocal)
                Energysd = sd(Elocal)
                write (*, '(*(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no')  &
                        beta, TAB, Energy, TAB, Energysd, N_LINE
                deallocate(Elocal)
                deallocate(delphibeta)
                deallocate(Edelphibeta)
        enddo

!        write (*, '(*(g0.7), A, *(g0.7), A)', advance='no') step, TAB, accept / (1. * N), N_LINE


!        write (*, '(*(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no') &
!                beta, TAB, Energy, TAB, accept / (1. * N), N_LINE

end program main
