program ising_model
        use ising
        implicit none
        character, parameter :: TAB = achar(9)
        character, parameter :: N_LINE = achar(10)
        integer, parameter :: NX = 20
        integer, parameter :: NY = 20
        type(smatrix) :: S
        real, dimension(10) :: r
        real, parameter :: magnet = 0
        real, parameter :: couple = 0
        integer :: i

        call random_seed()
        call intialize_random(S, NX, NY)
        call rlist(r, couple, magnet)

        do i = 1, 10000
                call sweep(S, r)
                write(*, '(A, *(g0.7), A, *(g0.7), A)', advance='no') &
                        "For Sweep no ", i, TAB, magnetization(S), N_LINE
        enddo

end program ising_model
