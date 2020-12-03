program ising_model
        use ising
        implicit none
        character, parameter :: TAB = achar(9)
        character, parameter :: N_LINE = achar(10)
        integer, parameter :: N = 20
        real, parameter :: magnet = 0
        real, parameter :: couple = -1
        real, parameter :: temprature = 10
        real :: t, avgmag

        type(smatrix) :: S
        real, dimension(10) :: r
        integer :: i

        call random_seed()
        call intialize_random(S, N)
        call rlist(r, couple, magnet, temprature)


        t = temprature

!        do i = 1, 1000
!                call sweep(S, r)
!                write(*, '(*(g0.7), A)', advance='no') abs(magnetization(S)), N_LINE
!        enddo

        avgmag = 0
        do
                if (t < 0)    then
                        stop
                endif

                do i = 1, 1000
                        call sweep(S, r)
                enddo

                do i = 1, 10000
                        call sweep(S, r)
                        avgmag = avgmag + abs(magnetization(S))
                enddo
                avgmag = avgmag / 10000
                write(*, '(*(g0.7), A, *(g0.7), A)', advance='no') &
                        t, TAB, avgmag, N_LINE
                
                t = t - 0.1
                call rlist(r, couple, magnet, t)
        enddo

end program ising_model
