program ising_model
        use ising
        implicit none
        character, parameter :: TAB = achar(9)
        character, parameter :: N_LINE = achar(10)
        real, parameter :: magnet = 10
        real, parameter :: couple = 1
        real, parameter :: temprature = 4
        real :: t, avgmag, avgsep
        integer :: i
        
        N = 20

        call random_seed()
        call intialize_random()
        call rlist(couple, magnet, temprature)

        t = temprature


        avgmag = 0
        avgsep = 0
        do
                if (t < 0)    then
                        stop
                endif

                do i = 1, 2000
                        call sweep()
                enddo

                do i = 1, 10000
                        call sweep()
                        avgmag = avgmag + Mag() / (N * N)
                        avgsep = avgsep + chi(t) / (N * N)
                enddo
                avgmag = avgmag / 10000
                avgsep = avgsep / 10000
!                write(*, '(*(g0.7), A, *(g0.7), A, *(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no') &
!                        t, TAB, avgmag, TAB, avgener, TAB, avgsep, TAB, avgcv, N_LINE
                write(*, '(*(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no') &
                        t, TAB, avgmag, TAB, avgsep, N_LINE
                
                t = t - 0.05
                call rlist(couple, magnet, t)
        enddo
!
        call cleanup()

end program ising_model
