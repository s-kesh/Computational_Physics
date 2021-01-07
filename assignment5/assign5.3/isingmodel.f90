program ising_model
        use ising
        implicit none
        character, parameter :: TAB = achar(9)
        character, parameter :: N_LINE = achar(10)
        real, parameter :: magnet = 0
        real, parameter :: couple = 1
        real, parameter :: temprature = 6
        real :: t, avgmag, avgsep, avgener, avgcv
        integer :: i
        !real :: E

        N = 20

        call random_seed()
        call intialize_random()
        call rlist(couple, magnet, temprature)

!        do i = 1, 5000
!                call sweep()
!                write(*, '(*(g0.7), A, *(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no') &
!                        Mag(), TAB, Ener(couple, magnet), TAB, &
!                        Chi(temprature), TAB, Cv(couple, magnet, temprature), N_LINE
!        enddo

        t = temprature


        avgmag = 0
        avgsep = 0
        avgener = 0
        avgcv = 0

        write(*, '(A, A, A, A, A, A, A, A, A, A)', advance='no') &
                "Temprature", TAB, "Magnetization", TAB, "Energy", TAB, &
                "Susceptibility", TAB, "Specific Heat", N_LINE
        do
                if (t < 0)    then
                        stop
                endif

                do i = 1, 2000
                        call sweep()
                enddo

                do i = 1, 10000
                        call sweep()
                        avgmag = avgmag + Mag()
                        avgsep = avgsep + Chi(t)
                        avgener = avgener + Ener(couple, magnet)
                        avgcv = avgcv + cv(couple, magnet, t)
                enddo
                avgmag = avgmag / 10000
                avgsep = avgsep / 10000
                avgener = avgener / 10000
                avgcv = avgcv / 10000
                write(*, '(*(g0.7), A, *(g0.7), A, *(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no') &
                        t, TAB, avgmag, TAB, avgener, TAB, avgsep, TAB, avgcv, N_LINE
!                write(*, '(*(g0.7), A, *(g0.7), A, *(g0.7), A)', advance='no') &
!                        t, TAB, avgmag, TAB, avgsep, N_LINE

                t = t - 0.01
                call rlist(couple, magnet, t)
        enddo
!
        call cleanup()

end program ising_model
