program mc_inte_weight
        implicit none
        character, parameter :: TAB = achar(9)
        character, parameter :: N_LINE = achar(10)
        real :: no, sumf, sumf2, fx, sumfw, sumf2w, fxw
        integer :: i, j, N
        real :: inte, sigma
        real :: intew, sigmaw

        write(*, '(A,A,A,A,A)', advance='no') "Without Weight Function", TAB, TAB, &
                "With Weight Function w(x) = (4 - 2x) / 3 ", N_LINE
        write(*, '(A,A)', advance='no') "------------------------------------&
                ---------------------------------------", N_LINE
        write(*, '(A,A,A,A,A,A,A,A,A,A,A)', advance='no') &
                "N", TAB, "integral", TAB, "Sigma", &
                TAB, TAB, "integral", TAB, "Sigma", N_LINE

        do j = 1, 9
                sumf = 0
                sumf2 = 0
                sumfw = 0
                sumf2w = 0
                N = 4**j

                do i = 1, N
                        call random(no)
                        fx = f(no)
                        fxw = f(x(no)) / w(x(no))
                        sumf = sumf + fx
                        sumf2 = sumf2 + fx**2
                        sumfw = sumfw + fxw
                        sumf2w = sumf2w + fxw**2
                end do

                inte = sumf / N
                sigma = sqrt((sumf2 / N - inte**2) / N )
                intew = sumfw / N
                sigmaw = sqrt((sumf2w / N - intew**2) / N )

                write(*, '(*(g0.7),A,*(g0.7),A,*(g0.7),A,*(g0.7),A,*(g0.7),A)', advance='no') &
                        N, TAB, inte, TAB, sigma, &
                        TAB, intew, TAB, sigmaw, N_LINE
        end do

        write(*, '(A,A)', advance='no') "------------------------------------&
                ---------------------------------------", N_LINE

        contains
                real function f(x)
                        implicit none
                        real, intent(in) :: x
                        f = 1 / (1 + x**2)
                end function

                real function w(x)
                        implicit none
                        real, intent(in) :: x
                        w = (4 - 2 * x) /3
                end function

                real function x(y)
                        implicit none
                        real, intent(in) :: y
                        x = 2 - sqrt(4 - 3 * y)
                end function

                subroutine random(x)
                        implicit none
                        real, intent(inout) :: x
                        call random_seed()
                        call random_number(x)
                end subroutine
end program
