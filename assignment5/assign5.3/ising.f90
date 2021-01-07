module ising
        implicit none
        public
                integer :: N
                real, dimension(10) :: list
                integer, allocatable, dimension(:) :: spins

        contains

                ! created list for calculating rejection
                subroutine rlist(j, b, t)
                        implicit none
                        real, intent(in) :: j, b, t
                        integer :: p, q
                        integer, dimension(5) :: f
                        integer, dimension(2) :: alpha
                        f = [-4, -2, 0, 2, 4]
                        alpha = [-1, 1]

                        do p = 1, 2
                                do q = 1, 5
                                        list((p -1) * 5 + q) &
                                                = exp((2 * alpha(p) &
                                                * (j * f(q) + b)) / t)
                                enddo
                        enddo
                end subroutine

                ! for calculation of metropils probability
                real function listmap(a, f)
                        implicit none
                        integer , intent(in) :: a, f
                        if ((a .eq. -1) .and. (f .eq. -4)) then
                                listmap = list(1)
                        else if ((a .eq. -1) .and. (f .eq. -2)) then
                                listmap = list(2)
                        else if ((a .eq. -1) .and. (f .eq. 0)) then
                                listmap = list(3)
                        else if ((a .eq. -1) .and. (f .eq. 2)) then
                                listmap = list(4)
                        else if ((a .eq. -1) .and. (f .eq. 4)) then
                                listmap = list(5)
                        else if ((a .eq. 1) .and. (f .eq. -4)) then
                                listmap = list(6)
                        else if ((a .eq. 1) .and. (f .eq. -2)) then
                                listmap = list(7)
                        else if ((a .eq. 1) .and. (f .eq. 0)) then
                                listmap = list(8)
                        else if ((a .eq. 1) .and. (f .eq. 2)) then
                                listmap = list(9)
                        else if ((a .eq. 1) .and. (f .eq. 4)) then
                                listmap = list(10)
                        endif
                end function

                ! Calculating f
                integer function fcal(x)
                        implicit none
                        integer, intent(in) :: x
                        integer, dimension(4) :: y

                        y(1) = ((x - 1) / N) * N + modu(x - 1, N)
                        y(2) = ((x - 1) / N) * N + modu(x + 1, N)
                        y(3) = modu(x - N, N * N)
                        y(4) = modu(x + N, N * N)

                        fcal = spins(y(1)) + spins(y(2)) &
                                + spins(y(3)) + spins(y(4))
                end function

                ! Custom made modulus function for calculating
                ! neighbours
                integer function modu(x, y)
                        implicit none
                        integer, intent(in) :: x, y
                        integer :: t1, t2
                        t1 = x / y
                        t2 = x - t1 * y
                        if (t2 .eq. 0 .or. t2 .lt. 0) then
                                modu = y + t2
                        else
                                modu = t2
                        endif
                end function

                ! intialize random spin system
                subroutine intialize_random()
                        implicit none
                        real, dimension(N * N) :: eta

                        allocate(spins(N * N))
                        call random_number(eta)
                        spins = floor(eta * 2) * 2  - 1
                end subroutine intialize_random

                ! cleaning up spins
                subroutine cleanup()
                        implicit none
                        deallocate(spins)
                end subroutine cleanup

                ! Metropils Algorithim
                subroutine mertopolis()
                        implicit none
                        integer :: salpha
                        integer :: x, f
                        real :: neta

                        call random_number(neta)
                        x = 1 + floor(N * N * neta)
                        salpha = -1 * spins(x)
                        f = fcal(x)

                        call random_number(neta)
                        if ( neta .lt. listmap(salpha, f)) then
                                spins(x) = salpha
                        endif
                end subroutine mertopolis

                ! A sweep. Runs Metropolis Algorithim N times
                subroutine sweep()
                        implicit none
                        integer :: i

                        do i = 1, N * N
                                call mertopolis()
                        enddo
                end subroutine sweep


                !! Thermal Observations

                real function Mag()
                        implicit none
                        Mag = abs(sum(spins)) / (N * N * 1.0)
                end function

                real function Ener(j, b)
                        implicit none
                        real, intent(in) :: j, b
                        integer :: p

                        Ener = 0.0
                        do p = 1, N*N
                                Ener = Ener + (-0.5 * j * spins(p) * fcal(p) - b * spins(p))
                        enddo
                        Ener =  Ener / (N * N * 1.0)
                end function

                real function Chi(T)
                        implicit none
                        real, intent(in) :: T
                        real :: msq

                        msq = (sum(spins**2)) / (N * N * 1.0)
                        chi = (1 / T) * ( msq - (Mag())**2)
                end function

                real function Cv(j , b, T)
                        implicit none
                        real, intent(in) :: j, b, T
                        real :: esqr
                        integer :: p

                        esqr = 0
                        do p = 1, N*N
                                esqr = esqr + (-0.5 * j * spins(p) * fcal(p) - b * spins(p))**2
                        enddo
                        esqr = esqr / (N * N * 1.0)
                        Cv = (1 / T**2) * (esqr - Ener(j, b)**2)
                end function
end module ising
