! Define function in this module
! Keep Same format. Just allocate function dimentions according to problem

module function
        implicit none
        contains
                ! ODE ẏ = -2y
                ! Solution y = C exp(-2t)
                function f(t,y)
                        implicit none
                        real, allocatable :: f(:)
                        real, intent(in) :: t, y(:)
                        allocate(f(1))
                        f = -2*y
                end function

                ! ODE ẏ = -ty
                ! Solution  y = C exp(-t² / 2)
                function g(t,y)
                        implicit none
                        real, allocatable :: g(:)
                        real, intent(in) :: t, y(:)
                        allocate(g(size(y)))
                        g = - t*y
                end function

                ! ODE ẏ = -2y - exp(3t)
                ! Solution y = ⅕ exp(3t) + C exp(-2t)
                function p(t,y)
                        implicit none
                        real, allocatable :: p(:)
                        real, intent(in) :: t, y(:)
                        allocate(p(size(y)))
                        p = -2 * y - exp(3*t)
                end function

                ! SHM
                ! ODE mÿ = -ky; k/m = 4
                function shm(x, y)
                        implicit none
                        real, allocatable :: shm(:)
                        real, intent(in) :: x, y(:)
                        allocate(shm(size(y)))
                        shm(1)  = y(2)
                        shm(2) = - 4 * y(1)
                end function
end module function
