! Module to Solve First Order ODE
! It impliment 4 methods
!       1) Eular Method
!       2) Average Method
!       3) Midpoint Method
!       4) Rounge Kutta Method

module ode1d
        implicit none

        private EUL, AVG, MID, RK4
        public odesolve

        contains
                ! Single Step for Eular Method
                real function EUL(fun, x, y, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        real, intent(in) :: x, y, xf, h
                        EUL = h * fun(x, y)
                end function

                ! Single Step for Average Method
                real function AVG(fun, x, y, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        real, intent(in) :: x, y, xf, h
                        AVG = (h/2) * (fun(x,y) + fun(x + h, y + h * fun(x,y)))
                end function

                ! Single Step for Midpoint Method
                real function MID(fun, x, y, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        real, intent(in) :: x, y, xf,  h
                        MID = h * fun(x + h/2, y + (h/2) * fun(x, y))
                end function

                ! Single Step for Rounge Kutta Method
                real function RK4(fun, x, y, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        real, intent(in) :: x, y, xf, h
                        real :: k1, k2, k3, k4
                        k1 = fun(x,y)
                        k2 = fun(x + h/2,  y + (h/2) * k1)
                        k3 = fun(x + h/2,  y + (h/2) * k2)
                        k4 = fun(x + h, y + h * k3)
                        RK4 = (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)
                end function

                ! Subroutine to solve ode.
                ! It will take 6 agruments.
                ! <function> <algorithim> <x₀> <y₀> <x_final> <step_size>
                ! It will print a table out with:
                ! 1st column => Method Used
                ! 2nd column => x
                ! 3rd column => y
                subroutine odesolve(fun, alg, x0, y0, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        character(3), intent(in) :: alg
                        real, intent(in) :: x0, y0, xf, h
                        real :: x, y, del_y
                        x = x0
                        y = y0
                        print *, "Method  ","x-axis             ","y-axis"
                        do
                                print *, alg, x, y
                                if (x >= xf)        then
                                        stop
                                endif
                                if (alg .eq. "EUL")      then
                                        y = y + EUL(fun, x, y, xf, h)
                                else if (alg .eq. "AVG")      then
                                        y = y + AVG(fun, x, y, xf, h)
                                else if (alg .eq. "MID")      then
                                        y = y + MID(fun, x, y, xf, h)
                                else if (alg .eq. "RK4")      then
                                        y = y + RK4(fun, x, y, xf, h)
                                endif
                                x = x + h
                        enddo
                end subroutine odesolve

end module ode1d
