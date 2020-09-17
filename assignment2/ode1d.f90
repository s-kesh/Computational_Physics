module ode1d
        implicit none

        private SE, ME, RK
        public odesolve

        contains
                real function SE(fun, x, y, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        real, intent(in) :: x, y, xf, h
                        SE = h * fun(x, y)
                end function

                real function ME(fun, x, y, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        real, intent(in) :: x, y, xf,  h
                        ME = fun(x + h/2, y + (h/2) * fun(x, y))
                end function

                real function RK(fun, x, y, xf, h)
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
                        RK = (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)
                end function

                subroutine odesolve(fun, alg, x0, y0, xf, h)
                        implicit none
                        interface
                                real function fun(x, y)
                                        real, intent(in) :: x, y
                                end function
                        end interface
                        character(2), intent(in) :: alg
                        real, intent(in) :: x0, y0, xf, h
                        real :: x, y, del_y
                        x = x0
                        y = y0
                        do
                                print *, alg, x, y
                                if (x >= xf)        then
                                        stop
                                endif
                                if (alg .eq. "SE")      then
                                        y = y + SE(fun, x, y, xf, h)
                                else if (alg .eq. "ME")      then
                                        y = y + ME(fun, x, y, xf, h)
                                else if (alg .eq. "RK")      then
                                        y = y + RK(fun, x, y, xf, h)
                                endif
                                x = x + h
                        enddo
                end subroutine odesolve

end module ode1d
