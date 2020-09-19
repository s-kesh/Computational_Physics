! Module to Solve ODE
! It impliment 4 methods
!       1) Eular Method
!       2) Average Method
!       3) Midpoint Method
!       4) Rounge Kutta Method

module ode
        implicit none
        public EUL, AVG, MID, RK4

        contains

                ! Classic Eular
                function EUL(fun, t0, y0, tf, h)
                        implicit none
                        real , allocatable :: EUL(:,:)
                        interface
                                function fun(t, y)
                                        real, allocatable :: fun(:)
                                        real, intent(in) :: t, y(:)
                                end function
                        end interface
                        real, intent(in) :: t0, y0(:), tf, h
                        integer :: i, N
                        real, allocatable :: yn(:)
                        real :: tn

                        tn = t0
                        N = int((tf - t0) / h)
                        allocate(EUL(N, size(y0)))
                        EUL(1, :) = y0

                        do i = 2, N
                                yn = EUL(i - 1, :)
                                EUL(i,:) = yn + h * fun(tn, yn)
                                tn = tn + h
                        enddo
                end function EUL

                ! Modified Eular for averaging
                function AVG(fun, t0, y0, tf, h)
                        implicit none
                        real , allocatable :: AVG(:,:)
                        interface
                                function fun(t, y)
                                        real, allocatable :: fun(:)
                                        real, intent(in) :: t, y(:)
                                end function
                        end interface
                        real, intent(in) :: t0, y0(:), tf, h
                        integer :: i, N
                        real :: tn
                        real, allocatable :: yn(:)

                        tn = t0
                        N = int((tf - t0) / h)
                        allocate(AVG(N, size(y0)))
                        AVG(1, :) = y0

                        do i = 2, N
                                yn = AVG(i - 1, :)
                                AVG(i, :) = yn + (h/2) * (fun(tn, yn) + fun(tn + h, yn + h * fun(tn, yn)))
                                tn = tn + h
                        enddo
                end function

                ! Midpoint Method
                function MID(fun, t0, y0, tf, h)
                        implicit none
                        real , allocatable :: MID(:,:)
                        interface
                                function fun(t, y)
                                        real, allocatable :: fun(:)
                                        real, intent(in) :: t, y(:)
                                end function
                        end interface
                        real, intent(in) :: t0, y0(:), tf, h
                        integer :: i, N
                        real, allocatable :: yn(:)
                        real :: tn

                        tn = t0
                        N = int((tf - t0) / h)
                        allocate(MID(N, size(y0)))
                        MID(1, :) = y0

                        do i = 2, N
                                yn = MID(i - 1, :)
                                MID(i, :) = yn + h * fun(tn + h/2, yn + (h/2) * fun(tn, yn))
                                tn = tn + h
                        enddo
                end function

                ! Classic Runge Kutta
                function RK4(fun, t0, y0, tf, h)
                        implicit none
                        real , allocatable :: RK4(:,:)
                        interface
                                function fun(t, y)
                                        real, allocatable :: fun(:)
                                        real, intent(in) :: t, y(:)
                                end function
                        end interface
                        real, intent(in) :: t0, y0(:), tf, h
                        integer :: i, N
                        real :: tn
                        real, allocatable :: k1(:), k2(:), k3(:), k4(:), yn(:)

                        tn = t0
                        N = int((tf - t0) / h)
                        allocate(RK4(N, size(y0)))
!                        allocate(k1(size(y0))); allocate(k2(size(y0))); allocate(k3(size(y0))); allocate(k4(size(y0)))
                        RK4(1, :) = y0

                        do i = 2, N
                                yn = RK4(i - 1, :)
                                k1 = fun(tn, yn)
                                k2 = fun(tn + h/2, yn + (h/2) * k1)
                                k3 = fun(tn + h/2, yn + (h/2) * k2)
                                k4 = fun(tn + h, yn + h * k3)
                                RK4(i, :) = yn + (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)
                                tn = tn + h
                        enddo
                end function
end module ode
