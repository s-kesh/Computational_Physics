! This program will solve ODE.
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204


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

module print_matrix
        implicit none
        contains
                subroutine pmatrix(t0, h, mat)
                        implicit none
                        real, intent(in) :: mat(:,:), t0, h
                        integer :: i, j, N
                        do i = 1, size(mat, 1)
                                print *, t0 + (i-1) * h, mat(i, :)
                        enddo
                end subroutine pmatrix
end module print_matrix

! Main Program
program diffeqnsolver
        use ode
        use function
        use print_matrix
        implicit none
        real :: t, t_final, step_size
        real, allocatable :: ysolve(:,:)
        real, allocatable :: y(:)
        character(50) :: filename
        integer :: i, choice
        integer, parameter :: d = 1

        allocate(y(d))

        print *, "This Program will solve given ode in function"
        print *, "Please enter intial time in sec"
        read (*,*) t
        print *, "Please Enter Initial Value at time =", t, "sec"
        do i = 1, size(y)
                read (*,*) y(i)
        enddo
        print *, "Please define step size"
        read (*,*) step_size
        print *, "Please enter final value of x"
        read (*,*) t_final

        print *, "Please Select the algorithim."
        print *, "Type '1' for 'Euler Method."
        print *, "Type '2' for 'Averging Method."
        print *, "Type '3' for 'Midpoint Method."
        print *, "Type '4' for 'Runge Kutta Method."
        read (*,*) choice

        print *, "===================================================="
        print *, "Report :"
        print *, "Intial Value at t =", t
        print *, (y(i), i = 1, size(y))
        print *, "step_size:", step_size
        print *, "Final time:", t_final

        if (choice .eq. 1)      then
                ysolve = EUL(g, t, y, t_final, step_size)
        else if (choice .eq. 2)      then
                ysolve = AVG(g, t, y, t_final, step_size)
        else if (choice .eq. 3)      then
                ysolve = MID(g, t, y, t_final, step_size)
        else if (choice .eq. 4)      then
                ysolve = RK4(g, t, y, t_final, step_size)
        else
                print *, "Error!!!!! You have selected wrong option."
                print *, "Exiting....."
                stop
        endif

        print *, "Method no:", choice
        print *, "===================================================="

        call pmatrix(t, step_size, ysolve)

end program diffeqnsolver
