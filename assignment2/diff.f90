! This program will solve ODE.
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204


! Define function in this module
! Keep Same format. Just allocate function dimentions according to problem

module function
        implicit none
        contains
                function f(x,y)
                        implicit none
                        real, allocatable :: f(:)
                        real, intent(in) :: x, y(:)
                        allocate(f(1))
                        f = -2*y
                end function

                function g(x,y)
                        implicit none
                        real, allocatable :: g(:)
                        real, intent(in) :: x, y(:)
                        allocate(g(size(y)))
                        g = - x*y
                end function

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
                subroutine pmatrix(mat)
                        implicit none
                        real, intent(in) :: mat(:,:)
                        integer :: i, j
                        do i = 1, size(mat, 1)
                                write (*, 10) (mat(i, j), achar(9), j = 1, size(mat, 2))
                        enddo
                        10      format(10F20.8)
                end subroutine pmatrix
end module print_matrix

! Main Program
program diffeqnsolver
!        use ode1d, only : odesolve ! Using only solving subroutine
        use ode
        use function
        use print_matrix
        implicit none
        real :: t, t_final, step_size
        real, allocatable :: ysolve(:,:)
        real, allocatable :: y(:)
        character(50) :: filename
        integer :: i, choice
        integer, parameter :: d = 2

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
!        print *, "Enter Output File Name"
!        read (*,*) filename

        print *, "Please Select the algorithim."
        print *, "Type '1' for 'Euler Method."
        print *, "Type '2' for 'Averging Method."
        print *, "Type '3' for 'Midpoint Method."
        print *, "Type '4' for 'Runge Kutta Method."
        read (*,*) choice

        if (choice .eq. 1)      then
                ysolve = EUL(shm, t, y, t_final, step_size)
        else if (choice .eq. 2)      then
                ysolve = AVG(shm, t, y, t_final, step_size)
        else if (choice .eq. 3)      then
                ysolve = MID(shm, t, y, t_final, step_size)
        else if (choice .eq. 4)      then
                ysolve = RK4(shm, t, y, t_final, step_size)
        else
                print *, "Error!!!!! You have selected wrong option."
                print *, "Exiting....."
                stop
        endif

        call pmatrix(ysolve)
end program diffeqnsolver
