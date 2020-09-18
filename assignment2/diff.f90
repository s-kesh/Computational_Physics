! This program will solve ODE.
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204


! Define function in this module
! You can define as much function as you want.
module function
        implicit none
        Public f
        contains
                real function f(x,y)
                        implicit none
                        real, intent(in) :: x, y
                        f = -2*y
                end function
end module function

! Main Program
program diffeqnsolver
        use ode1d, only : odesolve ! Using only solving subroutine
        use function
        implicit none
        real :: x, y, x_final, step_size
        character(100) :: xchar, ychar, xfchar, hchar
        integer :: choice

        if (command_argument_count() .ne. 4)    then
                print *, "Usage:"
                print *, "<program_name> <x_intial> <y_intial> <step_size> <x_final>"
                print *, "Example:"
                print *, "./diff 0 1 0.01 4"
            !    print *, "This Program will solve given ode in function"
            !    print *, "Please enter initial condition: x, y"
            !    read (*,*) x, y
            !    print *, "Please define step size"
            !    read (*,*) step_size
            !    print *, "Please enter final value of x"
            !    read (*,*) x_final
        else
                call get_command_argument(1, xchar)
                call get_command_argument(2, ychar)
                call get_command_argument(3, hchar)
                call get_command_argument(4, xfchar)
                read (xchar,*) x
                read (ychar,*) y
                read (hchar,*) step_size
                read (xfchar,*) x_final

                print *, "Please Select the algorithim."
                print *, "Type '1' for 'Euler Method."
                print *, "Type '2' for 'Averging Method."
                print *, "Type '3' for 'Midpoint Method."
                print *, "Type '4' for 'Runge Kutta Method."
                read (*,*) choice

                if (choice .eq. 1)      then
                        call odesolve(f, "EUL", x, y, x_final, step_size)
                else if (choice .eq. 2)      then
                        call odesolve(f, "AVG", x, y, x_final, step_size)
                else if (choice .eq. 3)      then
                        call odesolve(f, "MID", x, y, x_final, step_size)
                else if (choice .eq. 4)      then
                        call odesolve(f, "RK4", x, y, x_final, step_size)
                else
                        print *, "Error!!!!! You have selected wrong option."
                        print *, "Exiting....."
                endif
        endif

end program diffeqnsolver
