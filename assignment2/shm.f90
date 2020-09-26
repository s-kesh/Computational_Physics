! This program will SHM
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204

program simpleharmonicmotion
        use ode         ! ODE Solver Module
        ! shm(x, y(:)); ÿ = - 4 x
        use function, only : shm
        
        implicit none
        real, allocatable :: y(:), ysolve(:,:)
        real :: t, t_final, h
        integer :: i, choice
        integer, parameter :: d = 2
        character(50) :: hchar, cchar

        allocate(y(d))

        t = 0
        y(:) = [1, 0]
        t_final = 15

        if (command_argument_count() .ne. 2)    then
                print *, "Please run program again with command line arguments."
                print *, "Usage:"
                print *, "<program_name> <step_size> <choice>"
                print *, "use choice = 1 for Euler Method"
                print *, "use choice = 2 for Rounge-Kutta Method"
                stop
        else
                call get_command_argument(1, hchar)
                call get_command_argument(2, cchar)
                read (hchar,*)  h
                read (cchar,*)  choice
        endif

        if (choice .eq. 1)      then
                ysolve = EUL(shm, t, y, t_final, h)
        else if (choice .eq. 2) then
                ysolve = RK4(shm, t, y, t_final, h)
        else
                print *, "Please use choice either 1 or 2"
                stop
        endif

        do i = 1, size(ysolve,1)
                print *, t + (i-1) * h, ysolve(i, :)
        enddo

end program simpleharmonicmotion
