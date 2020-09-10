! This program will find roots of particle in a box problem.
! Give energy/potential in eV while half width in A⁰
! Hint to give a segment :
!       z should be between nπ/2 → (n - 1)π/2
!       keep n odd for even state and even for odd state
!       so, given segment should be calculated using
!       Energy = V₀ - (z / z₀a)², where z₀ = 0.511974
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204

real function f(x, z, flag)
        ! Function need to solve.
        ! It has two parts `cos` wala for odd and `sin` wala for even
        implicit none
        real, intent(in) :: x, z
        integer, intent(in) :: flag
        if (flag == 1)  then
                f = x * cos(x) + sqrt(z*z - x*x) * sin(x)
        else
                f = x * sin(x) - sqrt(z*z - x*x) * cos(x)
        endif
end function f

program bisect
        implicit none
        character(100) :: x0char, x1char, flchar, vchar, a0char
        integer :: iter, nmax
        real :: x0, x1, x, tol, a, v, z
        real :: f0, f1, fx
        real :: f
        integer :: flag

        ! scale of our problem defined as √((2m/hbar²) * (1eV) * (1A⁰)²)
        real :: z0
        z0 = 0.5119740204943241

        print *, "This Program will eigen value for particle in a box"

        if (command_argument_count() .ne. 5)  then
                ! User Input needed for problem
                print *, "Please enter half width of well in A⁰"
                read (*,*) a
                print *, "Please enter potential of well in eV"
                read (*,*) v

                ! Lower and upper bounds.
                ! Flag is used for even - odd functions.
                print *, "Please enter lower bound, upper bound (in eV) and tolerence."
                read (*,*) x0, x1, tol
                print *, "Enter maximum no of iteration you want to try."
                read (*,*) nmax
                print *, "Enter 'odd no' for odd function and 'even no' for even."
                read (*,*) flag
        else
                ! fixing few values.
                tol = 1e-5; nmax = 1000000000
                ! getting command line agruments
                call get_command_argument(1, vchar)
                call get_command_argument(2, a0char)
                call get_command_argument(3, x0char)
                call get_command_argument(4, x1char)
                call get_command_argument(5, flchar)
                read (vchar,*) v
                read (a0char,*) a
                read (x0char,*) x0
                read (x1char,*) x1
                read (flchar,*) flag
        endif

        z = z0 * a * sqrt(v)
        flag = mod(flag,2)

        ! converted bounds given in energy to dimensionless numbers
        ! with keeping our problem in mind
        x0 = z0 * a * sqrt(x0)
        x1 = z0 * a * sqrt(x1)
        tol = abs(tol)

        f0 = f(x0, z, flag); f1 = f(x1, z, flag)

        if (f0 * f1 < 0)        then
                do iter=0,nmax,1
                        ! Following block is whole logic
                        if (iter == nmax)       then
                                print *, "No of iteration has exceded &
                                        than your defined maximum."
                                exit
                        endif
                        x = (x0 + x1) / 2
                        fx = f(x, z, flag)
                        if (abs(fx) <= tol)    then
                                print *, "After iteration", iter
                                print *, "Eigen Value of equation is",&
                                        (x / (z0 * a))**2
                                exit
                        else if (f0 * fx < 0)   then
                                x1 = x
                                f1 = fx
                        else
                                x0 = x
                                f0 = fx
                        endif
                enddo
        else if (abs(f0) == 0) then
                print *, "Eigen Value of equation is",&
                        (x0 / (z0 * a))**2
        else if (abs(f1) == 0) then
                print *, "Eigen Value of equation is",&
                        v - (x1 / (z0 * a))**2
        else
                print *, "No bound state in this range. &
                        Please use different range or use different function (odd/even)."
        endif

end program bisect
