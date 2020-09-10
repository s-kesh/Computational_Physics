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

real function f(x, zx, flag)
        ! Function need to solve.
        ! It has two parts `cos` wala for odd and `sin` wala for even
        implicit none
        real, intent(in) :: x, zx
        integer, intent(in) :: flag
        if (flag == 1)  then
                f = x * cos(x) + sqrt(zx*zx - x*x) * sin(x)
        else
                f = x * sin(x) - sqrt(zx*zx - x*x) * cos(x)
        endif
end function f

real function df(x, zx, flag)
        implicit none
        real, intent(in) :: x, zx
        integer, intent(in) :: flag
        if (flag == 1)  then
                df = x * (-1) * sin(x) + cos(x)&
                        + sqrt(zx * zx - x*x) * cos(x)&
                        - (x / sqrt(zx*zx - x*x)) * sin(x)
        else
                df = x * cos(x) + sin(x)&
                        + sqrt(zx * zx - x*x) * sin(x)&
                        + (x / sqrt(zx*zx - x*x)) * cos(x)
        endif
end function df

program hybrid
        implicit none
        character(100) :: x0char, x1char, flchar
        integer :: iter, nmax
        real :: x0, x1, x, tol, a, v, zx
        real :: f0, f1, fx, t, ft, dft
        real :: f, df
        integer :: flag

        ! scale of our problem defined as √((2m/hbar²) * (1eV) * (1A⁰)²)
        real :: z0
        z0 = 0.5119740204943241

        print *, "This Program will eigen value for particle in a box"

        if (command_argument_count() .ne. 3)  then
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
                a = 10; v = 10
                tol = 1e-5; nmax = 1000000000
                call get_command_argument(1, x0char)
                call get_command_argument(2, x1char)
                call get_command_argument(3, flchar)
                read (x0char,*) x0
                read (x1char,*) x1
                read (flchar,*) flag
        endif

        zx = z0 * a * sqrt(v)
        flag = mod(flag,2)

        ! converted bounds given in energy to dimensionless numbers
        ! with keeping our problem in mind
        x0 = z0 * a * sqrt(v - x0)
        x1 = z0 * a * sqrt(v - x1)
        tol = abs(tol)

        f0 = f(x0, zx, flag); f1 = f(x1, zx, flag)

        if (f0 * f1 < 0)        then
                do iter=0,nmax,1
                        ! Following block is whole logic
                        if (iter == nmax)       then
                                print *, "No of iteration has exceded &
                                        than your defined maximum."
                                exit
                        endif
                        t = (x0 + x1) / 2
                        ft = f(t, zx, flag)
                        dft = df(t, zx, flag)
                        if (abs(dft) < 1e-14)   then
                                print *, "Hit local extreama"
                                exit
                        endif
                        x = t - ft / dft
                        fx = f(x, zx, flag)
                        if (x > x0 .and. x < x1)        then
                                continue
                        else
                                x = t; fx = ft
                        endif
                        if (abs(fx) <= tol)    then
                                print *, "After iteration", iter
                                print *, "Eigen Value of equation is",&
                                        v - (x / (z0 * a))**2
                                exit
                        else if (f0 * fx < 0)   then
                                x1 = x; f1 = fx
                        else
                                x0 = x; f0 = fx
                        endif
                enddo
        else if (abs(f0) == 0) then
                print *, "Eigen Value of equation is",&
                        v - (x0 / (z0 * a))**2
        else if (abs(f1) == 0) then
                print *, "Eigen Value of equation is",&
                        v - (x1 / (z0 * a))**2
        else
                print *, "No bound state in this range. &
                        Please use different range or use different function (odd/even)."
        endif

end program hybrid
