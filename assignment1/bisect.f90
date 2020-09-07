! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204
! You can use and modify this program as much you want.

real function f(x, n, m)
        implicit none
        real, intent(in) :: x, n, m
        f = sqrt(x) * sin(n * sqrt(x)) - sqrt(m - x) * cos(n * sqrt(x))
end function f

program bisect
        implicit none
        integer :: iter, nmax
        real :: x0, x1, x, tol, a, v
        real :: f0, f1, fx
        real :: f

        real :: a0, e0  !! atomic units for length and energy
        a0 = 0.529177210903        ! Taking value of width in atomic units. (a0) = bohr radius)
        e0 = 21.798723584403383    ! Taking value of potential in atomic units. (E0 = (hbar^2 / (2*m*a0^2)

        print *, "This Program will eigen value for particle in a box"
        print *, "Please enter width of well in A⁰"
        read (*,*) a
        print *, "Please enter potential of well in eV"
        a = a / a0      ! Taking value of width in atomic units.
        v = v / e0      ! Taking value of potential in atomic units

        print *, "Please enter lower bound, upper bound (in eV) and tolerence."
        read (*,*) x0, x1, tol
        print *, "Enter maximum no of iteration you want to try."
        read (*,*) nmax
        x0 = x0 / e0
        x1 = x1 / e0
        tol = abs(tol)

        f0 = f(x0, a, v); f1 = f(x1, a, v)

        if (f0 * f1 < 0)        then
                do iter=0,nmax,1
                        if (iter == nmax)       then
                                print *, "No if iteration has exceded than your defined maximum."
                                exit
                        endif
                        x = (x0 + x1) / 2
                        fx = f(x, a, v)
                        if (abs(fx) <= tol)    then
                                print *, "After iteration", iter
                                print *, "Eigen Value of equation is", e0 * x
                                exit
                        else if (f0 * fx < 0)   then
                                x1 = x
                        else
                                x0 = x
                                f0 = fx
                        endif
                enddo
        else if (abs(f0) < tol) then
                print *, "Eigen Value of equation is", e0 * x0
        else if (abs(f1) < tol) then
                print *, "Eigen Value of equation is", e0 * x1
        else
                print *, "Enter different segment"
        endif

end program bisect
