! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204
! You can use and modify this program as much you want.

real function f(x, zx)
        implicit none
        real, intent(in) :: x, zx
        f = x * sin(x) - sqrt(zx*zx - x*x) * cos(x)
end function f

program bisect
        implicit none
        integer :: iter, nmax
        real :: x0, x1, x, tol, a, v, zx
        real :: f0, f1, fx
        real :: f

        real :: z0  ! scale

        a = 10  ! half width 10 A⁰
        v = 10  ! potential 10eV

        z0 = 0.5119740204943241
        zx = a * sqrt(v) * z0

!        print *, "This Program will eigen value for particle in a box"
!        print *, "Please enter width of well in A⁰"
!        read (*,*) a
!        print *, "Please enter potential of well in eV"
!        read (*,*) v

!        print *, "Potential in units of (21.7987236eV)", v
!        print *, "half width of potential well in units of (0.529177A⁰)", a

        print *, "Please enter lower bound, upper bound (in eV) and tolerence."
        read (*,*) x0, x1
        tol = 1e-5
!        print *, "Enter maximum no of iteration you want to try."
!        read (*,*) nmax
        nmax = 100000
        x0 = z0 * sqrt(v - x0) * a
        x1 = z0 * sqrt(v - x1) * a
!        tol = abs(tol)

        print *, x0, x1

        f0 = f(x0, zx); f1 = f(x1, zx)
        print *, f0, f1

        if (f0 * f1 < 0)        then
                print *, "x0    x1      x       f1      f2      fx"
                do iter=0,nmax,1
                        if (iter == nmax)       then
                                print *, "No if iteration has exceded than your defined maximum."
                                exit
                        endif
                        x = (x0 + x1) / 2
                        fx = f(x, zx)
                        print *, x0, x1, x, f0, f1, fx
                        if (abs(fx) <= tol)    then
                                print *, "After iteration", iter, x
                                print *, "Eigen Value of equation is", v - (x / (z0 * a))**2
                                exit
                        else if (f0 * fx < 0)   then
                                x1 = x
                                f1 = fx
                        else
                                x0 = x
                                f0 = fx
                        endif
                enddo
        else if (abs(f0) == tol) then
                print *, "Eigen Value of equation is", v - (x0 / (z0 * a))**2
        else if (abs(f1) == tol) then
                print *, "Eigen Value of equation is", v - (x1 / (z0 * a))**2
        else
                print *, "Enter different segment"
        endif

end program bisect
