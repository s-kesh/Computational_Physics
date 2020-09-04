! This Program will find root of equation using hybrid method
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204
! You can use and modify this program as much you want.

real function f(x)
        implicit none
        real, intent(in) :: x
        f = exp(-1 * x) + cos(x)
!        f = 10*x*exp(-1 * x * x) - 1
end function f
program hybrid
        implicit none
        integer :: iter, nmax
        real :: x0, x1, tol, dx
        real :: z, dz, x, fx, f0

        real :: f

        print *, "This Program will find root of equation using hybrid method"
        print *, "Please enter lower bound, upper bound and tolerence."
        read (*,*) x0, x1, tol
        print *, "Enter maximum no of iteration you want to try."
        read (*,*) nmax
        tol = abs(tol)
        
        f0 = f(x0)
        if (f0 == 0)    then
                print *, "Root of equation is:", x0
        else if (f(x1) == 0)       then
                print *, "Root of equation is:", x1
        else
                do iter=0,nmax,1
                        if (iter == nmax)        then
                                print *, "Exceded maximum permissible iterations"
                        endif
                        z = (x0 + x1) /2
                        dz = (f(x1) * (x1 - x0)) / (f(x1) - f(x0))
                        if (abs(dz) < 1e-14)   then
                                print *, "Hit local extreama."
                                exit
                        endif
                        x = z - dz
                        if (x > x0 .and. x < x1) then
                                fx = f(x)
                                continue
                        else
                                x = z
                                fx = f(x)
                        endif
                        if (abs(fx) < tol)    then
                                print *, "After iteration ", iter
                                print *, "Root of equation is:", x
                                exit
                        endif
                        if (fx * f0 < 0)   then
                                x1 = x
                        else
                                x0 = x; f0 = f(x0)
                        endif
                enddo
        endif
end program hybrid
