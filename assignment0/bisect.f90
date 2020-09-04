! This Program will find root of equation using bisection method
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

program bisect
        implicit none
        integer :: iter, nmax
        real :: x0, x1, x, tol
        real :: f0, f1, fx

        real :: f

        print *, 'This Program will find root of equation using bisection method.'
        print *, 'Please enter lower bound, upper bound and tolerence.'
        read (*,*) x0, x1, tol
        print *, "Enter maximum no of iteration you want to try."
        read (*,*) nmax
        tol = abs(tol)
        
        f0 = f(x0); f1 = f(x1)

        if (f0 * f1 < 0)        then
                do iter=0,nmax,1
                        if (iter == nmax)       then
                                print *, "No if iteration has exceded than your defined maximum."
                        endif
                        x = (x0 + x1) / 2
                        fx = f(x)
                        if (abs(fx) <= tol)    then
                                print *, "After iteration", iter
                                print *, "Root of equation is :", x
                                exit
                        else if (f0 * fx < 0)   then
                                x1 = x
                        else
                                x0 = x
                                f0 = fx
                        endif
                enddo
        else if (abs(f0) < tol) then
                print *, "Root of equation is :", x0
        else if (abs(f1) < tol) then
                print *, "Root of equation is :", x1
        else
                print *, "Root does not lie in this interval"
        endif

end program bisect
