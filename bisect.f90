! This Program will find root of equation using bisection method
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204
! You can use and modify this program as much you want.

real function f(x)
        implicit none
        real, intent(in) :: x
        f = exp(-1 * x) + cos(x)
end function f

program bisect
        implicit none
        integer :: iter, nmax
        real :: x1, x2, x, tol
        real :: f1, f2, fx

        real :: f, average

        print *, 'This Program will find root of equation using bisection method.'
        print *, 'Please enter lower bound, upper bound and tolerence.'
        read (*,*) x1, x2, tol
        print *, "Enter maximum no of iteration you want to try."
        read (*,*) nmax
        tol = abs(tol)
        
        f1 = f(x1); f2 = f(x2)

        if (f1 * f2 < 0)        then
                do iter=0,nmax,1
                        if (iter > nmax)       then
                                print *, "No if iteration has exceded than your defined maximum."
                        endif
                        x = (x1 + x2) / 2
                        fx = f(x)
                        if (abs(fx) <= tol)    then
                                print *, "After iteration", iter
                                print *, "Root of equation is :", x
                                exit
                        else if (f1 * fx < 0)   then
                                x2 = x
                        else
                                x1 = x
                                f1 = fx
                        endif
                enddo
        else if (abs(f1) < tol) then
                print *, "Root of equation is :", x1
        else if (abs(f2) < tol) then
                print *, "Root of equation is :", x2
        else
                print *, "Root does not lie in this interval"
        endif

end program bisect
