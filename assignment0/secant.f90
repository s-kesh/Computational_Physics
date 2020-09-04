! This Program will find root of equation using secant method
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204
! You can use and modify this program as much you want.

real function f(x)
        implicit none
        real, intent(in) :: x
        f = 10*x*exp(-1 * x * x) - 1
end function f

program secant
        implicit none
        integer :: iter, nmax
        real :: x0, x1, tol, dx
        real :: f

        print *, 'This Program will find root of equation using secant method.'
        print *, 'Please enter lower bound, upper bound and tolerence.'
        read (*,*) x0, x1, tol
        print *, "Enter maximum no of iteration you want to try."
        read (*,*) nmax
        tol = abs(tol)

        if (x1 == x0)     then
                print *, "Keep x0 and x1 different"
        else if (f(x0) == 0) then
                print *, "Root of equation is:" , x0
        else if (f(x1) == 0) then
                print *, "Root of equation is:" , x1
        else
                do iter=0,nmax,1
                        if (iter == nmax)        then
                               print *, "Iteration exceded maximum limit"
                        endif
                        dx = (f(x1) * (x1 - x0)) / (f(x1) - f(x0))
                        if (abs(dx) < tol)      then
                                print *, "After iteration", iter
                                print *, "Root of equation is:", x1
                                exit
                        else
                                x0 = x1
                                x1 = x1 - dx
                        endif
                enddo

        endif
end program secant
