! This Program will find root of equation using secant method
! Langauage : Fortran
! Author : Keshav Sishodia
! Rollno : PH18D204
! You can use and modify this program as much you want.

real function f(x)
        implicit none
        real, intent(in) :: x
        f = exp(-1 * x) + cos(x)
end function f

program secant
        implicit none
        integer :: iter, nmax
        real :: x1, x2, tol, dx
        real :: f

        print *, 'This Program will find root of equation using secant method.'
        print *, 'Please enter lower bound, upper bound and tolerence.'
        read (*,*) x1, x2, tol
        print *, "Enter maximum no of iteration you want to try."
        read (*,*) nmax
        tol = abs(tol)

        if (x2 == x1)     then
                print *, "Keep x1 and x2 different"
        else if (f(x1) == 0) then
                print *, "Root of equation is:" , x1
        else if (f(x2) == 0) then
                print *, "Root of equation is:" , x2
        else
                do iter=0,nmax,1
                        if (iter > nmax)        then
                               print *, "Iteration exceded maximum limit"
                        endif
                        dx = (f(x2) * (x2 - x1)) / (f(x2) - f(x1))
                        if (abs(dx) < tol)      then
                                print *, "After iteration", iter
                                print *, "Root of equation is:", x2
                                exit
                        else
                                x1 = x2
                                x2 = x2 - dx
                        endif
                enddo

        endif
end program secant
