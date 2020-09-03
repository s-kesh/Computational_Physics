function fun(x) result(fx)
        implicit none
        real, intent(in) :: x
        real :: fx
        fx = exp(-1 * x) + cos(x)
end function fun

function average(a, b) result(x)
        implicit none
        real, intent(in) :: a
        real, intent(in) :: b
        real :: x

        x = (a + b) / 2
end function average

program bisect
        implicit none

        integer :: inter
        real :: x1, x2, x, tol
        real :: f1, f2, fx

        real :: fun, average

        print *, 'This Program will find root of equation using bisection method.'
        print *, 'Please enter lower bound, upper bound and tolerence.'
        read (*,*) x1, x2, tol
        tol = abs(tol)
        
        f1 = fun(x1); f2 = fun(x2)

        if (f1 * f2 < 0)        then
                do
                        x = average(x1, x2)
                        fx = fun(x)
                        if (abs(fx) <= tol)    then
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
