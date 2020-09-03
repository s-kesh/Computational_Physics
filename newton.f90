real function f(x)
        implicit none
        real, intent(in) :: x
        f = exp(-1 * x) + cos(x)
end function f

real function df(x)
        implicit none
        real, intent(in) :: x
        df = -1 * exp(-1 * x) - sin(x)
end function df

program newton
        implicit none
        integer :: iter, nmax
        real :: x, tol, dx, dfx

        real :: f, df

        print *, "This program with find root of given equation starting with user nput guess."
        print *, "Please enter your intial guess and tolerance"
        read (*,*) x, tol
        tol = abs(tol)
        print *, "Please enter maximum no of iteration allowed."
        read (*,*) nmax

        if (f(x) == 0)     then
                print *, "You got your root", x
        else
                iter = 0
                do
                        if (iter > nmax)        then
                                print *, "Iteration has exceded maximum allowed value"
                                exit
                        endif
                        dfx = df(x)
                        if (abs(dfx) < 1e-12)  then
                                print *, "Buh!! its local extreama."
                                exit
                        endif
                        dx = f(x) / dfx
                        if (abs(dx) < tol)      then
                                print *, "You got your root", x
                                exit
                        else
                                x = x - dx
                        endif
                        iter = iter + 1
                enddo
        endif

end program newton
