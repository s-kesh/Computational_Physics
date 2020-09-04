! This Program will find root of equation using newton rapsion method
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

real function df(x)
        implicit none
        real, intent(in) :: x
        df = -1 * exp(-1 * x) - sin(x)
!        df = 10 * (x * (-2) * x * exp(-1 * x * x) + exp(-1 * x * x))
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
                do iter=0,nmax,1
                        if (iter == nmax)        then
                                print *, "Iteration has exceded maximum allowed value"
                        endif
                        dfx = df(x)
                        if (abs(dfx) < 1e-12)  then
                                print *, "Buh!! its local extreama."
                                exit
                        endif
                        dx = f(x) / dfx
                        if (abs(dx) < tol)      then
                                print *, "After iteration", iter
                                print *, "You got your root", x
                                exit
                        else
                                x = x - dx
                        endif
                enddo
        endif

end program newton
