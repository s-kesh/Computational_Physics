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

real function df(x)
        implicit none
        real, intent(in) :: x
        df = -1 * exp(-1 * x) - sin(x)
!        df = 10 * (x * (-2) * x * exp(-1 * x * x) + exp(-1 * x * x))
end function df

program hybrid
        implicit none
        integer :: iter, nmax
        real :: x1, x2, tol
        real :: f1, f2, x, fx, z, fz, dfz

        real :: f, df

        print *, "This Program will find root of equation using hybrid method"
        print *, "Please enter lower bound, upper bound and tolerence."
        read (*,*) x1, x2, tol
        print *, "Enter maximum no of iteration you want to try."
        read (*,*) nmax
        tol = abs(tol)
        
        f1 = f(x1); f2 = f(x2)

        if (f1 == 0)    then
                print *, "Root of equation is:", x1
        else if (f2 == 0)       then
                print *, "Root of equation is:", x2
        else if (f1 * f2 < 0)   then
                do iter=0,nmax,1
                        if (iter == nmax)        then
                                print *, "Exceded maximum permissible iterations"
                        endif
                        z = (x1 + x2) /2
                        fz = f(z)
                        dfz = df(z)
                        if (abs(dfz) < 1e-14)   then
                                print *, "Hit local extreama."
                                exit
                        endif
                        x = z - fz / dfz
                        fx = f(x)
                        if (x > x1 .and. x < x2) then
                                continue
                        else
                                x = z
                                fx = fz
                        endif
                        if (abs(fx) < tol)    then
                                print *, "After iteration ", iter
                                print *, "Root of equation is:", x
                                exit
                        endif
                        if (fx * f1 < 0)   then
                                x2 = x; f2 = fx
                        else
                                x1 = x; f1 = fx
                        endif
                enddo
        else
                print *, "Choose another segment."
        endif
end program hybrid
