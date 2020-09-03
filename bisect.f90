program bisect
        implicit none

        integer :: inter
        real :: low
        real :: up
        real :: x
        real :: tol

        print *, 'This Program will find root of equation using bisection method.'
        print *, 'Please enter lower bound, upper bound and tolerence.'
        read (*,*) low
        read (*,*) up
        read (*,*) tol

        print *, low, up, tol

end program bisect
