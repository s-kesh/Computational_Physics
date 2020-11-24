real function f(x)
        implicit none
        real, intent(in) :: x
        integer, parameter :: pi =  4 * atan(1.0)
        f =  4 * x * (1 - x)
!        f = exp(-1 * x * x)
!        f = 2 - 1.5 * x + sin(pi * x)
end function

program diffusion1d
        use finitediff
        implicit none
        character, parameter :: TAB = achar(9)
        ! Total Space Points and Time Points
        integer :: pspace, ptime, i, j, ustat
        real :: li, lf, boundi, boundl, tfinal, al, D
        real, allocatable :: u(:,:)
        real :: f

!        D = 1.44
!        boundi = 2.0
!        boundl = 0.5
!        tfinal = 0.75
!        pspace = 10
!        ptime = 15
!        li = 0.0
!        lf = 1.0

        print *, "Please enter maximum no of space points."
        read (*,*) pspace
        print *, "Please enter maximum no of time steps."
        read (*,*) ptime

        pspace = pspace + 1
        allocate(u(ptime, pspace), stat = ustat)

        if (ustat .eq. 1)      then
                print *, "You have not entered grid. Program is exiting"
                stop
        endif

        print *, "Please Enter bounds at u(init,t) and u(last,t)"
        read (*,*) boundi
        read (*,*) boundl
        print *, "Please Enter initial distance"
        read (*,*) li
        print *, "Please Enter final time and distance"
        read (*,*) tfinal
        read (*,*) lf

        u(1, 1) = boundi
        do i = 2, pspace - 1
                u(1, i) = f( li + (i - 1) * (lf - li) / (pspace - 1))
        enddo
        u(1, pspace) = boundl

        al = alpha(pspace, ptime, li, lf, tfinal, D)
        
!        call ex(u, pspace, ptime, al)
        call im(u, pspace, ptime, al)
!        call cn(u, pspace, ptime, al)

        do i = 1, ptime
                write(*,'(*(g0.7), A)', advance='no') (u(i,j), TAB, j = 1, pspace - 1)
                write(*, '(*(g0.7))') u(i, pspace)
        enddo
        deallocate(u)
end program diffusion1d
