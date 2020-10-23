!**********************************
! Implimented the methods of finite difference
! to solve diffusion equation of form
!       dU                      d²U
!       --      =       D       ---
!       dt                      d²x
!**********************************
module finitediff
        implicit none
        public alpha, ex, im

        contains
                real function alpha(sp, tp, li, lf, tf, D)
                      implicit none
                      integer, intent(in) :: sp, tp
                      real , intent(in) :: li, lf, tf, D
                      alpha = ( (tf / tp) / ((lf - li) / sp)**2 ) * D
                end function

                subroutine ex(u, sp, tp, alp)
                        ! Operations are row wise
                        ! first row is initial value
                        implicit none
                        real, intent(in) :: alp
                        integer, intent(in) :: sp, tp
                        real, intent(inout) :: u(:,:)
                        integer :: sn, tn

                        do tn = 1, tp - 1
                                do sn = 2, sp - 1
                                        u(tn + 1 , sn) = (1 - 2 * alp) * u(tn, sn) &
                                                + alp * (u(tn, sn+1) + u(tn, sn -1))
                                enddo
                        enddo
                end subroutine ex

                subroutine im(u, sp, tp, alp)
                        ! Operations are row wise
                        ! first row is initial value
                        implicit none
                        real, intent(in) :: alp
                        integer, intent(in) :: sp, tp
                        real, intent(inout) :: u(:,:)
                        real, allocatable :: AU(:), AD(:), AL(:), uold(:)
                        real :: ib, lb
                        integer :: sn, tn, i, stat


                        allocate(AU(sp-1))
                        allocate(AD(sp))
                        allocate(AL(sp-1))
                        allocate(uold(sp))

                        do i = 1, sp -1
                                AU(i) = -1 * alp
                                AL(i) = AU(i)
                        end do
                        do i = 1, sp
                                AD(i) = 1 + 2 * alp
                        enddo

                        uold = u(1, :)
                        ib = uold(1)
                        lb = uold(sp)
                        do tn = 2, tp
                                call SGTSV(sp, 1, AL, AD, AU, uold, sp, stat)
                                uold(1) = ib
                                uold(sp) = lb
                                u(tn, :) = uold
                        enddo
                        deallocate(AU)
                        deallocate(AD)
                        deallocate(AL)
                        deallocate(uold)

                end subroutine im
end module finitediff
