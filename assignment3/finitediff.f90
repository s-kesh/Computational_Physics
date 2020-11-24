!**********************************
! Implimented the methods of finite difference
! to solve diffusion equation of form
!       dU                      d²U
!       --      =       D       ---
!       dt                      d²x
!**********************************
module finitediff
        implicit none
        public alpha, ex, im, cn

        contains
                real function alpha(sp, tp, li, lf, tf, D)
                      implicit none
                      integer, intent(in) :: sp, tp
                      real , intent(in) :: li, lf, tf, D
                      alpha = ( (tf / tp) / ((lf - li) / (sp - 1))**2 ) * D
                end function

                subroutine ex(u, sp, tp, alp)
                        ! Operations are row wise
                        ! first row is initial value
                        implicit none
                        real, intent(in) :: alp
                        integer, intent(in) :: sp, tp
                        real, intent(inout) :: u(:,:)
                        integer :: sn, tn
                        real :: ib, lb
                        ib = u(1,1)
                        lb = u(1,sp)

                        do tn = 1, tp - 1
                                do sn = 2, sp - 1
                                        u(tn + 1 , sn) = (1 - 2 * alp) * u(tn, sn) &
                                                + alp * (u(tn, sn+1) + u(tn, sn -1))
                                enddo
                                u(tn, 1) = ib
                                u(tn, sp) = lb
                        enddo
                end subroutine ex

                subroutine im(u, sp, tp, alp)
                        ! Operations are row wise
                        ! first row is initial value
                        implicit none
                        real, intent(in) :: alp
                        integer, intent(in) :: sp, tp
                        real, intent(inout) :: u(:,:)
                        ! Define A matrix.
                        ! AU is upper diagonal
                        ! AD is diagonal
                        ! AL is lower Diagonal
                        ! Eq Ax_new = x_old
                        real, dimension(sp) :: uold, AD, AD_orig
                        real, dimension(sp - 1) :: AU, AL, AU_orig, AL_orig
                        real :: ib, lb
                        integer :: tn, i, stat

                        do i = 1, sp - 1
                                AU_orig(i) = -1 * alp
                                AL_orig(i) = -1 * alp
                        end do
                        do i = 1, sp
                                AD_orig(i) = 1 + 2 * alp
                                uold(i) = u(1, i)
                        enddo

                        ib = u(1, 1)
                        lb = u(1, sp)

                        do tn = 2, tp
                                AL = AL_orig
                                AD = AD_orig
                                AU = AU_orig

                                ! Using subroutine from lapack,
                                ! It will save solution in uold and will modify AL, AD and AU
                                call SGTSV(sp, 1, AL, AD, AU, uold, sp, stat)

                                uold(1) = ib
                                uold(sp) = lb
                                u(tn, :) = uold
                        enddo
                end subroutine im

                subroutine cn(u, sp, tp, alp)
                        ! Operations are row wise
                        ! first row is initial value
                        implicit none
                        real, intent(in) :: alp
                        integer, intent(in) :: sp, tp
                        real, intent(inout) :: u(:,:)
                        ! Two Matrixes A and B.
                        ! Eq
                        ! Ax_new = Bx_old, we need to find x_new
                        real, dimension(sp) :: AD, AD_orig, Bxold
                        real, dimension(sp - 1) :: AU, AL, AU_orig, AL_orig
                        real :: ib, lb, alpha
                        integer :: tn, i, stat

                        alpha = alp / 2

                        do i = 1, sp - 1
                                AU_orig(i) = -1 * alp
                                AL_orig(i) = AU_orig(i)
                        enddo
                        do i = 1, sp
                                AD_orig(i) = 1 + 2 * alp
                        enddo
                        ib = u(1, 1)
                        lb = u(1, sp)

                        do tn = 2, tp

                                Bxold(1) = 2 * (1 - alpha) * u(tn - 1, 1) + alpha * u(tn - 1, 2)
                                do i = 2, sp - 1
                                        Bxold(i) = alpha * u(tn - 1, i - 1) &
                                                + (1 - 2 * alpha) * u(tn - 1, i) &
                                                + alpha * u(tn - 1, i + 1)
                                enddo
                                Bxold(sp) = alpha * u(tn - 1, sp - 1 )&
                                       + 2 * (1 - alpha) * u(tn - 1, sp)

                                AL = AL_orig
                                AD = AD_orig
                                AU = AU_orig

                                ! Using subroutine from lapack,
                                ! It will save solution in Bxold and will modify AL, AD and AU
                                call SGTSV(sp , 1, AL, AD, AU, Bxold, sp, stat)

                                Bxold(1) = ib
                                Bxold(sp) = lb
                                u(tn, :) = Bxold
                        enddo
                end subroutine cn
end module finitediff
