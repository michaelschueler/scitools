module scitools_gauss_lobatto
!! Contains the Gauss-Lobatto quadrature rules for the interval [-1,1].
!======================================================================================
   use scitools_debug
   use scitools_def,only: dp
   use scitools_utils,only: stop_error
   use scitools_legendre,only: legendre
   implicit none
   include "lobatto_inc.f90"
!--------------------------------------------------------------------------------------
   private
   public :: get_gauss_lobatto
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine get_gauss_lobatto(n, x, w, L_in, S_ni)
      integer, intent(in) :: n !! The number of points
      real(dp), allocatable, intent(out) :: x(:) !! The points
      real(dp), allocatable, intent(out) :: w(:) !! The weights
      real(dp), allocatable, intent(out), optional :: L_in(:,:) !! Vandermonde matrix
      real(dp), allocatable, intent(out), optional :: S_ni(:,:) !! inverse Vandermonde matrix
      integer :: i, j
      real(dp) :: Wn(n), L_in_tmp(n,n)

      if(n < 2 .or. n > 32) then
         call stop_error("get_gauss_lobatto: n must be in the range 2-32")
      end if 

      allocate(x(n), w(n))

      select case(n)
      case(2)
         x = x2(:)
         w = w2(:)
      case(3)
         x = x3(:)
         w = w3(:)
      case(4) 
         x = x4(:)
         w = w4(:)
      case(5)
         x = x5(:)
         w = w5(:)
      case(6)
         x = x6(:)
         w = w6(:)
      case(7)
         x = x7(:)
         w = w7(:)
      case(8)
         x = x8(:)
         w = w8(:)
      case(9)
         x = x9(:)
         w = w9(:)
      case(10)
         x = x10(:)
         w = w10(:)
      case(11)
         x = x11(:)
         w = w11(:)
      case(12)
         x = x12(:)
         w = w12(:)
      case(13)
         x = x13(:)
         w = w13(:)
      case(14)
         x = x14(:)
         w = w14(:)
      case(15)
         x = x15(:)
         w = w15(:)
      case(16)
         x = x16(:)
         w = w16(:)
      case(17)
         x = x17(:)
         w = w17(:)
      case(18)
         x = x18(:)
         w = w18(:)
      case(19)
         x = x19(:)
         w = w19(:)
      case(20)
         x = x20(:)
         w = w20(:)
      case(21)
         x = x21(:)
         w = w21(:)
      case(22)
         x = x22(:)
         w = w22(:)
      case(23)
         x = x23(:)
         w = w23(:)
      case(24)
         x = x24(:)
         w = w24(:)
      case(25)
         x = x25(:)
         w = w25(:)
      case(26)
         x = x26(:)
         w = w26(:)
      case(27)
         x = x27(:)
         w = w27(:)
      case(28)
         x = x28(:)
         w = w28(:)
      case(29)
         x = x29(:)
         w = w29(:)
      case(30)
         x = x30(:)
         w = w30(:)
      case(31)
         x = x31(:)
         w = w31(:)
      case(32)
         x = x32(:)
         w = w32(:)
      end select

      if(present(L_in) .or. present(S_ni)) then
         do i = 1, n
            L_in_tmp(:,i) = legendre(i-1, x)
         end do
      end if

      if(present(L_in)) then
         allocate(L_in(n,n))
         L_in = L_in_tmp
      end if

      if(present(S_ni)) then
         allocate(S_ni(n,n))

         do i = 1, n
            Wn(i) = 2.0_dp / (2.0_dp * i - 1.0_dp)
         end do
         Wn(n) = 2.0_dp / (n - 1.0_dp)

         do i = 1, n
            do j = 1, n
               S_ni(j, i) = L_in_tmp(i, j) * w(i) / Wn(j)
            end do
         end do

      end if

   end subroutine get_gauss_lobatto
!--------------------------------------------------------------------------------------

!======================================================================================
end module scitools_gauss_lobatto