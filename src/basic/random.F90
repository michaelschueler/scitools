module scitools_random
!======================================================================================
!! Routines to generate random numbers
   use scitools_debug
   use scitools_def,only: dp, iu
   implicit none
!-------------------------------------------------------------------------------------- 
   private
   public :: InitRandom, RandomUnitVector
!-------------------------------------------------------------------------------------- 
   interface RandomUnitVector
      module procedure DRandomUnitVector, ZRandomUnitVector
   end interface RandomUnitVector
!-------------------------------------------------------------------------------------- 
contains
!-------------------------------------------------------------------------------------- 
   subroutine InitRandom()
      integer, allocatable :: seed(:)
      integer :: nsize,count,i

      call random_seed(size=nsize)

      call system_clock(count)
      allocate(seed(nsize), source=count+37*[(i,i=0,nsize-1)])

      call random_seed(put=seed)

      deallocate(seed)

   end subroutine InitRandom
!-------------------------------------------------------------------------------------- 
   subroutine DRandomUnitVector(n,v)
      integer,intent(in) :: n
      real(dp),allocatable,intent(inout) :: v(:)
      real(dp) :: cn

      if(.not.allocated(v)) allocate(v(n))

      call assert(size(v) >= n)

      call random_number(v(1:n))
      v(1:n) = -1.0_dp + 2.0_dp * v(1:n)
      cn = norm2(v(1:n))

      v(1:n) = v(1:n) / cn

   end subroutine DRandomUnitVector
!-------------------------------------------------------------------------------------- 
   subroutine ZRandomUnitVector(n,v)
      integer,intent(in) :: n
      complex(dp),allocatable,intent(inout) :: v(:)
      real(dp) :: cn
      real(dp),allocatable :: x(:),y(:)

      if(.not.allocated(v)) allocate(v(n))

      call assert(size(v) >= n)

      allocate(x(n),y(n))
      call random_number(x(1:n))
      call random_number(y(1:n))
      x(1:n) = -1.0_dp + 2.0_dp * x(1:n)
      y(1:n) = -1.0_dp + 2.0_dp * y(1:n)
      v(1:n) = x(1:n) + iu * y(1:n)

      deallocate(x,y)

      cn = sqrt(abs(dot_product(v,v)))

      v(1:n) = v(1:n) / cn

   end subroutine ZRandomUnitVector
!-------------------------------------------------------------------------------------- 

!======================================================================================
end module scitools_random