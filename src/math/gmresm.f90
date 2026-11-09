module scitools_gmresm
!! Contains real and complex versions of the restarted GMRES(m) algorithm for
!! iteratively solving linear systems \(A x = b\).
!!
!! The implementation follows the standard algorithm with right preconditioning,
!! \(A M^{-1} u = b\), \(x = M^{-1} u\): the Arnoldi basis is built by modified
!! Gram-Schmidt orthogonalization and the least-squares problem is solved
!! progressively with Givens rotations. The residual norm \(\|b - A x\|\) is
!! therefore available in every iteration without extra matrix-vector products.
!!
!! References:
!!
!! * Y. Saad and M. H. Schultz, SIAM J. Sci. Stat. Comput. 7, 856 (1986).
!! * Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM (2003),
!!   Algorithms 6.9 and 9.5.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit
   use scitools_def,only: dp,zero
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: gmresm
!--------------------------------------------------------------------------------------
   interface gmresm
      !! Restarted GMRES(m) solver for real or complex linear systems \(A x = b\).
      !! The matrix and the preconditioner are only accessed through user-supplied
      !! routines.
      module procedure dgmresm, zgmresm
   end interface gmresm
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine dgmresm(m,n,x,b,matvec,psolve,dotprd,res,its,info)
      !! GMRES(m) for real systems.
      integer,  intent(in)    :: m !! dimension of the Krylov subspace (restart length)
      integer,  intent(in)    :: n !! dimension of the linear system
      real(dp), intent(inout) :: x(:) !! on input: initial guess; on output: solution
      real(dp), intent(in)    :: b(:) !! right-hand side \(b\)
      interface
      subroutine matvec(x,y)
         !! user-supplied matrix-vector product \(y = A x\)
         import :: dp
         real(dp),intent(in)    :: x(:) !! input vector
         real(dp),intent(inout) :: y(:) !! result \(A x\)
      end subroutine matvec

      subroutine psolve(x)
         !! user-supplied preconditioner, applied in place: \(x \rightarrow M^{-1} x\)
         import :: dp
         real(dp),intent(inout) :: x(:) !! vector to be preconditioned (overwritten)
      end subroutine psolve

      function dotprd(x,y) result(d)
         !! user-supplied inner product \((x|y)\)
         import :: dp
         real(dp),intent(in) :: x(:),y(:) !! first and second vector
         real(dp) :: d
      end function dotprd
      end interface
      real(dp), intent(inout) :: res !! on input: tolerance for \(\|b - A x\|\); on output: residual \(\|b - A x\|\) of the returned solution
      integer,  intent(inout) :: its !! on input: maximum number of iterations; on output: number of iterations performed
      integer,  intent(inout) :: info !! on input: 1 = print the residual in each iteration; on output: 0 = converged, 2 = maximum number of iterations reached
      real(dp), allocatable :: V(:,:), H(:,:), w(:), z(:)
      real(dp) :: c(m), s(m), g(m+1), y(m)
      real(dp) :: tol, beta, hnorm, r, temp
      integer :: imx, i, j, jlast
      logical :: verbose

      tol = res
      imx = its
      its = 0
      verbose = (info == 1)

      allocate(V(n,m+1), H(m+1,m), w(n), z(n))

      do ! restart cycles
         ! true residual of the current iterate
         call matvec(x, w)
         w = b - w
         beta = sqrt(abs(dotprd(w, w)))
         res = beta
         if(beta <= tol) then
            info = 0; exit
         end if
         if(its >= imx) then
            info = 2; exit
         end if

         V(:,1) = w / beta
         H = 0.0_dp
         g = 0.0_dp
         g(1) = beta
         jlast = 0

         do j = 1, m
            its = its + 1
            jlast = j

            ! Arnoldi step with modified Gram-Schmidt: w = A M^{-1} v_j
            z = V(:,j)
            call psolve(z)
            call matvec(z, w)
            do i = 1, j
               H(i,j) = dotprd(V(:,i), w)
               w = w - H(i,j) * V(:,i)
            end do
            hnorm = sqrt(abs(dotprd(w, w)))
            if(hnorm > 0.0_dp) V(:,j+1) = w / hnorm

            ! apply the previous Givens rotations to the new column
            do i = 1, j-1
               temp     =  c(i) * H(i,j) + s(i) * H(i+1,j)
               H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j)
               H(i,j)   = temp
            end do

            ! new rotation eliminating H(j+1,j) = hnorm
            r = hypot(H(j,j), hnorm)
            if(r > 0.0_dp) then
               c(j) = H(j,j) / r
               s(j) = hnorm / r
            else
               c(j) = 1.0_dp
               s(j) = 0.0_dp
            end if
            H(j,j) = r
            g(j+1) = -s(j) * g(j)
            g(j)   =  c(j) * g(j)

            res = abs(g(j+1))
            if(verbose) write(output_unit,'("gmresm: it = ",i6," res = ",E12.5)') its, res

            if(res <= tol .or. its >= imx .or. hnorm == 0.0_dp) exit
         end do

         ! solve the triangular system H y = g and update x = x + M^{-1} V y
         do i = jlast, 1, -1
            y(i) = (g(i) - sum(H(i,i+1:jlast) * y(i+1:jlast))) / H(i,i)
         end do
         z = matmul(V(:,1:jlast), y(1:jlast))
         call psolve(z)
         x = x + z
      end do

      deallocate(V, H, w, z)

   end subroutine dgmresm
!--------------------------------------------------------------------------------------
   subroutine zgmresm(m,n,x,b,matvec,psolve,dotprd,res,its,info)
      !! GMRES(m) for complex systems.
      !! The inner product must be conjugate-linear in its first argument, \((x|y) = x^\dagger y\).
      integer,     intent(in)    :: m !! dimension of the Krylov subspace (restart length)
      integer,     intent(in)    :: n !! dimension of the linear system
      complex(dp), intent(inout) :: x(:) !! on input: initial guess; on output: solution
      complex(dp), intent(in)    :: b(:) !! right-hand side \(b\)
      interface
      subroutine matvec(x,y)
         !! user-supplied matrix-vector product \(y = A x\)
         import :: dp
         complex(dp),intent(in)    :: x(:) !! input vector
         complex(dp),intent(inout) :: y(:) !! result \(A x\)
      end subroutine matvec

      subroutine psolve(x)
         !! user-supplied preconditioner, applied in place: \(x \rightarrow M^{-1} x\)
         import :: dp
         complex(dp),intent(inout) :: x(:) !! vector to be preconditioned (overwritten)
      end subroutine psolve

      function dotprd(x,y) result(d)
         !! user-supplied inner product \((x|y)\)
         import :: dp
         complex(dp),intent(in) :: x(:),y(:) !! first and second vector
         complex(dp) :: d
      end function dotprd
      end interface
      real(dp),    intent(inout) :: res !! on input: tolerance for \(\|b - A x\|\); on output: residual \(\|b - A x\|\) of the returned solution
      integer,     intent(inout) :: its !! on input: maximum number of iterations; on output: number of iterations performed
      integer,     intent(inout) :: info !! on input: 1 = print the residual in each iteration; on output: 0 = converged, 2 = maximum number of iterations reached
      complex(dp), allocatable :: V(:,:), H(:,:), w(:), z(:)
      complex(dp) :: s(m), g(m+1), y(m), temp, phase
      real(dp) :: c(m)
      real(dp) :: tol, beta, hnorm, r, habs
      integer :: imx, i, j, jlast
      logical :: verbose

      tol = res
      imx = its
      its = 0
      verbose = (info == 1)

      allocate(V(n,m+1), H(m+1,m), w(n), z(n))

      do ! restart cycles
         ! true residual of the current iterate
         call matvec(x, w)
         w = b - w
         beta = sqrt(abs(dotprd(w, w)))
         res = beta
         if(beta <= tol) then
            info = 0; exit
         end if
         if(its >= imx) then
            info = 2; exit
         end if

         V(:,1) = w / beta
         H = zero
         g = zero
         g(1) = beta
         jlast = 0

         do j = 1, m
            its = its + 1
            jlast = j

            ! Arnoldi step with modified Gram-Schmidt: w = A M^{-1} v_j
            z = V(:,j)
            call psolve(z)
            call matvec(z, w)
            do i = 1, j
               H(i,j) = dotprd(V(:,i), w)
               w = w - H(i,j) * V(:,i)
            end do
            hnorm = sqrt(abs(dotprd(w, w)))
            if(hnorm > 0.0_dp) V(:,j+1) = w / hnorm

            ! apply the previous Givens rotations to the new column
            do i = 1, j-1
               temp     =  c(i) * H(i,j) + s(i) * H(i+1,j)
               H(i+1,j) = -conjg(s(i)) * H(i,j) + c(i) * H(i+1,j)
               H(i,j)   = temp
            end do

            ! new rotation eliminating H(j+1,j) = hnorm (real)
            habs = abs(H(j,j))
            r = hypot(habs, hnorm)
            if(r == 0.0_dp) then
               c(j) = 1.0_dp
               s(j) = zero
            else if(habs == 0.0_dp) then
               c(j) = 0.0_dp
               s(j) = (1.0_dp, 0.0_dp)
               H(j,j) = r
            else
               phase = H(j,j) / habs
               c(j) = habs / r
               s(j) = phase * hnorm / r
               H(j,j) = phase * r
            end if
            g(j+1) = -conjg(s(j)) * g(j)
            g(j)   =  c(j) * g(j)

            res = abs(g(j+1))
            if(verbose) write(output_unit,'("gmresm: it = ",i6," res = ",E12.5)') its, res

            if(res <= tol .or. its >= imx .or. hnorm == 0.0_dp) exit
         end do

         ! solve the triangular system H y = g and update x = x + M^{-1} V y
         do i = jlast, 1, -1
            y(i) = (g(i) - sum(H(i,i+1:jlast) * y(i+1:jlast))) / H(i,i)
         end do
         z = matmul(V(:,1:jlast), y(1:jlast))
         call psolve(z)
         x = x + z
      end do

      deallocate(V, H, w, z)

   end subroutine zgmresm
!--------------------------------------------------------------------------------------

!======================================================================================
end module scitools_gmresm
