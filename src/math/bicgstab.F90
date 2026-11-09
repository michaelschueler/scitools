module scitools_bicgstab
!! Contains stabilized biconjugate-gradient linear solvers
!!
!! Reference: H. A. van der Vorst, SIAM J. Sci. Stat. Comput. 13, 631 (1992).
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use scitools_def,only: dp,zero
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: bicgstab
!--------------------------------------------------------------------------------------
   interface bicgstab
      !! Solves \(A x = b\) for real or complex \(A\) with the preconditioned BiCGSTAB method.
      !! The matrix is only accessed through the user-supplied product `matvec`.
      !! The solver works on the left-preconditioned system \(M^{-1} A x = M^{-1} b\);
      !! the residual is measured as \(\|M^{-1}(b - A x)\|\) with the user-supplied
      !! inner product.
      module procedure dbicgstab, zbicgstab
   end interface bicgstab
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine dbicgstab(ndim, x, b, matvec, psolve, dotp, info, maxiter, tol, res)
      !! BiCGSTAB solver for real linear systems \(A x = b\).
      integer,intent(in) :: ndim !! dimension of the linear system
      real(dp),dimension(ndim),intent(inout) :: x !! on input: initial guess; on output: solution
      real(dp),dimension(ndim),intent(in)    :: b !! right-hand side \(b\)
      interface
         subroutine matvec(x,y)
            !! user-supplied matrix-vector product \(y = A x\)
            import :: dp
            real(dp),intent(in) :: x(:) !! input vector
            real(dp),intent(inout) :: y(:) !! result \(A x\)
         end subroutine matvec

         subroutine psolve(x)
            !! user-supplied preconditioner, applied in place: \(x \rightarrow M^{-1} x\)
            import :: dp
            real(dp),intent(inout) :: x(:) !! vector to be preconditioned (overwritten)
         end subroutine psolve

         real(dp) function dotp(x,y)
            !! user-supplied inner product \((x|y)\)
            import :: dp
            real(dp),intent(in) :: x(:),y(:) !! first and second vector
         end function dotp
      end interface
      integer,intent(inout) :: info !! on output: 0 = converged, 1 = not converged within maxiter, 2 = breakdown
      integer,intent(in),optional :: maxiter !! maximum number of iterations (default: 1000)
      real(dp),intent(in),optional :: tol !! convergence threshold for the residual norm (default: 1e-6)
      real(dp),intent(out),optional :: res !! final residual norm \(\|M^{-1}(b - A x)\|\)
      integer :: iter_max, its
      real(dp) :: tol_, res_
      real(dp) :: rho, rho_new, alpha, beta, omega, r0v, t2
      real(dp),allocatable,dimension(:) :: pb, r, r0, p, v, s, t

      iter_max = 1000
      if(present(maxiter)) iter_max = maxiter
      tol_ = 1.0e-6_dp
      if(present(tol)) tol_ = tol

      allocate(pb(ndim), r(ndim), r0(ndim), p(ndim), v(ndim), s(ndim), t(ndim))

      ! left preconditioning: solve M^{-1} A x = M^{-1} b
      pb = b
      call psolve(pb)

      ! initial residual
      call matvec(x, t)
      call psolve(t)
      r = pb - t
      res_ = sqrt(abs(dotp(r, r)))

      info = 0
      if(res_ > tol_) then
         r0 = r
         rho = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp
         p = 0.0_dp; v = 0.0_dp
         info = 1

         do its = 1, iter_max
            rho_new = dotp(r0, r)
            if(rho_new == 0.0_dp) then
               info = 2; exit
            end if
            beta = (rho_new / rho) * (alpha / omega)
            p = r + beta * (p - omega * v)

            call matvec(p, v)
            call psolve(v)
            r0v = dotp(r0, v)
            if(r0v == 0.0_dp) then
               info = 2; exit
            end if
            alpha = rho_new / r0v

            s = r - alpha * v
            res_ = sqrt(abs(dotp(s, s)))
            if(res_ <= tol_) then
               x = x + alpha * p
               info = 0; exit
            end if

            call matvec(s, t)
            call psolve(t)
            t2 = dotp(t, t)
            if(t2 == 0.0_dp) then
               x = x + alpha * p
               info = 2; exit
            end if
            omega = dotp(t, s) / t2

            x = x + alpha * p + omega * s
            r = s - omega * t
            res_ = sqrt(abs(dotp(r, r)))
            if(res_ <= tol_) then
               info = 0; exit
            end if
            if(omega == 0.0_dp) then
               info = 2; exit
            end if
            rho = rho_new
         end do

         ! report the true residual of the returned solution
         call matvec(x, t)
         call psolve(t)
         r = pb - t
         res_ = sqrt(abs(dotp(r, r)))
         if(info == 0 .and. res_ > tol_) info = 1
      end if

      if(present(res)) res = res_

      deallocate(pb, r, r0, p, v, s, t)

   end subroutine dbicgstab
!--------------------------------------------------------------------------------------
   subroutine zbicgstab(ndim, x, b, matvec, psolve, dotp, iter, info, tol, res)
      !! BiCGSTAB solver for complex linear systems \(A x = b\).
      !! The inner product must be conjugate-linear in its first argument, \((x|y) = x^\dagger y\).
      integer,intent(in) :: ndim !! dimension of the linear system
      complex(dp),dimension(:),intent(inout) :: x !! on input: initial guess; on output: solution
      complex(dp),dimension(:),intent(in)    :: b !! right-hand side \(b\)
      interface
         subroutine matvec(x,y)
            !! user-supplied matrix-vector product \(y = A x\)
            import :: dp
            complex(dp),intent(in) :: x(:) !! input vector
            complex(dp),intent(inout) :: y(:) !! result \(A x\)
         end subroutine matvec

         subroutine psolve(x)
            !! user-supplied preconditioner, applied in place: \(x \rightarrow M^{-1} x\)
            import :: dp
            complex(dp),intent(inout) :: x(:) !! vector to be preconditioned (overwritten)
         end subroutine psolve

         complex(dp) function dotp(x,y)
            !! user-supplied inner product \((x|y)\)
            import :: dp
            complex(dp),intent(in) :: x(:),y(:) !! first and second vector
         end function dotp
      end interface
      integer,intent(inout) :: iter !! on input: maximum number of iterations; on output: number of iterations performed
      integer,intent(inout) :: info !! on input: 1 = print residual in each iteration; on output: 0 = converged, 1 = not converged within iter, 2 = breakdown
      real(dp),intent(in),optional :: tol !! convergence threshold for the residual norm (default: 1e-6)
      real(dp),intent(out),optional :: res !! final residual norm \(\|M^{-1}(b - A x)\|\)
      integer :: iter_max, its
      logical :: verbose
      real(dp) :: tol_, res_, t2
      complex(dp) :: rho, rho_new, alpha, beta, omega, r0v
      complex(dp),allocatable,dimension(:) :: pb, r, r0, p, v, s, t

      iter_max = iter
      verbose = (info == 1)
      tol_ = 1.0e-6_dp
      if(present(tol)) tol_ = tol

      allocate(pb(ndim), r(ndim), r0(ndim), p(ndim), v(ndim), s(ndim), t(ndim))

      ! left preconditioning: solve M^{-1} A x = M^{-1} b
      pb = b
      call psolve(pb)

      ! initial residual
      call matvec(x, t)
      call psolve(t)
      r = pb - t
      res_ = sqrt(abs(dotp(r, r)))

      info = 0
      iter = 0
      if(res_ > tol_) then
         r0 = r
         rho = (1.0_dp, 0.0_dp); alpha = (1.0_dp, 0.0_dp); omega = (1.0_dp, 0.0_dp)
         p = zero; v = zero
         info = 1

         do its = 1, iter_max
            iter = its
            rho_new = dotp(r0, r)
            if(rho_new == zero) then
               info = 2; exit
            end if
            beta = (rho_new / rho) * (alpha / omega)
            p = r + beta * (p - omega * v)

            call matvec(p, v)
            call psolve(v)
            r0v = dotp(r0, v)
            if(r0v == zero) then
               info = 2; exit
            end if
            alpha = rho_new / r0v

            s = r - alpha * v
            res_ = sqrt(abs(dotp(s, s)))
            if(res_ <= tol_) then
               x = x + alpha * p
               info = 0; exit
            end if

            call matvec(s, t)
            call psolve(t)
            t2 = real(dotp(t, t), kind=dp)
            if(t2 == 0.0_dp) then
               x = x + alpha * p
               info = 2; exit
            end if
            omega = dotp(t, s) / t2

            x = x + alpha * p + omega * s
            r = s - omega * t
            res_ = sqrt(abs(dotp(r, r)))

            if(verbose) then
               write(output_unit,'("bicgstab: it = ",i6," res = ",E12.5)') its, res_
            end if

            if(res_ <= tol_) then
               info = 0; exit
            end if
            if(omega == zero) then
               info = 2; exit
            end if
            rho = rho_new
         end do

         ! report the true residual of the returned solution
         call matvec(x, t)
         call psolve(t)
         r = pb - t
         res_ = sqrt(abs(dotp(r, r)))
         if(info == 0 .and. res_ > tol_) info = 1
      end if

      if(present(res)) res = res_

      deallocate(pb, r, r0, p, v, s, t)

   end subroutine zbicgstab
!--------------------------------------------------------------------------------------

!======================================================================================
end module scitools_bicgstab
