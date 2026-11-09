module scitools_MatrixExp
!! Provides tools to compute the matrix exponential \(e^{tH}\) of general
!! real or complex square matrices.
!!
!! The exponential is computed by scaling and squaring combined with diagonal
!! Padé approximants of degree 3, 5, 7, 9 or 13, chosen from the 1-norm of
!! \(tH\) such that the result is accurate to double precision:
!!
!! N. J. Higham, *The scaling and squaring method for the matrix exponential
!! revisited*, SIAM J. Matrix Anal. Appl. 26, 1179 (2005).
!!
!! All routines are thread safe.
!----------------------------------------------------------------------|
  use scitools_def,only: dp
  implicit none
!----------------------------------------------------------------------|
  private
  public :: &
       InitExpM,&
       CleanExpM,&
       ExpM,&
       ExpM_ThrSv
!----------------------------------------------------------------------|
  interface ExpM
     !! Matrix exponential \(e^{tH}\) of a general real or complex matrix \(H\).
     module procedure zexpm, dexpm
  end interface ExpM

  interface ExpM_ThrSv
     !! Matrix exponential \(e^{tH}\); 1x1 and 2x2 complex matrices are treated
     !! analytically. Kept for compatibility: [[ExpM]] is thread safe as well.
     module procedure zexpm_thrsv, dexpm
  end interface ExpM_ThrSv
!----------------------------------------------------------------------|
  ! maximal 1-norms for which the Padé approximants of degree 3, 5, 7, 9, 13
  ! give double-precision accuracy (Higham 2005, Table 2.3)
  integer,parameter :: pade_deg(5) = [3, 5, 7, 9, 13]
  real(dp),parameter :: theta(5) = [1.495585217958292e-2_dp, 2.539398330063230e-1_dp, &
       9.504178996162932e-1_dp, 2.097847961257068e0_dp, 5.371920351148152e0_dp]
  ! coefficients b_0 ... b_m of the Padé approximants
  real(dp),parameter :: b3(0:3) = [120.0_dp, 60.0_dp, 12.0_dp, 1.0_dp]
  real(dp),parameter :: b5(0:5) = [30240.0_dp, 15120.0_dp, 3360.0_dp, 420.0_dp, 30.0_dp, 1.0_dp]
  real(dp),parameter :: b7(0:7) = [17297280.0_dp, 8648640.0_dp, 1995840.0_dp, 277200.0_dp, &
       25200.0_dp, 1512.0_dp, 56.0_dp, 1.0_dp]
  real(dp),parameter :: b9(0:9) = [17643225600.0_dp, 8821612800.0_dp, 2075673600.0_dp, &
       302702400.0_dp, 30270240.0_dp, 2162160.0_dp, 110880.0_dp, 3960.0_dp, 90.0_dp, 1.0_dp]
  real(dp),parameter :: b13(0:13) = [64764752532480000.0_dp, 32382376266240000.0_dp, &
       7771770303897600.0_dp, 1187353796428800.0_dp, 129060195264000.0_dp, &
       10559470521600.0_dp, 670442572800.0_dp, 33522128640.0_dp, 1323241920.0_dp, &
       40840800.0_dp, 960960.0_dp, 16380.0_dp, 182.0_dp, 1.0_dp]
!----------------------------------------------------------------------|
contains
!----------------------------------------------------------------------|
  pure elemental complex(dp) function my_cosh(z)
  !! cosh for complex arguments as a workaround for older compilers
    complex(dp),intent(in) :: z

#ifdef COSHFIX
    my_cosh = 0.5_dp * (exp(z) + exp(-z))
#else
    my_cosh = cosh(z)
#endif

  end function my_cosh
!----------------------------------------------------------------------|
  pure elemental complex(dp) function my_sinh(z)
  !! sinh for complex arguments as a workaround for older compilers
    complex(dp),intent(in) :: z

#ifdef COSHFIX
    my_sinh = 0.5_dp * (exp(z) - exp(-z))
#else
    my_sinh = sinh(z)
#endif

  end function my_sinh
!----------------------------------------------------------------------|
  subroutine InitExpM(n)
    !! Does nothing. Kept for backward compatibility: no work space needs to be
    !! allocated before calling [[ExpM]].
    integer,intent(in)::n !! matrix dimension (unused)
  end subroutine InitExpM
!----------------------------------------------------------------------|
  subroutine CleanExpM
    !! Does nothing. Kept for backward compatibility.
  end subroutine CleanExpM
!----------------------------------------------------------------------|
  subroutine choose_degree(nrm, m, s)
    !! Chooses the Padé degree m and the number of squarings s for a matrix of 1-norm nrm.
    real(dp),intent(in) :: nrm !! 1-norm of \(tH\)
    integer,intent(out) :: m !! degree of the Padé approximant
    integer,intent(out) :: s !! number of squarings
    integer :: i

    s = 0
    do i = 1, 4
       if(nrm <= theta(i)) then
          m = pade_deg(i)
          return
       end if
    end do
    m = 13
    if(nrm > theta(5)) s = max(0, ceiling(log(nrm / theta(5)) / log(2.0_dp)))

  end subroutine choose_degree
!----------------------------------------------------------------------|
  function zexpm(t, H) result(expH)
    !! Computes \(e^{tH}\) for a general complex matrix \(H\).
    real(dp), intent(in) :: t !! time argument
    complex(dp), dimension(:,:), intent(in) :: H !! square matrix \(H\)
    complex(dp), dimension(size(H,1),size(H,2)) :: expH !! result \(e^{tH}\)
    complex(dp), allocatable, dimension(:,:) :: A, A2, A4, A6, U, V, Q
    integer :: n, m, s, i, info
    integer :: ipiv(size(H,1))
    real(dp) :: nrm

    n = size(H,1)
    A = t * H
    nrm = maxval(sum(abs(A), dim=1))
    call choose_degree(nrm, m, s)
    if(s > 0) A = A / 2.0_dp**s

    A2 = matmul(A, A)
    select case(m)
    case(3)
       U = b3(3) * A2
       V = b3(2) * A2
       call add_identity(U, b3(1)); call add_identity(V, b3(0))
    case(5)
       A4 = matmul(A2, A2)
       U = b5(5) * A4 + b5(3) * A2
       V = b5(4) * A4 + b5(2) * A2
       call add_identity(U, b5(1)); call add_identity(V, b5(0))
    case(7)
       A4 = matmul(A2, A2)
       A6 = matmul(A4, A2)
       U = b7(7) * A6 + b7(5) * A4 + b7(3) * A2
       V = b7(6) * A6 + b7(4) * A4 + b7(2) * A2
       call add_identity(U, b7(1)); call add_identity(V, b7(0))
    case(9)
       A4 = matmul(A2, A2)
       A6 = matmul(A4, A2)
       Q = matmul(A6, A2)  ! A^8
       U = b9(9) * Q + b9(7) * A6 + b9(5) * A4 + b9(3) * A2
       V = b9(8) * Q + b9(6) * A6 + b9(4) * A4 + b9(2) * A2
       call add_identity(U, b9(1)); call add_identity(V, b9(0))
    case default
       A4 = matmul(A2, A2)
       A6 = matmul(A4, A2)
       U = matmul(A6, b13(13) * A6 + b13(11) * A4 + b13(9) * A2) &
            + b13(7) * A6 + b13(5) * A4 + b13(3) * A2
       V = matmul(A6, b13(12) * A6 + b13(10) * A4 + b13(8) * A2) &
            + b13(6) * A6 + b13(4) * A4 + b13(2) * A2
       call add_identity(U, b13(1)); call add_identity(V, b13(0))
    end select
    U = matmul(A, U)

    ! r_m(A) = (V - U)^{-1} (V + U)
    Q = V - U
    expH = V + U
    call zgesv(n, n, Q, n, ipiv, expH, n, info)
    if(info /= 0) error stop "ExpM: singular Pade denominator (ZGESV failed)"

    do i = 1, s
       expH = matmul(expH, expH)
    end do

  contains

    subroutine add_identity(X, c)
      complex(dp), intent(inout) :: X(:,:)
      real(dp), intent(in) :: c
      integer :: k
      do k = 1, size(X,1)
         X(k,k) = X(k,k) + c
      end do
    end subroutine add_identity

  end function zexpm
!----------------------------------------------------------------------|
  function dexpm(t, H) result(expH)
    !! Computes \(e^{tH}\) for a general real matrix \(H\).
    real(dp), intent(in) :: t !! time argument
    real(dp), dimension(:,:), intent(in) :: H !! square matrix \(H\)
    real(dp), dimension(size(H,1),size(H,2)) :: expH !! result \(e^{tH}\)
    real(dp), allocatable, dimension(:,:) :: A, A2, A4, A6, U, V, Q
    integer :: n, m, s, i, info
    integer :: ipiv(size(H,1))
    real(dp) :: nrm

    n = size(H,1)
    A = t * H
    nrm = maxval(sum(abs(A), dim=1))
    call choose_degree(nrm, m, s)
    if(s > 0) A = A / 2.0_dp**s

    A2 = matmul(A, A)
    select case(m)
    case(3)
       U = b3(3) * A2
       V = b3(2) * A2
       call add_identity(U, b3(1)); call add_identity(V, b3(0))
    case(5)
       A4 = matmul(A2, A2)
       U = b5(5) * A4 + b5(3) * A2
       V = b5(4) * A4 + b5(2) * A2
       call add_identity(U, b5(1)); call add_identity(V, b5(0))
    case(7)
       A4 = matmul(A2, A2)
       A6 = matmul(A4, A2)
       U = b7(7) * A6 + b7(5) * A4 + b7(3) * A2
       V = b7(6) * A6 + b7(4) * A4 + b7(2) * A2
       call add_identity(U, b7(1)); call add_identity(V, b7(0))
    case(9)
       A4 = matmul(A2, A2)
       A6 = matmul(A4, A2)
       Q = matmul(A6, A2)  ! A^8
       U = b9(9) * Q + b9(7) * A6 + b9(5) * A4 + b9(3) * A2
       V = b9(8) * Q + b9(6) * A6 + b9(4) * A4 + b9(2) * A2
       call add_identity(U, b9(1)); call add_identity(V, b9(0))
    case default
       A4 = matmul(A2, A2)
       A6 = matmul(A4, A2)
       U = matmul(A6, b13(13) * A6 + b13(11) * A4 + b13(9) * A2) &
            + b13(7) * A6 + b13(5) * A4 + b13(3) * A2
       V = matmul(A6, b13(12) * A6 + b13(10) * A4 + b13(8) * A2) &
            + b13(6) * A6 + b13(4) * A4 + b13(2) * A2
       call add_identity(U, b13(1)); call add_identity(V, b13(0))
    end select
    U = matmul(A, U)

    ! r_m(A) = (V - U)^{-1} (V + U)
    Q = V - U
    expH = V + U
    call dgesv(n, n, Q, n, ipiv, expH, n, info)
    if(info /= 0) error stop "ExpM: singular Pade denominator (DGESV failed)"

    do i = 1, s
       expH = matmul(expH, expH)
    end do

  contains

    subroutine add_identity(X, c)
      real(dp), intent(inout) :: X(:,:)
      real(dp), intent(in) :: c
      integer :: k
      do k = 1, size(X,1)
         X(k,k) = X(k,k) + c
      end do
    end subroutine add_identity

  end function dexpm
!----------------------------------------------------------------------|
  function zexpm_thrsv(t, H) result(expH)
    !! Computes \(e^{tH}\) for a general complex matrix \(H\). For 1x1 and 2x2
    !! matrices the exponential is computed analytically, otherwise by [[zexpm]].
    real(dp),parameter :: eps=1.0e-6_dp
    real(dp), intent(in) :: t !! time argument
    complex(dp), dimension(:,:), intent(in) :: H !! square matrix \(H\)
    complex(dp), dimension(size(H,1),size(H,2)) :: expH !! result \(e^{tH}\)
    complex(dp) :: deth,ct,st,exad

    select case(size(H,1))
    case(1)
       expH = exp(t*H)
    case(2)
       deth = t * sqrt((H(1,1)-H(2,2))**2 + 4.0_dp * H(1,2)*H(2,1))
       exad = exp(t * 0.5_dp*(H(1,1)+H(2,2)))
       if(abs(deth) < eps) then
          ct = 1.0_dp + (0.5_dp*deth)**2/2.0_dp + (0.5_dp*deth)**4/24.0_dp
          st = 0.5_dp + deth**2/48.0_dp + deth**4/3840.0_dp
       else
          ct = my_cosh(0.5_dp*deth) ; st = my_sinh(0.5_dp*deth)/deth
       end if
       expH(1,1) = ct + (H(1,1)-H(2,2))*st*t
       expH(1,2) = 2.0_dp*H(1,2) * st * t
       expH(2,1) = 2.0_dp*H(2,1) * st * t
       expH(2,2) = ct - (H(1,1)-H(2,2))*st * t
       expH = exad * expH
    case default
       expH = zexpm(t, H)
    end select

  end function zexpm_thrsv
!----------------------------------------------------------------------|
end module scitools_MatrixExp
