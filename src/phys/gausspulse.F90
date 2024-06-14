module scitools_gausspulse
!======================================================================================
!!  Gauss pulse class
   use scitools_def,only: dp,iu
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: GaussPulse_t
!--------------------------------------------------------------------------------------
   type GaussPulse_t
      !..............................................
      real(dp)    :: Tdur,E0,omega,s0,t0,phi
      complex(dp) :: pol(3)
      !..............................................
      contains
      !..............................................
      procedure, public  :: Init
      procedure, public  :: ReadFromFile
      procedure, public  :: Afield
      procedure, public  :: Efield
      !..............................................
   end type GaussPulse_t 
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,E0,omega,Tdur,s0,phi,pol)
      class(GaussPulse_t) :: me
      real(dp),intent(in) :: E0
      real(dp),intent(in) :: omega
      real(dp),intent(in) :: Tdur
      real(dp),intent(in) :: s0  
      real(dp),intent(in) :: phi
      complex(dp),intent(in) :: pol(3) 

      me%E0 = E0
      me%omega = omega
      me%Tdur = Tdur
      me%s0 = s0
      me%phi = phi
      me%pol = pol

      me%t0 = me%s0 * me%Tdur

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine ReadFromFile(me,fname)
      class(GaussPulse_t) :: me
      character(len=*),intent(in) :: fname
      integer :: unit_inp
      real(dp) :: Tdur,E0,omega,s0,phi
      real(dp) :: pol_re(3),pol_im(3)
      complex(dp) :: pol(3)
      namelist/GAUSSPULSE/Tdur,E0,omega,s0,phi,pol_re,pol_im

      open(newunit=unit_inp, file=trim(fname), status='old', action='read')
      read(unit_inp, nml=GAUSSPULSE)
      close(unit_inp)

      pol = pol_re + iu * pol_im

      call me%Init(E0, omega, Tdur, s0, phi, pol)

   end subroutine ReadFromFile
!--------------------------------------------------------------------------------------
   function Afield(me,t) result(Af)
      class(GaussPulse_t),intent(in) :: me
      real(dp),intent(in)            :: t
      real(dp) :: fenv
      real(dp) :: Af(3)
      complex(dp) :: phase

      fenv = exp(-0.5_dp * (t - me%t0)**2 / me%Tdur**2 )
      phase = exp(-iu * (me%omega * (t - me%t0) - me%phi) )

      Af = me%E0/me%omega * fenv * real(iu * me%pol * phase, kind=dp )

   end function Afield
!--------------------------------------------------------------------------------------
   function Efield(me,t) result(Ef)
      class(GaussPulse_t),intent(in) :: me
      real(dp),intent(in)            :: t
      real(dp) :: fenv
      real(dp) :: Ef(3)
      complex(dp) :: phase

      fenv = exp(-0.5_dp * (t - me%t0)**2 / me%Tdur**2 )
      phase = exp(-iu * (me%omega * (t - me%t0) - me%phi) )

      Ef = me%E0 * fenv * real(me%pol * phase, kind=dp )

   end function Efield
!--------------------------------------------------------------------------------------
 
!======================================================================================    
end module scitools_gausspulse
