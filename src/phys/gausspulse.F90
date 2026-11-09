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
      !! Linearly or elliptically polarized laser pulse with Gaussian envelope,
      !! \( \mathbf{E}(t) \propto E_0\, e^{-(t-t_0)^2/(2 T_{dur}^2)} \mathrm{Re}[\mathbf{p}\, e^{-i(\omega (t-t_0) - \phi)}] \).
      !..............................................
      real(dp)    :: Tdur,E0,omega,s0,t0,phi !! duration, amplitude, frequency, center (in units of Tdur), center time, carrier-envelope phase
      real(dp)    :: tmin, tmax !! time window \(t_0 \pm 8 T_{dur}\)
      complex(dp) :: pol(3) !! complex polarization vector
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
      !! Initializes the pulse from its parameters.
      class(GaussPulse_t) :: me
      real(dp),intent(in) :: E0 !! field amplitude \(E_0\)
      real(dp),intent(in) :: omega !! carrier frequency \(\omega\)
      real(dp),intent(in) :: Tdur !! pulse duration \(T_{dur}\) (standard deviation of the envelope)
      real(dp),intent(in) :: s0 !! pulse center in units of Tdur, \(t_0 = s_0 T_{dur}\)
      real(dp),intent(in) :: phi !! carrier-envelope phase \(\phi\)
      complex(dp),intent(in) :: pol(3) !! complex polarization vector

      me%E0 = E0
      me%omega = omega
      me%Tdur = Tdur
      me%s0 = s0
      me%phi = phi
      me%pol = pol

      me%t0 = me%s0 * me%Tdur

      me%tmin = me%t0 - 8.0_dp * me%Tdur
      me%tmax = me%t0 + 8.0_dp * me%Tdur

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine ReadFromFile(me,fname)
      !! Reads the pulse parameters from the namelist `&GAUSSPULSE` (Tdur, E0, omega, s0, phi, pol_re, pol_im).
      class(GaussPulse_t) :: me
      character(len=*),intent(in) :: fname !! input file
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
      !! Vector potential \(\mathbf{A}(t)\).
      class(GaussPulse_t),intent(in) :: me
      real(dp),intent(in)            :: t !! time
      real(dp) :: fenv
      real(dp) :: Af(3)
      complex(dp) :: phase

      fenv = exp(-0.5_dp * (t - me%t0)**2 / me%Tdur**2 )
      phase = exp(-iu * (me%omega * (t - me%t0) - me%phi) )

      Af = me%E0/me%omega * fenv * real(iu * me%pol * phase, kind=dp )

   end function Afield
!--------------------------------------------------------------------------------------
   function Efield(me,t) result(Ef)
      !! Electric field \(\mathbf{E}(t)\).
      class(GaussPulse_t),intent(in) :: me
      real(dp),intent(in)            :: t !! time
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
