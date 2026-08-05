!==================================================================================
module random
!  An implementation of the Mersenne Twister algorithm for generating
!  psuedo-random sequences.
!
!  ORIGINAL ALGORITHM COPYRIGHT
!  ============================
!  Copyright (C) 1997,2002 Makoto Matsumoto and Takuji Nishimura.
!  Any feedback is very welcome. For any question, comments, see
!  http://www.math.keio.ac.jp/matumoto/emt.html or email
!  matumoto@math.keio.ac.jp
!========================================================================
!
!  For more information on this software package, please visit
!  Scott's web site, Coyote Gulch Productions, at:
!      http://www.coyotegulch.com
!========================================================================

  use consts, only: INT08, INT16, INT32, INT64, IEEE32, IEEE64, PR, PI, TWOPI

  implicit none

  !==============================================================================
  ! Everything is private unless explicitly made public
  private

  public :: SetRandomNumberSeed,RandomNumber,GaussianRandomNumber
  public :: GetRandomNumberState,SetRandomNumberState
  public :: mtprng_state
  public :: RandomUnitVector, RandomRotationMatrix, UniformRandomRotationMatrix

  !==============================================================================
  ! Constants
  integer(INT32), parameter :: N = 624_INT32
  integer(INT32), parameter :: M = 397_INT32

  !==============================================================================
  ! types
  type ::  mtprng_state
      integer(INT32)                   :: mti = -1
      integer(INT64), dimension(0:N-1) :: mt
  end type mtprng_state

  type(mtprng_state), save :: state

contains
  !==============================================================================
  !  Initializes the generator with "seed"
  !==============================================================================
  subroutine SetRandomNumberSeed(seed)
    ! arguments
    integer(INT32), intent(in) :: seed
    
    ! working storage
    integer :: i
    integer(INT64) :: s, b

    ! save seed        
    state%mt(0) = seed
    
    ! Set the seed using values suggested by Matsumoto & Nishimura, using
    !   a generator by Knuth. See original source for details.
    do i = 1, N - 1
      state%mt(i) = iand(4294967295_INT64,1812433253_INT64 * ieor(state%mt(i-1),ishft(state%mt(i-1),-30_INT64)) + i)
    end do
    
    state%mti = N

  end subroutine SetRandomNumberSeed
  
  !==============================================================================
  !   Obtain a psuedorandom real number in the range (0,1), i.e., a number
  !   greater than 0 and less than 1.
  !==============================================================================
  function RandomNumber() result(r)
  
    !return type
    real(IEEE64) :: r

    integer(INT64) :: ir

    ! internal constants
    integer(INT64), dimension(0:1), parameter :: mag01 = (/ 0_INT64, -1727483681_INT64 /)

    ! Period parameters
    integer(INT64), parameter :: UPPER_MASK =  2147483648_INT64
    integer(INT64), parameter :: LOWER_MASK =  2147483647_INT64

    ! Tempering parameters
    integer(INT64), parameter :: TEMPERING_B = -1658038656_INT64
    integer(INT64), parameter :: TEMPERING_C =  -272236544_INT64
    
    ! Note: variable names match those in original example
    integer(INT32) :: kk
    
    ! Local constant; precalculated to avoid division below
    real(IEEE64), parameter :: factor = 1.0_IEEE64 / 4294967296.0_IEEE64

    ! Generate N words at a time
    if (state%mti >= N) then
      ! The value -1 acts as a flag saying that the seed has not been set.
      if (state%mti == -1) call SetRandomNumberSeed(4357_INT32)
      
      ! Fill the mt array
      do kk = 0, N - M - 1
        ir = ior(iand(state%mt(kk),UPPER_MASK),iand(state%mt(kk+1),LOWER_MASK))
        state%mt(kk) = ieor(ieor(state%mt(kk + M),ishft(ir,-1_INT64)),mag01(iand(ir,1_INT64)))
      end do
      
      do kk = N - M, N - 2
        ir = ior(iand(state%mt(kk),UPPER_MASK),iand(state%mt(kk+1),LOWER_MASK))
        state%mt(kk) = ieor(ieor(state%mt(kk + (M - N)),ishft(ir,-1_INT64)),mag01(iand(ir,1_INT64)))
      end do
      
      ir = ior(iand(state%mt(N-1),UPPER_MASK),iand(state%mt(0),LOWER_MASK))
      state%mt(N-1) = ieor(ieor(state%mt(M-1),ishft(ir,-1)),mag01(iand(ir,1_INT64)))
      
      ! Start using the array from first element
      state%mti = 0
    end if
    
    ! Here is where we actually calculate the number with a series of
    !   transformations 
    ir = state%mt(state%mti)
    state%mti = state%mti + 1
    
    ir = ieor(ir,ishft(ir,-11))
    ir = iand(4294967295_INT64,ieor(ir,iand(ishft(ir, 7),TEMPERING_B)))
    ir = iand(4294967295_INT64,ieor(ir,iand(ishft(ir,15),TEMPERING_C)))
    ir = ieor(ir,ishft(ir,-18))
    
    r = (real(ir,IEEE64) + 0.5_IEEE64) * factor

  end function RandomNumber

  !============================================================================
  !  Generate gaussian distributed random numbers
  !============================================================================
  function GaussianRandomNumber()
    real(IEEE64) :: GaussianRandomNumber

    real(IEEE64) :: rsq,v1,v2
    real(IEEE64), save :: g
    logical, save :: gaus_stored=.false.

    if(gaus_stored)then
      GaussianRandomNumber=g
      gaus_stored=.false.
    else
      do
        v1=RandomNumber()
        v2=RandomNumber()
        v1=2._IEEE64*v1-1._IEEE64
        v2=2._IEEE64*v2-1._IEEE64
        rsq=v1**2+v2**2
        if(rsq > 0._IEEE64 .and. rsq < 1._IEEE64)exit
      end do
      rsq=sqrt(-2._IEEE64*log(rsq)/rsq)
      GaussianRandomNumber=v1*rsq
      g=v2*rsq
      gaus_stored=.true.
    end if
  end function GaussianRandomNumber
  !==============================================================================

  subroutine GetRandomNumberState(rstate)
    type(mtprng_state), intent(out) :: rstate

    rstate=state
  end subroutine GetRandomNumberState
  !==============================================================================

  subroutine SetRandomNumberState(rstate)
    type(mtprng_state), intent(in) :: rstate

    state=rstate
  end subroutine SetRandomNumberState
  !==============================================================================

  !==========================================================================
  !** Subroutine obtains a unit vector with random orientation
  !** Ref: Algorithm 40, page 410 'Understanding Molecular Simulation'
  !** Authors: Frenkel and Smit, 1st edition
  !==========================================================================
  subroutine RandomUnitVector(vec)
    real(PR), dimension(3), intent(OUT) :: vec

    real(PR) :: ran1,ran2,ransq,ranh

    ransq=2._PR
    do while(ransq >= 1._PR)
      ran1=1._PR-2._PR*RandomNumber()
      ran2=1._PR-2._PR*RandomNumber()
      ransq=ran1*ran1+ran2*ran2
    end do
    ranh=2._PR*sqrt(1._PR-ransq)
    vec(1)=ran1*ranh
    vec(2)=ran2*ranh
    vec(3)=(1._PR-2._PR*ransq)
  end subroutine RandomUnitVector

  !==========================================================================
  !** Subroutine obtains a Uniform Rotation Matrix for a random rotation using
  !** quarternions. Algorithm by Ken Shoemake
  !==========================================================================
  subroutine UniformRandomRotationMatrix(RotationMatrix)
    real(PR), dimension(3,3), intent(out) :: RotationMatrix

    real(PR) :: u1,u2,u3
    real(PR) :: q0,q1,q2,q3

    u1=RandomNumber()
    u2=RandomNumber()
    u3=RandomNumber()

    q0=sqrt(1._PR-u1)*sin(TWOPI*u2)
    q1=sqrt(1._PR-u1)*cos(TWOPI*u2)
    q2=sqrt(u1)*sin(TWOPI*u3)
    q3=sqrt(u1)*cos(TWOPI*u3)
    RotationMatrix(1,1)=q0**2+q1**2-q2**2-q3**2
    RotationMatrix(1,2)=2._pr*(q1*q2+q0*q3)
    RotationMatrix(1,3)=2._pr*(q1*q3-q0*q2)
    RotationMatrix(2,1)=2._pr*(q1*q2-q0*q3)
    RotationMatrix(2,2)=q0**2-q1**2+q2**2-q3**2
    RotationMatrix(2,3)=2._pr*(q2*q3+q0*q1)
    RotationMatrix(3,1)=2._pr*(q1*q3+q0*q2)
    RotationMatrix(3,2)=2._pr*(q2*q3-q0*q1)
    RotationMatrix(3,3)=q0**2-q1**2-q2**2+q3**2
  end subroutine UniformRandomRotationMatrix

  !==========================================================================
  !** Subroutine obtains a Rotation Matrix for a random rotation using
  !** Rodrigues Rotation Formula
  !==========================================================================
  subroutine RandomRotationMatrix(RotationMatrix,theta_max)
    real(PR), dimension(3,3), intent(out) :: RotationMatrix
    real(PR), intent(in)                  :: theta_max

    real(PR) :: theta,vec(3)
    real(PR) :: c,w,s

    !** Obtain a random axis of rotation
    call RandomUnitVector(vec)
    !** Obtain a random angle of rotation
    theta=(2._PR*RandomNumber()-1._PR)*theta_max
    !** Calculate Rotation Matrix using Rodrigues Rotation Formula
    c=cos(theta)
    w=1._PR-c
    s=sin(theta)
    RotationMatrix(1,1)=vec(1)*vec(1)*w+c
    RotationMatrix(1,2)=vec(1)*vec(2)*w-vec(3)*s
    RotationMatrix(1,3)=vec(1)*vec(3)*w+vec(2)*s
    RotationMatrix(2,1)=vec(1)*vec(2)*w+vec(3)*s
    RotationMatrix(2,2)=vec(2)*vec(2)*w+c
    RotationMatrix(2,3)=vec(2)*vec(3)*w-vec(1)*s
    RotationMatrix(3,1)=vec(1)*vec(3)*w-vec(2)*s
    RotationMatrix(3,2)=vec(2)*vec(3)*w+vec(1)*s
    RotationMatrix(3,3)=vec(3)*vec(3)*w+c
  end subroutine RandomRotationMatrix
end module random
!==================================================================================
