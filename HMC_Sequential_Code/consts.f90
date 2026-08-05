module consts
  implicit none
  public 
  save

  ! Kind types for 64-, 32-, 16-, and 8-bit signed integers
  integer, parameter :: INT64 = selected_int_kind(18)
  integer, parameter :: INT32 = selected_int_kind(9)
  integer, parameter :: INT16 = selected_int_kind(4)
  integer, parameter :: INT08 = selected_int_kind(2)

  ! Kind types for IEEE 754/IEC 60559 single- and double-precision reals
  integer, parameter :: IEEE32 = selected_real_kind(  6,  37 )
  integer, parameter :: IEEE64 = selected_real_kind( 15, 307 )


  !** Mathematical constants
  integer, parameter  :: PR=IEEE64
  real(PR), parameter :: PI=3.141592653589793238462643383279502884197_PR
  real(PR), parameter :: TWOPI=6.283185307179586476925286766559005768394_PR
  real(PR), parameter :: PIO2=1.57079632679489661923132169163975144209858_PR
  real(PR), parameter :: SQRT2=1.41421356237309504880168872420969807856967_PR

  
  !** Physical constants (Copied from NIST website)
  real(PR), parameter :: ATOMIC_MASS_UNIT=1.660538782e-27_PR
  real(PR), parameter :: ELECTRONIC_CHARGE=1.602176487e-19_PR
  real(PR), parameter :: BOLTZMANN_CONSTANT=1.3806504e-23_PR
  real(PR), parameter :: VACUUM_PERMITTIVITY=8.854187817e-12_PR
  real(PR), parameter :: AVOGADRO_NUMBER=1._PR/(ATOMIC_MASS_UNIT*1000._PR)
  real(PR), parameter :: MOLAR_GAS_CONSTANT=BOLTZMANN_CONSTANT*AVOGADRO_NUMBER

  !** Fundamental Units
  real(PR), parameter :: ANGSTROM=1.0e-10_PR
  real(PR), parameter :: PICO_SECOND=1.0e-12_PR
  real(PR), parameter :: UNIT_LENGTH=ANGSTROM
  real(PR), parameter :: UNIT_MASS=ATOMIC_MASS_UNIT
  real(PR), parameter :: UNIT_TIME=PICO_SECOND
  real(PR), parameter :: UNIT_CHARGE=ELECTRONIC_CHARGE

  !** Derived Units
  real(PR), parameter :: UNIT_VOLUME=UNIT_LENGTH**3
  real(PR), parameter :: UNIT_CONCENTRATION=1._PR/UNIT_VOLUME
  real(PR), parameter :: UNIT_DENSITY=UNIT_MASS/UNIT_LENGTH**3
  real(PR), parameter :: UNIT_ENERGY=UNIT_MASS*UNIT_LENGTH**2/UNIT_TIME**2
  real(PR), parameter :: UNIT_FORCE=UNIT_MASS*UNIT_LENGTH/UNIT_TIME**2
  real(PR), parameter :: UNIT_PRESSURE=UNIT_MASS/UNIT_LENGTH/UNIT_TIME**2

  !** Conversion Factors
  real(PR), parameter :: COULOMBIC_CONVERSION_FACTOR=UNIT_CHARGE**2/(4._PR*PI*VACUUM_PERMITTIVITY)/UNIT_LENGTH/UNIT_ENERGY
  real(PR), parameter :: K_B=BOLTZMANN_CONSTANT/UNIT_ENERGY

end module consts

