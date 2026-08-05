module datatypes
  use consts, only: PR
  implicit none
  public

  type :: SpeciesInfo
    real(PR), dimension(:), allocatable :: AtomicWeight
    real(PR) :: MolecularWeight
    integer :: NumberOfAtoms, NumberOfMolecules
  end type SpeciesInfo

  type :: BiasInfo 
    real(PR) :: k_n, k_m, k_s
    integer :: n0, m0, s0
  end type BiasInfo

  type :: ClusterInfo
    integer :: Size
    integer :: CrystalSize
    integer :: SolvationState

    real(PR), dimension(3) :: Centroid
    real(PR) :: Radius, MinimumSolventDistance
  end type ClusterInfo

end module datatypes

