module mcmoves
  use liblammps
  use, intrinsic :: ISO_C_binding, only : C_ptr, C_double
  use mpi_f08 , only: MPI_Bcast, MPI_comm, MPI_DOUBLE_PRECISION, MPI_LOGICAL
  use consts, only: PR, MOLAR_GAS_CONSTANT
  use nucleus, only: ComputeClusterProperties
  use random, only: RandomNumber, GaussianRandomNumber
  use datatypes, only: SpeciesInfo, ClusterInfo, BiasInfo
  implicit none
  private
  save

  real(PR), parameter :: deltabox_max=1._PR

  public :: MD_move

contains
  !====================================================================================================
  subroutine MD_move(lmp,MyRank,Communicator,Temperature,TrajectoryLength,Species,Bias,CutOffDistance,Cluster,naccpt)
    type(LAMMPS), intent(in) :: lmp
    integer, intent(in) :: MyRank
    type(MPI_comm), intent(in) :: Communicator
    real(PR), intent(in) :: Temperature
    integer, intent(in) :: TrajectoryLength
    type(SpeciesInfo), dimension(:), intent(in) :: Species
    type(BiasInfo), intent(in) :: Bias
    real(PR), dimension(2), intent(in) :: CutOffDistance
    type(ClusterInfo), intent(inout) :: Cluster
    integer, intent(inout) :: naccpt

    real(PR) :: beta, vf
    logical :: accept

    real(PR), dimension(:), allocatable :: x, xn, v
    real(PR), pointer :: boxlo_ptr => NULL(),boxhi_ptr=>NULL()
    real(PR), target :: box
    real(C_double), pointer :: PotentialEnergy=>NULL(), KineticEnergy=>NULL()
    real(PR) ::  Energy, NewEnergy
    real(PR) :: arg, deltaEnergy, deltaBias_m, deltaBias_s, deltaBias
    integer :: natoms3, ierror
    type(ClusterInfo) :: NewCluster
    character(len=100) :: command_str

    !** Temporary variables
    beta = 4184._8/MOLAR_GAS_CONSTANT/Temperature
    vf = 1.0e-7_8*MOLAR_GAS_CONSTANT*Temperature    ! LAMMPS units [real units, g A^2 * fs^-2 mol^-1]

    !** Determine box size **
    boxlo_ptr=lmp%extract_global('boxxlo')
    boxhi_ptr=lmp%extract_global('boxxhi')
    box=boxhi_ptr-boxlo_ptr

    !** Gather atoms positions
    call lmp%gather_atoms('x',3, x )
    natoms3=size(x)
    !** Allocate memory to temporary arrays
    allocate(xn(natoms3), v(natoms3))
    !** Generate fresh velocity distribution and scatter among procs
    call GenerateVelocities(x,v,Species,box,vf)
    call MPI_Bcast(v, natoms3 , MPI_DOUBLE_PRECISION, 0, Communicator, ierror)
    call lmp%scatter_atoms('v',v)

    ! Calculate the energies at the start of Trajectory
    call lmp%command('run 0')
    PotentialEnergy=lmp%extract_compute('thermo_pe', lmp%style%global, lmp%type%scalar)
    KineticEnergy=lmp%extract_compute('thermo_ke', lmp%style%global, lmp%type%scalar)
    Energy = PotentialEnergy + KineticEnergy

    ! Run HMC simulation for one trajectory
    ! write(command_str,'(a4,i10)')'run ',TrajectoryLength
    write(command_str,'(a4,i10,a16)')'run ',TrajectoryLength, ' pre no post no'
    call lmp%command(trim(command_str))

    PotentialEnergy=lmp%extract_compute('thermo_pe', lmp%style%global, lmp%type%scalar)
    KineticEnergy=lmp%extract_compute('thermo_ke', lmp%style%global, lmp%type%scalar)
    NewEnergy = PotentialEnergy + KineticEnergy

    ! Gather new positions
    call lmp%gather_atoms('x' , 3, xn )

    ! Accept or reject the trajectory
    deltaEnergy=NewEnergy-Energy
    if(MyRank == 0)then
      !** Compute Cluster Properties
      call ComputeClusterProperties(Species,xn,box,CutOffDistance,NewCluster)

      !** Compute Bias
      deltaBias_m=Bias%k_m*real((NewCluster%CrystalSize-Bias%m0)**2-(Cluster%CrystalSize-Bias%m0)**2,PR)
      deltaBias_s=Bias%k_s*real((NewCluster%SolvationState-Bias%s0)**2-(Cluster%SolvationState-Bias%s0)**2,PR)
      deltaBias=deltaBias_m+deltaBias_s

      arg=beta*deltaEnergy+deltaBias
      if(NewCluster%Size == Cluster%Size)then
        if(arg < 0._PR)then
          accept=.true.
        elseif(arg < 100._PR)then
          accept=(RandomNumber() < exp(-arg))
        else
          accept=.false.
        end if
      else
        accept=.false.
      end if
      if(accept)then
          naccpt=naccpt+1
          Cluster=NewCluster
      end if
    end if

    call MPI_Bcast(accept, 1 , MPI_LOGICAL, 0, Communicator, ierror)
    if(.not. accept)then
      !** Scatter positions
      call lmp%scatter_atoms('x', x)
    end if

    !** Deallocate memory from temporary arrays
    deallocate(x, xn, v)

  end subroutine MD_move

  !====================================================================================================
  subroutine GenerateVelocities(x,v,Species,box,vf)
    real(PR), dimension(:), intent(in) :: x
    real(PR), dimension(:), intent(inout) :: v
    type(SpeciesInfo), dimension(:), intent(in) :: Species
    real(PR), intent(in) :: box,vf

    integer :: spc, mol, atm, icount
    
    integer :: llimit,ulimit
    real(PR), dimension(:,:,:), allocatable :: AtomPosition
    real(PR), dimension(3) :: vec, CenterOfMass, Omega, TranslationalVelocity, RotationalVelocity
    real(PR), dimension(3) :: PrincipalAxis1, PrincipalAxis2, PrincipalAxis3
    real(PR) :: Omega1, Omega2, Omega3
    real(PR) :: sigma
    real(PR) :: mag

    ! Moments of Inertia of water molecule
    real(PR), parameter :: I11=1.3440_PR, I22=1.9407610796401196_PR,I33=0.5967610796401204_PR

    icount=0
    !** Generate velocities for cations and anions
    do spc=1,4
      do mol=1,Species(spc)%NumberOfMolecules
        do atm=1,Species(spc)%NumberOfAtoms
          sigma=sqrt(vf/Species(spc)%AtomicWeight(atm))
          icount=icount+1
          v(icount)=sigma*GaussianRandomNumber()
          icount=icount+1
          v(icount)=sigma*GaussianRandomNumber()
          icount=icount+1
          v(icount)=sigma*GaussianRandomNumber()
        end do
      end do
    end do

    !** Generate velocities for Rigid Water
    spc=5
    llimit=3*sum(Species(1:4)%NumberOfMolecules*Species(1:4)%NumberOfAtoms)+1
    ulimit=3*sum(Species(1:5)%NumberOfMolecules*Species(1:5)%NumberOfAtoms)
    allocate(AtomPosition(3,Species(spc)%NumberOfAtoms,Species(spc)%NumberOfMolecules))
    AtomPosition=reshape(x(llimit:ulimit),(/3, Species(spc)%NumberOfAtoms,Species(spc)%NumberOfMolecules/))

    do mol=1,Species(spc)%NumberOfMolecules
      !** Generate Translation Velocity from Maxwell Boltzmann Distribution
      TranslationalVelocity(1)=sqrt(vf/Species(spc)%MolecularWeight)*GaussianRandomNumber()
      TranslationalVelocity(2)=sqrt(vf/Species(spc)%MolecularWeight)*GaussianRandomNumber()
      TranslationalVelocity(3)=sqrt(vf/Species(spc)%MolecularWeight)*GaussianRandomNumber()
      !** Generate Angular Velocities from Maxwell Boltzmann Distribution
      Omega1 = sqrt(vf/I11)*GaussianRandomNumber()
      Omega2 = sqrt(vf/I22)*GaussianRandomNumber()
      Omega3 = sqrt(vf/I33)*GaussianRandomNumber()

      !** Modify Atom Positions according to PBC
      do atm=2,Species(spc)%NumberOfAtoms
        vec=AtomPosition(:,atm,mol)-AtomPosition(:,1,mol)
        vec=vec-box*anint(vec/box)
        AtomPosition(:,atm,mol)=AtomPosition(:,1,mol)+vec
      end do

      ! axis1 is COM with bisector of Hydrogens 
      PrincipalAxis1=(AtomPosition(:,2,mol)+AtomPosition(:,3,mol))/2._pr-AtomPosition(:,1,mol)
      ! axis2 is normal to the plane of molecule 
      vec=AtomPosition(:,2,mol)-AtomPosition(:,3,mol)
      PrincipalAxis2=CROSS_PRODUCT(PrincipalAxis1,vec)
      ! axis3 is cross product of 1 and 2 
      PrincipalAxis3=CROSS_PRODUCT(PrincipalAxis1,PrincipalAxis2)
      ! Normalize the Principal Axis
      mag=dot_product(PrincipalAxis1,PrincipalAxis1)
      PrincipalAxis1=PrincipalAxis1/sqrt(mag)
      mag=dot_product(PrincipalAxis2,PrincipalAxis2)
      PrincipalAxis2=PrincipalAxis2/sqrt(mag)
      mag=dot_product(PrincipalAxis3,PrincipalAxis3)
      PrincipalAxis3=PrincipalAxis3/sqrt(mag)
      !** Angular Velocity Vector
      Omega=Omega1*PrincipalAxis1+Omega2*PrincipalAxis2+Omega3*PrincipalAxis3

      ! Compute Center Of Mass
      CenterOfMass=0._PR
      do atm=1,Species(spc)%NumberOfAtoms
        CenterOfMass=CenterOfMass+AtomPosition(:,atm,mol)*Species(spc)%AtomicWeight(atm)
      end do
      CenterOfMass=CenterOfMass/Species(spc)%MolecularWeight

      do atm=1,Species(spc)%NumberOfAtoms
        vec=AtomPosition(:,atm,mol)-CenterOfMass
        RotationalVelocity=CROSS_PRODUCT(Omega,vec)
        icount=icount+1
        v(icount)=TranslationalVelocity(1)+RotationalVelocity(1)
        icount=icount+1
        v(icount)=TranslationalVelocity(2)+RotationalVelocity(2)
        icount=icount+1
        v(icount)=TranslationalVelocity(3)+RotationalVelocity(3)
      end do
    end do

  end subroutine GenerateVelocities

  !=========================================================================================
  ! Calculates the cross product of two vectors
  !=========================================================================================
  pure function CROSS_PRODUCT(array1,array2) result(array)
    real(PR), dimension(3), intent(In) :: array1,array2
    real(PR), dimension(3)             :: array

    array(1)=array1(2)*array2(3)- array1(3)*array2(2)
    array(2)=array1(3)*array2(1)- array1(1)*array2(3)
    array(3)=array1(1)*array2(2)- array1(2)*array2(1)
  end function CROSS_PRODUCT

  pure subroutine swap(a,b)
    real(PR), intent(inout) :: a,b
    real(PR) :: temp
    temp=a
    a=b
    b=temp
  end subroutine swap

end module mcmoves

