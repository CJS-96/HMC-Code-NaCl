module montecarlo
  use liblammps
  use, intrinsic :: ISO_C_binding, only : C_ptr, c_null_ptr
  use mpi_f08 , only: MPI_Bcast, MPI_comm, MPI_INTEGER
  use json_module, only: json_file
  use consts, only: PR
  use nucleus, only: ComputeClusterProperties
  use random, only: RandomNumber, SetRandomNumberSeed
  use datatypes, only: ClusterInfo, SpeciesInfo, BiasInfo
  use mcmoves, only: MD_move
  implicit none
  private

  public :: HMC_simulation
contains
  !====================================================================================================
  subroutine HMC_simulation(Communicator,MyRank,SimNo)
    type(MPI_comm), intent(in) :: Communicator
    integer, intent(in) :: MyRank, SimNo

    type(LAMMPS) :: lmp
    type(json_file) :: json
    type(SpeciesInfo), dimension(5) :: Species
    type(ClusterInfo) :: Cluster
    type(BiasInfo) :: Bias

    real(PR), dimension(:), allocatable :: x
    real(PR), dimension(:), allocatable :: CutOffDistance
    real(PR), pointer :: boxlo_ptr => NULL(),boxhi_ptr=>NULL()
    real(PR), target :: box
    real(PR) :: Temperature, Pressure, dt, Boxlo, Boxhi

    integer, dimension(:), allocatable :: atype, NumberOfMolecules
    integer, dimension(0:200,-20:200) :: Histogram
    integer :: NumberOfEquilibrationCycles, NumberOfProductionCycles
    integer :: TrajectoryLength, DisplayFrequency, WriteFrequency
    integer :: iseed, CycleNumber, NumberOfCycles, ierror, NumberOfIons
    integer :: log_unit, cls_unit
    integer :: naccpt

    character(len=100) :: command_str
    character(len=20) :: input_file,data_file,log_file, cls_file, open_args(5)
    logical :: is_found

    !** Configure filenames                               Filename change
    if(SimNo < 10)then
      write(input_file,'(a11,i1)')'input.json.',SimNo
      write(data_file,'(a9,i1)')'data.lmp.',SimNo
      write(log_file,'(a14,i1)')'SimulationLog.',SimNo
      write(cls_file,'(a9,i1)')'cls_data.',SimNo
    elseif(SimNo < 100)then
      write(input_file,'(a11,i2)')'input.json.',SimNo
      write(data_file,'(a9,i2)')'data.lmp.',SimNo
      write(log_file,'(a14,i2)')'SimulationLog.',SimNo
      write(cls_file,'(a9,i2)')'cls_data.',SimNo
    elseif(SimNo < 1000)then
      write(input_file,'(a11,i3)')'input.json.',SimNo
      write(data_file,'(a9,i3)')'data.lmp.',SimNo
      write(log_file,'(a14,i3)')'SimulationLog.',SimNo
      write(cls_file,'(a9,i3)')'cls_data.',SimNo
    else
      write(input_file,'(a11,i4)')'input.json.',SimNo
      write(data_file,'(a9,i4)')'data.lmp.',SimNo
      write(log_file,'(a14,i4)')'SimulationLog.',SimNo
      write(cls_file,'(a9,i4)')'cls_data.',SimNo
    end if

    !** Read Input json file
    call json%initialize()
    call json%load_file(trim(input_file))
    call json%get('Number of molecules', NumberOfMolecules, is_found)
    call json%get('Temperature', Temperature, is_found)
    call json%get('Number of equilibration cycles', NumberOfEquilibrationCycles, is_found)
    call json%get('Number of production cycles', NumberOfProductionCycles, is_found)
    call json%get('Trajectory length', TrajectoryLength, is_found)
    call json%get('Display frequency', DisplayFrequency, is_found)
    call json%get('Write frequency', WriteFrequency, is_found)
    call json%get('timestep',dt,is_found)
    call json%get('Cutoff distance', CutOffDistance, is_found)
    call json%get('iseed', iseed, is_found)
    call json%get('k_m', Bias%k_m, is_found)
    call json%get('m0', Bias%m0, is_found)
    call json%get('k_s', Bias%k_s, is_found)
    call json%get('s0', Bias%s0, is_found)
    call json%destroy()
    
    Species%NumberOfMolecules=NumberOfMolecules

    Species(1)%NumberOfAtoms=1         !** Cation
    allocate(Species(1)%AtomicWeight(Species(1)%NumberOfAtoms))
    Species(1)%AtomicWeight=(/22.99_8/)
    Species(1)%MolecularWeight=sum(Species(1)%AtomicWeight)

    Species(2)%NumberOfAtoms=1        !** Anion
    allocate(Species(2)%AtomicWeight(Species(2)%NumberOfAtoms))
    Species(2)%AtomicWeight=(/35.453_8/)
    Species(2)%MolecularWeight=sum(Species(2)%AtomicWeight)

    NumberOfIons=Species(1)%NumberOfMolecules+Species(2)%NumberOfMolecules

    Species(3)%NumberOfAtoms=1         !** Dissolved Cation
    allocate(Species(3)%AtomicWeight(Species(3)%NumberOfAtoms))
    Species(3)%AtomicWeight=(/22.99_8/)
    Species(3)%MolecularWeight=sum(Species(3)%AtomicWeight)

    Species(4)%NumberOfAtoms=1        !** Dissolved Anion
    allocate(Species(4)%AtomicWeight(Species(4)%NumberOfAtoms))
    Species(4)%AtomicWeight=(/35.453_8/)
    Species(4)%MolecularWeight=sum(Species(4)%AtomicWeight)

    Species(5)%NumberOfAtoms=3        !** Water
    allocate(Species(5)%AtomicWeight(Species(5)%NumberOfAtoms))
    Species(5)%AtomicWeight=(/15.99_8, 1.008_8, 1.008_8/)
    Species(5)%MolecularWeight=sum(Species(5)%AtomicWeight)

    !** Random Number
    call SetRandomNumberSeed(iseed)

    ! Open LAMMPS
    open_args(1) = 'liblammps'
    open_args(2) = '-log'
    open_args(3) = 'none'
    open_args(4) = '-screen'
    open_args(5) = 'none'
    lmp=lammps(open_args, Communicator%MPI_VAL)

    !** Read Input File
    call lmp%file('in.lmp.1')

    !** Read Data File
    write(command_str,'(2a)')'read_data ',trim(data_file)
    call lmp%command(trim(command_str))

    !** Read Input File
    call lmp%file('in.lmp.2')

    !** Set integrator
    call lmp%command('fix NVE all nve')

    !** Timestep
    write(command_str,'(a,f10.4)')'timestep',dt
    call lmp%command(trim(command_str))

    ! Create a Log file and cls_file   !! Append cls_file and create new one too
    if(MyRank == 0)then
      open(newunit=log_unit,file=trim(log_file),position='append',action='write')
      write(log_unit,'(a10,6(4x,a10))')'#Iter','n','m','s','radius','rW','naccpt'
      flush(unit=log_unit)
      open(newunit=cls_unit,file=trim(cls_file),form='unformatted',action='write')
    end if
    
    !===================================================================================================
    !** Preliminary steps before simulation
    !===================================================================================================
    !** Gather atoms positions
    call lmp%gather_atoms('x',3,x)

    !** Determine box size **
    boxlo_ptr=lmp%extract_global('boxxlo')
    boxhi_ptr=lmp%extract_global('boxxhi')
    box=boxhi_ptr-boxlo_ptr

    !** Initialize Accumulators
    Histogram=0
    naccpt=0
    CycleNumber=0

      !** Determine Cluster Properties
      if(MyRank == 0)then
         call ComputeClusterProperties(Species,x,box,CutOffDistance,Cluster)
         write(log_unit,'(4(i10,4x),2(f10.4,4x),i10)') &
              CycleNumber, Cluster%Size, Cluster%CrystalSize, Cluster%SolvationState, Cluster%Radius, Cluster%MinimumSolventDistance, naccpt
         flush(unit=log_unit)
      end if

    !===================================================================================================
    !** Simulation Cycle
    !===================================================================================================
    NumberOfCycles=NumberOfEquilibrationCycles+NumberOfProductionCycles
    do CycleNumber=1,NumberOfCycles
      call MD_move(lmp,MyRank,Communicator,Temperature,TrajectoryLength,Species,Bias,CutOffDistance,Cluster,naccpt)

      if(MyRank == 0)then
        ! Print to Log file
        if(mod(CycleNumber,DisplayFrequency) == 0)then
          write(log_unit,'(4(i10,4x),2(f10.4,4x),i10)') &
            CycleNumber, Cluster%Size, Cluster%CrystalSize, Cluster%SolvationState, Cluster%Radius, Cluster%MinimumSolventDistance, naccpt
          flush(unit=log_unit)
        end if

        !** Write cluster data
        write(cls_unit)Cluster%Size,Cluster%CrystalSize,Cluster%SolvationState

        !** Write to new cls_file

        !** Accumulate Histograms
        if(CycleNumber > NumberOfEquilibrationCycles)then
          Histogram(Cluster%CrystalSize,Cluster%SolvationState)=Histogram(Cluster%CrystalSize,Cluster%SolvationState)+1
        end if
      end if
      if(mod(CycleNumber,WriteFrequency) == 0)then
        if(SimNo<10)then
          write(command_str,'(a26,i1)')'write_data data.inter.lmp.',SimNo
        else if(SimNo<100)then
          write(command_str,'(a26,i2)')'write_data data.inter.lmp.',SimNo
        else if(SimNo<1000)then
          write(command_str,'(a26,i3)')'write_data data.inter.lmp.',SimNo
        else
          write(command_str,'(a26,i4)')'write_data data.inter.lmp.',SimNo
        endif
        call lmp%command(trim(command_str)) 
      end if
    end do

    !** Gather Information from LAMMPS
    call lmp%gather_atoms('type',1,atype)
    call lmp%gather_atoms('x',3,x)
    !** Determine box size **
    boxlo_ptr=lmp%extract_global('boxxlo')
    boxhi_ptr=lmp%extract_global('boxxhi')
    box=boxhi_ptr-boxlo_ptr
    if(MyRank == 0)then 
      close(unit=log_unit)
      call PrintResults(Bias,Histogram,SimNo)
      call PrintFinalConfig(x,atype,box,Species,NumberOfMolecules,CutOffDistance,SimNo)
    end if

    !** Save Final Configuration  ! Filename change
    if(SimNo<10)then
      write(command_str,'(a26,i1)')'write_data data.final.lmp.',SimNo
    else if(SimNo<100)then
      write(command_str,'(a26,i2)')'write_data data.final.lmp.',SimNo
    else if(SimNo<1000)then
      write(command_str,'(a26,i3)')'write_data data.final.lmp.',SimNo
    else
      write(command_str,'(a26,i4)')'write_data data.final.lmp.',SimNo
    endif
    call lmp%command(trim(command_str))

    !** Close LAMMPS
    call lmp%close(.TRUE.)

  end subroutine HMC_simulation

  !=========================================================================================
  subroutine PrintResults(Bias,Histogram,SimNo)
    type(BiasInfo), intent(in) :: Bias
    integer, dimension(0:200,-20:200), intent(in) :: Histogram
    integer, intent(in) :: SimNo


    character(len=100) :: bias_file, hist_file
    integer :: bias_unit, hist_unit
    integer :: m, s
    real(PR) :: bias_pot

     if(SimNo<10)then
       write(bias_file,'(a9,i1)')'BiasInfo.',SimNo
       write(hist_file,'(a10,i1)')'Histogram.',SimNo
     else if(SimNo< 100)then
       write(bias_file,'(a9,i2)')'BiasInfo.',SimNo
       write(hist_file,'(a10,i2)')'Histogram.',SimNo
     else if(SimNo< 1000)then
       write(bias_file,'(a9,i3)')'BiasInfo.',SimNo
       write(hist_file,'(a10,i3)')'Histogram.',SimNo
     else
       write(bias_file,'(a9,i4)')'BiasInfo.',SimNo
       write(hist_file,'(a10,i4)')'Histogram.',SimNo
     end if
       
    open(newunit=bias_unit,file=trim(bias_file),action='write')
    write(bias_unit,'(a5,2x,a1,f10.4)')'k_n',':',Bias%k_n
    write(bias_unit,'(a5,2x,a1,i4)')'n0',':',Bias%n0
    write(bias_unit,'(a5,2x,a1,f10.4)')'k_s',':',Bias%k_s
    write(bias_unit,'(a5,2x,a1,i4)')'s0',':',Bias%s0
    close(unit=bias_unit)

    !** Histogram   !! Filename change
    open(newunit=hist_unit,file=trim(hist_file),action='write')
    write(hist_unit,'(a15,4x,a15,4x,a15,4x,a15)')'#Crystal Size','Water Content','Count','Bias'
    do m=lbound(Histogram,1),ubound(Histogram,1)
      do s=lbound(Histogram,2),ubound(Histogram,2)
        bias_pot=Bias%k_m*real((m-Bias%m0)**2,PR)+Bias%k_s*real((s-Bias%s0)**2,PR)
        if(Histogram(m,s) /= 0)write(hist_unit,'(i15,4x,i15,4x,i15,4x,es15.4)')m,s,Histogram(m,s),bias_pot
      end do
    end do
    close(unit=hist_unit)

  end subroutine PrintResults

  !=========================================================================================
  subroutine PrintFinalConfig(x,atype,box,Species,NumberOfMolecules,CutOffDistance,SimNo)
    real(PR), dimension(:), intent(in) :: x
    integer, dimension(:), intent(in) :: atype
    real(PR), intent(in):: box
    type(SpeciesInfo), dimension(:), intent(in) :: Species
    integer, dimension(:), intent(in) :: NumberOfMolecules
    real(PR), dimension(:), intent(in) :: CutOffDistance
    integer, intent(in) :: SimNo

    character(len=100) :: salt_file, nacl_file, brine_file, water_file
    integer :: i, icount, n_ions
    integer :: salt_unit, nacl_unit, brine_unit, water_unit

    integer :: NumberOfIons, NumberOfDissolvedIons, NumberOfWaterAtoms
    type(ClusterInfo) :: Cluster
    logical, dimension(:), allocatable :: CrystalList
    real(PR), dimension(3) :: vec

    !** Allocate Memory to Working Arrays
    NumberOfIons=Species(1)%NumberOfMolecules+Species(2)%NumberOfMolecules
    NumberOfDissolvedIons=Species(3)%NumberOfMolecules+Species(4)%NumberOfMolecules
    NumberOfWaterAtoms=Species(5)%NumberOfMolecules*Species(5)%NumberOfAtoms
    allocate(CrystalList(NumberOfIons))

    !** Determine Cluster Properties
    call ComputeClusterProperties(Species,x,box,CutOffDistance,Cluster,CrystalList)

    !** Print Results to file   Filename change
    if(SimNo < 10)then
      write(salt_file,'(a5,i1,a4)')'Salt_',SimNo,'.xyz'
      write(nacl_file,'(a5,i1,a4)')'NaCl_',SimNo,'.xyz'
      write(brine_file,'(a6,i1,a4)')'Brine_',SimNo,'.xyz'
      write(water_file,'(a6,i1,a4)')'Water_',SimNo,'.xyz'
    else if(SimNo < 100)then
      write(salt_file,'(a5,i2,a4)')'Salt_',SimNo,'.xyz'
      write(nacl_file,'(a5,i2,a4)')'NaCl_',SimNo,'.xyz'
      write(brine_file,'(a6,i2,a4)')'Brine_',SimNo,'.xyz'
      write(water_file,'(a6,i2,a4)')'Water_',SimNo,'.xyz'
    else if(SimNo < 1000)then
      write(salt_file,'(a5,i3,a4)')'Salt_',SimNo,'.xyz'
      write(nacl_file,'(a5,i3,a4)')'NaCl_',SimNo,'.xyz'
      write(brine_file,'(a6,i3,a4)')'Brine_',SimNo,'.xyz'
      write(water_file,'(a6,i3,a4)')'Water_',SimNo,'.xyz'
    else
      write(salt_file,'(a5,i4,a4)')'Salt_',SimNo,'.xyz'
      write(nacl_file,'(a5,i4,a4)')'NaCl_',SimNo,'.xyz'
      write(brine_file,'(a6,i4,a4)')'Brine_',SimNo,'.xyz'
      write(water_file,'(a6,i4,a4)')'Water_',SimNo,'.xyz'
    end if

    !** Final Configuration
    n_ions=count(CrystalList,1)
    open(newunit=salt_unit,file=trim(salt_file),action='write')
    write(salt_unit,'(i10)')n_ions
    write(salt_unit,'(a10)')'Salt'
    open(newunit=nacl_unit,file=trim(nacl_file),action='write')
    write(nacl_unit,'(i10)')NumberOfIons-n_ions
    write(nacl_unit,'(a10)')'NaCl'
    open(newunit=brine_unit,file=trim(brine_file),action='write')
    write(brine_unit,'(i10)')NumberOfDissolvedIons
    write(brine_unit,'(a10)')'Brine'
    open(newunit=water_unit,file=trim(water_file),action='write')
    write(water_unit,'(i10)')NumberOfWaterAtoms
    write(water_unit,'(a10)')'Water'
    icount=1
    do i=1,size(atype)
      vec=x(icount:icount+2)-Cluster%Centroid
      vec=vec-box*anint(vec/box)
      select case(atype(i))
      case (1)
        if(CrystalList(i))then
          write(salt_unit,'(a3,3(4x,f10.4))')'Na',vec
        else
          write(nacl_unit,'(a3,3(4x,f10.4))')'Na',vec
        end if
      case (2)
        if(CrystalList(i))then
          write(salt_unit,'(a3,3(4x,f10.4))')'Cl',vec
        else
          write(nacl_unit,'(a3,3(4x,f10.4))')'Cl',vec
        end if
      case (3)
        write(brine_unit,'(a3,3(4x,f10.4))')'K',vec
      case (4)
        write(brine_unit,'(a3,3(4x,f10.4))')'F',vec
      case (5)
        write(water_unit,'(a3,3(4x,f10.4))')'O',vec
      case (6)
        write(water_unit,'(a3,3(4x,f10.4))')'H',vec
      end select
      icount=icount+3
    end do
    close(unit=salt_unit)
    close(unit=nacl_unit)
    close(unit=brine_unit)
    close(unit=water_unit)
    !** deAllocate Memory to Working Arrays
    deallocate(CrystalList)

  end subroutine PrintFinalConfig

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

end module montecarlo

