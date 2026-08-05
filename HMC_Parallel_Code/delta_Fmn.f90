program delta_Fmn
  use mpi_f08 , only: MPI_Init, MPI_Comm_rank, MPI_Comm_size, MPI_Finalize, MPI_comm, MPI_COMM_WORLD, MPI_Comm_Split
  use consts, only: PR
  use montecarlo, only: HMC_simulation

  implicit none

  type(MPI_comm) :: Communicator
  integer :: ierror, rank, nprocs
  integer :: NumberOfSimulations, NumberOfCores
  integer :: SimNo,MyRank
  integer :: color,key
  logical :: is_found
  character(len=20) :: value
  integer :: length, status

  call MPI_Init (ierror)
  call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierror)
  call MPI_Comm_size (MPI_COMM_WORLD, nprocs, ierror)


  !** Get number of simulations
  call get_command_argument(1,value,length,status)
  if(status /= 0)then
    if(rank == 0)write(6,*)'Argument specifying number of simulations missing'
    call MPI_Finalize(ierror)
    stop
  end if
  read(value,*)NumberOfSimulations
  
  !** Split the communicator for each simulation
  NumberOfCores=nprocs/NumberOfSimulations
  color=rank/NumberOfCores
  key=mod(rank,NumberOfCores)
  call MPI_Comm_Split(MPI_COMM_WORLD,color,key,Communicator,ierror)

  !** Perform the HMC simulation
  SimNo=color+1
  MyRank=key
  call HMC_simulation(Communicator,MyRank,SimNo)

  !** Clean up
  ! call MPI_Finalize (ierror)

end program delta_Fmn
