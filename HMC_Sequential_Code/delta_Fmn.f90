program delta_Fmn
  use mpi_f08 , only: MPI_Init, MPI_Comm_rank, MPI_Comm_size, MPI_Finalize, MPI_comm, MPI_COMM_WORLD, MPI_Comm_Split
  use consts, only: PR
  use montecarlo, only: HMC_simulation

  implicit none

  type(MPI_comm) :: Communicator
  integer :: ierror, rank, nprocs
  integer :: MyRank

  call MPI_Init (ierror)
  call MPI_Comm_rank (MPI_COMM_WORLD, MyRank, ierror)
  call MPI_Comm_size (MPI_COMM_WORLD, nprocs, ierror)

  !** Perform the HMC simulation
  Communicator=MPI_COMM_WORLD
  call HMC_simulation(Communicator,MyRank)

end program delta_Fmn
