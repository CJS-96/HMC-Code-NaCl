module nucleus
  use consts, only: PR, PI, INT64
  use spherical_harmonics, only: sp_hrmcs
  use datatypes, only: SpeciesInfo, ClusterInfo
  implicit none
  private
  save
  
  public :: ComputeClusterProperties

contains
  !====================================================================================================
  subroutine ComputeClusterProperties(Species,x,box,CutOffDistance,Cluster,CrystalList)
    type(SpeciesInfo), dimension(:), intent(in) :: Species
    real(PR), dimension(:), intent(in) :: x
    real(PR), intent(in) :: box
    real(PR), dimension(2), intent(in) :: CutOffDistance
    type(ClusterInfo), intent(out) :: Cluster
    logical, dimension(:), intent(out), optional :: CrystalList

    integer :: head
    integer :: NumberOfIons, NumberOfDissolvedIons, NumberOfWaterMolecules, TotalNumberOfIons, TotalNumberOfAtoms
    real(PR), dimension(2) :: CutOffDistanceSq

    logical, dimension(:,:,:), allocatable :: ConnectionMatrix
    logical, dimension(:,:), allocatable :: ClusterMatrix, CrystalMatrix
    logical, dimension(:), allocatable :: ClusterAtoms, ClusterHead, CrystalMolecules, CrystalHead
    integer, dimension(:), allocatable :: Neighbors
    real(PR), dimension(:,:), allocatable :: IonPosition, DissolvedIonPosition, OxygenPosition, AtomPosition
    integer, dimension(:), allocatable :: OrderParameterValue

    real(PR) :: vec(3), rijsq
    integer :: mol1, mol2

    NumberOfIons=Species(1)%NumberOfMolecules+Species(2)%NumberOfMolecules
    NumberOfDissolvedIons=Species(3)%NumberOfMolecules+Species(4)%NumberOfMolecules
    TotalNumberOfIons=NumberOfIons+NumberOfDissolvedIons
    NumberOfWaterMolecules=Species(5)%NumberOfMolecules
    TotalNumberOfAtoms=sum(Species%NumberOfMolecules*Species%NumberOfAtoms)
    CutOffDistanceSq=CutOffDistance**2

    !** Allocate Memory from Working Arrays
    Allocate(ConnectionMatrix(NumberOfIons, NumberOfIons, 2), &
             ClusterMatrix(NumberOfIons, NumberOfIons), &
             CrystalMatrix(NumberOfIons, NumberOfIons), &
             ClusterAtoms(NumberOfIons), &
             ClusterHead(NumberOfIons), &
             CrystalMolecules(NumberOfIons), &
             CrystalHead(NumberOfIons), &
             Neighbors(NumberOfIons), &
             IonPosition(3,NumberOfIons), &
             DissolvedIonPosition(3,NumberOfDissolvedIons), &
             OxygenPosition(3,NumberOfWaterMolecules), &
             AtomPosition(3,TotalNumberOfAtoms), &
             OrderParameterValue(NumberOfIons))

    !** Assign values to working arrays
    AtomPosition=reshape(x,(/3,TotalNumberOfAtoms/))
    IonPosition=AtomPosition(:,1:NumberOfIons)
    DissolvedIonPosition=AtomPosition(:,NumberOfIons+1:TotalNumberOfIons)
    OxygenPosition=AtomPosition(:,TotalNumberOfIons+1::3)

    !** Determine Connection Matrix
    Neighbors=0
    ConnectionMatrix=.false.
    do mol1=1,NumberOfIons
      ConnectionMatrix(mol1,mol1,:)=.true.
      do mol2=mol1+1,NumberOfIons
        vec=IonPosition(:,mol1)-IonPosition(:,mol2)
        vec=vec-box*anint(vec/box)
        rijsq=vec(1)**2+vec(2)**2+vec(3)**2
        if(rijsq < CutOffDistanceSq(1))then
          ConnectionMatrix(mol1,mol2,1)=.true.
          ConnectionMatrix(mol2,mol1,1)=.true.
          Neighbors(mol1)=Neighbors(mol1)+1
          Neighbors(mol2)=Neighbors(mol2)+1
        end if
        if(rijsq < CutOffDistanceSq(2))then
          ConnectionMatrix(mol1,mol2,2)=.true.
          ConnectionMatrix(mol2,mol1,2)=.true.
        end if
      end do
    end do

    !** Identify Clusters
    ClusterMatrix=ConnectionMatrix(:,:,1)
    call ClusterDeterminationAlgorithm(ClusterMatrix,ClusterHead)
    !** Identify Largest Cluster
    head=maxloc(count(ClusterMatrix,2),1)
    ClusterAtoms=ClusterMatrix(head,:)
    Cluster%Size=count(ClusterAtoms)
    if(Cluster%Size /= NumberOfIons)then
      return  !** Cluster criterion violated
    end if

    !** Check for association of solvated ions with cluster
    do mol1=1,NumberOfIons
      do mol2=NumberOfIons+1,TotalNumberOfIons
        vec=IonPosition(:,mol1)-DissolvedIonPosition(:,mol2-NumberOfIons)
        vec=vec-box*anint(vec/box)
        rijsq=vec(1)**2+vec(2)**2+vec(3)**2
        if(rijsq < CutOffDistanceSq(1))then
          Cluster%Size=Cluster%Size+1
          return                        !** Cluster criterion violated
        end if
      end do
    end do


    !** Compute and solvation state
    call ComputeClusterSolvationState(IonPosition,OxygenPosition,box,ClusterAtoms,Neighbors,Cluster)

    !** Determine crystallinity
    call ComputeOrderParameters(IonPosition,box,ConnectionMatrix(:,:,1),ClusterAtoms,OrderParameterValue)
    !** Identify CrystalMolecules
    where(OrderParameterValue >= 4)
      CrystalMolecules=.true.
    elsewhere
      CrystalMolecules=.false.
    end where

    !** Identify Crystals
    CrystalMatrix=ConnectionMatrix(:,:,2)
    call ClusterDeterminationAlgorithm(CrystalMatrix,CrystalHead,CrystalMolecules)
    !** Identify Largest Crystal
    head=maxloc(count(CrystalMatrix,2),1)
    Cluster%CrystalSize=maxval(count(CrystalMatrix,2))
    if(present(CrystalList))then
      CrystalList=CrystalMatrix(head,:)
    end if

    !** Deallocate Memory from Working Arrays
    deallocate(ConnectionMatrix, &
             ClusterMatrix, &
             CrystalMatrix, &
             ClusterAtoms, &
             ClusterHead, &
             CrystalMolecules, &
             CrystalHead, &
             Neighbors, &
             IonPosition, &
             DissolvedIonPosition, &
             OxygenPosition, &
             AtomPosition, &
             OrderParameterValue)

  end subroutine ComputeClusterProperties

  !====================================================================================================
  subroutine ComputeClusterSolvationState(IonPosition,OxygenPosition,box,ClusterAtoms,Neighbors,Cluster)
    real(PR), dimension(:,:), intent(in) :: IonPosition, OxygenPosition
    real(PR), intent(in) :: box
    logical, dimension(:), intent(in) :: ClusterAtoms
    integer, dimension(:), intent(in) :: Neighbors
    type(ClusterInfo), intent(inout) :: Cluster

    integer :: head, mol, icount, mol1, mol2
    integer :: NumberOfIons, NumberOfWaterMolecules
    real(PR) :: vec(3), rij, rijsq, SolvationState

    head=findloc(ClusterAtoms,.true.,1)
    NumberOfIons=size(IonPosition,2)
    NumberOfWaterMolecules=size(OxygenPosition,2)

    !Compute centroid of the Cluster
    Cluster%Centroid=IonPosition(:,head)*real(Neighbors(head)+1,PR)
    icount=Neighbors(head)+1
    do mol=head+1,NumberOfIons
      if(ClusterAtoms(mol))then
        vec=IonPosition(:,mol)-IonPosition(:,head) 
        vec=vec-box*anint(vec/box)
        Cluster%Centroid=Cluster%Centroid+(vec+IonPosition(:,head))*real(Neighbors(mol)+1,PR)
        icount=icount+Neighbors(mol)+1
      end if
    end do
    Cluster%Centroid=Cluster%Centroid/real(icount,PR)
     
    !Compute radius of cluster
    Cluster%Radius=0._pr
    do mol=1,NumberOfIons
      if(ClusterAtoms(mol))then
        vec=IonPosition(:,mol)-Cluster%Centroid
        vec=vec-box*anint(vec/box)
        rij=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
        if(rij > Cluster%Radius)Cluster%Radius=rij
      end if
    end do
  
    !Determine the solvation state of the nucleus
    Cluster%MinimumSolventDistance=100._pr
    do mol1=1,NumberOfWaterMolecules
      vec=OxygenPosition(:,mol1)-Cluster%Centroid
      vec=vec-box*anint(vec/box)
      rij=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
      Cluster%MinimumSolventDistance=min(Cluster%MinimumSolventDistance,rij)
    end do

    Cluster%SolvationState=int(10._PR*(Cluster%Radius-Cluster%MinimumSolventDistance))

  end subroutine ComputeClusterSolvationState

  !=========================================================================================================================
  subroutine ComputeOrderParameters(IonPosition,box,ConnectionMatrix,ClusterAtoms,OrderParameterValue)
    real(PR), dimension(:,:), intent(in) :: IonPosition
    real(PR), intent(in) :: box
    logical, dimension(:,:), intent(in) :: ConnectionMatrix
    logical, dimension(:), intent(in) :: ClusterAtoms
    integer, dimension(:), intent(out) :: OrderParameterValue

    integer :: mol1, mol2
    real(PR) :: rij(3), theta, phi
    integer :: N, m, NumberOfIons

    integer, parameter :: l=4
    complex(PR), dimension(:,:), allocatable :: qlm, qlm_tot,lqlm
    real(PR), dimension(:), allocatable :: ql,lql

    NumberOfIons=size(IonPosition,2)
    allocate(ql(NumberOfIons),qlm(-l:l,NumberOfIons))
    
    !** Initialize qlm and ql to zero
    ql=0._PR
    qlm=(0._PR, 0._PR)

    !** Unaveraged Order Parameter
    do mol1=1,NumberOfIons
      if(.not. ClusterAtoms(mol1))cycle
      N=0
      do mol2=1,NumberOfIons
        if(ConnectionMatrix(mol1,mol2) .and. mol1 /= mol2)then
          N=N+1
          rij=IonPosition(:,mol1)-IonPosition(:,mol2)
          rij=rij-box*anint(rij/box)
          theta=atan2(sqrt(rij(1)**2+rij(2)**2),rij(3))
          phi=atan2(rij(2),rij(1))
          do m=-l,l
            qlm(m,mol1)=qlm(m,mol1)+sp_hrmcs(l,m,theta,phi)
          end do
        end if
      end do
      qlm(:,mol1)=qlm(:,mol1)/real(N,PR)
      ql(mol1)=sqrt(real(dot_product(qlm(:,mol1),qlm(:,mol1))))
      qlm(:,mol1)=qlm(:,mol1)/ql(mol1)
    end do

    !** Averaged Order Parameter
    do mol1=1,NumberOfIons
      if(.not. ClusterAtoms(mol1))cycle
      N=0
      do mol2=1,NumberOfIons
        if(ConnectionMatrix(mol1,mol2) .and. mol1 /= mol2)then
          if(real(dot_product(qlm(:,mol1),qlm(:,mol2))) > 0.35_PR)N=N+1
        end if
      end do
      OrderParameterValue(mol1)=N
    end do

    deallocate(ql, qlm)
  end subroutine ComputeOrderParameters

  !====================================================================================================
  pure subroutine ClusterDeterminationAlgorithm(Matrix,Head,Subset_1)
    logical, dimension(:,:), intent(inout) :: Matrix
    logical, dimension(:), intent(out) :: Head
    logical, dimension(:), intent(in), optional :: Subset_1

    logical, dimension(:), allocatable :: Subset
    integer :: mol1, mol2, icount
    integer :: NumberOfAtoms

    NumberOfAtoms=size(Matrix,1)
    allocate(Subset(NumberOfAtoms))
    !** Remove non-subset atoms from Connection Matrix
    if(present(Subset_1))then
      Subset=Subset_1
      do concurrent (mol1=1:NumberOfAtoms)
        if(.not. Subset(mol1))then
          Matrix(mol1,:)=.false.
          Matrix(:,mol1)=.false.
        end if
      end do
    else
      Subset=.true.
    end if

    !** identify the clusters
    Head=.false.
    do mol1=1,NumberOfAtoms
      if (Subset(mol1)) then
        Subset(mol1)=.false.
        Head(mol1)=.true.
        !**keep combining rows until all a(mol1,mol2)=1 have solid(mol2)=false**
        do
          icount=0        
          do mol2=mol1+1,NumberOfAtoms
            if (Matrix(mol1,mol2) .and. Subset(mol2)) then
              icount=icount+1
              Subset(mol2)=.false.
              !**determine union of the two rows**
              Matrix(mol1,:)=Matrix(mol1,:) .or. Matrix(mol2,:)
            endif
          enddo
          if(icount == 0)exit
        end do
      endif
    enddo
    deallocate(Subset)
    
  end subroutine ClusterDeterminationAlgorithm

end module nucleus
