!gfortran G_6.10_DNA.f90 -fimplicit-none -fcheck=all -g -Wall -Werror -fmax-errors=5

PROGRAM First

! parameters
REAL, parameter :: pi=3.14159265
INTEGER, parameter :: NL(3,26) = reshape([ &
  -1, 0,  0, &
   1, 0,  0, &
   0, -1, 0, &
   0, 1,  0, &
   0, 0, -1, &
   0, 0,  1, &
  -1, -1, 0, &
  -1, 1,  0, &
   1, -1, 0, &
   1, 1,  0, &
  -1, 0, -1, &
  -1, 0,  1, &
   1, 0, -1, &
   1, 0,  1, &
   0, -1, -1, &
   0, -1, 1, &
   0, 1, -1, &
   0, 1,  1, &
  -1, -1, -1, &
  -1, -1, 1, &
  -1, 1, -1, &
  -1, 1,  1, &
   1, -1, -1, &
   1, -1, 1, &
   1, 1, -1, &
   1, 1,  1  &
], shape=[3,26])

! imported values from key file
INTEGER, DIMENSION(:,:), Allocatable :: Energy, Connection, LL
INTEGER, DIMENSION(:), Allocatable :: ConnectionN, StartDropletNk, ChemPot
INTEGER :: keyfileN, Tmax, Ti, Trelax, NkN, BoxSize(3), StartDropletR, StartSlabH
INTEGER, DIMENSION(:), Allocatable :: Nk, mNk, Vk, StartSlabNk
LOGICAL :: StartDroplet, StartImported, UseSlitherII, UseSlitherI, UseGrand, StartSlab
LOGICAL :: WriteLatticeEnd, WriteIndiDist_Cluster, WriteBoundTo
INTEGER :: E_Sanity, E_Lattice, E_Nk, E_LargestCluster, E_ClusterHist, E_BoundSites, E_BondTypes, E_Energy, E_PrintStep
INTEGER :: E_RadDist_Cluster, E_RadDist_All
INTEGER :: E_XYZDist_None, E_XYZDist_All, E_XYZDist_Cluster
INTEGER :: E_RG_All, E_RG_Cluster
LOGICAL :: HasDNA
! imported Joey values
REAL :: LSteric
REAL, DIMENSION(:,:), Allocatable :: SLEq
INTEGER, DIMENSION(:,:), Allocatable :: SLCon, SLPot
INTEGER, DIMENSION(:), Allocatable :: SLConN, Vak
INTEGER :: ESelf, LinkerType
INTEGER :: E_SpringLinkerVector, E_Network, E_SelfLoop
LOGICAL :: UseSpringLinkerEq, UseESelf, WriteCoreCluster
! Global Simulation parameters
INTEGER, DIMENSION(:,:,:), Allocatable :: Lattice
INTEGER, DIMENSION(:,:), Allocatable :: List
INTEGER, DIMENSION(:), Allocatable :: BoundTo
INTEGER :: ETC
! Global Book keeping
REAL :: MCMoves(8)
INTEGER :: Rej(1:8,1:4), Propose(1:8), MCMove
! Lookup tables
INTEGER, DIMENSION(:), Allocatable :: sNk, sVk, smNk, sNkVk
INTEGER, DIMENSION(:), Allocatable :: BeadType, BeadBounds
LOGICAL, Dimension(:), Allocatable :: CanRot
INTEGER :: UseSlitherI_ListN, UseSlitherII_ListN, UseGrand_ListN
INTEGER, DIMENSION(:), Allocatable :: UseSlitherI_List, UseSlitherII_List, UseGrand_List
INTEGER, DIMENSION(:), Allocatable :: BeadToProteinIndex, BeadToProteinType, BeadToMSI
INTEGER, DIMENSION(:), Allocatable :: ProteinToProteinType
! Global Joey Book keeping
INTEGER :: SelfLoopN


! Local Book Keeping
Character (len=100) :: FileName, Str, Str1, Str2, Str3, keyfile, Parm
Logical :: flag
INTEGER :: t,t2,m,i,i1,i2,i3,i4,i5,pot, CoreN, CoreHalfN, bonds
INTEGER :: j, j2, dx1(1:3),dx2(1:3),dx(1:3),dx_max
INTEGER :: EM, EC
INTEGER :: MS, MS2, MST, MST2, MSI, MSI2
INTEGER :: MP, MPT, MPT2
INTEGER :: x1(1:3),x2(1:3),x3(1:3),Move(1:3),R(1:3),ForwardRMN,BackRMN,RotMoves(1:27),RotMove, MCN
INTEGER :: AllHostNN, SlitherTimes, SlitherForward,NB_Cluster
INTEGER, DIMENSION(:,:), Allocatable :: Xall,XProt,AllHostN, SlitherList
INTEGER, DIMENSION(:), Allocatable :: Lall,Neighbor,AllHost,FreeHost, SlitherRot, Fisher(:), DNAHost
INTEGER :: MS_front, MS_back, FisherN, FreeHostN, DNAHostN
REAL :: logForwardRMN, logBackRMN, y1(3), y2(3), y3(3), GT(3,3), b, c, d, p, q, theta, phi, rr
REAL :: rr1, rr3(3), COM_All(3), COM_Cluster(3), leq
INTEGER :: tclock1,tclock2,tclock3,tclock4,tclock_sanity,iclock
LOGICAL, DIMENSION(:), Allocatable :: Analysis_log, IsHost
INTEGER, DIMENSION(:,:,:), Allocatable :: RadDist_ij,XDist_ij,YDist_ij,ZDist_ij
INTEGER, DIMENSION(:,:), Allocatable :: BondTypes,RadDist,XDist,YDist,ZDist
INTEGER, DIMENSION(:), Allocatable :: BoundSites,ClusterHist,LargestCluster
REAL, DIMENSION(:), Allocatable :: RG_All,RG_Cluster
REAL, DIMENSION(:,:), Allocatable :: COMi

! local Joey variables
INTEGER, DIMENSION(:,:), Allocatable :: LinkerCon
INTEGER, DIMENSION(:), Allocatable :: Core, CoreHalf
LOGICAL, DIMENSION(:), Allocatable :: IsFull,IsHalfFull
REAL, DIMENSION(:), Allocatable :: RG_Core,RG_CoreHalf,Occupancy,SLVec

INTEGER :: DebugInt(1:10)
REAL :: DebugReal(1:10)

REAL :: time1,time2

!!!
!!!

call cpu_time(time1)

t=0
DebugInt=0
DebugReal=0
call system_clock(tclock1)
call ReadKey
write(*,*) 'Tmax:'
write(*,*) Tmax
call BasicSetup


write(*,*) 'MCMoves:'
write(*,*) MCMoves/sum(MCMoves)
do i1=2,8
  MCMoves(i1)=MCMoves(i1-1)+MCMoves(i1)
enddo
MCMoves=MCMoves/MCMoves(8)

call InitialConditions

!write(*,*) 'Sanity check before starting'
!call Sanity

write(*,*) 'Analysis of initial conditions'
call Analysis

write(*,*) 'All looks good!'
call system_clock(tclock2)
write(*,'(A, I17)')   'Time Spent Preparing Simulation (milliseconds):', tclock2-tclock1
write(*,*) 'Starting Simulation'

call system_clock(tclock1)

do t=1,Tmax
  !Nucleation process, if there is a nucleation time set as Trelax
  do t2=1,Ti
    if (t .le. Trelax) then
      call random_number(rr1)
      rr1=rr1*MCMoves(5)
    elseif (sNk(NkN+1)==0) then
      call MoveGrand
      cycle
    else
      call random_number(rr1)
    endif
    !write(*,*) MCmove
    
    
    do MCMove=1,8
    if (rr1<MCMoves(MCMove)) then
      call ChooseMove
      exit
    endif
    enddo
  enddo
  call Analysis
enddo

call cpu_time(time2)


call EndofProgram

write(*,*) ETC

!!!
!!!
!!!
!!!
Contains
subroutine FindParm(TargetParm)
character (len=*) :: TargetParm
!write(*,*) Targetparm
open(2,file=keyfile)
do i=1,1000
read(2,*) Str
if (Str == TargetParm) then
read(2,'(A)') Parm
exit
elseif (Str == 'END') then
write(*,*) 'Missing Key word: ', TargetParm, ', Crashing the program'
stop
endif
enddo
close(2)
end subroutine FindParm
!!!!
!!!!
!!!!
!!!!
subroutine ReadKey
integer :: nargs

nargs=Command_Argument_Count()
if (nargs/=1) then
write(*,*) 'The program requires a single input, the key file'
stop
endif
call Get_Command_Argument(1,keyfile,keyfileN,m)
if (m/=0) then
write(*,*) 'I forgot what the "m" flag tells me but it should be zero'
write(*,*) 'I think I should program it to crash just in case'

stop
endif
!!
write(*,*) 'Reading Key File'
call FindParm('Tmax')
read(Parm,*) Tmax
call FindParm('Ti')
read(Parm,*) Ti
call FindParm('Trelax')!Tunc means there is no cluster move in first few Tunc steps
read(Parm,*) Trelax
call FindParm('BoxSize')
read(Parm,*) BoxSize
call Findparm('NkN')
read(Parm,*) NkN
Allocate( Nk(1:NkN) )
Allocate( mNk(1:NkN) )
Allocate( Vk(1:NkN) )
Allocate( Vak(1:NkN) )
Allocate( StartDropletNk(1:NkN) )
Allocate( StartSlabNk(1:NkN) )
call Findparm('Nk')
read(Parm,*) Nk
call Findparm('Vk')
read(Parm,*) Vk
call Findparm('Vak')
read(Parm,*) Vak
!! DNA start
call Findparm('HasDNA')
read(Parm,*) HasDNA
if (HasDNA .and. Vk(1) /= BoxSize(3)) then
  write(*,*) 'DNA must match box size! Crashing program'
  stop
endif
if (HasDNA .and. Nk(1)/=1) then
  write(*,*) 'Only supports having 1 stretched DNA strand in the box and you have ',Nk(1)
  stop
endif
!!
call Findparm('UseESelf')
read(Parm,*) UseESelf
if (UseESelf) then
  call Findparm('LSteric')
  read(Parm,*) LSteric
  call Findparm('ESelf')
  read(Parm,*) ESelf
  if (ESelf==0) then ! if not using, then not using
    UseESelf=.false.
  endif
endif

call Findparm('StartImported')
read(Parm,*) StartImported
call Findparm('StartDroplet')
read(Parm,*) StartDroplet
call Findparm('StartSlab')
read(Parm,*) StartSlab

if (HasDNA .and. StartDroplet) then
  write(*,*) 'DNA is incompatible with starting with a droplet! Crashing program'
  stop
endif
if (HasDNA .and. StartSlab) then
  write(*,*) 'DNA is incompatible with starting with a slab! Crashing program'
  stop
endif
if (StartDroplet .and. StartSlab) then
  write(*,*) 'You need to pick between starting with a slab and starting with a droplet'
  stop
endif
if (StartDroplet) then
  call Findparm('StartDropletR')
  read(Parm,*) StartDropletR
  call Findparm('StartDropletNk')
  read(Parm,*) StartDropletNk
else
  StartDropletNk=0
endif
if (StartSlab) then
  call Findparm('StartSlabH')
  read(Parm,*) StartSlabH
  call Findparm('StartSlabNk')
  read(Parm,*) StartSlabNk
else
  StartSlabNk=0
endif

MCMoves=0
call Findparm('iMoveRot')
read(Parm,*) i1
MCMoves(1)=real(i1)
call Findparm('iMoveTransI')
read(Parm,*) i1
MCMoves(2)=real(i1)
call Findparm('iMoveTransII')
read(Parm,*) i1
MCMoves(3)=real(i1)
call Findparm('iMoveSlitherI')
read(Parm,*) i1
MCMoves(4)=real(i1)
call Findparm('iMoveSlitherII')
read(Parm,*) i1
MCMoves(5)=real(i1)
call Findparm('iMoveClusterI')
read(Parm,*) i1
MCMoves(6)=real(i1)
call Findparm('iMoveClusterII')
read(Parm,*) i1
MCMoves(7)=real(i1)
call Findparm('UseGrand')
read(Parm,*) UseGrand
if (UseGrand) then
  UseGrand_ListN=0
  Allocate( UseGrand_List(1:NkN) )
  Allocate( ChemPot(1:NkN) )
  call Findparm('ChemPot')
  read(Parm,*) ChemPot
  do i1=merge(2,1,HasDNA),NkN !never use stretched DNA
  if (ChemPot(i1)/=0) then
    UseGrand_ListN=UseGrand_ListN+1
    UseGrand_List(UseGrand_ListN)=i1
  endif
  enddo
  if (UseGrand_ListN==0) then
    write(*,*) 'Nothing is using the grand canonical moves, turning them off'
    UseGrand=.false.
  endif
  call Findparm('iMoveGrand')
  read(Parm,*) i1
  MCMoves(8)=real(i1)
  if (i1==0) then
    UseGrand=.false.
  endif
endif
if (UseGrand) then
  call Findparm('MaxNk')
  read(Parm,*) mNk
else
  MCMoves(8)=0.
  mNk=Nk
endif

! lookup table for common sums
ALLOCATE ( sVk(1:NkN+1) )
ALLOCATE ( sNk(1:NkN+1) )
ALLOCATE ( smNk(1:NkN+1) )
ALLOCATE ( sNkVk(1:NkN+1) )
ALLOCATE ( BeadBounds(1:NkN+1) )
sVk=0
sNk=0
smNk=0
sNkVk=0
BeadBounds=0
do i1=1,NkN
  sVk(i1+1)=sum(Vk(1:i1))
  sNk(i1+1)=sum(Nk(1:i1))
  smNk(i1+1)=sum(mNk(1:i1))
  sNkVk(i1+1)=sum(Nk(1:i1)*Vk(1:i1))
  BeadBounds(i1+1) = sum(mNk(1:i1)*Vk(1:i1))
enddo


call Findparm('E_Sanity')
read(Parm,*) E_Sanity
if (E_Sanity>Tmax) E_Sanity=0
call Findparm('E_Lattice')
read(Parm,*) E_Lattice
if (E_Lattice>Tmax) E_Lattice=0
if (E_Lattice>0) then
  call Findparm('WriteBoundTo')
  read(Parm,*) WriteBoundTo
endif
if (UseGrand) then
  call Findparm('E_Nk')
  read(Parm,*) E_Nk
  if (E_Lattice/=0 .and. E_Nk>E_Lattice) then
    write(*,*) 'must print Nk as often as lattice if using GrandCannonical moves'
    write(*,*) 'Adjusting E_Nk from',E_Nk,'to',E_Lattice
    E_Nk=E_Lattice
  endif
  if (mod(E_Lattice,E_Nk)/=0) then
    write(*,*) 'Im trying to make sure that every lattice has an Nk printed with this check.'
    write(*,*) 'It looks like it fails but I may have done this math wrong'
    write(*,*) 'Crashing until you fix the code'
    stop
  endif
else
  E_Nk=0
endif


call Findparm('E_LargestCluster')
read(Parm,*) E_LargestCluster
if (E_LargestCluster>Tmax) E_LargestCluster=0
call Findparm('E_ClusterHist')
read(Parm,*) E_ClusterHist
if (E_ClusterHist>Tmax) E_ClusterHist=0
call Findparm('E_RG_All')
read(Parm,*) E_RG_All
if (E_RG_All>Tmax) E_RG_All=0
call Findparm('E_RG_Cluster')
read(Parm,*) E_RG_Cluster
if (E_RG_Cluster>Tmax) E_RG_Cluster=0
call Findparm('E_BoundSites')
read(Parm,*) E_BoundSites
if (E_BoundSites>Tmax) E_BoundSites=0
call Findparm('E_BondTypes')
read(Parm,*) E_BondTypes
if (E_BondTypes>Tmax) E_BondTypes=0
call Findparm('E_Energy')
read(Parm,*) E_Energy
if (E_Energy>Tmax) E_Energy=0
call Findparm('E_PrintStep')
read(Parm,*) E_PrintStep
if (E_PrintStep>Tmax) E_PrintStep=0   
call Findparm('E_RadDist_All')
read(Parm,*) E_RadDist_All
if (E_RadDist_All>Tmax) E_RadDist_All=0
call Findparm('E_XYZDist_None')
read(Parm,*) E_XYZDist_None
if (E_XYZDist_None>Tmax) E_XYZDist_None=0
call Findparm('E_XYZDist_All')
read(Parm,*) E_XYZDist_All
if (E_XYZDist_All>Tmax) E_XYZDist_All=0
call Findparm('E_RadDist_Cluster')
read(Parm,*) E_RadDist_Cluster
if (E_RadDist_Cluster>Tmax) E_RadDist_Cluster=0
call Findparm('E_XYZDist_Cluster')
read(Parm,*) E_XYZDist_Cluster
if (E_XYZDist_Cluster>Tmax) E_XYZDist_Cluster=0
call Findparm('WriteIndiDist_Cluster')
read(Parm,*) WriteIndiDist_Cluster
If(WriteIndiDist_Cluster) then  ! Allocate memory
  ALLOCATE(COMi(3,1:NkN))
  COMi=0
  write(*,*) COMi
endif
call Findparm('E_Network')
read(Parm,*) E_Network
call Findparm('E_SelfLoop')
read(Parm,*) E_SelfLoop
call Findparm('WriteCoreCluster')
read(Parm,*) WriteCoreCluster
if(WriteCoreCluster .and. StartSlab) then
  WriteCoreCluster = .false.
  write(*,*) 'Currently core cluster is only defined in terms of sperical droplets'
  stop
endif
If(WriteCoreCluster) then  ! Allocate memory
  ALLOCATE(Core(1:smNk(NkN+1)))
  ALLOCATE(CoreHalf(1:smNk(NkN+1))) ! CoreHalf is for half-occupied protein
  ALLOCATE(IsFull(1:smNk(NkN+1)))
  ALLOCATE(IsHalfFull(1:smNk(NkN+1)))
  ALLOCATE(Occupancy(1:smNk(NkN+1)))
endif
call Findparm('WriteLatticeEnd')
read(Parm,*) WriteLatticeEnd
call Findparm('E_SpringLinkerVector')
read(Parm,*) E_SpringLinkerVector
call Findparm('Seed')
read(Parm,*) iclock
if (iclock==0) then
  call system_clock(tclock2)
  i1=getpid()
  tclock2 = ieor(tclock2, int(i1,kind(tclock2)))
  iclock=mod(tclock2,1000000000)
  write(*,'(A, I12)') ' Assigning a random seed:', iclock
  ! call random_seed(put=iclock*[1,1,1,1,1,1,1,1,1,1,1,1])
  call random_seed(put=iclock*[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])

else
  write(*,'(A, I12)') ' Assigning the keyfile seed:', iclock
  ! call random_seed(put=iclock*[1,1,1,1,1,1,1,1,1,1,1,1])
  call random_seed(put=iclock*[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
endif
open(1,file='B_seed.txt')
write(1,*) iclock
close(1)

!! Import Beads
write(*,*) 'Importing bead types'
call FindParm('BeadTypeFile')
FileName=Parm
Allocate( BeadType(1:sVk(NkN+1)) )
open(1,file=FileName)
do i=1,sVk(NkN+1)
  read(1,*) BeadType(i)
enddo
close(1)


!! Import Energy
write(*,*) 'Importing Energy'
call Findparm('EnergyFile')
FileName=Parm
Allocate(Energy(1:maxval(BeadType),1:maxval(BeadType)))
Energy=0
open(1,file=FileName)
do i=1,size(Energy,1)
read(1,*) Energy(i,:)
end do
close(1)

!sanity for energy
do i1=1,size(Energy,1)
do i2=1,size(Energy,1)
if (Energy(i1,i2)/=Energy(i2,i1)) then
write(*,*) 'Energy Not Symmetric'
write(*,*) 'I don`t know what to do with myself without a symmetric energy table!'
write(*,*) i1,i2
write(*,*) Energy(i1,i2),Energy(i2,i1)
stop
endif
enddo
enddo


write(*,*) 'Importing ConnectionN, Connection, and Linker Length'
Allocate(ConnectionN(1:sVk(NkN+1)))! ConnectionN[i] is the number of connections of i
ConnectionN=0
call Findparm('ConnectionFile')
FileName=Parm
open(1,file=FileName)
do
  i1=0
  i3=0
  read(1,*,iostat=m) i1,i2,i3
  if (m < 0) exit ! end of file
  if (i1==0 .and. i3==0) cycle ! blank line
  if (i1==i2) then
    write(*,*) 'Connection file has beads bonded to themselves? crashing'
    stop
  elseif (i3==0) then
    write(*,*) 'Connection file has incomplete lines, crashing'
    stop
  endif
  ConnectionN(i1)=ConnectionN(i1)+1
  ConnectionN(i2)=ConnectionN(i2)+1
end do
Allocate( Connection(1:max(maxval(ConnectionN),1),1:sVk(NkN+1)) )
Allocate( LL(1:maxval(ConnectionN),1:sVk(NkN+1)) )
rewind(1)
Connection=0
ConnectionN=0
do
  i1=0
  i3=0
  read(1,*,iostat=m) i1,i2,i3
  if (m < 0) exit ! end of file
  if (i1==0 .and. i3==0) cycle ! blank line
  ConnectionN(i1)=ConnectionN(i1)+1
  Connection(ConnectionN(i1),i1)=i2
  LL(ConnectionN(i1),i1)=i3
  ConnectionN(i2)=ConnectionN(i2)+1
  Connection(ConnectionN(i2),i2)=i1
  LL(ConnectionN(i2),i2)=i3
enddo
close(1)
! sanity

do i1=1,sVk(NkN+1)
if (ConnectionN(i1)>0) then
  if (count(Connection(1:ConnectionN(i1),i1)<i1)>1) then ! loop?
    write(*,*) 'You either have a loop or have ordered your beads in a weird and unexpected way'
    write(*,*) 'crashing'
    stop
  endif
  if (Connection(1,i1)/=minval(Connection(1:ConnectionN(i1),i1))) then
    write(*,*) 'First connection isnt with the smallest indexed bead it could interact with'
    write(*,*) 'Something funny and unexpected with how your ordered the beads'
    write(*,*) 'crashing'
    stop 
  endif
endif
enddo
do i1=1,NkN
if (Vk(i1) - sum(ConnectionN(1+sVk(i1):sVk(i1+1)))/2 /=1) then
  write(*,*) 'You dont have a loop but the number of bonds dont add up correctly'
  write(*,*) 'crashing'
  stop
endif
enddo
do i1=1,NkN
do i2=1+sVk(i1),sVk(i1+1)
if (any(Connection(1:ConnectionN(i2),i2)<=sVk(i1)) &
    .or. any(Connection(1:ConnectionN(i2),i2)>sVk(i1+1)) ) then
  write(*,*) 'Protein connection file says that a protein is covalently bound to another protein!'
  write(*,*) 'You messed up the connection file'
  write(*,*) 'Crashing'
  stop
endif
enddo
enddo

if (UseGrand) then ! adjust the chempot for the approximate entropy of the chain
do i1=1,UseGrand_ListN
  MPT=UseGrand_List(i1)
  MSI=sVk(MPT)
  do i2=2,Vk(MPT)
    ChemPot(MPT) = ChemPot(MPT) - ceiling(log(real((2*LL(1,MSI+i2)+1)**3-1))*1000)
  enddo
enddo
endif

!Import SLCon, SLEq and SLPot
call Findparm('UseSpringLinkerEq')
read(Parm,*) UseSpringLinkerEq
if (UseSpringLinkerEq) then
  write(*,*) 'Using Equilibrium Linker Length Potential'
  call Findparm('SpringLinkerEqFile')
  FileName=Parm
  open(1,file=FileName)
  ! Import Connection
  Allocate( SLConN(1:sVk(NkN+1)) )! SLConN[i] is the number of SL connections of i
  write(*,*) 'Importing Equilibrium Spring Linker ConnectionN, ConnectionL, Length and Potential'
  SLConN=0
  do
    read(1,*,iostat=m) i,j,leq,pot
    if (m<0) exit
    if (m > 0) then
      write(*,*) 'The Spring linker file was not formatted correctly.  Crashing'
      stop
    endif
    if(i >= j) then
      write(*,*) "Error in Spring Linker Connection File, index i must be smaller than j"
      stop
    endif
    SLConN(i)=SLConN(i)+1
    SLConN(j)=SLConN(j)+1
  end do
  Allocate( SLCon(maxval(SLConN),sVk(NkN+1)) )
  Allocate(SLEq( maxval(SLConN),sVk(NkN+1)))
  Allocate(SLPot(maxval(SLConN),sVk(NkN+1)))
  ! Read in the data
  
  
  SLCon=0
  SLConN(i)=0
  SLConN(j)=0
  rewind(1)
  do
    read(1,*,iostat=m) i,j,leq,pot
    if (m<0) exit
    SLConN(i)=SLConN(i)+1
    SLConN(j)=SLConN(j)+1
    SLCon(SLConN(i),i)=j
    SLCon(SLConN(j),j)=i
    SLEq( SLConN(i),i)=leq
    SLEq( SLConN(j),j)=leq
    SLPot(SLConN(i),i)=pot
    SLPot(SLConN(j),j)=pot
  end do
  close(1)
  !sanity check
  do i1=1,sVk(NkN+1)
  do i2=1,SLConN(i1)-1
    if (any(SLCon(i2+1:SLConN(i1),i1)==SLCon(i2,i1))) then
      write(*,*) 'Spring linker file doesnt have unique linkers, crashing'
      stop
    endif
  enddo
  enddo
  do i1=1,NkN
  do i2=1,Vk(i1)
    MSI=i2+sVk(i1)
    do i3=1,SLConN(MSI)
      if (SLCon(i3,MSI)>sVk(i1+1)) then
        write(*,*) 'Spring linkers crosses between protein types????, crashing!'
        stop
      endif
    enddo
  enddo
  enddo
  
  print *,   "Finished Equilibrium Spring Linker Importation"
  print *, "Spring Linker ConnectionN:"
  do i1=1,sVk(NkN+1)
    print *, SLConN(i1)
  enddo
  print *, "Spring Linker Connection:"
  do i1=1,sVk(NkN+1)
    print *, SLCon(:,i1)
  enddo
  print *, "Spring Linker Equilibrium Length:"
  do i1=1,sVk(NkN+1)
    print *, SLEq(:,i1)
  enddo
  print *, "Spring Linker Length Potential:"
  do i1=1,sVk(NkN+1)
    print *, SLPot(:,i1)
  enddo

  !Use analysis for Spring linker or not
  if (E_SpringLinkerVector /= 0) then
    write(*,*) 'Write Spring Linker Vector Information'
  endif
else
    write(*,*) 'Not Using Equilibrium Spring Linker Length Potential'
endif


end subroutine ReadKey
!!!!
!!!!
!!!!
!!!!
subroutine BasicSetup
write(*,*) 'Starting Basic Setup'
write(*,*) 'Box Size is', BoxSize
!Check if rotmove can work, turns the move type off if there is nothing that can rotate

! Can a bead rotate
ALLOCATE ( CanRot(maxval(BeadType)) )
CanRot=.false.
do i1=1,maxval(BeadType)
if (any(Energy(i1,:)/=0)) then
 CanRot(i1)=.true.
endif
enddo
if (count(CanRot)==0 .and. MCMoves(1)>0.) then
  write(*,*) 'You are not allowed to have rotation moves if nothing has rotation states.'
  write(*,*) 'Readjusting the move set from: ', MCMoves(1)
  MCMoves(1)=0.
  write(*,*) 'To: ', MCMoves(1)
endif

! lookup tables
ALLOCATE ( BeadToProteinIndex(1:BeadBounds(NkN+1)) )
ALLOCATE ( BeadToProteinType(1:BeadBounds(NkN+1)) )
ALLOCATE ( BeadToMSI(1:BeadBounds(NkN+1)) )
ALLOCATE ( ProteinToProteinType(1:smNk(NkN+1)) )
do i1=1,NkN
  BeadToProteinType(1 + BeadBounds(i1):BeadBounds(i1+1)) = i1
  ProteinToProteinType( 1 + smNk(i1):smNk(i1+1)) = i1
  do i2=1,mNk(i1)
    BeadToProteinIndex(1+(i2-1)*Vk(i1) + BeadBounds(i1) : i2*Vk(i1) + BeadBounds(i1)) = smNk(i1) + i2
    BeadToMSI(1+(i2-1)*Vk(i1) + BeadBounds(i1) : i2*Vk(i1) + BeadBounds(i1)) = [( i3,i3=sVk(i1)+1,sVk(i1+1) )]
  enddo
enddo

!! SlitherI Key words  (remove some from one end and add them to the other end) (Must be a homopolymer)
if (MCMoves(4)==0.) then
  UseSlitherI=.false.
else
  UseSlitherI=.true.
  Allocate ( UseSlitherI_List(1:NkN) )
  UseSlitherI_List=0
  UseSlitherI_ListN=0
  do i1=merge(2,1,HasDNA),NkN !check each protein type individually if it can slither 1, !!! DNA
    if (Vk(i1)==1) cycle !must be a polymer, single beads cant slither
    MS=sVk(i1)
    if (any(Connection(1,MS+2:MS+Vk(i1)) /= [(j, j=MS+1, MS+Vk(i1)-1)])) cycle !linear and correctly ordered
    if (any(LL(1,MS+2:MS+Vk(i1)) /= LL(1,MS+1))) cycle ! not all LLs are the same
    if (any(BeadType(MS+2:MS+Vk(i1)) /= BeadType(MS+1))) cycle ! not all bead interaction types are the same
!    if (Energy(BeadType(MS+1),BeadType(MS+1))/=0) cycle ! !Checks that there aren't in-cis interactions This is a requirement because I was lazy in coding the slither move for how it deals with reasigning the interactions after a move
    UseSlitherI_ListN=UseSlitherI_ListN+1
    UseSlitherI_List(UseSlitherI_ListN)=i1
  end do !finished checking how many slither proteins there are
  if (UseSlitherI_ListN==0) then
    UseSlitherI=.false.
    write(*,*) 'SlithermovesI arent appropriate for this system, turning them off'
    write(*,*) 'These moves are appropriate only for homopolymers'
    write(*,*) 'Move Probability old:',MCMoves(4)
    MCMoves(4)=0.
    write(*,*) 'Move Probability new:',MCMoves(4)
  endif
end if

!! SlitherII Key Words
if (MCMoves(5)==0) then
  UseSlitherII=.false.
else
  UseSlitherII=.true.
  Allocate ( UseSlitherII_List(1:NkN) )
  UseSlitherII_List=0
  UseSlitherII_ListN=0
  bead: do i1=merge(2,1,HasDNA),NkN !check each protein individually  !! DNA
    if (Vk(i1)==1) cycle !must be a polymer
    if (any(LL(1,sVk(i1)+2:sVk(i1+1))/=LL(1,sVk(i1)+1))) cycle ! all LLs must be the same
    do i2=sVk(i1)+2,sVk(i1+1)
      if (Connection(1,i2)/=i2-1) cycle bead
    enddo
    UseSlitherII_ListN=UseSlitherII_ListN+1
    UseSlitherII_List(UseSlitherII_ListN)=i1
  enddo bead
  if (UseSlitherII_ListN==0) then
    UseSlitherII=.false.
    write(*,*) 'SlithermovesII arent appropriate for this system, turning them off'
    write(*,*) 'Move Probability old:',MCMoves(5)
    MCMoves(5)=0
    write(*,*) 'Move Probability new:',MCMoves(5)
  endif
end if

if (UseGrand) then
  write(*,*) 'cant remember why I started this block'
endif

if (UseSlitherI) then
  i1=maxval(Vk(UseSlitherI_List(1:UseSlitherI_ListN)))
else
  i1=0
endif
if (UseSlitherII) then
  i2=maxval(Vk(UseSlitherII_List(1:UseSlitherII_ListN)))
else
  i2=0
endif
if (UseGrand) then
  i3=maxval(Vk(UseGrand_List(1:UseGrand_ListN)))
else
  i3=0
endif

if ( UseSlitherI .or. UseSlitherII ) then
  Allocate ( SlitherList(3,1:max(i1,i2)) )
endif
if ( UseSlitherII .or. UseGrand ) then
  Allocate ( SlitherRot(1:max(i2,i3)) )
endif



if (Trelax .ne. 0) then 
  write(*,*) 'You have set up time for nucleation as',Trelax
  write(*,*) 'Cluster Move I&II and grand move probability are temporaroly set as 0'
endif


if (StartDroplet) then
  if (sum(StartDropletNk*Vk)>3*StartDropletR**3) then
    write(*,*) 'not enough room for all the proteins you want in the dense phase!'
    write(*,*) 'It isnt allowed to have a starting density over 71 percent (9/4/pi)'
    stop
  endif
  do i1=1,NkN
  if (StartDropletNk(i1)>Nk(i1)) then
    write(*,*) 'You asked for more inside the starting droplet than there is material for component ',i1
    write(*,*) 'Decreasing the amount so the starting droplet has all of component ',i1,' in the droplet'
    StartDropletNk(i1)=Nk(i1)
  endif
  enddo
endif

if (any(Nk>mNk)) then
  write(*,*) 'Not enough room to fit all the proteins youve asked for in Nk (MaxNk smaller?)!'
  write(*,*) 'Crashing'
  stop
endif


if (StartSlab .and. StartImported) then
  write(*,*) 'You can not both start imported and slab, turn off one of them'
  stop
endif

if (StartSlab) then
!  if (BoxSize(3)<2*BoxSize(1) .or. BoxSize(3)<2*BoxSize(2)) then
!    write(*,*) 'Your size of slab direction is not long enough, make it longer'
!   stop
!  endif

  if (sum(StartSlabNk*Vk)>0.7*BoxSize(1)*BoxSize(2) * StartSlabH) then
    write(*,*) 'not enough room for all the proteins you want in the slab dense phase!'
    write(*,*) 'It isnt allowed to have a starting density over 70 percent'
    stop
  endif
  do i1=1,NkN
  if (StartSlabNk(i1)>Nk(i1)) then
    write(*,*) 'You asked for more inside the starting slab than there is material for component ',i1
    write(*,*) 'Decreasing the amount so the starting slab has all of component ',i1,' in the droplet'
    StartSlabNk(i1)=Nk(i1)
    stop
  endif
  enddo
else
  if(BoxSize(1)/=BoxSize(2) .or. BoxSize(1)/=BoxSize(3)) then
    write(*,*) 'Your box is not a cube, you better be starting Slab simulation.'
  endif
endif

!allocate stuff
ALLOCATE ( XProt(1:3,1:maxval(Vk(merge(2,1,HasDNA):NkN))) )
ALLOCATE ( Xall(1:maxval(ConnectionN)*2,1:3) )
ALLOCATE ( Lall(1:maxval(ConnectionN)*2) )
ALLOCATE ( Neighbor(1:maxval(ConnectionN)) )
ALLOCATE ( Lattice(1:BoxSize(1),1:BoxSize(2),1:BoxSize(3)) )
ALLOCATE ( List(1:3,1:BeadBounds(NkN+1)) )
ALLOCATE ( BoundTo(1:BeadBounds(NkN+1)) )
ALLOCATE ( Analysis_log(1:smNk(NkN+1)) )
ALLOCATE ( AllHostN(1:smNk(NkN+1),1:2) )
ALLOCATE ( AllHost(1:smNk(NkN+1)) )
if (HasDNA) then
  ALLOCATE ( DNAHost(1:smNk(NkN+1)) )
endif
ALLOCATE ( FreeHost(1:smNk(NkN+1)) )
ALLOCATE ( Fisher(1:smNk(NkN+1)) )
ALLOCATE ( IsHost(1:smNk(NkN+1)) )


tclock_sanity=0
Lattice=0
List=0
BoundTo=0
!neighbor list, optimized for cache

Propose=0
Rej=0
RotMoves=0
ETC=0
end subroutine BasicSetup
!!!!
!!!!
!!!!
!!!!
subroutine InitialConditions !!!!!! This is massively memory inefficienct.  I should un-build all invalid proteins instead of using temp lattices and lists
INTEGER :: r2,r1a
INTEGER, DIMENSION(:,:,:), Allocatable :: Lattice_temp
INTEGER, DIMENSION(:,:), Allocatable :: List_temp
INTEGER, DIMENSION(:), Allocatable :: Vec_NkN, Vec_ConN


ALLOCATE ( Lattice_temp(1:BoxSize(1),1:BoxSize(2),1:BoxSize(3)) )
ALLOCATE ( List_temp(1:3,1:BeadBounds(NkN)) )
r2=0


if (StartImported) then
  write(*,*) 'Importing Initial Conditions'
  write(*,*) 'Be careful!'
  ALLOCATE ( Vec_NkN(1:NkN) )
  ALLOCATE ( Vec_ConN(1:maxval(ConnectionN)) ) 
  
  if (UseGrand) then
    write(*,*) 'I need to debug the grand stuff with imported initial conditions, crashing!'
    stop
  endif
  
  open(unit=1,file='StartParms.txt')
  read(1,*) x1
  if (any(x1 /= BoxSize)) then
    write(*,*) 'New BoxSize (',BoxSize,') doesnt match the old BoxSize (',x1,')'
    stop
  endif
  read(1,*) i1
  if (i1.ne.NkN) then
    write(*,*) 'New NkN (',NkN,') doesnt match the old NkN (',i1,')'
    stop
  endif
  read(1,*) Vec_NkN
  if (UseGrand) then
    Nk=Vec_NkN
  elseif (any(Vec_NkN.ne.Nk)) then
    write(*,*) 'New Nk (',Nk,') doesnt match the old Nk (',Vec_NkN,')'
    stop
  endif
  read(1,*) Vec_NkN
  if (any(Vec_NkN.ne.Vk)) then
    write(*,*) 'New Vk (',Vk,') doesnt match the old Vk (',Vec_NkN,')'
    stop
  endif

  do i1=1,sVk(NkN+1)
    read(1,*) i2
    if (i2.ne.ConnectionN(i1)) then
      write(*,*) 'New ConnectionN (',ConnectionN(i1),') doesnt match the old ConnectionN (',i2,')'
      stop
    endif
    read(1,*) Vec_ConN(1:ConnectionN(i1))
    if (any(Vec_ConN(1:ConnectionN(i1)).ne.Connection(1:ConnectionN(i1),i1))) then
      write(*,*) 'New Connection (',Connection(1:ConnectionN(i1),i1),&
        ') doesnt match the old Connection (',Vec_ConN(1:ConnectionN(i1)),')'
      stop
    endif
    read(1,*) Vec_ConN(1:ConnectionN(i1))
    if (any(Vec_ConN(1:ConnectionN(i1)).ne.LL(1:ConnectionN(i1),i1))) then
      write(*,*) 'New LL (',LL(1:ConnectionN(i1),i1), &
        ') doesnt match the old LL (',Vec_ConN(1:ConnectionN(i1)),')'
      stop
    endif
  enddo
  close(1)
  
  open(unit=1,file='Start_ListForm.txt')
  open(unit=2,file='Start_BoundTo.txt')
  do i1=1,NkN
  do i2=1+BeadBounds(i1),Nk(i1)*Vk(i1)+BeadBounds(i1)
    read(2,*) BoundTo(i2)
    read(1,*) List(:,i2)
    Lattice(List(1,i2),List(2,i2),List(3,i2))=i2
    if (BoundTo(MS)>MS) then
      MS2=BoundTo(MS)
      MST=BeadType(i3+sVk(i1))
      MST2=BeadType(BeadToMSI(MS2))
      ETC=ETC+Energy(MST,MST2)
      if (UseESelf) then
      if (BeadToProteinIndex(MS2)==i2 + smNk(i1)) then
        ETC=ETC+floor(ESelf*exp(real(-(abs(MS-MS2)-1)/LSteric)))
      endif
      endif
    endif
  enddo
  enddo
  close(1)
  close(2)
  write(*,*) 'finished importing?'

else 
  if (StartSlab) then
    write(*,*) 'Inventing initial conditions (slab), assuming the third dimension is slab direction'
    do i1=1,NkN
    do i2=1,Nk(i1)
      a1: do r1a=1,1000  !number of attempts to place a specific protein
        if (r1a==1000) then 
          if (i2<=StartSlabNk(i1)) then
            write(*,*) 'Crashed on forming initial conditions (slab).'
          else
            write(*,*) 'Crashed on forming initial conditions (dilute).'
          endif
          write(*,*) 'Consider adding more room or trying with a different seed'
          write(*,*) 'Added ',i2+sNk(i1),' proteins out of ',sNk(NkN+1)
          stop
        endif
        
        x1=randiBS(BoxSize)
        if (i2<=StartSlabNk(i1)) then !Put it into the dense phase
          x1(3)=randi(StartSlabH)+(BoxSize(3)-StartSlabH)/2
        endif
        
        XProt=AddProtein_robust(i1,x1)!Protein coordinates
        if (XProt(1,1)==0) then
          write(*,*) 'Couldnt place a protein without having a steric clash with itself'
          write(*,*) 'Crashing on bead number',XProt(2,1)
          stop
        endif
        
        do i3=1,Vk(i1)
          XProt(:,i3) = PB3(XProt(:,i3))
          if (Lattice(XProt(1,i3),XProt(2,i3),XProt(3,i3)) /= 0) then
            cycle a1
          endif
        enddo
        
        List(:,1+(i2-1)*Vk(i1)+BeadBounds(i1):(i2)*Vk(i1)+BeadBounds(i1)) = XProt(:,1:Vk(i1))
        do i3=1,Vk(i1)
          Lattice(XProt(1,i3),XProt(2,i3),XProt(3,i3))=i3+(i2-1)*Vk(i1)+BeadBounds(i1)
        enddo
        exit a1
      enddo a1
    enddo
    enddo
    
  elseif (HasDNA) then  !!! DNA (start)
    write(*,*) 'Inventing initial conditions (dilute + DNA)'
    do i3=1,Vk(1) ! adding DNA molecule
      Lattice(BoxSize(1)/2,BoxSize(2)/2,i3)=i3
      List(:,i3)=[BoxSize(1)/2,BoxSize(2)/2,i3]
    enddo
    do i1=2,NkN
    do i2=1,Nk(i1)
      a2: do r1a=1,1000 !number of attempts to place a specific protein
      
        if (r1a==1000) then
          write(*,*) 'Crashed on forming initial conditions (dilute + DNA).'
          write(*,*) 'Consider adding more room or trying with a different seed'
          write(*,*) 'Added ',i2+sNk(i1),' proteins out of ',sNk(NkN+1)
          stop
        endif
        
        x1=randiBS(BoxSize)
        XProt=AddProtein_robust(i1,x1)!Protein coordinates
        if (XProt(1,1)==0) then
          write(*,*) 'Couldnt place a protein without having a steric clash with itself'
          write(*,*) 'Crashing on bead number',XProt(2,1)
          stop
        endif
        do i3=1,Vk(i1)
          XProt(:,i3) = PB3(XProt(:,i3))
          if (Lattice(XProt(1,i3),XProt(2,i3),XProt(3,i3)) /= 0) then
            cycle a2 ! protein overl
          endif
        enddo
        
        List(:,1+(i2-1)*Vk(i1)+BeadBounds(i1):(i2)*Vk(i1)+BeadBounds(i1)) = XProt(:,1:Vk(i1))
        do i3=1,Vk(i1)
          Lattice(XProt(1,i3),XProt(2,i3),XProt(3,i3))=i3+(i2-1)*Vk(i1)+BeadBounds(i1)
        enddo
        exit a2
      enddo a2
    enddo
    enddo !!! DNA  (end)
      
  else
    if (all(StartDropletNk==0)) then
      write(*,*) 'Inventing initial conditions (dilute)'
    else
      write(*,*) 'Inventing initial conditions (dense)'
    endif
    do i1=1,NkN
    do i2=1,Nk(i1)
      a3: do r1a=1,1000 !number of attempts to place a specific protein
      
        if (r1a==1000) then
          if (i2<=StartDropletNk(i1)) then
            write(*,*) 'Crashed on forming initial conditions (dense).'
          else
            write(*,*) 'Crashed on forming initial conditions (dilute).'
          endif
          write(*,*) 'Consider adding more room or trying with a different seed'
          write(*,*) 'Added ',i2+sNk(i1),' proteins out of ',sNk(NkN+1)
          write(*,*) ''
          stop
        endif
        
        if (i2<=StartDropletNk(i1)) then !Put it into the dense phase
          call random_number(rr3)
          theta=rr3(1)*2.0*pi
          phi = acos(rr3(2)*2.0 - 1.0)
          rr=rr3(3)**(1.0/3.0)*StartDropletR
          x1=Floor([rr*sin(phi)*cos(theta),rr*sin(phi)*sin(theta),rr*cos(phi)]+real(BoxSize)/2.0+.5)
        else ! else puts it in the bulk phase
          x1=randiBS(BoxSize)
        endif
        
        XProt=AddProtein_robust(i1,x1)!Protein coordinates
        if (XProt(1,1)==0) then
          write(*,*) 'Couldnt place a protein without having a steric clash with itself'
          write(*,*) 'Crashing on bead number',XProt(2,1)
          stop
        endif
        do i3=1,Vk(i1)
          XProt(:,i3) = PB3(XProt(:,i3))
          if (Lattice(XProt(1,i3),XProt(2,i3),XProt(3,i3)) /= 0) then
            cycle a3
          endif
        enddo
        
        List(:,1+(i2-1)*Vk(i1)+BeadBounds(i1):(i2)*Vk(i1)+BeadBounds(i1)) = XProt(:,1:Vk(i1))
        do i3=1,Vk(i1)
          Lattice(XProt(1,i3),XProt(2,i3),XProt(3,i3))=i3+(i2-1)*Vk(i1)+BeadBounds(i1)
        enddo
        exit a3
      enddo a3
    enddo
    enddo
  endif
  
  write(*,*) 'Forming random bonds to stabilize dense phase'
  do i1=1,NkN
  do i2=1,Nk(i1)
  do i3=1,Vk(i1)
    MS = i3 + (i2-1)*Vk(i1) + BeadBounds(i1)
    MST= BeadType(i3+sVk(i1))
    if (BoundTo(MS)==0) then
      ForwardRMN=1
      do i4=1,26 !lets look around to see who it could interacting with
        R=PB3(List(:,MS)+NL(:,i4))
        MS2=Lattice(R(1),R(2),R(3))
        if (MS2/=0) then
          MST2=BeadType(BeadToMSI(MS2))
          if (Energy(MST,MST2)/=0  & ! is sticky
                .and. BoundTo(MS2)==0) then ! is single and ready to mingle
            ForwardRMN=ForwardRMN+1
            RotMoves(ForwardRMN)=MS2
          endif
        endif
      enddo
      if (ForwardRMN>1) then
        RotMove=RotMoves(randi(ForwardRMN))
        if (RotMove > 0) then
          BoundTo(MS)=RotMove
          BoundTo(RotMove)=MS
          MST2=BeadType(BeadToMSI(RotMove))
          ETC=ETC+Energy(MST,MST2)
          if (UseESelf) then
          if (i2 + smNk(i1)==BeadToProteinIndex(RotMove)) then
            ETC=ETC+floor(ESelf*exp(real(-(abs(MS-RotMove)-1)/LSteric)))
          endif
          endif
        endif
      endif
    endif
  enddo
  enddo
  enddo
endif !end of choice between compact and dilute

if (UseSpringLinkerEQ) then
do i1=1,NkN
do i2=1,Nk(i1)
do i3=1,Vk(i1)
  MS =i3 + (i2-1)*Vk(i1) + BeadBounds(i1)
  MSI=i3 + sVk(i1)
  do i4=1, SLconN(MSI)
  if (SLcon(i4,MSI)>MSI) then
    x1=List(:,MS)
    x2=List(:,MS-MSI+SLcon(i4,MSI))
    dx=DistPB(x1,x2)
    dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI))!3d distance differnce to Eq Len for MS1
    ETC=ETC + SLPot(i4,MSI)*dx_max**2
  endif
  enddo
enddo
enddo
enddo
endif

write(*,*) 'Starting Energy:', ETC

! Open all the export files
if (E_Lattice/=0) then
!  open(unit=8,file='AListForm.r4',form='unformatted',access='direct',recl=4*3*sum(Vk*Nk))
!  open(unit=8,file='ABoundTo.r4',form='unformatted',access='direct',recl=4*sum(Vk*Nk))
  open(unit=1,file='A_ListForm.xyz')
  if (WriteBoundTo) open(unit=2,file='A_BoundTo.txt')
  
  if (.not.UseGrand) then
    open(4,file='A_Bonds.vmd')
    write(4,'(A)')'topo clearbonds'
    i4=0
    do i1=1,NkN
      write(4,'(A,i6)') 'set NkN', i1
      write(4,'(A,i6)') 'set Vk', Vk(i1)
      write(4,'(A,i6)') 'set end', BeadBounds(i1+1)
      do i2=sVk(i1)+1,sVk(i1+1)
        do i3=1,ConnectionN(i2)
          if(Connection(i3,i2)>i2) then
            write(4,'(A,i6)') 'set current', BeadBounds(i1) +i2-1
            write(4,'(A,i6)') 'set partner', Connection(i3, i2)-i2 !defined as distance to current
            write(4,*) 'for {set j $current} {$j < $end} {incr j $Vk} {'
            write(4,*) '         topo addbond $j [expr {$j + $partner}]'
            write(4,*) '}'
          endif
        enddo 
      enddo
    enddo
    close(4)
  endif
endif
if (E_Nk/=0) then
  open(unit=3,file='A_Nk.txt')
endif

if (E_SelfLoop/=0) then
  do i1=1, NkN! assuming only one linker
    if (Vak(i1)==2) then
      LinkerType=i1
    endif
  enddo
  !If linker left == linker right: SelfLoopN+=1
  write(*,*) 'Self loop analysis is on, assuming only one Linker which is protein type', LinkerType
  write(*,*) 'There are', Nk(LinkerType), 'Linkers'
  write(*,*) 'Be very careful that your Linker are only allowed to have two interacting beads at both ends!'
  open(unit=97,file='A_SelfLoop.txt')
endif


if (E_Network/=0) then
  do i1=1, NkN! assuming only one linker
    if (Vak(i1)==2) then
      LinkerType=i1
    endif
  enddo
  write(*,*) 'Network loop analysis is on, assuming only one Linker which is protein type', LinkerType
  write(*,*) 'There are', Nk(LinkerType), 'Linkers'
  open(unit=98,file='A_LinkerConnection.txt')
  ALLOCATE ( LinkerCon(1:Nk(LinkerType),3) )! first one for the linker index, then each linker could have two connections
  LinkerCon=0
endif


if(UseSpringLinkerEq) then
if (E_SpringLinkerVector/=0) then
  open(unit=99,file='A_SpringLinkerVector.txt')
  ALLOCATE ( SLVec(1:7) )
  SLVec=0
endif
endif

if (E_Energy/=0) then
  open(unit=101,file='A_Energy.txt')
endif
if (E_LargestCluster/=0) then
  open(unit=102,file='A_LargestCluster.txt')
  if (HasDNA) then
    open(unit=110,file='A_DNA_Cluster')
  endif
  if (WriteCoreCluster) then  !(NkN+1)*4+NkN
    ALLOCATE ( LargestCluster(1:5*(NkN+1)) )
  else
    ALLOCATE ( LargestCluster(1:3*(NkN+1)) )
  endif
  LargestCluster=0
endif
if (E_RG_All/=0) then
  open(unit=103,file='A_RG_All.txt')
  ALLOCATE ( RG_All(1:4) )
  RG_All=0
endif
if (E_RG_Cluster/=0) then
  open(unit=107,file='A_RG_Cluster.txt')
  ALLOCATE ( RG_Cluster(1:8) )
  RG_Cluster=0
  if(WriteCoreCluster) then
    open(unit=108,file='A_RG_Core.txt')
  ALLOCATE ( RG_Core(1:8) )
  RG_Core=0
    open(unit=109,file='A_RG_CoreHalf.txt')
  ALLOCATE ( RG_CoreHalf(1:8) )
  RG_CoreHalf=0
  endif
endif
if (E_ClusterHist/=0) then
  open(unit=104,file='A_ClusterHist.txt')
  ALLOCATE ( ClusterHist(1:21) )
  ClusterHist=0
endif
if (E_BoundSites/=0) then
  open(unit=105,file='A_BoundSites.txt')
  ALLOCATE ( BoundSites(1:sVk(NkN+1)) )
  BoundSites=0
endif
if (E_BondTypes/=0) then
  open(unit=106,file='A_BondTypes.txt')
  ALLOCATE ( BondTypes(1:sVk(NkN+1),1:sVk(NkN+1)) )
  BondTypes=0
endif

if (E_XYZDist_None/=0 .or. E_XYZDist_All/=0 .or. E_XYZDist_Cluster/=0) then
  Allocate (XDist(1:BoxSize(1),NkN))
  Allocate (YDist(1:BoxSize(2),NkN))
  Allocate (ZDist(1:BoxSize(3),NkN))
endif
  
if (E_XYZDist_None/=0) then
  do i1=1,NkN
    write(Str1,*) i1
    do i2=1,3
      Str3='A_XYZDist_None_Protein'//trim(adjustl(Str1))//'_'//CHAR(64+23+i2)//'.txt'
      open(unit=200+3*i1+i2,file=Str3)
    enddo
  enddo
endif

if (E_XYZDist_All/=0) then
  do i1=1,NkN
    write(Str1,*) i1
    do i2=1,3
      Str3='A_XYZDist_All_Protein'//trim(adjustl(Str1))//'_'//CHAR(64+23+i2)//'.txt'
      open(unit=300+3*i1+i2,file=Str3)
    enddo
  enddo
endif

if (E_XYZDist_Cluster/=0) then
  do i1=1,NkN
    write(Str1,*) i1
    do i2=1,3
      Str3='A_XYZDist_Cluster_Protein'//trim(adjustl(Str1))//'_'//CHAR(64+23+i2)//'.txt'
      open(unit=400+3*i1+i2,file=Str3)
    enddo
  enddo
  if(WriteIndiDist_Cluster) then
    do i1=1,NkN
    write(Str1,*) i1
    do i2=1,NkN
    write(Str2,*) i2
    do i3=1,3!XYZ axis
      Str3='A_XYZDist_Cluster_P'//trim(adjustl(Str1))//'CenteredOnP'//trim(adjustl(Str2))//'_'//CHAR(64+23+i3)//'.txt'
      open(unit=10000+ i1*NkN*3 + i2*3 + i3,file=Str3)
    enddo
    enddo
    enddo
    Allocate (XDist_ij(1:BoxSize(1),NkN,NkN))
    Allocate (YDist_ij(1:BoxSize(2),NkN,NkN))
    Allocate (ZDist_ij(1:BoxSize(3),NkN,NkN))
  endif
endif

if (E_RadDist_Cluster/=0 .or. E_RadDist_All/=0) then
  if (WriteIndiDist_Cluster) then
    Allocate ( RadDist_ij(1:floor(sqrt(3.)/2*maxval(BoxSize)+1),NkN,NkN) )
  endif
  Allocate ( RadDist(1:floor(sqrt(3.)/2*maxval(BoxSize)+1),NkN) )
endif
   

if (E_RadDist_Cluster/=0) then
  do i1=1,NkN
    write(Str1,*) i1
    if(WriteIndiDist_Cluster) then
    do i2=1,NkN
      write(Str2,*) i2
      Str3='A_ClusterRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOn'//trim(adjustl(Str2))//'.txt'
      open(unit=500 + i1*NkN + i2,file=Str3)
    enddo
    endif
    Str3='A_ClusterRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOnLargestCluster.txt'
    open(unit=500+i1,file=Str3)
  enddo
  if(WriteCoreCluster) then
    !Open Core Radial Hist
    do i1=1,NkN
      write(Str1,*) i1
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(Str2,*) i2
        Str3='A_CoreRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOn'//trim(adjustl(Str2))//'.txt'
        open(unit=20000+2*i1*NkN+i2*2+1,file=Str3)
      enddo
      endif
      Str3='A_CoreRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOnCore.txt'
      open(unit=20000+2*NkN*NkN+2*NkN+3,file=Str3)
    enddo
    !Open CoreHalf Radial Hist
    do i1=1,NkN
      write(Str1,*) i1
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(Str2,*) i2
        Str3='A_CoreHalfRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOn'//trim(adjustl(Str2))//'.txt'
        open(unit=20000+2*i1*NkN+i2*2+2,file=Str3)
      enddo
      endif
      Str3='A_CoreHalfRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOnCoreHalf.txt'
      open(unit=20000+2*NkN*NkN+2*NkN+4,file=Str3)
    enddo
  endif
endif

end subroutine InitialConditions
!!!!
!!!!
!!!!
!!!!
function randi(j)
INTEGER :: j,randi
call random_number(rr1)
!randi=mod(floor(rr1*j),j)+1
randi=mod(int(rr1*j),j)+sign(1,j)
end function randi
!
function randi3(j)
INTEGER :: j,randi3(3)
call random_number(rr3)
!randi3=mod(floor(rr3*j),j)+1
randi3=mod(int(rr3*j),j)+sign(1,j)
end function randi3
!
function randi33(j)
INTEGER :: j(3),randi33(3)
call random_number(rr3)
!randi3=mod(floor(rr3*j),j)+1
randi33=mod(int(rr3*j),j)+sign(1,j)
end function randi33
!
function randiBS(j)
INTEGER :: j(3),randiBS(3)
call random_number(rr3)
!randiBS=mod(floor(rr3*j),j)+1
randiBS=mod(int(rr3*j),j)+sign(1,j)
end function randiBS
!
function PB1(j1,i1)
INTEGER :: j1,i1,PB1!j1 is the coordinates, i1 is the dimension for PB
!PB1=mod(j1+BoxSize(i1)-1,BoxSize(i1))+1
PB1=modulo(j1-1,BoxSize(i1))+1
end function PB1
!
function PB3(j3)
INTEGER :: j3(3),PB3(3)
!PB3=mod(j3+BoxSize-1,BoxSize)+1
PB3=modulo(j3-1,BoxSize)+1
end function PB3
!
function DistPB(x,y)
INTEGER :: x(3),y(3),DistPB(3)
!DistPB=mod(x-y+BoxSize/2+BoxSize,BoxSize)-BoxSize/2
DistPB=modulo(x-y+BoxSize/2,BoxSize)-BoxSize/2
end function DistPB
!
function rDistPB(x,y)
REAL :: x(3),y(3),rDistPB(3)
!rDistPB=mod(x-y+BoxSize/2+BoxSize,real(BoxSize))-BoxSize/2
rDistPB=modulo(x-y+BoxSize/2,real(BoxSize))-BoxSize/2
end function rDistPB
!
!
!pure function ProteinToProteinType(j1, mNk, NkN)
!integer :: i, ProteinToProteinType
!integer, intent(in) :: NkN, j1
!integer, intent(in) :: mNk(NkN)
!ProteinToProteinType = -1
!do i=1,NkN
!if (j1 <= sum(mNk(1:i))) then
!  ProteinToProteinType=i
!  return
!endif
!enddo
!end function ProteinToProteinType
!
!
!pure function ProteinToBead0(j1, mNk, Vk, BeadBounds, NkN)
!integer :: i, ProteinToBead0
!integer, intent(in) :: NkN, j1
!integer, intent(in) :: mNk(NkN), Vk(NkN), BeadBounds(NkN+1)
!ProteinToBead0 = -1
!do i=1,NkN
!if (j1 <= sum(mNk(1:i))) then
!  ProteinToBead0= (j1-sum(mNk(1:i-1))-1) *Vk(i) + BeadBounds(i)
!  return
!endif
!enddo
!end function ProteinToBead0
!
!
!pure function BeadToProteinType(j1, BeadBounds, NkN)
!integer :: i, BeadToProteinType
!integer, intent(in) :: NkN, j1
!integer, intent(in) :: BeadBounds(NkN+1)
!BeadToProteinType = -1
!do i = 1, NkN
!if (j1 <= BeadBounds(i+1)) then
!  BeadToProteinType = i
!  return
!endif
!enddo
!end function BeadToProteinType
!
!
!pure function BeadToProteinIndex(j1, BeadBounds, NkN, Vk)  ! needed for self interactions
!integer :: i, BeadToProteinIndex
!integer, intent(in) :: NkN, j1
!integer, intent(in) :: BeadBounds(NkN+1), Vk(NkN)
!BeadToProteinIndex = -1
!do i = 1, NkN
!if (j1 <= BeadBounds(i+1)) then
!  BeadToProteinIndex = (j1 - BeadBounds(i) -1)/Vk(i)+1 + sum(mNk(1:i-1))
!  return
!endif
!enddo
!end function BeadToProteinIndex
!
!
!pure function BeadToBeadIndex(j1, BeadBounds, NkN, Vk)  ! needed for linker type calculation
!integer :: i, BeadToBeadIndex
!integer, intent(in) :: NkN, j1
!integer, intent(in) :: BeadBounds(NkN+1), Vk(NkN)
!BeadToBeadIndex = -1  ! error number
!do i = 1, NkN
!if (j1 <= BeadBounds(i+1)) then
!  BeadToBeadIndex = mod(j1 - BeadBounds(i)-1,Vk(i))+1 + sum(Vk(1:i-1))
!  return
!endif
!enddo
!end function BeadToBeadIndex
!
!
!pure function BeadToBeadType(j1, BeadBounds, NkN, Vk) ! needed for energy type calculation
!integer :: i, BeadToBeadType
!integer, intent(in) :: NkN, j1
!integer, intent(in) :: BeadBounds(NkN+1), Vk(NkN)
!BeadToBeadType = -1
!do i = 1, NkN
!if (j1 <= BeadBounds(i+1)) then
!  BeadToBeadType = beadtype(mod(j1 - BeadBounds(i)-1,Vk(i))+1 + sum(Vk(1:i-1)))
!  return
!endif
!enddo
!end function BeadToBeadType
!
!
function AddProtein_robust(j1,x) !which protein number, 860
INTEGER :: AddProtein_robust(1:3,1:maxval(Vk(merge(2,1,HasDNA):NkN))) !protein beads, coordintes
INTEGER, intent(in) :: j1, x(3)
INTEGER :: j3,j4
! first bead is special
AddProtein_robust(:,1) = x
MSI = 1 + sVk(j1)
! rest of beads are the same

bead: do j2=2,Vk(j1) !bead number
  MSI=MSI+1
  try: do j3=1,100 !Attempt a bunch until something doesn't overlap with self
    x1 = AddProtein_robust(:,Connection(1,MSI)-sVk(j1)) + randi3(2*LL(1,MSI)+1)-LL(1,MSI)-1
    do j4=1,j2-1  ! check all previous beads
    if (all(x1 == AddProtein_robust(:,j4))) then
      cycle try
    endif
    enddo
    AddProtein_robust(:,j2)=x1
    cycle bead
  enddo try
  AddProtein_robust(1:2,1)=[0, j2]
  return
enddo bead

end function AddProtein_robust
!
!

subroutine ChooseMove

if (MCMove==8) then
  call MoveGrand
elseif (sNk(NkN+1)==0) then
  return
elseif (MCMove==1) then
  call MoveRot
elseif (MCMove==2) then
  call MoveTransI
elseif (MCMove==3) then
  call MoveTransII
elseif (MCMove==4) then
  call MoveSlitherI
elseif (MCMove==5) then
  call MoveSlitherII
elseif (MCMove==6) then
  call MoveClusterI
elseif (MCMove==7) then
  call MoveClusterII
endif

end subroutine ChooseMove
!!!!
!!!!
!!!!
!!!!
subroutine MoveRot

Propose(3)=Propose(3)+1
MS=randi(sNkVk(NkN+1))
! correct the indexing if there are voids for grand moves
if (UseGrand) then
do i1=merge(2,1,HasDNA),NkN
if (MS>BeadBounds(i1)+Nk(i1)*Vk(i1)) then
  MS=MS+(mNk(i1)-Nk(i1))*Vk(i1)
else
  return
endif
enddo
endif
x1=List(:,MS)
MST=BeadType(BeadToMSI(MS))

EC=0
EM=0
! old position
MS2=BoundTo(MS)
if (MS2>0) then
  MST2=BeadType(BeadToMSI(MS2))
  EC=EC+Energy(MST,MST2)
  if (UseESelf) then
  if (BeadToProteinIndex(MS)==BeadToProteinIndex(MS2)) then
    EC=EC+floor(ESelf*exp(real(-(abs(MS-MS2)-1)/LSteric)))
  endif
  endif
endif
! new position
ForwardRMN=1
do i=1,26
  R=PB3(x1+NL(:,i))
  MS2=Lattice(R(1),R(2),R(3))
  if (MS2/=0) then
    MST2=BeadType(BeadToMSI(MS2))
    if ((Energy(MST,MST2)/=0) .and. &
        (BoundTo(MS2)==0 .or. BoundTo(MS2)==MS)) then
      ForwardRMN=ForwardRMN+1
      RotMoves(ForwardRMN)=MS2
    endif
  endif
enddo
! choose a new position
RotMove=RotMoves(randi(ForwardRMN))
if (RotMove/=0) then
  MST2=BeadType(BeadToMSI(RotMove))
  EM=EM+Energy(MST,MST2)
  if (UseESelf) then
  if (BeadToProteinIndex(MS)==BeadToProteinIndex(RotMove)) then
    EM=EM+floor(ESelf*exp(real(-(abs(MS-RotMove)-1)/LSteric)))
  endif
  endif
endif

! Accept or Reject the move
call random_number(rr1)
if (rr1<exp(real(EC-EM)/1000)) then
  ETC=ETC+EM-EC
  if (BoundTo(MS)>0) then
    BoundTo(BoundTo(MS))=0
  endif
  BoundTo(MS)=RotMove
  if (RotMove/=0) then
    BoundTo(RotMove)=MS
  endif
else
  Rej(3,1)=Rej(3,1)+1
endif

end subroutine MoveRot
!!!!
!!!!
!!!!
!!!!
subroutine MoveTransI


Propose(1)=Propose(1)+1
if (HasDNA) then  !! DNA
  MS=randi(sNkVk(NkN+1)-Vk(1)) +Vk(1)
else
  MS=randi(sNkVk(NkN+1))
endif
! correct the indexing if there are voids for grand moves
if (UseGrand) then
do i1=merge(2,1,HasDNA),NkN
if (MS>BeadBounds(i1)+Nk(i1)*Vk(i1)) then
  MS=MS+(mNk(i1)-Nk(i1))*Vk(i1)
else
  return
endif
enddo
endif


MSI=BeadToMSI(MS)
MCN=ConnectionN(MSI) !number of linker connections (this is ConnectionN(MST2) or ConnectionN(ListToProtein(MS,3)))
! choose new coordinate
if (MCN==0) then !if it is a monomer domain that isn't tethered
  Move=PB3(List(:,MS)+randi3(5)-3)   ! can move from -2 to 2
else
  x1=List(:,MS)
  do i1=1,MCN !run through all of it's fixed bonds
    x2=List(:,MS-MSI+Connection(i1,MSI))
    Xall(i1,:)=DistPB(x2,x1) !this is the list of all displacements
  end do
!  do i2=1,3
!    i3=maxval(Xall(1:MCN,i2)-LL(1:MCN,MSI)) !minimum displacement
!    i4=minval(Xall(1:MCN,i2)+LL(1:MCN,MSI)) !maximum displacement
!    Move(i2)=PB1(randi(i4-i3+1)+i3+List(i2,MS)-1,i2) !random final position
!  end do

  x1=maxval(Xall(1:MCN,:)-spread(LL(1:MCN,MSI),2,3),1) !minimum displacement
  x2=minval(Xall(1:MCN,:)+spread(LL(1:MCN,MSI),2,3),1) !maximum displacement
  Move = PB3(randi33( x2-x1+1 ) + x1 + List(:,MS) -1)
endif
if (Lattice(Move(1),Move(2),Move(3))/=0 .and. Lattice(Move(1),Move(2),Move(3))/=MS) then !checks for steric clash
  ! reject the move now if there is a steric clash
  Rej(1,2)=Rej(1,2)+1
  return
endif

! choose new rotation
EC=0
EM=0
BackRMN=1
ForwardRMN=1
MST=BeadType(BeadToMSI(MS))
if (CanRot(MST)) then !this is asking if it can interact with anything at all
  if (BoundTo(MS)>0) then !if it was already bound we need to calculate the energy of that bond
    !EC=EC+Energy(MST1,BeadToType(BoundTo(MS)))
    MS2=BoundTo(MS)
    MST2=BeadType(BeadToMSI(MS2))
    EC=Energy(MST,MST2)
    if (UseESelf) then
    if (BeadToProteinIndex(MS)==BeadToProteinIndex(MS2)) then
      EC=EC+floor(ESelf*exp(real(-(abs(MS-MS2)-1)/LSteric)))
    endif
    endif
  endif
  do i=1,26 
    !lets look around to see who it could have been interacting with to count the reverse move probability
    R=PB3(List(:,MS)+NL(:,i))
    MS2=Lattice(R(1),R(2),R(3))
    if (MS2/=0) then
      MST2=BeadType(BeadToMSI(MS2))
      if (Energy(MST,MST2)/=0 .and. &
          (BoundTo(MS2)==0 .or. BoundTo(MS2)==MS)) then
        BackRMN=BackRMN+1
      endif
    endif
    !lets look around to see who it can interact with to find the forward move options
    R=PB3(Move+NL(:,i))
    MS2=Lattice(R(1),R(2),R(3))
    if (MS2/=0 .and. MS2/=MS) then
      MST2=BeadType(BeadToMSI(MS2))
      if (Energy(MST,MST2)/=0 .and. &
          (BoundTo(MS2)==0 .or. BoundTo(MS2)==MS)) then
        ForwardRMN=ForwardRMN+1
        RotMoves(ForwardRMN)=MS2
      endif
    endif
  enddo
endif
RotMove=RotMoves(randi(ForwardRMN))
if (RotMove/=0) then
  MST2=BeadType(BeadToMSI(RotMove))
  EM=Energy(MST,MST2)
  if (UseESelf) then
  if (BeadToProteinIndex(MS)==BeadToProteinIndex(RotMove)) then
    EM=EM+floor(ESelf*exp(real(-(abs(MS-RotMove)-1)/LSteric)))
  endif
  endif
endif
if (UseSpringLinkerEq) then
do i4=1,SLConN(MSI)
  x1=List(:,MS)
  x2=List(:,MS-MSI+SLCon(i4,MSI))
  dx=DistPB(x1,x2)
  dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI))!3d distance differnce to Eq Len for MS1
  EC=EC + SLPot(i4,MSI)*dx_max**2
  dx=DistPB(Move,x2)
  dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI))
  EM=EM + SLPot(i4,MSI)*dx_max**2
enddo
endif

! Accept or Reject the move
call random_number(rr1)
if (rr1<exp(real(EC-EM)/1000)*ForwardRMN/BackRMN) then
  ETC=ETC+EM-EC
  Lattice(List(1,MS),List(2,MS),List(3,MS))=0
  Lattice(Move(1),Move(2),Move(3))=MS
  List(:,MS)=Move
  if (BoundTo(MS) /= RotMove) then
    if (BoundTo(MS)>0) then
      BoundTo(BoundTo(MS))=0
    endif
    BoundTo(MS)=RotMove
    if (RotMove/=0) then
      BoundTo(RotMove)=MS
    endif
  endif
else
  Rej(1,1)=Rej(1,1)+1
endif

end subroutine MoveTransI
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
subroutine MoveTransII  !! this section needs to have the functions ported in still

Propose(2)=Propose(2)+1
if (HasDNA) then  !! DNA
  MS=randi(sNkVk(NkN+1)-Vk(1))+Vk(1)
else
  MS=randi(sNkVk(NkN+1))
endif


! correct the indexing if there are voids for grand moves
if (UseGrand) then
do i1=merge(2,1,HasDNA),NkN
if (MS>BeadBounds(i1)+Nk(i1)*Vk(i1)) then
  MS=MS+(mNk(i1)-Nk(i1))*Vk(i1)
else
  return
endif
enddo
endif

MS2=BoundTo(MS)
if (HasDNA) then
if (MS2/=0 .and. MS2<=Vk(1)) then
  ! can't move the DNA
  Rej(2,2)=Rej(2,2)+1
  return
endif
endif

MSI=BeadToMSI(MS)
!MCN=ConnectionN(MSI1) !number of linker connections (this is ConnectionN(MST2) or ConnectionN(ListToProtein(MS,3)))

MCN=0 !Counter for the number of linkers that must be preserved
x1=List(:,MS)
! find the distances that are tethered to MS
do i1=1,ConnectionN(MSI) !run through all of it's fixed bonds
if (MS2-MS /= Connection(i1,MSI)-MSI) then ! self fixed bonds don't matter
  MCN=MCN+1
  x3=List(:,MS-MSI+Connection(i1,MSI))
  Xall(MCN,:)=DistPB(x3,x1) !this is the list of all displacements
  Lall(MCN)=LL(i1,MSI)
endif
enddo
! find the distances that are tethered to boundto(MS)
if (MS2/=0) then
  MSI2=BeadToMSI(MS2)
  x2=List(:,MS2)
  do i1=1,ConnectionN(MSI2) !run through all of it's fixed bonds
  if (MS-MS2 /= Connection(i1,MSI2)-MSI2) then ! self fixed bonds don't matter
    MCN=MCN+1
    x3=List(:,MS2-MSI2+Connection(i1,MSI2))
    Xall(MCN,:)=DistPB(x3,x2) !this is the list of all displacements
    Lall(MCN)=LL(i1,MSI2)
  endif
  enddo
endif
! find the proposed move
if (MCN==0) then !if there are no covalent bonds
  Move=randi3(5)-3! goes from -2 to 2
else
  do i2=1,3
    i3=maxval(Xall(1:MCN,i2)-Lall(1:MCN)) !minimum displacement
    i4=minval(Xall(1:MCN,i2)+Lall(1:MCN)) !maximum displacement
    Move(i2)=randi(i4-i3+1)+i3-1 !random final position
  end do
endif
x1=PB3(x1+Move)
if (Lattice(x1(1),x1(2),x1(3))/=0 .and. Lattice(x1(1),x1(2),x1(3))/=MS2) then !checks if MS clashes in new position
  Rej(2,2)=Rej(2,2)+1
  return !clash, can't proceed.
endif
if (MS2/=0) then
  x2=PB3(x2+Move)
  if (Lattice(x2(1),x2(2),x2(3))/=0 .and. Lattice(x2(1),x2(2),x2(3))/=MS) then !checks if boundto(MS) clashes in new position  
    Rej(2,2)=Rej(2,2)+1
    return
  endif
endif

if (UseSpringLinkerEq) then
  EM=0
  EC=0
  do i4=1,SLConN(MSI)
    dx1=List(:,MS)
    dx2=List(:,MS-MSI+SLCon(i4,MSI))
    dx=DistPB(dx1,dx2)
    dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI))!differnce to Eq Leng for MS1
    EC=EC + SLPot(i4,MSI)*dx_max**2
    dx=DistPB(x1,dx2)
    dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI))!differnce to Eq Leng for MS1
    EM=EM + SLPot(i4,MSI)*dx_max**2
  enddo
  if (MS2/=0) then
  do i4=1,SLConN(MSI2)
    dx1=List(:,MS2)
    dx2=List(:,MS2-MSI2+SLCon(i4,MSI2))
    dx=DistPB(dx1,dx2)
    dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI2))!differnce to Eq Leng for MS2
    EC=EC + SLPot(i4,MSI2)*dx_max**2
    dx=DistPB(x2,dx2)
    dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI2))!differnce to Eq Leng for MS2
    EM=EM + SLPot(i4,MSI2)*dx_max**2
  enddo
  endif
  call random_number(rr1)
  if (rr1<exp(real(EC-EM)/1000)) then
    ETC=ETC+EM-EC !applies the move below
  else
    Rej(2,1)=Rej(2,1)+1
    return ! move was rejected
  endif
endif

Lattice(List(1,MS),List(2,MS),List(3,MS))=0
if (MS2/=0) then
  Lattice(List(1,MS2),List(2,MS2),List(3,MS2))=0
endif
List(:,MS)=x1
Lattice(x1(1),x1(2),x1(3))=MS
if (MS2/=0) then
  List(:,MS2)=x2
  Lattice(x2(1),x2(2),x2(3))=MS2
endif

end subroutine MoveTransII
!!!!
!!!!
!!!!
!!!!
subroutine MoveSlitherI  ! homopolymer slither, moves exactly one bead

Propose(4)=Propose(4)+1

EM=0
EC=0
SlitherList=0

MP = randi(sum(Nk(UseSlitherI_List(1:UseSlitherI_ListN))))
do i1=1,UseSlitherI_ListN
if (MP<=sum(Nk(UseSlitherI_List(1:i1)))) then
  MPT=UseSlitherI_List(i1)
  MP = MP - sum(Nk(UseSlitherI_List(1:i1-1)))
  exit
endif
enddo

MS = (MP-1)*Vk(MPT) + BeadBounds(MPT)
MP = MP + smNk(MPT)
MSI = BeadToMSI(MS+1)-1

R=randi3(2*LL(1,MSI+1)+1)-LL(1,MSI+1)-1! goes from -LL to LL
call random_number(rr1)
if (rr1<.5) then !lets go forward!
  SlitherForward=1
  MS_back=MS+1  ! removed bead
  MS_front=MS+Vk(MPT)
  Move=PB3(List(:,MS_front)+R)
  SlitherList(:,1:Vk(MPT)-1)=List(:,MS+2:MS+Vk(MPT))
  SlitherList(:,Vk(MPT))=Move
else!lets go Backwards!
  SlitherForward=-1
  MS_back=MS+Vk(MPT)
  MS_front=MS+1
  Move=PB3(List(:,MS_front)+R)
  SlitherList(:,2:Vk(MPT))=List(:,MS+1:MS+Vk(MPT)-1)
  SlitherList(:,1)=Move
endif
MST=BeadType(BeadToMSI(MS_back))

if (Lattice(Move(1),Move(2),Move(3))/=0 .and. Lattice(Move(1),Move(2),Move(3))/=MS_back) then
  ! steric clash, rejecting now
  Rej(4,2)=Rej(4,2)+1
  return
endif

! calculate the acceptance probability
BackRMN=1
ForwardRMN=1
if (.not. CanRot(MST)) then
  RotMove=0
else !check if it can bind
  ! evaluate current state
  MS2=BoundTo(MS_back)
  if (MS2/=0) then ! analyzes what it was bound to at the back
    BackRMN=2
    MST2 = BeadType(BeadToMSI(MS2))
    EC=EC + Energy(MST,MST2)
    if (UseESelf) then
    if (MP == BeadToProteinIndex(MS2)) then
      EC=EC+floor(ESelf*exp(real(-(abs(MS_back-MS2)-1)/LSteric)))
    endif
    endif
  endif
  do i1=1,26 !analyzes what it could bind to in the front
    R=PB3(List(:,MS_back)+NL(:,i1))
    MS2=Lattice(R(1),R(2),R(3))
    if (MS2/=0) then
    if (BoundTo(MS2)==0) then
      MST2=BeadType(BeadToMSI(MS2))
      if (Energy(MST,MST2)/=0) then
        BackRMN=BackRMN+1
      endif
    endif
    endif
    
    R=PB3(Move+NL(:,i1))
    MS2=Lattice(R(1),R(2),R(3))
    if (MS2/=0) then
      if (Energy(MST,MST)/=0) then  !self interaction
      if (MS2==MS_front .or. (MS2>MS+1 .and. MS2<MS+Vk(MPT)) ) then
      if (BoundTo(MS2)==0) then
        ForwardRMN=ForwardRMN+1
        RotMoves(ForwardRMN)=MS2-SlitherForward
      endif
      endif
      endif
      if (MS2<=MS .or. MS2>MS+Vk(MPT)) then !interaction with not self
        MST2=BeadType(BeadToMSI(MS2))
        if (Energy(MST,MST2)/=0 &
            .and. (BoundTo(MS2)==0 .or. BoundTo(MS2)==MS_back)) then
          ForwardRMN=ForwardRMN+1
          RotMoves(ForwardRMN)=MS2
        endif
      endif
    endif
  enddo
  RotMove=RotMoves(randi(ForwardRMN))
  if (RotMove/=0) then
    MST2=BeadType(BeadToMSI(RotMove))
    EM=EM+Energy(MST,MST2)
    if (UseESelf) then
    if (MP == BeadToProteinIndex(RotMove)) then
      EM=EM+floor(ESelf*exp(real(-(abs(MS_front-RotMove)-1)/LSteric)))
    endif
    endif
  endif
endif

if (UseSpringLinkerEq) then
do i1=1,Vk(MPT)
  MSI=i1+sVk(MPT)
  do i2=1,SLConN(MSI)
    if (MSI<SLCon(MSI,i2)) then
      x1=List(:,MS+i1)
      x2=List(:,MS+SLCon(MSI,i2)-sVk(MPT))
      dx=DistPB(x1,x2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i2,MSI))!maximal differnce to Eq Leng 
      EC=EC + SLPot(i2,MSI)*dx_max**2
      x1=SlitherList(:,i1)
      x2=SlitherList(:,SLCon(i2,MSI)-sVk(MPT))
      dx=DistPB(x1,x2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i2,MSI))!maximal differnce to Eq Leng 
      EM=EM + SLPot(i2,MSI)*dx_max**2
    endif
  enddo
enddo
endif

!Accept or Reject the move
call random_number(rr1)
if (rr1<exp(real(EC-EM)/1000)*ForwardRMN/BackRMN) then
  ETC=ETC+EM-EC
  
  !clear out the end that is removed
  if (BoundTo(MS_back)>0) then
    BoundTo(BoundTo(MS_back))=0
  endif
  Lattice(List(1,MS_back),List(2,MS_back),List(3,MS_back))=0
  ! update the middle beads
  do i1=MS_back,MS_front-SlitherForward,SlitherForward
    BoundTo(i1)=BoundTo(i1+SlitherForward)
    if (BoundTo(i1)>0) then
      BoundTo(BoundTo(i1))=i1
    endif
    Lattice(List(1,i1+SlitherForward),List(2,i1+SlitherForward),List(3,i1+SlitherForward))=i1
    List(:,i1)=List(:,i1+SlitherForward)
  enddo
  BoundTo(MS_front)=RotMove
  if (BoundTo(MS_front)/=0) then
    BoundTo(BoundTo(MS_front))=MS_front
  endif
  List(:,MS_front)=Move
  Lattice(Move(1),Move(2),Move(3))=MS_front
else
  Rej(4,1)=Rej(4,1)+1
endif

end subroutine MoveSlitherI
!!!!
!!!!
!!!!
!!!!
subroutine MoveSlitherII

Propose(5)=Propose(5)+1

MP = randi(sum(Nk(UseSlitherII_List(1:UseSlitherII_ListN))))
do i1=1,UseSlitherII_ListN
if (MP<=sum(Nk(UseSlitherII_List(1:i1)))) then
  MPT=UseSlitherII_List(i1)
  MP = MP - sum(Nk(UseSlitherII_List(1:i1-1)))
  exit
endif
enddo

MS_back = (MP-1)*Vk(MPT) + BeadBounds(MPT)
MP = MP + smNk(MPT)
MSI = BeadToMSI(MS_back+1)-1

EM=0
EC=0
logForwardRMN=0
logBackRMN=0
SlitherRot=0
SlitherList=0
if (Vk(MPT)>5) then
  Slithertimes=randi(5)
else
  SlitherTimes=randi(Vk(MPT)-1)
endif

!find the potential move first
call random_number(rr1)
if (rr1<.5) then !lets go forward!
  SlitherForward=1
  do i1=1,Vk(MPT)!-SlitherTimes
    if (i1<=Vk(MPT)-SlitherTimes) then
      SlitherList(:,i1)=List(:,MS_back+SlitherTimes+i1)
    else
      R=randi3(2*LL(1,MSI+i1)+1)-LL(1,MSI+i1)-1 ! goes from -LL to LL
      x1=PB3(SlitherList(:,i1-1)+R)
      MS2=Lattice(x1(1),x1(2),x1(3)) !proposed move
      if (MS2/=0 .and. (MS2<=MS_back .or. MS2>MS_back+SlitherTimes) ) then
        rej(5,2)=rej(5,2)+1 !isn't open
        return
      endif
      do i2=Vk(MPT)-SlitherTimes+1,i1-1
        if (all(x1==SlitherList(:,i2))) then
          rej(5,2)=rej(5,2)+1 !overlaps with other new beads
          return
        endif
      enddo
      SlitherList(:,i1)=x1
    endif
  enddo
else  !lets go backwards!
  SlitherForward=-1
  do i1=Vk(MPT),1,-1
    if (i1>SlitherTimes) then
      SlitherList(:,i1)=List(:,MS_back-SlitherTimes+i1)
    else
      R=randi3(2*LL(1,MSI+i1)+1)-LL(1,MSI+i1)-1 ! goes from -LL to LL
      x1=PB3(SlitherList(:,i1+1)+R)
      MS2=Lattice(x1(1),x1(2),x1(3))
      if (MS2/=0 .and. (MS2<=MS_back+Vk(MPT)-SlitherTimes .or. MS2>MS_back+Vk(MPT)) ) then
        rej(5,2)=rej(5,2)+1 !isn't open
        return
      endif
      do i2=i1+1,SlitherTimes
        if (all(x1==SlitherList(:,i2))) then
          rej(5,2)=rej(5,2)+1 !overlaps with other new beads
          return
        endif
      enddo
      SlitherList(:,i1)=x1
    endif
  enddo
endif

! lets find all the interaction terms:
! we break every bond and reform every bond
do i1=1,Vk(MPT)!lets find all the interaction energies, this index is the segment number in the moving protein
MST=BeadType(MSI+i1)
if (CanRot(MST)) then !checks that the selected bead can interact with something
  i3=1 !index for 'could have been bound'
  i4=1 !index for 'could now bind to'
  MS=BoundTo(MS_back+i1)
  if (MS/=0 .and. (MS<=MS_back .or. MS>MS_back+i1)) then !if it is currently bound to a module and that module is (not on the same protein and a lower bead number) to prevent double counting
    MST2=BeadType(BeadToMSI(MS))
    EC=EC+Energy(MST,MST2)
    if (UseESelf) then
    if (MP == BeadToProteinIndex(BoundTo(MS_back+i1))) then
      EC=EC+floor(ESelf*exp(real(-(abs(MS_back+i1-MS)-1)/LSteric)))
    endif
    endif
  endif
!!! Maybe it's now reversed properly for detailed balance
! reverse probabilities - i3
  if (MS<=MS_back .or. MS>MS_back+i1) then ! if it is not 'preclaimed'
  do i2=1,26
    R=PB3(List(:,MS_back+i1)+NL(:,i2))  ! check the current number of rotations
    MS=Lattice(R(1),R(2),R(3))
    if (MS/=0) then ! is there a bead
      MST2=BeadType(BeadToMSI(MS))
      if (Energy(MST,MST2)/=0) then  ! is that bead sticky for our bead
        MS2=BoundTo(MS)
        if (MS<=MS_back .or. MS>MS_back+Vk(MPT)) then !checks our bead of interest is interacting with another chain
          if (MS2==0 .or. & !and is free
            (MS2>=MS_back+i1.and.MS2<=MS_back+Vk(MPT))) then !or bound to MP, and >current
              i3=i3+1
          endif
        elseif (MS>MS_back+i1 .and. &  !checks if it's on same chain and higher index
            (MS2==0 .or. & !and is free
            (MS2<=MS_back .or. MS2>=MS_back+i1))) then !or would have been unbound
              i3=i3+1
        endif
      endif
    endif
  enddo
  endif
  ! calculate the forward probability - i4 
  if (SlitherRot(i1)==0) then ! if it is not 'preclaimed'
    x1=SlitherList(:,i1)
    do i2=1,26  ! all binding sites from beads that will remain where they are
      R=PB3(x1+NL(:,i2)) ! Check all possible future rotations
      MS=Lattice(R(1),R(2),R(3))
      if (MS/=0 .and. (MS<=MS_back.or.MS>MS_back+Vk(MPT))) then ! not self
        MST2=BeadType(BeadToMSI(MS))
        MS2=BoundTo(Lattice(R(1),R(2),R(3))) !is bead of interest bound to something else?
        if (Energy(MST,MST2)/=0 .and. &  !is sticky
            (MS2==0 .or. & !and is free
            (MS2>MS_back.and.MS2<=MS_back+Vk(MPT)))) then !or bound to MP
          if (all(SlitherRot(1:i1-1)/=MS)) then !not bound to a new slither position
            i4=i4+1
            RotMoves(i4)=MS
          endif
        endif
      endif
    enddo
    do i2=i1+1,Vk(MPT) ! Checks to see if it will bind to it's own chain, but only chains of higher index
      MST2=BeadType(MSI+i2)
      if (Energy(MST,MST2)/=0) then ! is sticky for self
        x2=SlitherList(:,i2)
        if (SlitherRot(i2)==0 .and. all(abs(DistPB(x1,x2))<=1)) then !is free and a neighbor
          i4=i4+1
          RotMoves(i4)=MS_back+i2
        endif
      endif
    enddo
    SlitherRot(i1)=RotMoves(randi(i4))
    if (SlitherRot(i1)/=0) then !is bound to a bead
      if (SlitherRot(i1)>MS_back+i1 .and. SlitherRot(i1)<=MS_back+Vk(MPT)) then !is bound to self then we need to assign the other bead
        SlitherRot(SlitherRot(i1)-MS_back)=MS_back+i1
      elseif (SlitherRot(i1)>MS_back .and. SlitherRot(i1)<MS_back+i1) then !this should never happen, eventually remove it
        write(*,*) 'Slither is binding to something it should not be able to!'
        stop
      elseif (SlitherRot(i1)==MS_back+i1) then
        write(*,*) 'Slither is binding to its self!'
        write(*,*) t
        write(*,*) i1
        write(*,*) SlitherRot(1:i1)
        write(*,*)
        write(*,*) RotMoves
        write(*,*) i4
        stop
      endif
      MST2=BeadType(BeadToMSI(SlitherRot(i1)))
      EM=EM+Energy(MST,MST2)
      if (UseESelf) then
      if (MP == BeadToProteinIndex(SlitherRot(i1))) then
        EM=EM+floor(ESelf*exp(real(-(abs(MS_back+i1-SlitherRot(i1))-1)/LSteric)))
      endif
      endif
    endif
  endif
  !lets do a few detailed balance calculations
  logBackRMN=logBackRMN+log(real(i3))
  logForwardRMN=logForwardRMN+log(real(i4))
endif
enddo
!!!If use spring linker
if (UseSpringLinkerEq) then
do i1=1,Vk(MPT)
do i2=1,SLConN(MSI+i1)
if (SLCon(i2,MSI+i1)>MSI+i1) then
  x1=List(:,MS_back+i1)
  x2=List(:,MS_back-MSI+SLCon(i2,MSI+i1))
  dx=DistPB(x1,x2)
  dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i2,MSI+i1))
  EC=EC + SLPot(i2,MSI+i1)*dx_max**2
  x1=SlitherList(:,i1)
  x2=SlitherList(:,SLCon(i2,MSI+i1)-MSI)
  dx=DistPB(x1,x2)
  dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i2,MSI+i1))
  EM=EM + SLPot(i2,MSI+i1)*dx_max**2
endif
enddo
enddo
endif

!!! Accept or reject the move
call random_number(rr1)
if (rr1<exp(real(EC-EM)/1000 + logForwardRMN - logBackRMN)) then
  ETC=ETC+EM-EC
  do i1=1,Vk(MPT) !clear old bonds
    if (BoundTo(MS_back+i1)/=0) then
      BoundTo(BoundTo(MS_back+i1))=0
    endif
    Lattice(List(1,MS_back+i1),List(2,MS_back+i1),List(3,MS_back+i1))=0 !clear lattice
  enddo
  do i1=1,Vk(MPT) !build new bonds
    BoundTo(MS_back+i1)=SlitherRot(i1) ! set what I'm bound to
    if (SlitherRot(i1)/=0) then
      BoundTo(SlitherRot(i1))=MS_back+i1 !set what I'm bound to to be bound to me
    endif
    Lattice(SlitherList(1,i1),SlitherList(2,i1),SlitherList(3,i1))=MS_back+i1 !build new lattice
    List(:,MS_back+i1)=SlitherList(:,i1)
  enddo
else
  rej(5,1)=rej(5,1)+1
endif

end subroutine MoveSlitherII
!!!!
!!!!
!!!!
!!!!
subroutine MoveClusterI

Propose(6)=Propose(6)+1
if (HasDNA) then
  MP=randi(sNk(NkN+1)-1)+1
else
  MP=randi(sNK(NkN+1))
endif
do i1=1,NkN
if (MP<=sNk(i1+1)) then
  MPT=i1
  MP = MP - sNk(i1)+ smNk(i1)
  exit
endif
enddo

AllHostNN=1  !number of proteins in this cluster
AllHost(1)=MP !the proteins in this cluster
IsHost=.true.  !has the protein been seen by this cluster yet?
IsHost(MP)=.false.  ! it sees itself
i2=1 !current scan number
flag=.true. !flag for if move has failed from the cluster being too large for this type of move
do while (i2<=AllHostNN .and. flag)
  MPT=ProteinToProteinType(AllHost(i2))
  MS= (AllHost(i2) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
  do i3=MS+1,MS+Vk(MPT)  !move each segment in protein
    if (BoundTo(i3)/=0) then
      i4=BeadToProteinIndex(BoundTo(i3))

      if (HasDNA .and. i4==1) then
        Rej(7,2)=Rej(7,2)+1
        return
      endif

      if (Vk(ProteinToProteinType(i4))>100) then
        Rej(7,2)=Rej(7,2)+1
        return
      endif
      
      if (IsHost(i4) .eqv. .true.) then
        AllHostNN=AllHostNN+1
        if (AllHostNN>2) then
          Rej(6,2)=Rej(6,2)+1  ! too big
          return
        endif
        AllHost(AllHostNN)=i4
        IsHost(i4)=.false.
      endif
    endif
  enddo
  i2=i2+1 !start analyzing the next protein
enddo

Move=randiBS(BoxSize)-1
if (all(Move==0)) then !Does it move the protein a total distance of zero?  Lets just ignore this move
  Rej(6,2)=Rej(6,2)+1
  return
endif
do i2=1,AllHostNN  !run thru all proteins in host
  MPT=ProteinToProteinType(AllHost(i2))
  MS= (AllHost(i2) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
  do i3=MS+1,MS+Vk(MPT)  !move each segment in protein
    x1=PB3(List(:,i3)+Move) !new loc
    if (Lattice(x1(1),x1(2),x1(3))/=0) then !if new loc is occupied
      Rej(6,1)=Rej(6,1)+1
      return
    endif
  enddo
enddo  !end of examining if the move is acceptable

do i2=1,AllHostNN  !run thru all proteins in host and delete them from the lattice
  MPT=ProteinToProteinType(AllHost(i2))
  MS= (AllHost(i2) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
  do i3=MS+1,MS+Vk(MPT)  !move each segment in protein
    Lattice(List(1,i3),List(2,i3),List(3,i3))=0
  enddo
enddo
do i2=1,AllHostNN  !run thru all proteins in host and rebuild them on the lattice and list
  MPT=ProteinToProteinType(AllHost(i2))
  MS= (AllHost(i2) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
  do i3=MS+1,MS+Vk(MPT)  !move each segment in protein
    x1=PB3(List(:,i3)+Move)
    List(:,i3)=x1
    Lattice(x1(1),x1(2),x1(3))=i3
  enddo
enddo

end subroutine MoveClusterI
!!!!
!!!!
!!!!
!!!!
subroutine MoveClusterII

call NetworkAnalysis

! AllHostNN = number of clusters
! AllHostN = start and end of each cluster


i5=min(AllHostNN-1,10) !number of moves
Propose(7)=Propose(7)+i5

FisherN=AllHostNN-1 !number of available (set up fisher-yates algorithm)
j2=minval(maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1))) ! the largest to remove
!builds the fisher list of available clusters  !skips the largest cluster
Fisher(1:AllHostNN-1) = [(i, i = 1, j2 - 1), (i, i = j2 + 1, AllHostNN)]

ClusterMoveNumber: do i1=1,i5  !!!!! Does up to 10 network moves each time MoveCluster is called
  MPT=randi(FisherN) !random element from fisher list
  MP=Fisher(MPT)    !cluster number pulled from the fisher list
  Fisher(MPT)=Fisher(FisherN) !take the last element from the fisher list and move it to where it will be available in the future
  FisherN=FisherN-1 !shrink the population of fish

  Move=randiBS(BoxSize)
  if (all(Move==0)) then !Does it move the protein a total distance of zero?  Lets just ignore this move
    Rej(7,2)=Rej(7,2)+1
    cycle
  endif
  do i2=AllHostN(MP,1),AllHostN(MP,2)  !run thru all proteins in host
    MPT=ProteinToProteinType(AllHost(i2))
    MS= (AllHost(i2) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
    do i3=MS+1,MS+Vk(MPT)  !move each segment in protein
      x1=PB3(List(:,i3)+Move) !new loc
      if (Lattice(x1(1),x1(2),x1(3))/=0) then !if new loc is occupied
        Rej(7,1)=Rej(7,1)+1
        cycle ClusterMoveNumber
      endif
    enddo
  enddo  !end of examining if the move is acceptable
  do i2=AllHostN(MP,1),AllHostN(MP,2)  !run thru all proteins in host
    MPT=ProteinToProteinType(AllHost(i2))
    MS=(AllHost(i2)-smNk(MPT)-1)*Vk(MPT) + BeadBounds(MPT) 
    do i3=MS+1,MS+Vk(MPT)  !remove each segment in protein
      Lattice(List(1,i3),List(2,i3),List(3,i3))=0
    enddo
  enddo
  do i2=AllHostN(MP,1),AllHostN(MP,2)  !run thru all proteins in host
    MPT=ProteinToProteinType(AllHost(i2))
    MS=(AllHost(i2)-smNk(MPT)-1)*Vk(MPT) + BeadBounds(MPT) 
    do i3=MS+1,MS+Vk(MPT)  !move each segment in protein
      x1=PB3(List(:,i3)+Move)
      List(:,i3)=x1
      Lattice(x1(1),x1(2),X1(3))=i3
    enddo
  enddo
enddo ClusterMoveNumber

end subroutine MoveClusterII
!!!!
!!!!
!!!!
!!!!
subroutine MoveGrand

Propose(8)=Propose(8)+1


MPT=UseGrand_List(randi(UseGrand_ListN))
call random_number(rr3)
if (t==3989) then
  write(*,*) 'here34'
  write(*,*) rr3(1)
endif
if (rr3(1)<.5) then !Add a protein
  if (Nk(MPT)==mNk(MPT)) then ! must fit in the list form
    write(*,*) 'warning, Nk=mNk=',mNk,' and Im trying to add a protein.  t=',t,'breaks detailed balance'
    Rej(8,2)=Rej(8,2)+1
    return
  endif
  
  MP=1 + Nk(MPT) + smNk(MPT)
  XProt=AddProtein_Robust(MPT,randi33(BoxSize))
  
  if (XProt(1,1)==0) then !checks if the protein was placed w/ no self clashes  Note: This breaks detailed balance
    write(*,*) 'warning, protein placement has self-clash.  This breaks detailed balance.'
    Rej(8,2)=Rej(8,2)+1
    return
  endif
  
  do i1=1,Vk(MPT) ! check if the lattice is empty
    x1=PB3(XProt(:,i1))
    SlitherList(:,i1)=x1
    if (Lattice(x1(1),x1(2),x1(3))/=0) then
      Rej(8,2)=Rej(8,2)+1
      return
    endif
  enddo
  
  EM=0
  MS = Nk(MPT)*Vk(MPT) + BeadBounds(MPT) ! start of the new protein
  MSI= sVk(MPT)
  logForwardRMN=0
  SlitherRot=0
  do i1=1,Vk(MPT)
    ! check the entropy of chain placement
    if (i1>1) then
      ForwardRMN = (2*LL(1,MSI+i1)+1)**3-1  
      do i2=1,i1-2 ! remove all the double counting
      if (all( XProt(:,i2)>=XProt(:,i1-1)-LL(1,MSI+i1) .or. XProt(:,i2)<=XProt(:,i1-1)+LL(1,MSI+i1) ) ) then
        ForwardRMN=ForwardRMN-1
      endif
      enddo
      logForwardRMN = logForwardRMN + log(real(ForwardRMN))
    endif
    ! spring linker energy
    if (UseSpringLinkerEq) then
    do i4=1,SLConN(MSI+i1)
      if (SLCon(i4,MSI+i1)>MSI+i1) cycle 
      x1=XProt(:,i1)
      x2=XProt(:,SLCon(i4,MSI+i1)-MSI)
      dx=x1-x2 !! this doesn't need PB
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI+i1))!3d distance differnce to Eq Len for MS1
      EM=EM + SLPot(i4,MSI+i1)*dx_max**2
    enddo
    endif
    ! choose interactions
    MST= BeadType(MSI+i1)
    if (.not.CanRot(MST)) cycle
    RotMoves=0
    ForwardRMN=1
    do i2=1,26
      x1=PB3(SlitherList(:,i1)+NL(:,i2))  ! not self interactions
      MS2=Lattice(x1(1),x1(2),x1(3))
      if (MS2==0) cycle
      if (BoundTo(MS2)/=0) cycle    
      if (any(SlitherRot(1:i1-1)==MS2)) cycle
      MST2=BeadType(BeadToMSI(MS2))
      if (Energy(MST,MST2)/=0) then
        ForwardRMN=ForwardRMN+1
        RotMoves(ForwardRMN)=MS2
      endif
    enddo
    do i2=1,26
      x1=XProt(:,i1)+NL(:,i2)  ! self interactions
      do i3=1,i1-1
        if (SlitherRot(i3)==0) then
        if (all(x1==XProt(:,i3))) then
          ForwardRMN=ForwardRMN+1
          RotMoves(ForwardRMN)=MS+i3
        endif
        endif
      enddo
    enddo
    
    logForwardRMN=logForwardRMN+log(real(ForwardRMN))
    ! choose a rotation
    SlitherRot(i1)=RotMoves(randi(ForwardRMN))
    if (SlitherRot(i1)/=0) then
      MST2=BeadType(BeadToMSI(SlitherRot(i1)))
      EM=EM+Energy(MST,MST2)
      if (MP==BeadToProteinIndex(SlitherRot(i1))) then
        SlitherRot(SlitherRot(i1)-MS)=MS+i1
        if (UseESelf) then
          EM=EM+floor(ESelf*exp(real(-(abs(MS+i1-SlitherRot(i1))-1)/LSteric)))
        endif
      endif
    endif
  enddo
 
    
  ! accept or reject the move
  rr1=exp(real(ChemPot(MPT)-EM)/1000 + logForwardRMN)*product(BoxSize)/(Nk(MPT)+1) 
  if (rr3(2)<rr1) then
    Nk(MPT)=Nk(MPT)+1
    BoundTo(MS+1:MS+Vk(MPT))=SlitherRot(1:Vk(MPT))
    List(:,MS+1:MS+Vk(MPT))=SlitherList(:,1:Vk(MPT))
    do i1=1,Vk(MPT)
      x1=SlitherList(:,i1)
      if (SlitherRot(i1)/=0) BoundTo(SlitherRot(i1))=MS+i1
      Lattice(x1(1),x1(2),x1(3))=MS+i1
    enddo
    ETC=ETC+EM
    if (t==3989) then
      write(*,*) 'herklsdf'
      write(*,*) SlitherRot
    endif
    
  else
    Rej(8,1)=Rej(8,1)+1
  endif
else !remove a protein, then shift what ever protein is at the end to fill it's spot
  if (Nk(MPT)==0) then !can't go below Nk=0
    Rej(8,2)=Rej(8,2)+1
    return
  endif
  MP=randi(Nk(MPT))
  MS = (MP-1)*Vk(MPT) + BeadBounds(MPT)
  MP = MP + smNk(MPT)
  MSI= sVk(MPT)

  EC=0
  logBackRMN=0
  do i1=1,Vk(MPT)
    ! check the entropy of chain placement
    if (i1>1) then
      BackRMN = (2*LL(1,MSI+i1)+1)**3-1  
      do i2=1,i1-2 ! remove all the double counting
      if (all( List(:,MS+i2)>=List(:,MS+i1-1)-LL(1,MSI+i1) .or. List(:,MS+i2)<=List(:,MS+i1-1)+LL(1,MSI+i1) ) ) then
        BackRMN=BackRMN-1
      endif
      enddo
      logBackRMN = logBackRMN + log(real(BackRMN))
    endif
    ! spring linker energy
    if (UseSpringLinkerEq) then
    do i4=1,SLConN(MSI+i1)
      if (SLCon(i4,MSI+i1)>MSI+i1) cycle 
      x1=List(:,MS+i1)
      x2=List(:,MS+SLCon(i4,MSI+i1)-MSI)
      dx=DistPB(x1,x2) !! this DOES need PB
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI+i1))!3d distance differnce to Eq Len for MS1
      EC=EC + SLPot(i4,MSI+i1)*dx_max**2
    enddo
    endif
    ! choose interactions
    MST= BeadType(MSI+i1)
    if (.not.CanRot(MST)) cycle
!    RotMoves=0
    BackRMN=1
    SlitherRot(1:Vk(MPT))=[(MS + i2, i2 = 1, Vk(MPT))]
    do i2=1,26
      x1=PB3(List(:,MS+i1)+NL(:,i2))  ! is filled? 
      MS2=Lattice(x1(1),x1(2),x1(3))
      if (MS2==0) cycle !nothing to bind to
      if (MS2>MS+i1 .and. MS2<=MS+Vk(MPT)) cycle !it wouldn't have been placed yet
      if ( (BoundTo(MS2)/=0 .and. BoundTo(MS2)<MS+i1) .or. BoundTo(MS2)>MS+Vk(MPT) )  cycle !that bead would have been bound to something
      BackRMN=BackRMN+1
    enddo
    logBackRMN = logBackRMN + log(real(BackRMN))
    ! choose a rotation
    MS2=BoundTo(MS+i1)
    if (MS2/=0) then
      if (MP/=BeadToProteinIndex(MS2)) then ! different protein
        MST2=BeadType(BeadToMSI(MS2))
        EC=EC + Energy(MST,MST2)
      elseif (MS2<MS+i1) then
        MST2=BeadType(BeadToMSI(MS2))
        EC=EC + Energy(MST,MST2)
        if (UseESelf) then
          EC=EC+floor(ESelf*exp(real(-(abs(MS+i1-MS2)-1)/LSteric)))
        endif
      endif
    endif
  enddo
  
  rr1 = exp(real(EC-ChemPot(MPT))/1000 - logBackRMN)*(Nk(MPT))/product(BoxSize)
  if (rr3(2)<rr1) then
    do i1=MS+1,MS+Vk(MPT)  ! remove the protein
      Lattice(List(1,i1),List(2,i1),List(3,i1))=0
    enddo
    if (MP<Nk(MPT)) then !shift the protein at the end to fill what was deleted
      MS2=(Nk(MPT)-1)*Vk(MPT) + BeadBounds(MPT)
      do i1=1,Vk(MPT)  ! Run through all beads
        BoundTo(MS+i1)=BoundTo(MS2+i1)
        if (BoundTo(MS+i1)/=0) then
          BoundTo(BoundTo(MS+i1))=MS+i1
        endif
        List(:,MS+i1)=List(:,MS2+i1)
        Lattice(List(1,MS+i1),List(2,MS+i1),List(3,MS+i1))=MS+i1
      enddo
    endif
    Nk(MPT)=Nk(MPT)-1
    ETC=ETC - EC
  else
    Rej(8,1)=Rej(8,1)+1
  endif
endif

end subroutine MoveGrand
!!!!
!!!!
!!!!
!!!!
subroutine NetworkAnalysis
! constructs how all the proteins are networked
! AllHost gives a list of the proteins in the order they are bound to complexes
! AllHostN gives the list of the start and end position of each cluster
! AllHostNN gives the total list length of AllHostN

IsHost=.true.  !the name is legacy.  It's actually a logical for if this protein has been seen in the script yet
AllHostNN=0  !number of clusters
flag=.true.  !marks if it is the first cluster on the list because that has to be treated differently
! new code


if (HasDNA) then  ! separates the DNA cluster from the rest of the clusters
  FreeHostN=1  !number of polymers in this cluster
  FreeHost(1)=1  !list of polymers in this cluster, just added self to the list
  IsHost(1)=.false.  !prevents it from being picked again in its own cluster
  i3=1  !current number of proteins bound and counted
  do while (i3<=FreeHostN) !indexing list of proteins connected and counted
    MPT=ProteinToProteinType(FreeHost(i3))
    MS= (FreeHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
    do i4=MS+1,MS+Vk(MPT)  !segment in protein
      if (BoundTo(i3)/=0) then
        MP=BeadToProteinIndex(BoundTo(i3))
        if (IsHost(MP)) then !checks if the protein needs to be connected to this cluster or is already counted
          FreeHostN=FreeHostN+1  !add it to the cluster
          FreeHost(FreeHostN)=MP
          IsHost(MP)=.false.   !prevent it from being counted again
        endif
      endif
    enddo
    i3=i3+1
  enddo
  DNAHostN=FreeHostN
  DNAHost(1:DNAHostN)=FreeHost(1:DNAHostN)
endif

do i1=1,NkN
do i2=1 + smNk(i1) ,Nk(i1) + smNk(i1) 
if (IsHost(i2)) then ! if it is a host, collect up what it's bound to, otherwise it's in someone else's cluster
  AllHostNN=AllHostNN+1 !number of clusters, we just found another cluster
  FreeHostN=1  !number of polymers in this cluster
  FreeHost(1)=i2  !list of polymers in this cluster, just added self to the list
  IsHost(i2)=.false.  !prevents it from being picked again in its own cluster
  i3=1  !current number of proteins bound and counted
  do while (i3<=FreeHostN) !indexing list of proteins connected and counted
    MPT=ProteinToProteinType(FreeHost(i3))
    MS= (FreeHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
    do i4=MS+1,MS+Vk(MPT)  !segment in protein
    if (BoundTo(i4)/=0) then
      MP=BeadToProteinIndex(BoundTo(i4))
      if (IsHost(MP)) then !checks if the protein needs to be connected to this cluster or is already counted
        FreeHostN=FreeHostN+1  !add it to the cluster
        FreeHost(FreeHostN)=MP
        IsHost(MP)=.false.   !prevent it from being counted again
      endif
    endif
    enddo
    i3=i3+1
  enddo
  if (flag) then
    AllHostN(AllHostNN,:)=[1,FreeHostN]
    flag=.false.
  else
    AllHostN(AllHostNN,:)=[AllHostN(AllHostNN-1,2)+1,AllHostN(AllHostNN-1,2)+FreeHostN]
  endif
  AllHost(AllHostN(AllHostNN,1):AllHostN(AllHostNN,2))=FreeHost(1:FreeHostN)
endif
enddo
enddo

if (WriteCoreCluster) then
  Core = 0 !Core means fully occupied polymers
  CoreHalf = 0!CoreHalf means more than half occupied polymers
  CoreN = 0
  CoreHalfN = 0
  i1 = maxloc(AllHostN(1:AllHostNN, 2) - AllHostN(1:AllHostNN, 1), 1)
  do i2 = AllHostN(i1, 1), AllHostN(i1, 2) ! Runs through all proteins
    MPT = ProteintoProteinType(i2)
    MS = (i2-smNk(MPT)-1)*Vk(MPT) + BeadBounds(MPT)
    bonds = count( BoundTo(MS+1:MS+Vk(MPT))/=0 )
    if (bonds==Vak(MPT)) then
      CoreN = CoreN + 1
      CoreHalfN = CoreHalfN + 1
      Core(CoreN) = AllHost(i2)
      CoreHalf(CoreHalfN) = AllHost(i2)
    elseif (bonds>=Vak(MPT)/2) then
      CoreHalfN = CoreHalfN + 1
      CoreHalf(CoreHalfN) = AllHost(i2)
    endif
  enddo
endif

end subroutine NetworkAnalysis
!!!!
!!!!
!!!!
!!!!
subroutine Analysis

if (E_PrintStep/=0 .and. mod(t,E_PrintStep)==0) then
  write(*,*) t
endif
if (E_Sanity/=0 .and. mod(t,E_Sanity)==0 .and. t>0)  then
  write(*,*) 'Periodic sanity check'
  call Sanity
endif
if (E_Energy/=0 .and. mod(t,E_Energy)==0) then
  write(101,*) ETC
endif
!Call Network Analysis
if ( (E_LargestCluster/=0 .and. mod(t,E_LargestCluster)==0) .or. &
     (E_ClusterHist/=0 .and. mod(t,E_ClusterHist)==0) .or. &
     (E_RadDist_Cluster/=0 .and. mod(t,E_RadDist_Cluster)==0) .or. & 
     (E_XYZDist_Cluster/=0 .and. mod(t,E_XYZDist_Cluster)==0) .or. & 
     (E_RG_Cluster/=0 .and. mod(t,E_RG_Cluster)==0)) then
  call NetworkAnalysis
endif
if ( (E_RadDist_Cluster/=0 .and. mod(t,E_RadDist_Cluster)==0) .or. & 
     (E_XYZDist_Cluster/=0 .and. mod(t,E_XYZDist_Cluster)==0) .or. & 
     (E_RG_Cluster/=0 .and. mod(t,E_RG_Cluster)==0)) then
  !Calculate Center of Mass of the largest Cluster
  y1=0
  y2=0
  i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
  do i3=AllHostN(i1,1),AllHostN(i1,2)
    MPT=ProteinToProteinType(AllHost(i3))
    MS= (AllHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
    do i4=MS+1,MS+Vk(MPT)  !segment in protein
      y1=y1+sin(pi*2/BoxSize*List(:,i4))
      y2=y2+cos(pi*2/BoxSize*List(:,i4))
    enddo
  enddo
  COM_Cluster=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the cluster
    !Calculate Center of Mass of individual protein from the Cluster
  if(WriteIndiDist_Cluster) then
    COMi=0
    i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
    do i2=1,NkN !Which specific component is used for the COM?
      y1=0
      y2=0
      do i3=AllHostN(i1,1),AllHostN(i1,2) !runs through all proteins
        MPT=ProteinToProteinType(AllHost(i3))
        MS= (AllHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
        if (MPT==i2) then !if right protein type, run through all beads
        do i4=MS+1,MS+Vk(MPT)
          y1=y1+sin(pi*2/BoxSize*List(:,i4)) ! uses fourier series to find COM of protein type i2
          y2=y2+cos(pi*2/BoxSize*List(:,i4))
        enddo
        endif
      enddo
      COMi(:,i2)=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize)) 
    enddo
  endif
endif

if ( (E_RadDist_All/=0 .and. mod(t,E_RadDist_All)==0) .or. &
     (E_XYZDist_All/=0 .and. mod(t,E_XYZDist_All)==0) .or. &
     (E_RG_All/=0 .and. mod(t,E_RG_All)==0)) then
  ! calc Center of Mass of all polymers
  y1=0
  y2=0
  do i1=1,NkN
  do i2=1,Nk(i1)
  do i3=1 + (i2-1)*Vk(i1) + BeadBounds(i1), (i2)*Vk(i1) + BeadBounds(i1) 
    y1=y1+sin(pi*2/BoxSize*List(:,i3))
    y2=y2+cos(pi*2/BoxSize*List(:,i3))
  enddo
  enddo
  enddo
  COM_All=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the system
  !write(*,*) 'Center of Mass is', COM_All
endif

! cluster size of the three largest clusters and the polymers that are in them
if (E_LargestCluster/=0 .and. mod(t,E_LargestCluster)==0)  then
  Analysis_log=.true.
  LargestCluster=0
  do i1=1,min(3,AllHostNN) !three largest clusters
    i2=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1,Analysis_log(1:AllHostNN))
    LargestCluster(1+(NkN+1)*(i1-1))=AllHostN(i2,2)-AllHostN(i2,1)+1
    do i3=AllHostN(i2,1),AllHostN(i2,2) !all proteins in the cluster
      MPT=ProteinToProteinType(AllHost(i3))
      LargestCluster(1+(NkN+1)*(i1-1)+MPT) = &
          LargestCluster(1+(NkN+1)*(i1-1)+MPT)+1
    enddo
    Analysis_log(i2)=.false.
  enddo
  if(WriteCoreCluster) then
    LargestCluster(1+(NkN+1)*3)=CoreN !Core Cluster detail in second last columns
    do i3=1,CoreN
      MPT=ProteinToProteinType(Core(i3))
      write(*,*) 'allocate this properly!'
      stop
      LargestCluster(1+(NkN+1)*3+MPT) = &
          LargestCluster(1+(NkN+1)*3+MPT)+1
    enddo
        LargestCluster(1+(NkN+1)*4)=CoreHalfN
    do i3=1,CoreHalfN !Core Cluster detail in last columns
      MPT=ProteinToProteinType(CoreHalf(i3))
      LargestCluster(1+(NkN+1)*4+MPT) = &
          LargestCluster(1+(NkN+1)*4+MPT)+1
    enddo
  endif
  write(102,*) LargestCluster
  
  if (HasDNA) then
    LargestCluster=0
    LargestCluster(1)=DNAHostN
    do i3=1,DNAHostN
      MPT=ProteinToProteinType(DNAHost(i3))
      LargestCluster(1+MPT) = LargestCluster(1+MPT)+1
    enddo
    write(110,*) LargestCluster(1:NkN+1)
  endif
endif

! histogram of the cluster sizes (up to 20, all else gets put into the 21 bin)
if (E_ClusterHist/=0 .and. mod(t,E_ClusterHist)==0)  then
  ClusterHist=0
  do i1=1,AllHostNN
    if (AllHostN(i1,2)-AllHostN(i1,1)+1>20) then
      ClusterHist(21)=ClusterHist(21)+1
    else
      ClusterHist(AllHostN(i1,2)-AllHostN(i1,1)+1)= &
          ClusterHist(AllHostN(i1,2)-AllHostN(i1,1)+1)+1
    endif
  enddo
  write(104,*) ClusterHist(:)
endif



!XYZDistribution for all proteins, NOT RECENTERED!
if (E_XYZDist_None/=0 .and. mod(t,E_XYZDist_None)==0)  then
  XDist=0
  YDist=0
  ZDist=0
  do i1=1,NkN
  do i2=1+BeadBounds(i1),Nk(i1)*Vk(i1)+BeadBounds(i1)
    !now calculate the histogram
    XDist(List(1,i2),i1)= XDist(List(1,i2),i1) + 1
    YDist(List(2,i2),i1)= YDist(List(2,i2),i1) + 1
    ZDist(List(3,i2),i1)= ZDist(List(3,i2),i1) + 1
  enddo
  enddo
  do i1=1,NkN
    write(200+3*i1+1,*) XDist(:,i1)
    write(200+3*i1+2,*) YDist(:,i1)
    write(200+3*i1+3,*) ZDist(:,i1)
  enddo
endif

!XYZDistribution for all proteins, recentered on COM
if (E_XYZDist_All/=0 .and. mod(t,E_XYZDist_All)==0)  then
  XDist=0
  YDist=0
  ZDist=0
  do i1=1,NkN
  do i2=1+BeadBounds(i1),Nk(i1)*Vk(i1)+BeadBounds(i1)
    !now calculate the histogram
    x1=PB3(List(:,i2)-floor(COM_All)+BoxSize/2) !new coordinates re-centered on COM_All
    XDist(x1(1),i1)= XDist(x1(1),i1) + 1
    YDist(x1(2),i1)= YDist(x1(2),i1) + 1
    ZDist(x1(3),i1)= ZDist(x1(3),i1) + 1
  enddo
  enddo
  do i1=1,NkN
    write(300+3*i1+1,*) XDist(:,i1)
    write(300+3*i1+2,*) YDist(:,i1)
    write(300+3*i1+3,*) ZDist(:,i1)
  enddo
endif

!! !XYZDistribution (spacial organization) for the largest cluster, recentered on COM: 
if (E_XYZDist_Cluster/=0 .and. mod(t,E_XYZDist_Cluster)==0)  then
  i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
  XDist=0
  YDist=0
  ZDist=0
  if (WriteIndiDist_Cluster) then
    XDist_ij=0
    YDist_ij=0
    ZDist_ij=0
  endif
  if (AllHostN(i1,2)-AllHostN(i1,1)+1>=50) then ! no histogram if the cluster isn't at least 50 big
    !now calculate the histogram
    do i3=AllHostN(i1,1),AllHostN(i1,2)
      MPT=ProteinToProteinType(AllHost(i3))
      MS= (AllHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)
        x1=PB3(List(:,i4)-floor(COM_Cluster)+BoxSize/2) !new coordinates re-centered on COM
        XDist(x1(1),MPT)= XDist(x1(1),MPT) + 1
        YDist(x1(2),MPT)= YDist(x1(2),MPT) + 1
        ZDist(x1(3),MPT)= ZDist(x1(3),MPT) + 1
        if(WriteIndiDist_Cluster) then
          do i2=1,NkN
            x1=PB3(List(:,i4)-floor(COMi(:,i2))+BoxSize/2) !MPT centered on COM_i2
            XDist_ij(x1(1),MPT,i2)= XDist_ij(x1(1),MPT,i2) + 1
            YDist_ij(x1(2),MPT,i2)= YDist_ij(x1(2),MPT,i2) + 1
            ZDist_ij(x1(3),MPT,i2)= ZDist_ij(x1(3),MPT,i2) + 1
          enddo
        endif
      enddo
    enddo
  endif
  do i1=1,NkN
    write(400+3*i1+1,*) XDist(:,i1)
    write(400+3*i1+2,*) YDist(:,i1)
    write(400+3*i1+3,*) ZDist(:,i1)
  enddo

  if(WriteIndiDist_Cluster) then
  do i1=1,NkN
    do i2=1,NkN
      write(10000+ i1*NkN*3 + i2*3 + 1,*) XDist_ij(:,i1,i2)
      write(10000+ i1*NkN*3 + i2*3 + 2,*) YDist_ij(:,i1,i2)
      write(10000+ i1*NkN*3 + i2*3 + 3,*) ZDist_ij(:,i1,i2)
    enddo
  enddo
  endif
endif

!! Radial Histogram: (spacial organization) for the largest cluster
if (E_RadDist_Cluster/=0 .and. mod(t,E_RadDist_Cluster)==0)  then
  i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
  RadDist=0
  if (AllHostN(i1,2)-AllHostN(i1,1)+1>=50) then ! no histogram if the cluster isn't at least 50 big
    if(WriteIndiDist_Cluster) then
    RadDist_ij=0
    do i2=1,NkN !COM of protein i2
      !now calculate the histogram
      do i3=AllHostN(i1,1),AllHostN(i1,2)
        MPT=ProteinToProteinType(AllHost(i3))
        MS= (AllHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
        do i4=MS+1,MS+Vk(MPT)  !segment in protein
          y2=real(List(:,i4))
          i5=floor(sqrt(sum(rDistPB(y2,COMi(:,i2))**2))+1)
          MPT2=ProteinToProteinType(AllHost(i3))
          RadDist_ij(i5,MPT2,i2)= RadDist_ij(i5,MPT2,i2) + 1
        enddo
      enddo
    enddo
    endif
    ! COM_Cluster is based on all proteins in the largest cluster
    i2=NkN+1!I did not change the way of assign memory, which is not most effcient in terms of memory usage, but current format is easier in programming.
    do i3=AllHostN(i1,1),AllHostN(i1,2)
      MPT=ProteinToProteinType(AllHost(i3))
      MS= (AllHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)  !segment in protein
        y2=real(List(:,i4))
        i5=floor(sqrt(sum(rDistPB(y2,COM_Cluster)**2))+1)
        MPT2=ProteinToProteinType(AllHost(i3))
        RadDist(i5,MPT2)= RadDist(i5,MPT2) + 1
      enddo
    enddo
  endif
  do i1=1,NkN
    if (WriteIndiDist_Cluster) then
    do i2=1,NkN
      write(500 + i1*NkN + i2,*) RadDist_ij(:,i1,i2)
    enddo
    endif
    write(500+i1,*) RadDist(:,i1)
  enddo

  ! Radial Histgram for Core_Cluster, very expensive!
  If(WriteCoreCluster) then
    if(WriteIndiDist_Cluster) then
    RadDist_ij=0
    do i2=1,NkN !Which specific component is used for the com?
      y1=0
      y2=0
      do i3=1,CoreN !runs through all proteins
        MPT=ProteinToProteinType(Core(i3))
        MS= (Core(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
        if (MPT==i2) then !if right protein type, run through all beads
        do i4=MS+1, MS+Vk(MPT)
          y1=y1+sin(pi*2/BoxSize*List(:,i4)) ! uses fourier series to find COM
          y2=y2+cos(pi*2/BoxSize*List(:,i4))
        enddo
        endif
      enddo
      y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the protein type i2 in Core_cluster
      do i3=1,CoreN
        MPT=ProteinToProteinType(Core(i3))
        MS= (Core(i3) - sNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
        do i4=MS+1,MS+Vk(MPT)  !segment in protein
          y2=real(List(:,i4))
          i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
          RadDist_ij(i5,MPT,i2)= RadDist_ij(i5,MPT,i2) + 1
        enddo
      enddo
    enddo
    endif

    RadDist=0
    i2=NkN+1
    y1=0
    y2=0
    do i3=1,CoreN
      MPT=ProteinToProteinType(Core(i3))
      MS= (Core(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)  !segment in protein
        y1=y1+sin(pi*2/BoxSize*List(:,i4)) ! uses fourier series to find COM
        y2=y2+cos(pi*2/BoxSize*List(:,i4))
      enddo
    enddo
    y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the Core_cluster
    do i3=1,CoreN
      MPT=ProteinToProteinType(Core(i3))
      MS= (Core(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)  !segment in protein
        y2=real(List(:,i4))
        i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
        RadDist(i5,MPT)= RadDist(i5,MPT) + 1
      enddo
    enddo
    do i1=1,NkN
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(20000+2*i1*NkN+i2*2+1,*) RadDist_ij(:,i1,i2)
      enddo
      endif
      write(20000+2*NkN*NkN+2*NkN+3,*) RadDist(:,i1)
    enddo

    if(WriteIndiDist_Cluster) then
    RadDist_ij=0
    do i2=1,NkN !Which specific component is used for the com?
      y1=0
      y2=0
      do i3=1,CoreHalfN !runs through all proteins
        MPT=ProteinToProteinType(CoreHalf(i3))
        MS= (CoreHalf(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
        if (MPT==i2) then !if right protein type, run through all beads
        do i4=MS+1, MS+Vk(MPT)
          y1=y1+sin(pi*2/BoxSize*List(:,i4)) ! uses fourier series to find COM
          y2=y2+cos(pi*2/BoxSize*List(:,i4))
        enddo
        endif
      enddo
      y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the protein type i2 in CoreHalf_cluster
      do i3=1,CoreHalfN
        MPT=ProteinToProteinType(CoreHalf(i3))
        MS= (CoreHalf(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
        do i4=MS+1,MS+Vk(MPT)
          y2=real(List(:,i4))
          i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
          RadDist_ij(i5,MPT,i2)= RadDist_ij(i5,MPT,i2) + 1
        enddo
      enddo
    enddo
    endif

    RadDist=0
    y1=0
    y2=0
    do i3=1,CoreHalfN
      MPT=ProteinToProteinType(CoreHalf(i3))
      MS= (CoreHalf(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)
        y1=y1+sin(pi*2/BoxSize*List(:,i4)) ! uses fourier series to find COM
        y2=y2+cos(pi*2/BoxSize*List(:,i4))
      enddo
    enddo
    y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the CoreHalf_cluster
    do i3=1,CoreHalfN
      MPT=ProteinToProteinType(CoreHalf(i3))
      MS= (CoreHalf(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)
        y2=real(List(:,i4))
        i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
        RadDist(i5,MPT)= RadDist(i5,MPT) + 1
      enddo
    enddo
    do i1=1,NkN
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(Str1,*) i1
        write(Str2,*) i2
        write(20000+2*i1*NkN+i2*2+2,*) RadDist_ij(:,i1,i2)
      enddo
      endif
      write(20000+2*NkN*NkN+2*NkN+4,*) RadDist(:,i1)
    enddo
  endif
endif

if (E_RG_All/=0 .and. mod(t,E_RG_All)==0)  then
  !calc Gyration Tensor
  GT=0
  do i1=1,NkN
  do i2=1+BeadBounds(i1),Nk(i1)*Vk(i1)+BeadBounds(i1)
    y3=rDistPB(real(List(:,i2)),COM_All)
    do i3=1,3
    do i4=i3,3
      GT(i3,i4)=GT(i3,i4)+y3(i3)*y3(i4)/sNkVk(NkN+1)
    enddo
    enddo
  enddo
  enddo
  GT(2,1)=GT(1,2)
  GT(3,1)=GT(1,3)
  GT(3,2)=GT(2,3)

  !diagonalize Gyration Tensor
  b=GT(1,1)+GT(2,2)+GT(3,3)
  c=GT(1,1)*GT(2,2)+GT(1,1)*GT(3,3)+GT(2,2)*GT(3,3)-GT(1,2)**2-GT(1,3)**2-GT(2,3)**2
  d=GT(1,1)*GT(2,3)**2+GT(2,2)*GT(1,3)**2+GT(3,3)*GT(1,2)**2-GT(1,1)*GT(2,2)*GT(3,3)-2*GT(1,2)*GT(1,3)*GT(2,3)
  p=b**2-3*c !this can be negative?
  q=2*b**3-9*b*c-27*d
  d=acos(q/2/sqrt(p**3))
  RG_All(1) = 1./3*(b+2.*sqrt(p)*cos(d/3))
  RG_All(2) = 1./3*(b+2.*sqrt(p)*cos(d/3+pi*2/3))
  RG_All(3) = 1./3*(b+2.*sqrt(p)*cos(d/3-pi*2/3))
  RG_All(4) = sqrt(RG_All(1)+RG_All(2)+RG_All(3))! the real RG
  write(103,*) RG_All
  RG_All=0
end if

if (E_RG_Cluster/=0 .and. mod(t,E_RG_Cluster)==0)  then
  !calc Center of Mass
  i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
  NB_Cluster=0
  do i3=AllHostN(i1,1),AllHostN(i1,2)
    NB_Cluster=NB_Cluster+Vk(ProteinToProteinType(AllHost(i3)))
  enddo
  !calc Gyration Tensor
  GT=0
  do i3=AllHostN(i1,1),AllHostN(i1,2)
    MPT=ProteinToProteinType(AllHost(i3))
    MS= (AllHost(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
    do i4=MS+1,MS+Vk(MPT)
      y3=rDistPB(real(List(:,i4)),COM_Cluster)
      do i2=1,3
      do i5=i2,3
        GT(i2,i5)=GT(i2,i5)+y3(i2)*y3(i5)/ NB_Cluster
        GT(i5,i2)=GT(i2,i5)
      enddo
      enddo
    enddo
  enddo
  !diagonalize Gyration Tensor
  b=GT(1,1)+GT(2,2)+GT(3,3)
  c=GT(1,1)*GT(2,2)+GT(1,1)*GT(3,3)+GT(2,2)*GT(3,3)-GT(1,2)**2-GT(1,3)**2-GT(2,3)**2
  d=GT(1,1)*GT(2,3)**2+GT(2,2)*GT(1,3)**2+GT(3,3)*GT(1,2)**2-GT(1,1)*GT(2,2)*GT(3,3)-2*GT(1,2)*GT(1,3)*GT(2,3)
  p=b**2-3*c !this can be negative?
  q=2*b**3-9*b*c-27*d
  d=acos(q/2/sqrt(p**3))
  RG_Cluster(1) = 1./3*(b+2.*sqrt(p)*cos(d/3)) !First Eigven value
  RG_Cluster(2) = 1./3*(b+2.*sqrt(p)*cos(d/3+pi*2/3))
  RG_Cluster(3) = 1./3*(b+2.*sqrt(p)*cos(d/3-pi*2/3))
  RG_Cluster(4) = sqrt(RG_Cluster(1)+RG_Cluster(2)+RG_Cluster(3))! the real RG
  RG_Cluster(5) = sqrt(5.*RG_Cluster(1)) !Actually semiaxis size
  RG_Cluster(6) = sqrt(5.*RG_Cluster(2))
  RG_Cluster(7) = sqrt(5.*RG_Cluster(3))
  RG_Cluster(8) = 4*pi*RG_Cluster(5)*RG_Cluster(6)*RG_Cluster(7)/3 
  write(107,*) RG_Cluster
  RG_Cluster=0

  if(WriteCoreCluster) then
    !calc Center of Mass of Core
    y1=0
    y2=0
    NB_Cluster=0
    do i3=1,CoreN
      MPT=ProteinToProteinType(Core(i3))
      MS= (Core(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)
        y1=y1+sin(pi*2/BoxSize*List(:,i4))
        y2=y2+cos(pi*2/BoxSize*List(:,i4))
        NB_Cluster=NB_Cluster+1
      enddo
    enddo
    y1=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the cluster
    !write(*,*) 'Center of Mass of Largest Cluster is', y1

    
    !calc Gyration Tensor
    GT=0
    do i3=1,CoreN
      MPT=ProteinToProteinType(Core(i3))
      MS= (Core(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)
        y3=rDistPB(real(List(:,i4)),y1)
        do i2=1,3
        do i5=i2,3
          GT(i2,i5)=GT(i2,i5)+y3(i2)*y3(i5)/ NB_Cluster
          GT(i5,i2)=GT(i2,i5)
        enddo
        enddo
      enddo
    enddo

    !diagonalize Gyration Tensor
    b=GT(1,1)+GT(2,2)+GT(3,3)
    c=GT(1,1)*GT(2,2)+GT(1,1)*GT(3,3)+GT(2,2)*GT(3,3)-GT(1,2)**2-GT(1,3)**2-GT(2,3)**2
    d=GT(1,1)*GT(2,3)**2+GT(2,2)*GT(1,3)**2+GT(3,3)*GT(1,2)**2-GT(1,1)*GT(2,2)*GT(3,3)-2*GT(1,2)*GT(1,3)*GT(2,3)
    p=b**2-3*c !this can be negative?
    q=2*b**3-9*b*c-27*d
    d=acos(q/2/sqrt(p**3))
    RG_Core(1) = 1./3*(b+2.*sqrt(p)*cos(d/3)) !First Eigven value
    RG_Core(2) = 1./3*(b+2.*sqrt(p)*cos(d/3+pi*2/3))
    RG_Core(3) = 1./3*(b+2.*sqrt(p)*cos(d/3-pi*2/3))
    RG_Core(4) = sqrt(RG_Core(1)+RG_Core(2)+RG_Core(3))! the real RG
    RG_Core(5) = sqrt(5.*RG_Core(1)) !Actually semiaxis size
    RG_Core(6) = sqrt(5.*RG_Core(2))
    RG_Core(7) = sqrt(5.*RG_Core(3))
    RG_Core(8) = 4*pi*RG_Core(5)*RG_Core(6)*RG_Core(7)/3 
    write(108,*) RG_Core
    RG_Core=0

    !calc Center of Mass of CoreHalf
    y1=0
    y2=0
    NB_Cluster=0
    do i3=1,CoreHalfN
      MPT=ProteinToProteinType(CoreHalf(i3))
      MS= (CoreHalf(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)
        y1=y1+sin(pi*2/BoxSize*List(:,i4))
        y2=y2+cos(pi*2/BoxSize*List(:,i4))
        NB_Cluster=NB_Cluster+1
      enddo
    enddo
    y1=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the cluster
    !write(*,*) 'Center of Mass of Largest Cluster is', y1
    !calc Gyration Tensor
    GT=0
    do i3=1,CoreHalfN
      MPT=ProteinToProteinType(CoreHalf(i3))
      MS= (CoreHalf(i3) - smNk(MPT) -1)*Vk(MPT) + BeadBounds(MPT)
      do i4=MS+1,MS+Vk(MPT)
        y3=rDistPB(real(List(:,i4)),y1)
        do i2=1,3
        do i5=i2,3
          GT(i2,i5)=GT(i2,i5)+y3(i2)*y3(i5)/ NB_Cluster
          GT(i5,i2)=GT(i2,i5)
        enddo
        enddo
      enddo
    enddo
    !diagonalize Gyration Tensor
    b=GT(1,1)+GT(2,2)+GT(3,3)
    c=GT(1,1)*GT(2,2)+GT(1,1)*GT(3,3)+GT(2,2)*GT(3,3)-GT(1,2)**2-GT(1,3)**2-GT(2,3)**2
    d=GT(1,1)*GT(2,3)**2+GT(2,2)*GT(1,3)**2+GT(3,3)*GT(1,2)**2-GT(1,1)*GT(2,2)*GT(3,3)-2*GT(1,2)*GT(1,3)*GT(2,3)
    p=b**2-3*c !this can be negative?
    q=2*b**3-9*b*c-27*d
    d=acos(q/2/sqrt(p**3))
    RG_CoreHalf(1) = 1./3*(b+2.*sqrt(p)*cos(d/3)) !First Eigven value
    RG_CoreHalf(2) = 1./3*(b+2.*sqrt(p)*cos(d/3+pi*2/3))
    RG_CoreHalf(3) = 1./3*(b+2.*sqrt(p)*cos(d/3-pi*2/3))
    RG_CoreHalf(4) = sqrt(RG_CoreHalf(1)+RG_CoreHalf(2)+RG_CoreHalf(3))! the real RG
    RG_CoreHalf(5) = sqrt(5.*RG_CoreHalf(1)) !Actually semiaxis size
    RG_CoreHalf(6) = sqrt(5.*RG_CoreHalf(2))
    RG_CoreHalf(7) = sqrt(5.*RG_CoreHalf(3))
    RG_CoreHalf(8) = 4*pi*RG_CoreHalf(5)*RG_CoreHalf(6)*RG_CoreHalf(7)/3 
    write(109,*) RG_CoreHalf
    RG_CoreHalf=0
  endif
end if

if (E_BoundSites/=0 .and. mod(t,E_BoundSites)==0)  then
  BoundSites=0
  do i1=1,NkN
  do i2=1,Nk(i1)
  do i3=1,Vk(i1)
    MS=i3+(i2-1)*Vk(i1)+BeadBounds(i1)
    if (BoundTo(MS) /=0) then
      BoundSites(i3+sVk(i1))=BoundSites(i3+sVk(i1))+1
    endif
  enddo
  enddo
  enddo
  write(105,*) BoundSites
endif

if (E_BondTypes/=0 .and. mod(t,E_BondTypes)==0) then
  BondTypes=0
  do i1=1,NkN
  do i2=1,Nk(i1)
  do i3=1,Vk(i1)
    MS=i3+(i2-1)*Vk(i1)+BeadBounds(i1)
    if (BoundTo(MS) /=0) then
      BondTypes(i3+sVk(i1),BeadToMSI(BoundTo(i2)))= &
        BondTypes(i3+sVk(i1),BeadToMSI(BoundTo(i2)))+1
    endif
  enddo
  enddo
  enddo
  do i1=1,sVk(NkN+1)
    write(106,*) BondTypes(i1,:)
  enddo
endif

if(UseSpringLinkerEq)then
if(E_SpringLinkerVector/=0 .and. mod(t,E_SpringLinkerVector)==0 ) then
  do i1=1,NkN
  do i2=1,Vk(i1)
    MSI=i2 + sVk(i1)
    do i3=1,SLConN(MSI)
    if (SLCon(i3,MSI)>MSI) then
      SLVec=0
      do i4=1,Nk(i1)
        MS= i2 + (i4-1)*Vk(i1) + BeadBounds(i1)
        MS2= SLCon(i3,MSI) - sVk(i1) + (i4-1)*Vk(i1) + BeadBounds(i1)
        x1=List(:,MS)
        x2=List(:,MS2)
        dx=DistPB(x1,x2)
        SLVec(1:3)=SLVec(1:3)+dx! average values of x y z
        SLVec(4)=SLVec(4)+sqrt(real(sum(dx**2))) ! average length
        SLVec(5)=SLVec(5)+acos(SLVec(1)/SLVec(4)) ! average alpha (angle between x axis)
        SLVec(6)=SLVec(6)+acos(SLVec(2)/SLVec(4))
        SLVec(7)=SLVec(7)+acos(SLVec(3)/SLVec(4))
      enddo
      SLVec=SLVec/Nk(i1)
      write(99,*) SLVec
    endif
    enddo
  enddo
  enddo
endif
endif

if(E_SelfLoop/=0 .and. mod(t,E_SelfLoop)==0 ) then
  SelfLoopN=0
  do i1=1, Nk(LinkerType)
    MS = BoundTo(1+(i1-1)*Vk(LinkerType)+BeadBounds(LinkerType))
    MS2= BoundTo((i1)*Vk(LinkerType)+BeadBounds(LinkerType))
    if (MS/=0 .and. MS2/=0) then !both ends bound
      if (BeadToProteinIndex(MS) ==BeadToProteinIndex(MS2)) then !both ends bound to same protein
        SelfLoopN=SelfLoopN+1
      endif
    endif
  enddo
  write(97,*) SelfLoopN
endif

if(E_Network/=0 .and. mod(t,E_Network)==0 ) then
  LinkerCon=0
  do i1=1, Nk(LinkerType)
    LinkerCon(i1,1)=i1
    i3=1
    MS=(i1-1)*Vk(LinkerType)+BeadBounds(LinkerType)
    do i2=MS+1,MS+Vk(LinkerType)
    if (BoundTo(i2)/=0) then
      i3=i3+1
      LinkerCon(i1,i3)=BeadToProteinIndex(BoundTo(i2))
    endif
    enddo
  enddo
  do i1=1, Nk(LinkerType)
    write(98,*) LinkerCon(i1,:)
  enddo
endif

if (E_Lattice/=0 .and. mod(t,E_Lattice)==0)  then
  write(1,*) BeadBounds(NkN+1)
  write(1,*)
  i3=0
  str2='   '
  do i1=1,NkN
  do i2=1,Nk(i1)
  do i3=1,Vk(i1)
    MS=i3+(i2-1)*Vk(i1)+BeadBounds(i1)
    write(str1,'(A,I1,A)') '(A,I',floor(log10(real(i2))+1),',A, I5, I5, I5)'
    write(1, str1) Char(64 + i1),i2,str2(1:4-floor(log10(real(i2))+1)), List(:,MS)
    if (WriteBoundTo)  write(2,*) BoundTo(MS)
  enddo
  enddo
  enddo
endif

if (E_Nk/=0 .and. mod(t,E_Nk)==0) then
  write(3,*) Nk
endif
  

end subroutine Analysis
!!!!
!!!!
!!!!
!!!!
subroutine Sanity
call system_clock(tclock3)


! check list->lattice
do i1=1,NkN
do i2=1,Nk(i1)*Vk(i1)
  MS=i2+BeadBounds(i1)
  if (any(List(:,MS)<=0 .or. List(:,MS)>BoxSize)) then
    write(*,*) 'List coordinate is outside the box??'
    write(*,*) 'MoveType:',MCMove
    write(*,*) 'Step Number',t
    write(*,*) 'Protein:',BeadToProteinIndex(MS)
    write(*,*) 'Bead:',MS
    write(*,*) 'Listed as at:', List(:,MS)
    stop
  endif
  if (Lattice(List(1,MS),List(2,MS),List(3,MS))/=MS) then
    write(*,*) 'Lattice doesn`t have a segment listed in List'
    write(*,*) 'MoveType:',MCMove
    write(*,*) 'Step Number',t
    write(*,*) 'Protein:',BeadToProteinIndex(MS)
    write(*,*) 'Bead:',MS
    write(*,*) 'Listed as at:', List(:,MS)
    write(*,*) 'This Lattice site has:', Lattice(List(1,MS),List(2,MS),List(3,MS))
    write(*,*) BeadToProteinIndex(Lattice(List(1,MS),List(2,MS),List(3,MS)))
    if (Lattice(List(1,MS),List(2,MS),List(3,MS))/=0) then
      write(*,*) List(:,Lattice(List(1,MS),List(2,MS),List(3,MS)))
    endif
    
    do i3=1,BoxSize(1)
    do i4=1,BoxSize(2)
    do i5=1,BoxSize(3)
    if (Lattice(i3,i4,i5)==MS) then
      write(*,*) 'Found segment at a different location:',MS
      write(*,*) i3,i4,i5
      stop
    endif
    enddo
    enddo
    enddo
    write(*,*) 'Segment not found anywhere!',MS
    stop
  endif
enddo
enddo

! check Lattice->List
do i1=1,BoxSize(1)
do i2=1,BoxSize(2)
do i3=1,BoxSize(3)
  MS=Lattice(i1,i2,i3)
  if (MS/=0) then
  if (any(List(:,MS)/=[i1,i2,i3])) then
    write(*,*) 'List doesn`t match the segment placement on lattice'
    write(*,*) 'Step Number',t
    write(*,*) MS
    write(*,*) List(:,MS)
    write(*,*) [i1,i2,i3]
    write(*,*) 'here1'
    write(*,*) sNk(NkN+1)
    write(*,*) rr3(1)
    write(*,*) 
    write(*,*) MCMove
    write(*,*) MCMoves
    stop
  endif
  endif
enddo
enddo
enddo

! Run through proteins
EM=0
EC=0
do i1=1,NkN
do i2=1,Nk(i1)*Vk(i1)
  MS=i2+BeadBounds(i1)
  MSI=BeadToMSI(MS)
  MS2=BoundTo(MS)
  if (MS2/=0) then
    MST=BeadType(MSI)
    ! BoundTo->BoundTo(BoundTo)
    if (BoundTo(MS2)/=MS) then
      write(*,*) 'Something isn`t bound to what is bound to it'
      write(*,*) 'Step Number:',t
      write(*,*) 'Last move',MCMove
      write(*,*)
      write(*,*) 'bead',MS
      write(*,*) 'bound to',MS2
      write(*,*) 'which is instead bound to', BoundTo(MS2)
      if (BoundTo(MS2)/=0) write(*,*) 'That is bound to',BoundTo(BoundTo(MS2))
      stop
    endif
    ! Binding length =1
    x1=List(:,MS)
    x2=List(:,MS2)
    if (maxval(abs(DistPB(x1,x2)))>1) then
      write(*,*) 'Bound Modules are too far apart to actually bind'
      write(*,*) 'failed on timestep:', t
      write(*,*) MS
      write(*,*) x1
      write(*,*) MS2
      write(*,*) x2
      write(*,*) 'Should be bound to something bound to itself:'
      write(*,*) BoundTo(MS2)
      write(*,*) 'Scaled to closest neighbors'
      write(*,*) DistPB(x1,x2)
      write(*,*) 'Last move type:'
      write(*,*) MCMove
      write(*,*) MCMoves
      write(*,*) 'Bead ',MS,' is on protein ',BeadToProteinIndex(MS)
      write(*,*) 'Bead ',MS2,' is on protein ',BeadToProteinIndex(MS2)
      stop
    endif
    if (MS2.gt.MS) then
      MST2=BeadType(BeadToMSI(MS2))
      EM=EM+Energy(MST,MST2)
      if (UseESelf) then
      if (BeadToProteinIndex(MS)==BeadToProteinIndex(MS2)) then
        EM=EM+floor(ESelf*exp(real(-(abs(MS-MS2)-1)/LSteric)))
      endif
      endif
      if (Energy(MST,MST2)==0) then
        write(*,*) 'An interaction has formed between domains that arent allowed to interact'
        write(*,*) t
        write(*,*) MCMoves
        write(*,*) MCMove
        write(*,*) i1
        write(*,*) MS
        write(*,*) MST
        write(*,*) MS2
        write(*,*) MST2
        stop
      endif
    endif
    
    ! Linker length
    do i4=1,ConnectionN(MSI)
    if (Connection(i4,MSI).gt.MSI) then
      x2=List(:,MS-MSI+Connection(i4,MSI))
      x3=DistPB(x1,x2)
      if (maxval(abs(x3)).gt.LL(i4,MSI)) then
        write(*,*) 'Linkers are too far apart, bead and its tethered bead:'
        write(*,*) MS !bead number
        write(*,*) MS-MSI+Connection(i4,MSI)
        write(*,*) x1 !connected bead
        write(*,*) x2
        write(*,*) 'Move',MCmove
        write(*,*) 'time', t
        write(*,*) 'Distances'
        write(*,*) x1
        write(*,*) x2
        write(*,*) x3
        write(*,*) LL(i4,MSI)
        write(*,*)
        write(*,*) MPT,MS,i4
        stop
      endif
    endif
    enddo
  endif
  
  if (UseSpringLinkerEq) then
  do i4=1,SLConN(MSI)
  if (SLCon(i4,MSI)>MSI) then
    x1=List(:,MS)
    x2=List(:,MS-MSI+SLCon(i4,MSI))
    dx=DistPB(x1,x2)
    dx_max=floor(sqrt(real(sum(dx**2)))-SLEq(i4,MSI))!3d distance differnce to Eq Leng for MS1
    EM=EM + SLPot(i4,MSI)*dx_max**2
  endif
  enddo
  endif
enddo
enddo

if (ETC/=EM) then
  write(*,*) 
  write(*,*) 'Energy is wrong.  This is a symptom of a greater problem.'
  write(*,*) 'Move type:',MCMove
  write(*,*) 'Current Step: ',t
  write(*,*) 'Current Energy: ',EM
  write(*,*) 'Energy Saved From Simulation ETC: ',ETC
  write(*,*) 'Difference in energy from expected and seen',ETC-EM
  stop
endif

if (HasDNA) then
  if (any(List(1,1:Vk(1))/=BoxSize(1)/2)) then
    write(*,*) 'DNA has moved in X!  Crashing'
    write(*,*) t
    write(*,*) MCMove
    write(*,*) BoxSize/2
    write(*,*)
    write(*,*) List(1,1:Vk(1))
    write(*,*) (List(1,1:Vk(1))/=BoxSize(1)/2)
    write(*,*) any(List(1,1:Vk(1))/=BoxSize(1)/2)
    write(*,*)
    do i1=1,Vk(1)
      write(*,*) List(:,i1)
    enddo
    stop
  elseif (any(List(2,1:Vk(1))/=BoxSize(2)/2)) then
    write(*,*)
    write(*,*) 'DNA has moved in Y!  Crashing'
    do i1=1,Vk(1)
      write(*,*) List(:,i1)
    enddo
    stop
  elseif (any(List(3,1:Vk(1)) /= [(i, i = 1, Vk(1))])) then
    write(*,*)
    write(*,*) 'DNA has moved in Z!  Crashing'
    do i1=1,Vk(1)
      write(*,*) List(:,i1)
    enddo
    stop
  endif
endif

call system_clock(tclock4)
tclock_sanity=tclock_sanity+tclock4-tclock3
end subroutine Sanity
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
subroutine EndofProgram
if (E_SelfLoop/=0) then
  close(97)
endif

if (E_Network/=0) then
  close(98)
endif

if (E_SpringLinkerVector/=0) then
  close(99)
endif

if (E_Energy/=0) then
  close(101)
endif

if (E_LargestCluster/=0)  then
  close(102)
  if (HasDNA) then
    close(110)
  endif
endif

if (E_RG_All/=0)  then
  close(103)
endif

if (E_RG_Cluster/=0)  then
  close(107)
  if(WriteCoreCluster) then
    close(108)
    close(109)
  endif
endif

if (E_ClusterHist/=0)  then
  close(104)
endif

if (E_BoundSites/=0)  then
  close(105)
endif

if (E_BondTypes/=0) then
  close(106)
endif

!! Slab Histogram: (spacial organization)
if (E_XYZDist_None/=0)  then
  do i1=1,NkN
    do i2=1,3
      close(200+3*i1+i2)
    enddo
  enddo
endif

if (E_XYZDist_All/=0)  then
  do i1=1,NkN
    do i3=1,3
      close(300+3*i1+i3)
    enddo
  enddo
endif

if (E_XYZDist_Cluster/=0)  then
  do i1=1,NkN
  do i3=1,3
    close(400+3*i1+i3)
  enddo
  enddo
  if(WriteIndiDist_Cluster) then
  do i1=1,NkN
  do i2=1,NkN
  do i3=1,3
    close(10000+ i1*NkN*3 + i2*3 + i3)
  enddo
  enddo
  enddo
  endif
endif

!! Radial Histogram: (spacial organization)
if (E_RadDist_Cluster/=0)  then
  do i1=1,NkN
    if(WriteIndiDist_Cluster) then
    do i2=1,NkN
      close(500 + i1*NkN + i2)
    enddo
    endif
    close(500+i1)
  enddo
  if(WriteCoreCluster) then
    !Open Core Radial Hist
    do i1=1,NkN
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        close(20000+2*i1*NkN+i2*2+1)
        close(20000+2*i1*NkN+i2*2+2)
      enddo
      endif
      close(20000+2*NkN*NkN+2*NkN+3)
      close(20000+2*NkN*NkN+2*NkN+4)
    enddo
  endif
endif


if (E_Lattice/=0) then
  close(1)
  if (WriteBoundTo) close(2)
endif
if (E_Nk/=0) then
  close(3)
endif


if (WriteLatticeEnd) then
  open(unit=1,file='End_Parms.txt')
  write(1,*) BoxSize
  write(1,*) NkN
  write(1,*) Nk
  write(1,*) Vk
  do i1=1,sVk(NkN+1)
    write(1,*) ConnectionN(i1)
    write(1,*) Connection(1:ConnectionN(i1),i1)
    write(1,*) LL(1:ConnectionN(i1),i1)
  enddo
  close(1)
  
  open(unit=1,file='End_ListForm.txt')
  open(unit=2,file='End_BoundTo.txt')
  do i1=1,NkN
  do i2=1+BeadBounds(i1),Nk(i1)*Vk(i1)+BeadBounds(i1)
    write(1,*) List(:,i2)
    write(2,*) BoundTo(i2)
  enddo
  enddo
  close(1)
  close(2)
endif


call system_clock(tclock2)
write(*,*)   'Time Spent During Simulation :'
write(*,'(I17,A,I10,A)') tclock2-tclock1, ' milliseconds, ~', (tclock2-tclock1)/3600000, ' hours'
write(*,*)   'Time Spent During Sanity :'
write(*,'(I17,A,I10,A)') tclock_sanity, ' milliseconds, ~', tclock_sanity/3600000, ' hours'
write(*,*)
write(*,*) '         No Moves,   Rejected,   Proposed,   Accepted'
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'TransI ', real(Rej(1,2)), real(Rej(1,1)), real(Propose(1)), &
    real(Propose(1) - sum(Rej(1,:)))
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'TransII', real(Rej(2,2)), real(Rej(2,1)), real(Propose(2)), &
    real(Propose(2) - sum(Rej(2,:)))
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'Rot    ', real(Rej(3,2)), real(Rej(3,1)), real(Propose(3)), &
    real(Propose(3) - sum(Rej(3,:)))
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'SltI   ', real(Rej(4,2)), real(Rej(4,1)), real(Propose(4)), &
    real(Propose(4) - sum(Rej(4,:)))
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'SltII  ', real(Rej(5,2)), real(Rej(5,1)), real(Propose(5)), &
    real(Propose(5) - sum(Rej(5,:)))
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'ClstI  ', real(Rej(6,2)), real(Rej(6,1)), real(Propose(6)), &
    real(Propose(6) - sum(Rej(6,:)))
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'ClstII ', real(Rej(7,2)), real(Rej(7,1)), real(Propose(7)), &
    real(Propose(7) - sum(Rej(7,:)))
write(*,'(A, ES12.2, ES12.2, ES12.2, ES12.2)') 'Grand  ', real(Rej(8,2)), real(Rej(8,1)), real(Propose(8)), &
    real(Propose(8) - sum(Rej(8,:)))

write(*,*)
write(*,'(A, I12,I12,I12,I12)') 'TransII', Rej(1,2), Rej(1,1), Propose(1), Propose(1)-sum(Rej(1,:))
write(*,'(A, I12,I12,I12,I12)') 'TransII', Rej(2,2), Rej(2,1), Propose(2), Propose(2)-sum(Rej(2,:))
write(*,'(A, I12,I12,I12,I12)') 'Rot    ', Rej(3,2), Rej(3,1), Propose(3), Propose(3)-sum(Rej(3,:))
write(*,'(A, I12,I12,I12,I12)') 'SltI   ', Rej(4,2), Rej(4,1), Propose(4), Propose(4)-sum(Rej(4,:))
write(*,'(A, I12,I12,I12,I12)') 'SltII  ', Rej(5,2), Rej(5,1), Propose(5), Propose(5)-sum(Rej(5,:))
write(*,'(A, I12,I12,I12,I12)') 'ClstI  ', Rej(6,2), Rej(6,1), Propose(6), Propose(6)-sum(Rej(6,:))
write(*,'(A, I12,I12,I12,I12)') 'ClstII ', Rej(7,2), Rej(7,1), Propose(7), Propose(7)-sum(Rej(7,:))
write(*,'(A, I12,I12,I12,I12)') 'Grand  ', Rej(8,2), Rej(8,1), Propose(8), Propose(8)-sum(Rej(8,:))
end subroutine EndofProgram

END PROGRAM First
