PROGRAM First



! imported values from key file
Character (len=100) :: PotFile, Str, Str1, Str2, Str3, keyfile, Parm
Character (len=100) :: ConFile, SLEqFile
Character (len=26) :: Alph
INTEGER, DIMENSION(:,:), Allocatable :: Energy, Connection, LL_temp, SLCon, SLPot_temp
INTEGER, DIMENSION(:), Allocatable :: ConnectionN, Connectioni, Connectionj, SLConN, SLConi, SLConj, StartDropletNk, ChemPot
INTEGER :: keyfileN, nargs, Tmax, Ti,Trelax, NkN,mNk, BoxSize(3), StartDropletR, StartSlabH, ESelf
INTEGER, DIMENSION(:), Allocatable :: Nk, Vk, Vak, StartSlabNk
Logical :: StartDroplet, StartImported, UseSlitherII, UseSlitherI, StartSlab
Logical :: WriteLatticeEnd, UseGrand, UseSpringLinkerEq, WriteCoreCluster, WriteIndiDist_Cluster
INTEGER, DIMENSION(:), Allocatable :: UseSlitherI_List, UseSlitherII_List
INTEGER :: UseSlitherI_ListN, UseSlitherII_ListN, LinkerType, SelfLoopN
INTEGER :: E_Sanity, E_Lattice, E_LargestCluster, E_ClusterHist, E_BoundSites, E_BondTypes, E_Energy, E_PrintStep
INTEGER :: E_RadDist_Cluster, E_SpringLinkerVector, E_Network, E_SelfLoop
INTEGER :: E_XYZDist_All, E_XYZDist_Cluster, E_RG_All, E_RG_Cluster
INTEGER :: iMoveRot, iMoveTransI,iMoveTransII,iMoveSlitherI,iMoveSlitherII,iMoveClusterI,iMoveClusterII,iMoveGrand
REAL :: LSteric
! Global Book keeping
INTEGER, DIMENSION(:,:,:), Allocatable :: Lattice
INTEGER, DIMENSION(:,:), Allocatable :: List, LinkerCon
INTEGER, DIMENSION(:), Allocatable :: BoundTo, ProteinToType
Logical, Dimension(:), Allocatable :: CanRot
INTEGER :: NL(1:26,3), Rej(1:8,1:4), Propose(1:8),EM,EC,ETC,VkMax, MCmove
INTEGER :: mVk
REAL :: pi=3.14159265


! Local Book Keeping
Logical :: flag
INTEGER :: t,t2,m,i,i1,i2,i3,i4,i5,len,pot, CoreN, CoreHalfN
INTEGER :: j, j2, j3, dx1(1:3),dx2(1:3),dx(1:3),dx_max
INTEGER :: MPT1, MCN, MS1, MS2, MS,MST1,MST2,MP,MPT!, MS_local, MS_global
INTEGER :: x1(1:3),x2(1:3),x3(1:3),Move(1:3),R(1:3),ForwardRMN,BackRMN,RotMoves(1:27),RotMove
INTEGER :: AllHostNN, SlitherTimes, SlitherForward,NB_Cluster
INTEGER, DIMENSION(:,:), Allocatable :: Xall,XProt,AllHostN, SlitherList
INTEGER, DIMENSION(:), Allocatable :: Lall,Neighbor,AllHost,FreeHost, SlitherRot, Fisher(:), Core, CoreHalf
integer :: MS_front, MS_back, FisherN, FreeHostN
REAL :: logForwardRMN, logBackRMN,y1(3), y2(3), y3(3), GT(3,3), b, c, d, p, q, theta, phi, rr
REAL :: rr1, rr3(3), COM_All(3), COM_Cluster(3),leq
INTEGER :: tclock1,tclock2,tclock3,tclock4,tclock_sanity,iclock
LOGICAL, DIMENSION(:), Allocatable :: Analysis_log, IsHost,IsFull,IsHalfFull
INTEGER :: iMCMoves(1:8)
INTEGER, DIMENSION(:,:,:), Allocatable :: RadDistCluster,RadDistCore,RadDistCoreHalf
INTEGER, DIMENSION(:,:,:), Allocatable :: XDist_Clusterij,YDist_Clusterij,ZDist_Clusterij
INTEGER, DIMENSION(:,:), Allocatable :: BondTypes,XDist_All,YDist_All,ZDist_All,XDist_Cluster,YDist_Cluster,ZDist_Cluster
INTEGER, DIMENSION(:), Allocatable :: BoundSites,ClusterHist,LargestCluster
REAL, DIMENSION(:), Allocatable :: RG_All,RG_Cluster,RG_Core,RG_CoreHalf,Occupancy,SLVec
INTEGER, DIMENSION(:), Allocatable :: ReadImport
REAL, DIMENSION(:,:), Allocatable :: COMi,SLEq_temp

INTEGER :: DebugInt(1:10)
REAL :: DebugReal(1:10)

!!!
!!!

t=0
Alph='abcdefghijklmnopqrstuvwxyz'
DebugInt=0
DebugReal=0
call system_clock(tclock1)
call ReadKey
write(*,*) 'Tmax:'
write(*,*) Tmax
call BasicSetup
iMCMoves=[iMoveTransI,iMoveTransII,iMoveRot,iMoveSlitherI,iMoveSlitherII,iMoveClusterI,iMoveClusterII,iMoveGrand]
do i=2,8
iMCMoves(i)=iMCMoves(i-1)+iMCMoves(i)
enddo
write(*,*) 'iMCMoves:'
write(*,*) iMCMoves

call InitialConditions

write(*,*) 'Sanity check before starting'
call Sanity

write(*,*) 'Analysis of initial conditions'
call Analysis

write(*,*) 'All looks good!'
call system_clock(tclock2)
write(*,'(A, I17)')   'Time Spent Preparing Simulation (milliseconds):', tclock2-tclock1
write(*,*) 'Starting Simulation'


call system_clock(tclock1)

do t=1,Tmax
!Nucleartion process, if there is a nucleation time set as Trelax
do t2=1,Ti
  if (t .le. Trelax) then
  MCmove=randi(iMCMoves(5))-1
  else
  MCmove=randi(iMCMoves(8))-1
  endif
  !write(*,*) MCmove
  if (sum(Nk)>0 .or. MCmove>=iMCMoves(7)) then !
  if (MCmove<iMCmoves(1)) then
   call MoveTransI
  elseif (MCmove<iMCmoves(2)) then
   call MoveTransII
  elseif (MCmove<iMCmoves(3)) then
   call MoveRot
  elseif (MCmove<iMCmoves(4)) then
   call MoveSlitherI
  elseif (MCmove<iMCMoves(5)) then
    call MoveSlitherII
  elseif (MCmove<iMCMoves(6)) then
    call MoveClusterI
  elseif (MCmove<iMCMoves(7)) then
    call MoveClusterII
  elseif (MCmove<iMCMoves(8)) then
    call MoveGrand
  endif
  endif
enddo


call Analysis


enddo
call EndofProgram

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
endsubroutine FindParm
!!!!
!!!!
!!!!
!!!!
subroutine ReadKey
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
Allocate( Vk(1:NkN) )
Allocate( Vak(1:NkN) )
Allocate( ReadImport(1:NkN) )
Allocate( StartDropletNk(1:NkN) )
Allocate( StartSlabNk(1:NkN) )
call Findparm('Nk')
read(Parm,*) Nk
call Findparm('Vk')
read(Parm,*) Vk
VkMax=maxval(Vk)
call Findparm('Vak')
read(Parm,*) Vak
call Findparm('LSteric')
read(Parm,*) LSteric
call Findparm('ESelf')
read(Parm,*) ESelf

call Findparm('EnergyFile')
PotFile=Parm
call Findparm('ConnectionFile')
ConFile=Parm


call Findparm('StartImported')
read(Parm,*) StartImported
call Findparm('StartDroplet')
read(Parm,*) StartDroplet
call Findparm('StartSlab')
read(Parm,*) StartSlab
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

call Findparm('iMoveRot')
read(Parm,*) iMoveRot
call Findparm('iMoveTransI')
read(Parm,*) iMoveTransI
call Findparm('iMoveTransII')
read(Parm,*) iMoveTransII
call Findparm('iMoveSlitherI')
read(Parm,*) iMoveSlitherI
call Findparm('iMoveSlitherII')
read(Parm,*) iMoveSlitherII
call Findparm('iMoveClusterI')
read(Parm,*) iMoveClusterI
call Findparm('iMoveClusterII')
read(Parm,*) iMoveClusterII
call Findparm('UseGrand')
read(Parm,*) UseGrand
if (UseGrand) then
    call Findparm('iMoveGrand')
    read(Parm,*) iMoveGrand
    allocate( ChemPot(1:NkN) )
    call Findparm('ChemPot')
    read(Parm,*) ChemPot
    call Findparm('MaxSumNk')
    read(Parm,*) mNk
else
    iMoveGrand=0
    mNk=sum(Nk)
endif

call Findparm('E_Sanity')
read(Parm,*) E_Sanity
if (E_Sanity>Tmax) E_Sanity=0
call Findparm('E_Lattice')
read(Parm,*) E_Lattice
if (E_Lattice>Tmax) E_Lattice=0
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
call Findparm('E_Energy')
read(Parm,*) E_Energy
if (E_Energy>Tmax) E_Energy=0
call Findparm('E_PrintStep')
read(Parm,*) E_PrintStep
if (E_PrintStep>Tmax) E_PrintStep=0
call Findparm('E_XYZDist_All')
read(Parm,*) E_XYZDist_All
if (E_XYZDist_All>Tmax) E_XYZDist_All=0
call Findparm('E_RadDist_Cluster')
read(Parm,*) E_RadDist_Cluster
if (E_RadDist_Cluster>Tmax) E_RadDist_Cluster=0
if (StartDroplet .eqv. .false. ) E_RadDist_Cluster=0
call Findparm('E_XYZDist_Cluster')
read(Parm,*) E_XYZDist_Cluster
if (E_XYZDist_Cluster>Tmax) E_XYZDist_Cluster=0
call Findparm('WriteIndiDist_Cluster')
read(Parm,*) WriteIndiDist_Cluster
If(WriteIndiDist_Cluster) then  ! Allocate memory
  ALLOCATE(COMi(1:NkN,3))
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
  ALLOCATE(Core(1:mNk))
  ALLOCATE(CoreHalf(1:mNk)) ! CoreHalf is for half-occupied protein
  ALLOCATE(IsFull(1:mNk))
  ALLOCATE(IsHalfFull(1:mNk))
  ALLOCATE(Occupancy(1:mNk))
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
!! Import Energy
write(*,*) 'Importing Energy'
Allocate(Energy(1:sum(Vk),1:sum(Vk)))
Energy=0
open(1,file=Potfile)
do i=1,sum(Vk)
read(1,*) Energy(i,:)
end do
close(1)
!sanity for energy
do i1=1,sum(Vk)
do i2=1,sum(Vk)
if (Energy(i1,i2)/=Energy(i2,i1)) then
write(*,*) 'Energy Not Symmetric'
write(*,*) 'I don`t know what to do with myself without a symmetric energy table!'
write(*,*) i1,i2
write(*,*) Energy(i1,i2),Energy(i2,i1)
stop
endif
enddo
enddo
open(1,file=ConFile)

! Import Connection
Allocate(ConnectionN(1:sum(Vk)))! ConnectionN[i] is the number of connections of i
Allocate(Connectioni(1:sum(Vk)))! Connectioni[i] is the number of connections of i with i<j
Allocate(Connectionj(1:sum(Vk)))! Connectionj[i] is used to count how many connections has been considered for bead i
write(*,*) 'Importing ConnectionN, Connection, and Linker Length'
ConnectionN=0
Connectioni=0
do
  read(1,*,iostat=m) i,j,len
  if (m /= 0) exit
  if(i .ge. j) then
    write(*,*) "Error in Connection File, index i must be smaller than j"
    stop
  endif
  Connectioni(i)=Connectioni(i)+1
  ConnectionN(i)=ConnectionN(i)+1
  ConnectionN(j)=ConnectionN(j)+1
end do
Allocate(Connection(sum(Vk),maxval(ConnectionN)))
Allocate(LL_temp(sum(Vk),maxval(ConnectionN)))
rewind(1)
Connection=0
LL_temp=0
connectionj=0
! Read in the data
do i1= 1, sum(Vk)
do i2=1,Connectioni(i1)
  read(1, *) i,j,len
  Connectionj(i)=Connectionj(i)+1
  Connection(i,Connectionj(i))=j
  LL_temp(i,Connectionj(i))=len
  Connectionj(j)=Connectionj(j)+1
  Connection(j,Connectionj(j))=i
  LL_temp(j,Connectionj(j))=len
end do
end do
!sanity check
do i=1,sum(Vk)
if (Connectionj(i) /= ConnectionN(i)) then
  write(*,*) "Check your Connection file again, there must be error in ", i, "connection"
  stop
endif
end do
print *, "Finished Connection Importation"
! print *, "ConnectionN:"
! do i1=1,sum(Vk)
!   print *, ConnectionN(i1)
! enddo
!   print *, "Connection:"
! do i1=1,sum(Vk)
!   print *, Connection(i1,:)
! enddo
!   print *, "Linker Length:"
! do i1=1,sum(Vk)
!   print *, LL_temp(i1,:)
!enddo


!Import SLCon, SLEq and SLPot
call Findparm('UseSpringLinkerEq')
read(Parm,*) UseSpringLinkerEq
if (UseSpringLinkerEq) then
  write(*,*) 'Using Equilibrium Linker Length Potential'
  call Findparm('SpringLinkerEqFile')
  SLEqFile=Parm
  open(1,file=SLEqFile)
  ! Import Connection
  Allocate(SLConN(1:sum(Vk)))! SLConN[i] is the number of SL connections of i
  Allocate(SLConi(1:sum(Vk)))! SLConi[i] is the number of SL connections of i with i<j
  Allocate(SLConj(1:sum(Vk)))! SLConj[i] is used to count how many SL connections has been considered for bead i
  write(*,*) 'Importing Equilibrium Spring Linker ConnectionN, ConnectionL, Length and Potential'
  SLConN=0
  SLConi=0
  do
    read(1,*,iostat=m) i,j,leq,pot
    if (m /= 0) exit
    if(i .ge. j) then
      write(*,*) "Error in Spring Linker Connection File, index i must be smaller than j"
      stop
    endif
    SLConi(i)=SLConi(i)+1
    SLConN(i)=SLConN(i)+1
    SLConN(j)=SLConN(j)+1
  end do
  Allocate(SLCon(sum(Vk),maxval(SLConN)))
  Allocate(SLEq_temp(sum(Vk),maxval(SLConN)))
  Allocate(SLPot_temp(sum(Vk),maxval(SLConN)))
  rewind(1)
  SLCon=0
  SLEq_temp=0.0
  SLConj=0
  ! Read in the data
  do i1= 1, sum(Vk)
  do i2=1,SLConi(i1)
    read(1, *) i,j,leq,pot
    SLConj(i)=SLConj(i)+1
    SLCon(i,SLConj(i))=j
    SLEq_temp(i,SLConj(i))=leq
    SLPot_temp(i,SLConj(i))=pot
    SLConj(j)=SLConj(j)+1
    SLCon(j,SLConj(j))=i
    SLEq_temp(j,SLConj(j))=leq
    SLPot_temp(j,SLConj(j))=pot
  end do
  end do
  !sanity check
  do i=1,sum(Vk)
  if (SLConj(i) /= SLConN(i)) then
    write(*,*) "Check your SpringLinkerEqFile again, there must be error in ", i, "Spring Linker"
    stop
  endif
  end do
  close(1)
  print *,   "Finished Equilibrium Spring Linker Importation"
  print *, "Spring Linker ConnectionN:"
  do i1=1,sum(Vk)
  print *, SLConN(i1)
  enddo
  print *, "Spring Linker Connection:"
  do i1=1,sum(Vk)
    print *, SLCon(i1,:)
  enddo
  print *, "Spring Linker Equilibrium Length:"
  do i1=1,sum(Vk)
    print *, SLEq_temp(i1,:)
  enddo
  print *, "Spring Linker Length Potential:"
  do i1=1,sum(Vk)
    print *, SLPot_temp(i1,:)
  enddo

  !Use analysis for Spring linker or not
  if (E_SpringLinkerVector /= 0) then
    write(*,*) 'Write Spring Linker Vector Information'
  endif
else
    write(*,*) 'Not Using Equilibrium Spring Linker Length Potential'
endif


!check all the covalent bonds
do i1=1,sum(Vk)
  do i2=1,ConnectionN(i1)
    flag=.false.
    do i3=1,ConnectionN(Connection(i1,i2))
      if (Connection(Connection(i1,i2),i3).eq.i1 .and. LL_temp(Connection(i1,i2),i3).eq.LL_temp(i1,i2)) then
        flag=.true.
      end if
    end do
    if (flag.eqv..false.) then
      write(*,*) 'You fucked up one of the connection or linker files or listed the wrong valence dumbass.'
      write(*,*) ConnectionN
      write(*,*) Vk
      stop
    end if
  end do
end do
end subroutine ReadKey
!!!!
!!!!
!!!!
!!!!
subroutine BasicSetup
write(*,*) 'Starting Basic Setup'
write(*,*) 'Box Size is', BoxSize
!Check if rotmove can work, turns the move type off if there is nothing that can rotate
flag=.false.
do i1=1,sum(Vk)
  do i2=i1,sum(Vk)
    if (Energy(i1,i2)/=0) then
      flag=.true.
      exit
    endif
  enddo
  if (flag) then
    exit
  endif
enddo
if (flag.eqv..false. .and. iMoveRot>0) then
  write(*,*) 'You are not allowed to have rotation moves if nothing has rotation states.'
  write(*,*) 'Readjusting the move set from: ', iMoveRot
  iMoveRot=0
  write(*,*) 'To: ', iMoveRot
endif
!! SlitherI Key words  (remove some from one end and add them to the other end) (Must be a homopolymer)
if (iMoveSlitherI>0) then
  UseSlitherI=.true.
  Allocate ( UseSlitherI_List(1:NkN) )
  UseSlitherI_List=0
  UseSlitherI_ListN=0
  do i1=1,NkN !check each protein type individually if it can slither 1
    if (Vk(i1)>1) then !must be a polymer
      flag=.true.
      if (ConnectionN(1+sum(Vk(1:i1-1))).gt.1) flag=.false.  !is the first bead in the list bound to more than one?
      if (ConnectionN(sum(Vk(1:i1))).gt.1) flag=.false.      !is the last bead in the list bound to more than one?
      do i2=2,Vk(i1)-1 !run through all valences checking the linkers lengths and number
        if (ConnectionN(i2+sum(Vk(1:i1-1))).ne.2) flag=.false. !Cant do slither moves on branched polymers
        if (Connection(i2+sum(Vk(1:i1-1)),1).ne.sum(Vk(1:i1-1))+i2-1) flag=.false.
        if (Connection(i2+sum(Vk(1:i1-1)),2).ne.sum(Vk(1:i1-1))+i2+1) flag=.false.
        if (LL_temp(i2+sum(Vk(1:i1-1)),1)/=LL_temp(i2+sum(Vk(1:i1-1)),2)) then !check that the linkers are identical.  
          flag=.false.
        endif
      end do
      do i2=1,Vk(i1)-1
        if (any(Energy(sum(Vk(1:i1-1))+i2,:)/=Energy(sum(Vk(1:i1-1))+i2+1,:))) then
          flag=.false.
        endif
        do i3=sum(Vk(1:i1-1))+1,sum(Vk(1:i1)) !Checks that there aren't in-cis interactions This is a requirement because I was lazy in coding the slither move for how it deals with reasigning the interactions after a move
            if (Energy(sum(Vk(1:i1-1))+i2,i3)/=0) then
                flag=.false.
            end if
        end do
      end do
      if (flag) then
        UseSlitherI_ListN=UseSlitherI_ListN+1
        UseSlitherI_List(UseSlitherI_ListN)=i1
      endif
    end if
  end do !finished checking how many slither proteins there are
  if (UseSlitherI_ListN==0) then
    UseSlitherI=.false.
    write(*,*) 'SlithermovesI arent appropriate for this system, turning them off'
    write(*,*) 'These moves are appropriate only for homopolymers'
    write(*,*) 'Because of how I coded it, I cant deal with homopolymers that can interact with themselves =('
    write(*,*) 'Some day Ill fix this.'
    write(*,*) 'Move Probability old:',iMoveSlitherI
    iMoveSlitherI=0
    write(*,*) 'Move Probability new:',iMoveSlitherI
  endif
else
  UseSlitherI=.false.
end if

if ((iMoveSlitherI>0) .or. (iMoveSlitherII>0)) then
  Allocate ( SlitherList(1:maxval(Vk),3) )
endif

!! SlitherII Key Words
if (iMoveSlitherII>0) then
  UseSlitherII=.true.
  Allocate ( UseSlitherII_List(1:NkN) )
  Allocate ( SlitherRot(1:maxval(Vk)) )
  UseSlitherII_List=0
  UseSlitherII_ListN=0
  do i1=1,NkN !check each protein individually
    if (Vk(i1)>1) then !must be a linear polymer with identical linker lengths
      flag=.true.
      if (ConnectionN(1+sum(Vk(1:i1-1))).gt.1) flag=.false.
      if (ConnectionN(sum(Vk(1:i1))).gt.1) flag=.false.
      do i2=2,Vk(i1)-1 !run through all valences
        if (ConnectionN(i2+sum(Vk(1:i1-1))).ne.2) flag=.false. !Cant do slither moves on branched polymers
        if (LL_temp(i2+sum(Vk(1:i1-1)),1)/=LL_temp(i2+sum(Vk(1:i1-1)),2)) then !check that the linkers are identical.  
        flag=.false.
        endif
      end do
      if (flag) then
        UseSlitherII_ListN=UseSlitherII_ListN+1
        UseSlitherII_List(UseSlitherII_ListN)=i1
      endif
    endif
  enddo
  if (UseSlitherII_ListN==0) then
    UseSlitherII=.false.
    write(*,*) 'SlithermovesII arent appropriate for this system, turning them off'
    write(*,*) 'Move Probability old:',iMoveSlitherII
    iMoveSlitherII=0
    write(*,*) 'Move Probability new:',iMoveSLitherII
  endif
else
  UseSlitherII=.false.
end if
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

if (sum(Nk)>mNk) then
  write(*,*) 'Not enough room to fit all the proteins youve asked for in MaxSumNk!'
  write(*,*) 'Crashing'
  stop
endif
mVk=maxval(Vk)
! write(*,*) 'Max Vk:',mVk,' Min Vk:',minval(Vk)


if (StartSlab) then
  if(StartImported)then
    write(*,*) 'You can not both start imported and slab, turn off one of them'
    stop
  endif
  if (BoxSize(3)<2*BoxSize(1) .or. BoxSize(3)<2*BoxSize(1)) then
    write(*,*) 'Your size of slab direction is not long enough, make it longer'
    stop
  endif

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
    write(*,*) 'Your box is not a cube, you better are starting Slab.'
  endif
endif

!allocate stuff
ALLOCATE ( XProt(1:mVk,1:3) )
ALLOCATE ( Xall(1:maxval(ConnectionN)*2,1:3) )
ALLOCATE ( Lall(1:maxval(ConnectionN)*2) )
ALLOCATE ( Neighbor(1:maxval(ConnectionN)) )
ALLOCATE ( CanRot(sum(Vk)) )
ALLOCATE ( Lattice(1:BoxSize(1),1:BoxSize(2),1:BoxSize(3)) )
ALLOCATE ( List(1:mVk*mNk,1:3) )
ALLOCATE ( BoundTo(1:mVk*mNk) )
ALLOCATE ( ProteinToType(1:mNk) )
ALLOCATE ( Analysis_log(1:mNk) )
ALLOCATE ( AllHostN(1:mNk,1:2) )
ALLOCATE ( AllHost(1:mNk) )
ALLOCATE ( FreeHost(1:mNk) )
ALLOCATE ( Fisher(1:mNk) )
ALLOCATE ( IsHost(1:mNk) )



tclock_sanity=0
Lattice=0
List=0
BoundTo=0
!neighbor list
i=1
do i1=-1,1
do i2=-1,1
do i3=-1,1
if ((i1==0 .and. i2==0 .and. i3==0)) then
  cycle
endif
  NL(i,:)=[i1,i2,i3]
  i=i+1
enddo
enddo
enddo
! Can a bead rotate
CanRot=.false.
do i1=1,sum(Vk)
 do i2=1,sum(Vk)
  if (Energy(i1,i2)/=0) then
   CanRot(i1)=.true.
  endif
 enddo
enddo

Propose=0
Rej=0
RotMoves=0
ETC=0
endsubroutine BasicSetup
!!!!
!!!!
!!!!
!!!!
subroutine InitialConditions !!!!!! This is massively memory inefficienct.  I should un-build all invalid proteins instead of using temp lattices and lists
INTEGER :: r2,r1a
INTEGER, DIMENSION(:,:,:), Allocatable :: Lattice_temp
INTEGER, DIMENSION(:,:), Allocatable :: List_temp

ALLOCATE ( Lattice_temp(1:BoxSize(1),1:BoxSize(2),1:BoxSize(3)) )
ALLOCATE ( List_temp(1:mVk*mNk,1:3) )
r2=0

if (StartImported) then
  write(*,*) 'Importing Initial Conditions'
  write(*,*) 'Be careful!'
  
  open(unit=1,file='StartParms.txt')
  read(1,*) x1
  do i1=1,3
    if (x1(i1).ne.BoxSize(i1)) then
      write(*,*) 'New BoxSize (',BoxSize,') doesnt match the old BoxSize (',x1,')'
      stop
    endif
  enddo
  read(1,*) i1
  if (i1.ne.NkN) then
    write(*,*) 'New NkN (',NkN,') doesnt match the old NkN (',i1,')'
    stop
  endif
  read(1,*) ReadImport
  if (any(ReadImport.ne.Nk)) then
    write(*,*) 'New Nk (',Nk,') doesnt match the old Nk (',ReadImport,')'
    stop
  endif
  read(1,*) ReadImport
  if (any(ReadImport.ne.Vk)) then
    write(*,*) 'New Nk (',Nk,') doesnt match the old Nk (',ReadImport,')'
    stop
  endif

  do i1=1,sum(Vk)
    read(1,*) i2
    if (i2.ne.ConnectionN(i1)) then
      write(*,*) 'New ConnectionN (',ConnectionN(i1),') doesnt match the old ConnectionN (',i2,')'
      stop
    endif
    read(1,*) Lall(1:ConnectionN(i1))
    if (any(Lall(1:ConnectionN(i1)).ne.Connection(i1,1:ConnectionN(i1)))) then
      write(*,*) 'New Connection (',Connection(i1,1:ConnectionN(i1)),&
        ') doesnt match the old Connection (',Lall(1:ConnectionN(i1)),')'
      stop
    endif
    read(1,*) Lall(1:ConnectionN(i1))
    if (any(Lall(1:ConnectionN(i1)).ne.LL_temp(i1,1:ConnectionN(i1)))) then
      write(*,*) 'New LL_temp (',LL_temp(i1,1:ConnectionN(i1)), &
        ') doesnt match the old LL_temp (',Lall(1:ConnectionN(i1)),')'
      stop
    endif
  enddo
  close(1)
  
  open(unit=1,file='StartListForm.txt')
  open(unit=2,file='StartBoundTo.txt')
  open(unit=3,file='StartProteinToType.txt')
  do i1=1,sum(Nk)
    read(3,*) ProteinToType(i1)
    do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
      read(2,*) BoundTo(i2)
      read(1,*) List(i2,:)
    enddo
  enddo
  close(1)
  close(2)
  write(*,*) 'finished importing?'
  do i1=1,sum(Nk)
  do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
    Lattice(List(i2,1),List(i2,2),List(i2,3))=i2
    if (BoundTo(i2)>i2) then
      i=BeadToType(i2)
      j=BeadToType(BoundTo(i2))
      ETC=ETC+Energy(i,j)
      if (BeadToProtein(i2)==BeadToProtein(BoundTo(i2))) then
        ETC=ETC+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
      endif
    endif
  enddo
  enddo
  write(*,*) 'finished importing2?'

else 
  if (StartSlab) then
    write(*,*) 'Inventing initial conditions (slab), assuming the third dimension is slab direction'
    do i1=1,NkN
    do i2=1,Nk(i1)
      do r1a=1,10000  !number of attempts to place a specific protein
        x1=randiBS(BoxSize)
        if (i2<=StartSlabNk(i1)) then !Put it into the dense phase
            x1(3)=randi(StartSlabH)+(BoxSize(3)-StartSlabH)/2
        endif
        XProt=AddProtein(i1)!Protein coordinates

        flag=.true.
        
        do i3=1,Vk(i1)
          XProt(i3,:)=PB3(XProt(i3,:)+x1)
          if (Lattice(XProt(i3,1),XProt(i3,2),XProt(i3,3)) /= 0) then
            flag=.false. !this checks if it overlaps with something else
          endif
          do i4=1,i3-1  
          if (all(XProt(i3,:)==XProt(i4,:))) then
            flag=.false.!this checks if it overlaps with itself
          endif
          enddo
        enddo
        if (flag) then
          exit  !successfully built, add to the lattice.
        endif
      enddo !end of loop to try a protein until it fits
      if (flag) then
        do i3=1,Vk(i1)
          Lattice(XProt(i3,1),XProt(i3,2),XProt(i3,3))=(i2+sum(Nk(1:i1-1))-1)*mVk + i3
          List((i2+sum(Nk(1:i1-1))-1)*mVk + i3,:)=XProt(i3,:)
          ProteinToType(i2+sum(Nk(1:i1-1)))=i1
        enddo
      else
        if (i2<=StartSlabNk(i1)) then
          write(*,*) 'Crashed on forming initial conditions (slab).'
        else
          write(*,*) 'Crashed on forming initial conditions (dilute).'
        endif
        write(*,*) 'Consider adding more room or trying with a different seed'
        write(*,*) 'Added ',i2+sum(Nk(1:i1-1)),' proteins out of ',sum(Nk)
        stop
      endif
    enddo
    enddo

  else
    if (all(StartDropletNk==0)) then
      write(*,*) 'Inventing initial conditions (dilute)'
    else
      write(*,*) 'Inventing initial conditions (dense)'
    endif
    do i1=1,NkN
    do i2=1,Nk(i1)
      do r1a=1,10000 !number of attempts to place a specific protein
        if (i2<=StartDropletNk(i1)) then !Put it into the dense phase
            call random_number(rr3)
            theta=rr3(1)*2.0*pi
            phi = acos(rr3(2)*2.0 - 1.0)
            rr=rr3(3)**(1.0/3.0)*StartDropletR
            x1=Floor([rr*sin(phi)*cos(theta),rr*sin(phi)*sin(theta),rr*cos(phi)]+real(BoxSize)/2.0+.5)
        else ! else puts it in the bulk phase
            x1=randiBS(BoxSize)
        endif
        XProt=AddProtein(i1)
        flag=.true.
        
        do i3=1,Vk(i1)
          XProt(i3,:)=PB3(XProt(i3,:)+x1)
          if (Lattice(XProt(i3,1),XProt(i3,2),XProt(i3,3)) /= 0) then
            flag=.false. !this checks if it overlaps with something else
          endif
          do i4=1,i3-1  
          if (all(XProt(i3,:)==XProt(i4,:))) then
            flag=.false.!this checks if it overlaps with itself
          endif
          ! write(*,*) ' Crashed at adding the',i2+sum(Nk(1:i1-1)),' protein out of ',sum(Nk)
          ! write(*,*) 'Crashed at adding the',i3,' bead out of ',Vk(i1)
          enddo
        enddo
        if (flag) then
          exit  !successfully built, add to the lattice.
        endif
      enddo !end of loop to try a protein until it fits
      if (flag) then
        do i3=1,Vk(i1)
          Lattice(XProt(i3,1),XProt(i3,2),XProt(i3,3))=(i2+sum(Nk(1:i1-1))-1)*mVk + i3
          List((i2+sum(Nk(1:i1-1))-1)*mVk + i3,:)=XProt(i3,:)
          ProteinToType(i2+sum(Nk(1:i1-1)))=i1
        enddo
      else
        if (i2<=StartDropletNk(i1)) then
          write(*,*) 'Crashed on forming initial conditions (dense).'
        else
          write(*,*) 'Crashed on forming initial conditions (dilute).'
        endif
        write(*,*) 'Consider adding more room or trying with a different seed'
        write(*,*) 'Added ',i2+sum(Nk(1:i1-1)),' proteins out of ',sum(Nk)
        write(*,*) ''
        stop
      endif
    enddo
    enddo
  endif
  
  write(*,*) 'Forming random bonds to stabilize dense phase'
  MS=0
  do i1=1,NkN
  do i2=1,Nk(i1)
  do i3=1,Vk(i1)
  MS=MS+1
  MST1= BeadToType(MS)
  if (BoundTo(MS)==0) then
    ForwardRMN=1
    do i4=1,26 !lets look around to see who it could interacting with
      R=PB3(List(MS,:)+NL(i4,:))
      if (Lattice(R(1),R(2),R(3))/=0) then
      if (Energy(MST1,BeadToType(Lattice(R(1),R(2),R(3))))/=0 .and. & ! is sticky
          BoundTo(Lattice(R(1),R(2),R(3)))==0) then ! is single and ready to mingle
          ForwardRMN=ForwardRMN+1
          RotMoves(ForwardRMN)=Lattice(R(1),R(2),R(3))
      endif
      endif
    enddo
    if (ForwardRMN>1) then
      RotMove=RotMoves(randi(ForwardRMN))
      if (RotMove.gt.0) then
        BoundTo(MS)=RotMove
        BoundTo(RotMove)=MS
        i=BeadToType(MS)
        j=BeadToType(RotMove)
        ETC=ETC+Energy(i,j)
        if (BeadToProtein(MS)==BeadToProtein(RotMove)) then
          ETC=ETC+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
        endif
      endif
    endif
  endif
  enddo
  enddo
  enddo
endif !end of choice between compact and dilute

if (UseSpringLinkerEq) then
do i1=1,sum(Nk) ! Spring linker length potential calculation
  do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
    MST1=BeadToType(i2)
    do i4=1, SLconN(MST1)
      if (SLcon(MST1,i4)>MST1) then
        x1=List(i2,:)
        x2=List(i2-MST1+SLcon(MST1,i4),:)
        dx=DistPB(x1,x2)
        dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1,i4))!3d distance differnce to Eq Len for MS1
        ETC=ETC + SLPot_temp(MST1,i4)*dx_max**2
      endif
    enddo
  enddo
enddo
endif
write(*,*) 'Starting Energy:', ETC

! Open all the export files
if (E_Lattice/=0) then
!  open(unit=8,file='AListForm.r4',form='unformatted',access='direct',recl=4*3*sum(Vk*Nk))
!  open(unit=8,file='ABoundTo.r4',form='unformatted',access='direct',recl=4*sum(Vk*Nk))
  open(unit=1,file='AListForm.xyz')
  open(unit=2,file='ABoundTo.txt')
  if (UseGrand) then
    open(unit=3,file='AType.txt')
  endif
  
  open(4,file='ABonds.vmd')
  write(4,'(A)')'topo clearbonds'
  i4=0
  do i1=1,NkN
    write(4,'(A,i6)') 'set NkN', i1
    write(4,'(A,i6)') 'set Vk', Vk(i1)
    write(4,'(A,i6)') 'set end', sum(Nk(1:i1-1)*Vk(1:i1-1))+Nk(i1)*Vk(i1)
    do i2=sum(Vk(1:i1-1))+1,sum(Vk(1:i1-1))+Vk(i1)
      do i3=1,ConnectionN(i2)
        if(Connection(i2,i3)>i2) then
          write(4,'(A,i6)') 'set current', sum(Nk(1:i1-1)*Vk(1:i1-1))+i2-1
          write(4,'(A,i6)') 'set partner', Connection(i2, i3)-i2 !defined as distance to current
          write(*,*) Connection(i2, i3), i2, Connection(i2, i3)-i2
          write(4,*) 'for {set j $current} {$j < $end} {incr j $Vk} {'
          write(4,*) '         topo addbond $j [expr {$j + $partner}]'
          write(4,*) '}'
        endif
      enddo 
    enddo
  enddo
  close(4)
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
  open(unit=97,file='ASelfLoop.txt')
endif


if (E_Network/=0) then
  do i1=1, NkN! assuming only one linker
    if (Vak(i1)==2) then
      LinkerType=i1
    endif
  enddo
  write(*,*) 'Network loop analysis is on, assuming only one Linker which is protein type', LinkerType
  write(*,*) 'There are', Nk(LinkerType), 'Linkers'
  open(unit=98,file='ALinkerConnection.txt')
  ALLOCATE ( LinkerCon(1:Nk(LinkerType),3) )! first one for the linker index, then each linker could have two connections
  LinkerCon=0
endif

if(UseSpringLinkerEq) then
if (E_SpringLinkerVector/=0) then
  open(unit=99,file='ASpringLinkerVector.txt')
  ALLOCATE ( SLVec(1:7) )
  SLVec=0
endif
endif
if (E_Energy/=0) then
  open(unit=101,file='AEnergy.txt')
endif
if (E_LargestCluster/=0) then
  open(unit=102,file='ALargestCluster.txt')
  ALLOCATE ( LargestCluster(1:5*(NkN+1)) )
  LargestCluster=0
endif
if (E_RG_All/=0) then
  open(unit=103,file='ARG_All.txt')
  ALLOCATE ( RG_All(1:4) )
  RG_All=0
endif
if (E_RG_Cluster/=0) then
  open(unit=107,file='ARG_Cluster.txt')
  ALLOCATE ( RG_Cluster(1:8) )
  RG_Cluster=0
  if(WriteCoreCluster) then
    open(unit=108,file='ARG_Core.txt')
  ALLOCATE ( RG_Core(1:8) )
  RG_Core=0
    open(unit=109,file='ARG_CoreHalf.txt')
  ALLOCATE ( RG_CoreHalf(1:8) )
  RG_CoreHalf=0
  endif
endif
if (E_ClusterHist/=0) then
  open(unit=104,file='AClusterHist.txt')
  ALLOCATE ( ClusterHist(1:21) )
  ClusterHist=0
endif
if (E_BoundSites/=0) then
  open(unit=105,file='ABoundSites.txt')
  ALLOCATE ( BoundSites(1:sum(Vk)) )
  BoundSites=0
endif
if (E_BondTypes/=0) then
  open(unit=106,file='ABondTypes.txt')
  ALLOCATE ( BondTypes(1:sum(Vk),1:sum(Vk)) )
  BondTypes=0
endif

if (E_XYZDist_All/=0) then
  do i1=1,NkN
    write(Str1,*) i1
    do i2=1,3
    Str3='AXYZDist_All_Protein'//trim(adjustl(Str1))//'_'//CHAR(64+23+i2)//'.txt'
    open(unit=200+i1+10*i2,file=Str3)
    enddo
  enddo
  Allocate (XDist_All(1:BoxSize(1),NkN))
  Allocate (YDist_All(1:BoxSize(2),NkN))
  Allocate (ZDist_All(1:BoxSize(3),NkN))
  XDist_All=0
  YDist_All=0
  ZDist_All=0
endif

if (E_XYZDist_Cluster/=0) then
  do i1=1,NkN
    write(Str1,*) i1
    do i3=1,3
    Str3='AXYZDist_Cluster_Protein'//trim(adjustl(Str1))//'_'//CHAR(64+23+i3)//'.txt'
    open(unit=250+i1+10*i3,file=Str3)
    enddo
  enddo
  Allocate (XDist_Cluster(1:BoxSize(1),NkN))
  Allocate (YDist_Cluster(1:BoxSize(2),NkN))
  Allocate (ZDist_Cluster(1:BoxSize(3),NkN))
  XDist_Cluster=0
  YDist_Cluster=0
  ZDist_Cluster=0
  if(WriteIndiDist_Cluster) then
    do i1=1,NkN
    write(Str1,*) i1
    do i2=1,NkN
    write(Str2,*) i2
    do i3=1,3!XYZ axis
    Str3='AXYZDist_Cluster_P'//trim(adjustl(Str1))//'CenteredOnP'//trim(adjustl(Str2))//'_'//CHAR(64+23+i3)//'.txt'
    open(unit=100*i1+10000*i2+i3,file=Str3)
    enddo
    enddo
    enddo
    Allocate (XDist_Clusterij(1:BoxSize(1),NkN,NkN))
    Allocate (YDist_Clusterij(1:BoxSize(2),NkN,NkN))
    Allocate (ZDist_Clusterij(1:BoxSize(3),NkN,NkN))
    YDist_Clusterij=0
    XDist_Clusterij=0
    ZDist_Clusterij=0
  endif
endif

if (E_RadDist_Cluster/=0) then
  do i1=1,NkN
    write(Str1,*) i1
    if(WriteIndiDist_Cluster) then
    do i2=1,NkN
      write(Str2,*) i2
      Str3='AClusterRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOn'//trim(adjustl(Str2))//'.txt'
      open(unit=100*i1+10000*i2,file=Str3)
    enddo
    endif
    Str3='AClusterRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOnLargestCluster.txt'
    open(unit=200+i1,file=Str3)
  enddo
  Allocate ( RadDistCluster(1:floor(sqrt(3.)/2*BoxSize(1))+1,NkN,NkN+1))
  RadDistCluster=0
  if(WriteCoreCluster) then
    !Open Core Radial Hist
    do i1=1,NkN
      write(Str1,*) i1
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(Str2,*) i2
        Str3='ACoreRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOn'//trim(adjustl(Str2))//'.txt'
        open(unit=100*i1+10000*i2+1,file=Str3)
      enddo
      endif
      Str3='ACoreRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOnCore.txt'
      open(unit=200+i1+30,file=Str3)
    enddo
    Allocate ( RadDistCore(1:floor(sqrt(3.)/2*BoxSize(1))+1,NkN,NkN+1))
    RadDistCore=0
    !Open CoreHalf Radial Hist
    do i1=1,NkN
      write(Str1,*) i1
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(Str2,*) i2
        Str3='ACoreHalfRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOn'//trim(adjustl(Str2))//'.txt'
        open(unit=100*i1+10000*i2+2,file=Str3)
      enddo
      endif
      Str3='ACoreHalfRadDist_Protein'//trim(adjustl(Str1))//'_CenteredOnCoreHalf.txt'
      open(unit=200+i1+60,file=Str3)
    enddo
    Allocate ( RadDistCoreHalf(1:floor(sqrt(3.)/2*BoxSize(1))+1,NkN,NkN+1))
    RadDistCore=0
  endif
endif

endsubroutine InitialConditions
!!!!
!!!!
!!!!
!!!!
function randi(j)
INTEGER :: j,randi
call random_number(rr1)
!randi=mod(floor(rr1*j),j)+1
randi=mod(int(rr1*j),j)+sign(1,j)
endfunction randi
!
function randi3(j)
INTEGER :: j,randi3(3)
call random_number(rr3)
!randi3=mod(floor(rr3*j),j)+1
randi3=mod(int(rr3*j),j)+sign(1,j)
endfunction randi3
!
function randiBS(j)
INTEGER :: j(3),randiBS(3)
call random_number(rr3)
!randiBS=mod(floor(rr3*j),j)+1
randiBS=mod(int(rr3*j),j)+sign(1,j)
endfunction randiBS
!
function PB1(j1,i1)
INTEGER :: j1,i1,PB1!j1 is the coordinates, i1 is the dimension for PB
!PB1=mod(j1+BoxSize(i1)-1,BoxSize(i1))+1
PB1=modulo(j1-1,BoxSize(i1))+1
endfunction PB1
!
function PB3(j3)
INTEGER :: j3(3),PB3(3)
!PB3=mod(j3+BoxSize-1,BoxSize)+1
PB3=modulo(j3-1,BoxSize)+1
endfunction PB3
!
function DistPB(x,y)
INTEGER :: x(3),y(3),DistPB(3)
!DistPB=mod(x-y+BoxSize/2+BoxSize,BoxSize)-BoxSize/2
DistPB=modulo(x-y+BoxSize/2,BoxSize)-BoxSize/2
endfunction DistPB
!
function rDistPB(x,y)
REAL :: x(3),y(3),rDistPB(3)
!rDistPB=mod(x-y+BoxSize/2+BoxSize,real(BoxSize))-BoxSize/2
rDistPB=modulo(x-y+BoxSize/2,real(BoxSize))-BoxSize/2
endfunction rDistPB
!
function BeadToProtein(j1)
INTEGER :: j1,BeadToProtein
BeadToProtein=(j1-1)/mVk+1
endfunction BeadToProtein
!
function BeadToType(j1)
INTEGER :: j1, BeadToType !bead type with respect to other beads
j2=ProteinToType(BeadToProtein(j1))
BeadToType = sum(Vk(1:j2-1)) + mod(j1-1,mVk)+1
endfunction BeadToType
!
!!!!
!!!!
!!!!
!!!!

! function AddProtein(MPT) !which protein number, 860
! INTEGER :: AddProtein(1:mVk,1:3)
! INTEGER :: MPT
! AddProtein(1,:)=[0,0,0]
! do j2=2,Vk(MPT)
!   MST1=sum(Vk(1:MPT-1))+j2
!   MCN=ConnectionN(MST1)
!   flag=.true. 
!   do j3=1,MCN
!   if (Connection(MST1,j3)<MST1) then
!     if (flag) then
!       x2=AddProtein(Connection(MST1,j3)-sum(Vk(1:MPT-1)),:)-LL_temp(MST1,j3)
!       x3=AddProtein(Connection(MST1,j3)-sum(Vk(1:MPT-1)),:)+LL_temp(MST1,j3)
!       flag=.false.
!     else!this part is for loop polymers
!       x2=max(x2,AddProtein(MST1-Connection(MST1,j3),:)-LL_temp(MST1,j3))
!       x3=min(x3,AddProtein(MST1-Connection(MST1,j3),:)+LL_temp(MST1,j3))
!     endif
!   endif
!   enddo
!   do j3=1,3
!     AddProtein(j2,j3)=randi(x3(j3)-x2(j3)+1)-1+x2(j3)
!   enddo
! enddo
! endfunction AddProtein

function AddProtein(MPT) !which protein number, 860
INTEGER :: AddProtein(1:mVk,1:3) !protein beads, coordintes
INTEGER :: MPT
AddProtein(1,:)=[0,0,0]
do j2=2,Vk(MPT) !bead number
do j4=1,100 !Attempts
  MST1=sum(Vk(1:MPT-1))+j2
  MCN=ConnectionN(MST1)
  flag=.true.
  do j3=1,MCN
  if (Connection(MST1,j3)<MST1) then
    if (flag) then
      x2=AddProtein(Connection(MST1,j3)-sum(Vk(1:MPT-1)),:)-LL_temp(MST1,j3)
      x3=AddProtein(Connection(MST1,j3)-sum(Vk(1:MPT-1)),:)+LL_temp(MST1,j3)
      flag=.false.
    else
      write(*,*) 'this has to be debuuged'
      stop
      x1=max(x2,AddProtein(MST1-Connection(MST1,j3),:)-LL_temp(MST1,j3))
      x2=min(x3,AddProtein(MST1-Connection(MST1,j3),:)+LL_temp(MST1,j3))
    endif
  endif
  enddo
  ! do j3=1,3
  !   AddProtein(j2,j3)=randi(x3(j3)-x2(j3)+1)-1+x2(j3)
  ! enddo

  AddProtein(j2,:)=PB3(randiBS(x3-x2+1)-1+x2)

  flag=.true.
  do j3=1,j2-1
  if (all(AddProtein(j3,:)==AddProtein(j2,:))) then
    flag=.false.
  endif
  enddo
  if (flag) then
    exit
  endif
enddo
enddo
endfunction AddProtein

!!!!
!!!!
!!!!
!!!!
subroutine MoveTransI
EM=0
EC=0
Propose(1)=Propose(1)+1
MP=randi(sum(Nk))
MPT=ProteinToType(MP)
MS=(MP-1)*mVk+randi(Vk(MPT)) !bead number to be moved
MST1=BeadToType(MS)
MCN=ConnectionN(MST1) !number of linker connections (this is ConnectionN(MST2) or ConnectionN(ListToProtein(MS,3)))

if (MCN==0) then !if it is a monomer domain that isn't tethered
  Move=PB3(List(MS,:)+[randi(5),randi(5),randi(5)]-3)   ! can move from -2 to 2
else
  x1=List(MS,:)
  do i1=1,MCN !run through all of it's fixed bonds
    x2=List(MS-MST1+Connection(MST1,i1),:)
    Xall(i1,:)=DistPB(x2,x1) !this is the list of all displacements
  end do
  do i2=1,3
    i3=maxval(Xall(1:MCN,i2)-LL_temp(MST1,1:MCN)) !minimum displacement
    i4=minval(Xall(1:MCN,i2)+LL_temp(MST1,1:MCN)) !maximum displacement
    Move(i2)=PB1(randi(i4-i3+1)+i3+List(MS,i2)-1,i2) !random final position
  end do
endif
BackRMN=1
ForwardRMN=1
if (Lattice(Move(1),Move(2),Move(3))==0 .or. Lattice(Move(1),Move(2),Move(3))==MS) then !checks for steric clash
  if (CanRot(MST1)) then !this is asking if it can interact with anything at all
    if (BoundTo(MS)>0) then !if it was already bound we need to calculate the energy of that bond
      !EC=EC+Energy(MST1,BeadToType(BoundTo(MS)))
      i=MST1
      j=BeadToType(BoundTo(MS))
      EC=EC+Energy(i,j)
      if (BeadToProtein(MS)==BeadToProtein(BoundTo(MS))) then
        EC=EC+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
      endif

    endif
    do i=1,26 !lets look around to see who it could have been interacting with to count the reverse move probability
      R=PB3(List(MS,:)+NL(i,:))
      if (Lattice(R(1),R(2),R(3))/=0) then
      if (Energy(MST1,BeadToType(Lattice(R(1),R(2),R(3))))/=0 .and. &
          (BoundTo(Lattice(R(1),R(2),R(3)))==0 .or. BoundTo(Lattice(R(1),R(2),R(3)))==MS)) then
        BackRMN=BackRMN+1
      endif
      endif
    enddo
    do i=1,26
      R=PB3(Move+NL(i,:))
      if (Lattice(R(1),R(2),R(3))/=0 .and. Lattice(R(1),R(2),R(3))/=MS) then
      if (Energy(MST1,BeadToType(Lattice(R(1),R(2),R(3))))/=0 .and. &
          (BoundTo(Lattice(R(1),R(2),R(3)))==0 .or. BoundTo(Lattice(R(1),R(2),R(3)))==MS)) then
        ForwardRMN=ForwardRMN+1
        RotMoves(ForwardRMN)=Lattice(R(1),R(2),R(3))
      endif
      endif
    enddo
  endif
  flag=.true. !reports if there is no clash
else
  flag=.false.
endif

! Accept or Reject the move
if (flag) then
  RotMove=RotMoves(randi(ForwardRMN))
  if (RotMove/=0) then
      i=MST1
      j=BeadToType(RotMove)
      EM=EM+Energy(i,j)
      if (BeadToProtein(MS)==BeadToProtein(RotMove)) then
        EM=EM+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
      endif
  endif

  if (UseSpringLinkerEq) then
    do i4=1,SLConN(MST1)
      x1=List(MS,:)
      x2=List(MS-MST1+SLCon(MST1,i4),:)
      dx=DistPB(x1,x2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1,i4))!3d distance differnce to Eq Len for MS1
      EC=EC + SLPot_temp(MST1,i4)*dx_max**2
      dx=DistPB(Move,x2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1,i4))
      EM=EM + SLPot_temp(MST1,i4)*dx_max**2
    enddo
  endif

  call random_number(rr1)
  if (rr1<exp(real(EC-EM)/1000)*ForwardRMN/BackRMN) then
    ETC=ETC+EM-EC
    Lattice(List(MS,1),List(MS,2),List(MS,3))=0
    Lattice(Move(1),Move(2),Move(3))=MS
    List(MS,:)=Move
    if (BoundTo(MS)>0) then
      BoundTo(BoundTo(MS))=0
    endif
    BoundTo(MS)=RotMove
    if (RotMove/=0) then
      BoundTo(RotMove)=MS
    endif
  else
    Rej(1,1)=Rej(1,1)+1
  endif
else
  Rej(1,2)=Rej(1,2)+1
endif
endsubroutine MoveTransI
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
subroutine MoveTransII  !! this section needs to have the functions ported in still
EM=0
EC=0
Propose(2)=Propose(2)+1
MP=randi(sum(Nk))
MPT=ProteinToType(MP)
MS1=(MP-1)*mVk+randi(Vk(MPT)) !bead number to be moved
MS2=BoundTo(MS1)
MST1=BeadToType(MS1)
MCN=0 !Counter for the number of linkers that must be preserved
flag=.true. !if the move should be proposed
x1=List(MS1,:)
! find the distances that are tethered to MS
do i1=1,ConnectionN(MST1) !run through all of it's fixed bonds
if (MS2-MS1 /= Connection(MST1,i1)-MST1) then
  MCN=MCN+1
  x2=List(MS1-MST1+Connection(MST1,i1),:)
  Xall(MCN,:)=DistPB(x2,x1) !this is the list of all displacements
  Lall(MCN)=LL_temp(MST1,i1)
endif
enddo
! find the distances that are tethered to boundto(MS)
if (MS2/=0) then
  MST2=BeadToType(MS2)
  x2=List(MS2,:)
  do i1=1,ConnectionN(MST2) !run through all of it's fixed bonds
  if (MS1-MS2 /= Connection(MST2,i1)-MST2) then
    MCN=MCN+1
    x3=List(MS2-MST2+Connection(MST2,i1),:)
    Xall(MCN,:)=DistPB(x3,x2) !this is the list of all displacements
    Lall(MCN)=LL_temp(MST2,i1)
  endif
  enddo
endif
! find the proposed move
if (MCN==0) then !if it is a monomer domain that isn't tethered
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
  flag=.false. !clash, can't proceed.
endif
if (MS2/=0) then
  x2=PB3(x2+Move)
  if (Lattice(x2(1),x2(2),x2(3))/=0 .and. Lattice(x2(1),x2(2),x2(3))/=MS1) then !checks if boundto(MS) clashes in new position
    flag=.false.
  endif
endif

! When UseSpringLinkerEq==FALSE: Allways Accept the move if flag==true
if (flag) then
  MST2=BeadToType(MS2)
  if (UseSpringLinkerEq) then
    do i4=1,SLConN(MST1)
      dx1=List(MS1,:)
      dx2=List(MS1-MST1+SLCon(MST1,i4),:)
      dx=DistPB(dx1,dx2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1,i4))!differnce to Eq Leng for MS1
      EC=EC + SLPot_temp(MST1,i4)*dx_max**2
      dx=DistPB(x1,dx2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1,i4))!differnce to Eq Leng for MS1
      EM=EM + SLPot_temp(MST1,i4)*dx_max**2

    enddo
    do i4=1,SLConN(MST2)

      dx1=List(MS2,:)
      dx2=List(MS2-MST2+SLCon(MST2,i4),:)
      dx=DistPB(dx1,dx2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST2,i4))!differnce to Eq Leng for MS2
      EC=EC + SLPot_temp(MST2,i4)*dx_max**2
      dx=DistPB(x2,dx2)
      dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST2,i4))!differnce to Eq Leng for MS2
      EM=EM + SLPot_temp(MST2,i4)*dx_max**2
    enddo
    call random_number(rr1)
    if (rr1<exp(real(EC-EM)/1000)) then
      ETC=ETC+EM-EC
      Lattice(List(MS1,1),List(MS1,2),List(MS1,3))=0
      if (MS2/=0) then
        Lattice(List(MS2,1),List(MS2,2),List(MS2,3))=0
      endif
      List(MS1,:)=x1
      Lattice(x1(1),x1(2),x1(3))=MS1
      if (MS2/=0) then
        List(MS2,:)=x2
        Lattice(x2(1),x2(2),x2(3))=MS2
      endif
    endif

  else
    Lattice(List(MS1,1),List(MS1,2),List(MS1,3))=0
    if (MS2/=0) then
      Lattice(List(MS2,1),List(MS2,2),List(MS2,3))=0
    endif
    List(MS1,:)=x1
    Lattice(x1(1),x1(2),x1(3))=MS1
    if (MS2/=0) then
      List(MS2,:)=x2
      Lattice(x2(1),x2(2),x2(3))=MS2
    endif
  endif
else
    Rej(2,2)=Rej(2,2)+1
endif
endsubroutine MoveTransII
!!!!
!!!!
!!!!
!!!!
subroutine MoveRot

EC=0
EM=0
Propose(3)=Propose(3)+1
MP=randi(sum(Nk))
MPT=ProteinToType(MP)
MS=(MP-1)*mVk+randi(Vk(MPT)) !bead number to be moved
MST1=BeadToType(MS)
x1=List(MS,:)
ForwardRMN=1
if (BoundTo(MS)>0) then
  i=MST1
  j=BeadToType(BoundTo(MS))
  EC=EC+Energy(i,j)
  if (BeadToProtein(MS)==BeadToProtein(BoundTo(MS))) then
    EC=EC+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
  endif
endif
do i=1,26
  R=PB3(x1+NL(i,:))
  if (Lattice(R(1),R(2),R(3))/=0) then
    i2=Lattice(R(1),R(2),R(3))
    if ((Energy(MST1,BeadToType( i2 ))/=0) .and. &
        (BoundTo(i2)==0 .or. BoundTo(i2)==MS)) then
      ForwardRMN=ForwardRMN+1
      RotMoves(ForwardRMN)=i2
    endif
  endif
enddo
! Accept or Reject the move
RotMove=RotMoves(randi(ForwardRMN))
if (RotMove/=0) then
  i=MST1
  j=BeadToType(RotMove)
  EM=EM+Energy(i,j)
  if (BeadToProtein(MS)==BeadToProtein(RotMove)) then
    EM=EM+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
  endif
endif

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

endsubroutine MoveRot
!!!!
!!!!
!!!!
!!!!
subroutine MoveSlitherI  ! homopolymer slither
EM=0
EC=0
SlitherList=0
Propose(4)=Propose(4)+1
MPT=UseSlitherI_List(randi(UseSlitherI_ListN))


! logic: We chose which type, now we must find a protein of that type
! We need to find where the 'i2'th protein of type MPT is.
! We increment forward in protein number (MP)
! If the MP'th protein is of type MPT, we count it with i1
MP=0
i1=0
i2=randi(Nk(MPT)) !number of type MPT that we need to find
do while (i1<i2)
  MP=MP+1
  if (ProteinToType(MP)==MPT) then
    i1=i1+1
  endif
enddo
MS=(MP-1)*mVk

R=randi3(2*LL_temp(MST1,1)+1)-LL_temp(MST1,1)-1! goes from -LL to LL
call random_number(rr1)
if (rr1<.5) then !lets go forward!
 SlitherForward=1
 MS_back=MS+1  ! removed bead
 MS_front=MS+Vk(MPT)
 Move=PB3(List(MS_front,:)+R)
 SlitherList(1:Vk(MPT)-1,:)=List(MS+2:MS+Vk(MPT),:)
 SlitherList(Vk(MPT),:)=Move
else!lets go Backwards!
 SlitherForward=-1
 MS_back=MS+Vk(MPT)
 MS_front=MS+1
 Move=PB3(List(MS_front,:)+R)
 SlitherList(2:Vk(MPT),:)=List(MS+1:MS+Vk(MPT)-1,:)
 SlitherList(1,:)=Move
endif
MST1=BeadToType(MS_back)
MST2=BeadToType(MS_front)


BackRMN=1
ForwardRMN=1
if (Lattice(Move(1),Move(2),Move(3))==0 .or. Lattice(Move(1),Move(2),Move(3))==MS_back) then
  ! calculate the acceptance probability
  if (CanRot(MST1)) then !check if it can bind
    if (BoundTo(MS_back)>0) then
      EC=EC+Energy(MST1,BeadToType(BoundTo(MS_back)))
    endif
    do i=1,26
      R=PB3(List(MS_back,:)+NL(i,:))
      if (Lattice(R(1),R(2),R(3))/=0) then
      if (Energy(MST1,BeadToType( Lattice(R(1),R(2),R(3))))/=0 &
          .and. (BoundTo(Lattice(R(1),R(2),R(3)))==0 .or. BoundTo(Lattice(R(1),R(2),R(3)))==MS_back)) then
        BackRMN=BackRMN+1
      endif
      endif
      R=PB3(Move+NL(i,:))
      if (Lattice(R(1),R(2),R(3))/=0) then
      if (Energy(MST1,BeadToType(Lattice(R(1),R(2),R(3))))/=0 &
          .and. (BoundTo(Lattice(R(1),R(2),R(3)))==0 .or. BoundTo(Lattice(R(1),R(2),R(3)))==MS_back)) then
        ForwardRMN=ForwardRMN+1
        RotMoves(ForwardRMN)=Lattice(R(1),R(2),R(3))
      endif
      endif
    enddo
  endif
  !Accept or Reject the move
  RotMove=RotMoves(randi(ForwardRMN))
  if (RotMove/=0) then
    i=MST2
    j=BeadToType(RotMove)
    EM=EM+Energy(i,j)
    if (BeadToProtein(MS_front)==BeadToProtein(RotMove)) then
      EM=EM+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
    endif
  endif

  if (UseSpringLinkerEq) then
      do i1=sum(Vk(1:MPT-1))+1,sum(Vk(1:MPT-1))+Vk(MPT)
      do i2=1,SLConN(i1)
        if (i1<SLCon(i1,i2)) then
          x1=List(MS+i1-sum(Vk(1:MPT-1)),:)
          x2=List(MS+SLCon(i1,i2)-sum(Vk(1:MPT-1)),:)
          dx=DistPB(x1,x2)
          dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(i1,i2))!maximal differnce to Eq Leng 
          EC=EC + SLPot_temp(i1,i2)*dx_max**2
        endif
      enddo
      enddo
      do i1=sum(Vk(1:MPT-1))+1,sum(Vk(1:MPT-1))+Vk(MPT)
      do i2=1,SLConN(i1)
        if (i1<SLCon(i1,i2)) then
          x1=SlitherList(i1-sum(Vk(1:MPT-1)),:)
          x2=SlitherList(SLCon(i1,i2)-sum(Vk(1:MPT-1)),:)
          dx=DistPB(x1,x2)
          dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(i1,i2))!maximal differnce to Eq Leng 
          EM=EM + SLPot_temp(i1,i2)*dx_max**2
        endif
      enddo
      enddo
  endif

  call random_number(rr1)
  if (rr1<exp(real(EC-EM)/1000)*ForwardRMN/BackRMN) then
    ETC=ETC+EM-EC
    !clear out the end that is removed
    if (BoundTo(MS_back)>0) then
      BoundTo(BoundTo(MS_back))=0
    endif
    Lattice(List(MS_back,1),List(MS_back,2),List(MS_back,3))=0
    do i=MS_back,MS_front-SlitherForward,SlitherForward
      BoundTo(i)=BoundTo(i+SlitherForward)
      if (BoundTo(i)>0) then
        BoundTo(BoundTo(i))=i
      endif
      Lattice(List(i+SlitherForward,1),List(i+SlitherForward,2),List(i+SlitherForward,3))=i
      List(i,:)=List(i+SlitherForward,:)
    enddo
    BoundTo(MS_front)=RotMove
    if (BoundTo(MS_front)/=0) then
      BoundTo(BoundTo(MS_front))=MS_front
    endif
    List(MS_front,:)=Move
    Lattice(Move(1),Move(2),Move(3))=MS_front
  else
    Rej(4,1)=Rej(4,1)+1
  endif
else
  Rej(4,2)=Rej(4,2)+1
endif
endsubroutine MoveSlitherI
!!!!
!!!!
!!!!
!!!!
subroutine MoveSlitherII
EM=0
EC=0
logForwardRMN=0
logBackRMN=0
SlitherRot=0
SlitherList=0
Propose(5)=Propose(5)+1
MPT=UseSlitherII_List(randi(UseSlitherII_ListN))

! logic: We chose which type, now we must find a protein of that type
! We need to find where the 'i2'th protein of type MPT is.
! We increment forward in protein number (MP)
! If the MP'th protein is of type MPT, we count it with i1
MP=0
i1=0
i2=randi(Nk(MPT)) !number of type MPT that we need to find
do while (i1<i2)
  MP=MP+1
  if (ProteinToType(MP)==MPT) then
    i1=i1+1
  endif
enddo
MS=(MP-1)*mVk
MST1=sum(Vk(1:MPT-1))

flag=.false.  ! says if the move was placed in a BAD spot
if (Vk(MPT)>5) then
  Slithertimes=randi(5)
else
  SlitherTimes=randi(Vk(MPT)-1)
endif
call random_number(rr1)
if (rr1<.5) then !lets go forward!
  SlitherForward=1
  MS_back=MS  ! Bead before the start of polymer
  !find the potential move first
  do i1=1,Vk(MPT)-SlitherTimes
    SlitherList(i1,:)=List(MS_back+SlitherTimes+i1,:)
  enddo
  do i1=Vk(MPT)-SlitherTimes+1,Vk(MPT)
    R=randi3(2*LL_temp(MST1+i1,1)+1)-LL_temp(MST1+i1,1)-1 ! goes from -LL to LL
    x1=PB3(SlitherList(i1-1,:)+R)
    if (Lattice(x1(1),x1(2),x1(3))==0 .or. &
        &(Lattice(x1(1),x1(2),x1(3))>MS .and. Lattice(x1(1),x1(2),x1(3))<=MS+SlitherTimes)) then
      do i2=1,i1-1
        if (x1(1)==SlitherList(i1,1).and.x1(2)==SlitherList(i1,2).and.x1(3)==SlitherList(i1,3)) then
          flag=.true.
          exit
        endif
      enddo
      SlitherList(i1,:)=x1
    else
      flag=.true.
      exit
    endif
  enddo
else  !lets go backwards!
  SlitherForward=-1
  MS_back=MS  ! Bead before the start of polymer
  !find the potential move first
  do i1=Vk(MPT),SlitherTimes+1,-1
    SlitherList(i1,:)=List(MS_back-SlitherTimes+i1,:)
  enddo
  do i1=SlitherTimes,1,-1
    R=randi3(2*LL_temp(MST1+i1,1)+1)-LL_temp(MST1+i1,1)-1 ! goes from -LL to LL
    x1=PB3(SlitherList(i1+1,:)+R)
    if (Lattice(x1(1),x1(2),x1(3))==0 .or. &
        &(Lattice(x1(1),x1(2),x1(3))>MS+Vk(MPT)-SlitherTimes .and. Lattice(x1(1),x1(2),x1(3))<=MS+Vk(MPT))) then
      do i2=i1+1,Vk(MPT)
        if (x1(1)==SlitherList(i1,1).and.x1(2)==SlitherList(i1,2).and.x1(3)==SlitherList(i1,3)) then
          flag=.true.
          exit
        endif
      enddo
      SlitherList(i1,:)=x1
    else
      flag=.true.
      exit
    endif
  enddo
endif
if (.not.flag) then  !checks that the proposed move doesn't sterically clash with itself.
do i1=1,Vk(MPT)-1
do i2=i1+1,Vk(MPT)
if (all(SlitherList(i1,:)==SlitherList(i2,:))) then
  flag=.true.
  exit
endif
enddo
enddo
endif

if (flag) then
  rej(5,2)=rej(5,2)+1
else
  ! lets find all the interaction terms:
  do i1=1,Vk(MPT)!lets find all the interaction energies, this index is the segment number in the moving protein
  if (CanRot(MST1+i1)) then !checks that the selected bead can interact with something
    i3=1 !index for 'could have been bound'
    i4=1 !index for 'could now bind to'
    if (BoundTo(MS+i1)/=0 .and. (BoundTo(MS+i1)<=MS .or. BoundTo(MS+i1)>MS+i1)) then !if it is currently bound to a module and that module is (not on the same protein and a lower bead number) to prevent double counting
      !EC=EC+Energy(MST1+i1,BeadToType(BoundTo(MS+i1)))
      i=MST1+i1
      j=BeadToType(BoundTo(MS+i1))
      EC=EC+Energy(i,j)
      if (BeadToProtein(MS+i1)==BeadToProtein(BoundTo(MS+i1))) then
        EC=EC+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
      endif
    endif
!!! Mabye it's now reversed properly for detailed balance
    if (BoundTo(MS+i1)<=MS .or. BoundTo(MS+i1)>MS+i1) then ! if it is not 'preclaimed'
    do i2=1,26
      R=PB3(List(MS+i1,:)+NL(i2,:))  ! check the current number of rotations
      if (Lattice(R(1),R(2),R(3))/=0) then ! is there a bead
      if (Energy(MST1+i1,BeadToType(Lattice(R(1),R(2),R(3))))/=0) then  ! is that bead sticky for our bead
      if (Lattice(R(1),R(2),R(3))<=MS .or. Lattice(R(1),R(2),R(3))>MS+Vk(MPT)) then !checks if it's on another chain
        if (BoundTo(Lattice(R(1),R(2),R(3)))==0 .or. & !and is free
          (BoundTo(Lattice(R(1),R(2),R(3)))>=MS+i1.and.BoundTo(Lattice(R(1),R(2),R(3)))<=MS+Vk(MPT))) then !or bound to MP, and >current
            i3=i3+1
        endif
      elseif (Lattice(R(1),R(2),R(3))>MS+i1 .and. &  !checks if it's on same chain and higher index
          (BoundTo(Lattice(R(1),R(2),R(3)))==0 .or. & !and is free
          (BoundTo(Lattice(R(1),R(2),R(3)))<=MS.or.BoundTo(Lattice(R(1),R(2),R(3)))>=MS+i1))) then !or would have been unbound
            i3=i3+1
      endif
      endif
      endif
    enddo
    endif
    if (SlitherRot(i1)==0) then ! if it is not 'preclaimed'
      x1=SlitherList(i1,:)
      do i2=1,26  ! all binding sites from beads that will remain where they are
        R=PB3(x1+NL(i2,:)) ! Check all possible future rotations
        if (Lattice(R(1),R(2),R(3))/=0 .and. (Lattice(R(1),R(2),R(3))<=MS.or.Lattice(R(1),R(2),R(3))>MS+Vk(MPT))) then
        if (Energy(MST1+i1,BeadToType(Lattice(R(1),R(2),R(3))))/=0 .and. &  !is sticky
            (BoundTo(Lattice(R(1),R(2),R(3)))==0 .or. & !and is free
            (BoundTo(Lattice(R(1),R(2),R(3)))>MS.and.BoundTo(Lattice(R(1),R(2),R(3)))<=MS+Vk(MPT)))) then !or bound to MP
          if (all(SlitherRot(1:i1-1)/=Lattice(R(1),R(2),R(3)))) then !not bound to a new slither position
            i4=i4+1
            RotMoves(i4)=Lattice(R(1),R(2),R(3))
          endif
        endif
        endif
      enddo
      do i2=i1+1,Vk(MPT) ! Checks to see if it will bind to it's own chain, but only chains of higher index
        if (Energy(MST1+i1,MST1+i2)/=0) then ! is sticky for self
        x2=SlitherList(i2,:)
        if (SlitherRot(i2)==0 .and. all(abs(DistPB(x1,x2))<=1)) then !is free and a neighbor
          i4=i4+1
          RotMoves(i4)=MS+i2
        endif
        endif
      enddo
      SlitherRot(i1)=RotMoves(randi(i4))
      if (SlitherRot(i1)/=0) then !is bound to a bead
        if (SlitherRot(i1)>MS+i1 .and. SlitherRot(i1)<=MS+Vk(MPT)) then !is bound to self then we need to assign the other bead
          SlitherRot(SlitherRot(i1)-MS)=MS+i1
        elseif (SlitherRot(i1)>MS .and. SlitherRot(i1)<MS+i1) then !this should never happen, eventually remove it
          write(*,*) 'Slither is binding to something it should not be able to!'
          stop
        elseif (SlitherRot(i1)==MS+i1) then
          write(*,*) 'Slither is binding to its self!'
          write(*,*) t
          write(*,*) i1
          write(*,*) SlitherRot(1:i1)
          write(*,*)
          write(*,*) RotMoves
          write(*,*) i4
          stop
        endif
        i=MST1+i1
        j=BeadToType(SlitherRot(i1))
        EM=EM+Energy(i,j)
        if (BeadToProtein(MS+i1)==BeadToProtein(SlitherRot(i1))) then
          EM=EM+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
        endif
      endif
    endif
    !lets do a few detailed balance calculations
    logBackRMN=logBackRMN+log(real(i3))
    logForwardRMN=logForwardRMN+log(real(i4))
  endif
  enddo
  !!!If use LLEq or not
  if (UseSpringLinkerEq) then
      do i1=1,Vk(MPT)
        x1=List(MS+i1,:)
        do i2=1,SLConN(MST1+i1)
          if (SLCon(MST1+i1,i2)>MST1+i1) then
            x2=List(MS-MST1+SLCon((MST1+i1),i2),:)
            dx=DistPB(x1,x2)
            dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1+i1,i2))
            EC=EC + SLPot_temp(MST1+i1,i2)*dx_max**2
          endif
        enddo
      enddo
      !do i1=Vk(MPT)-SlitherTimes+1,Vk(MPT)
      do i1=1,Vk(MPT)
        x1=SlitherList(i1,:)
!        print *,"SlitherList=",SlitherList(i1,:)
        do i2=1,SLConN(MST1+i1)
          if (SLCon(MST1+i1,i2)>MST1+i1) then
            x2=SlitherList(SLCon((MST1+i1),i2)-MST1,:)
            dx=DistPB(x1,x2)
            dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1+i1,i2))
            EM=EM + SLPot_temp(MST1+i1,i2)*dx_max**2
          endif
        enddo
      enddo
  endif

  !!! Accept or reject the move
  call random_number(rr1)
  if (rr1<exp(real(EC-EM)/1000 + logForwardRMN - logBackRMN)) then
    ETC=ETC+EM-EC
    do i1=1,Vk(MPT) !clear old bonds
      if (BoundTo(MS+i1)/=0) then
        BoundTo(BoundTo(MS+i1))=0
      endif
      Lattice(List(MS+i1,1),List(MS+i1,2),List(MS+i1,3))=0 !clear lattice
    enddo
    do i1=1,Vk(MPT) !build new bonds
      BoundTo(MS+i1)=SlitherRot(i1) ! set what I'm bound to
      if (SlitherRot(i1)/=0) then
        BoundTo(SlitherRot(i1))=MS+i1 !set what I'm bound to to be bound to me
      endif
      Lattice(SlitherList(i1,1),SlitherList(i1,2),SlitherList(i1,3))=MS+i1 !build new lattice
      List(MS+i1,:)=SlitherList(i1,:)
    enddo
  else
    rej(5,1)=rej(5,1)+1
  endif
endif
endsubroutine MoveSlitherII
!!!!
!!!!
!!!!
!!!!
subroutine MoveClusterI


call NetworkAnalysis
!i'd like this to do min(AllHostNN-1,10) times where it doesn't do the same cluster twice... that's for a future day to program.

! AllHostNN = number of clusters
! AllHostN = start and end of each cluster
if (AllHostNN>1) then
  i5=min(AllHostNN-1,10) !number of moves
  Propose(6)=Propose(6)+i5

  FisherN=AllHostNN-1 !number of available (set up fisher-yates algorithm)
  j2=minval(maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1))) ! the largest to remove
  do i1=1,j2-1  !builds the fisher list of available clusters
    Fisher(i1)=i1
  enddo
  do i1=j2+1,AllHostNN !skips the largest cluster
    Fisher(i1-1)=i1
  enddo
  do i1=1,i5  !!!!! Does up to 10 network moves each time MoveCluster is called
    MPT=randi(FisherN) !random element from fisher list
    MP=Fisher(MPT)    !cluster number pulled from the fisher list
    Fisher(MPT)=Fisher(FisherN) !take the last element from the fisher list and move it to where it will be available in the future
    FisherN=FisherN-1 !shrink the population of fish
    flag=.false. !says if the move sterically clashed
    
    Move=[randiBS(BoxSize)]
    if (Move(1)==0 .and. Move(2)==0 .and. Move(3)==0) then !Does it move the protein a total distance of zero?  Lets just ignore this move
      Rej(6,2)=Rej(6,2)+1
      cycle
    endif
    do i2=AllHostN(MP,1),AllHostN(MP,2)  !run thru all proteins in host
      MPT1=ProteinToType(AllHost(i2))
      do i3=(AllHost(i2)-1)*mVk+1,(AllHost(i2)-1)*mVk+Vk(MPT1)
        x1=PB3(List(i3,:)+Move) !new loc
        if (Lattice(x1(1),x1(2),x1(3))/=0) then !if new loc is occupied
          flag=.true.
          exit
        endif
      enddo
      if (flag) then
        exit
      endif
    enddo  !end of examining if the move is acceptable
    if (.not.flag) then
      do i2=AllHostN(MP,1),AllHostN(MP,2)  !run thru all proteins in host
        MPT1=ProteinToType(AllHost(i2))
        do i3=(AllHost(i2)-1)*mVk+1,(AllHost(i2)-1)*mVk+Vk(MPT1)  !move each segment in protein
          Lattice(List(i3,1),List(i3,2),List(i3,3))=0
        enddo
      enddo
      do i2=AllHostN(MP,1),AllHostN(MP,2)  !run thru all proteins in host
        MPT1=ProteinToType(AllHost(i2))
        do i3=(AllHost(i2)-1)*mVk+1,(AllHost(i2)-1)*mVk+Vk(MPT1)  !move each segment in protein
          x1=PB3(List(i3,:)+Move)
          List(i3,:)=x1
          Lattice(x1(1),x1(2),X1(3))=i3
        enddo
      enddo
    else
      Rej(6,1)=Rej(6,1)+1
    endif
  enddo
else
  Rej(6,2)=Rej(6,2)+10
endif
endsubroutine MoveClusterI
!!!!
!!!!
!!!!
!!!!
subroutine MoveClusterII
Propose(7)=Propose(7)+1
MP=randi(sum(Nk))
MPT=ProteinToType(MP)
AllHostNN=1  !number of proteins in this cluster
AllHost(1)=MP !the proteins in this cluster
IsHost(1:sum(Nk))=.true.  !has the protein been seen by this cluster yet?
IsHost(MP)=.false.  ! it sees itself
i2=1 !current scan number
flag=.true. !flag for if move has failed from the cluster being too large for this type of move
do while (i2<=AllHostNN .and. flag)
  MPT1=ProteinToType(AllHost(i2))
  do i3=(AllHost(i2)-1)*mVk+1,(AllHost(i2)-1)*mVk+Vk(MPT1)  !move each segment in protein
    if (BoundTo(i3)/=0) then
      i4=BeadToProtein(BoundTo(i3))
      if (IsHost(i4).eqv..true.) then
        AllHostNN=AllHostNN+1
        AllHost(AllHostNN)=i4
        IsHost(i4)=.false.
      endif
    endif
  enddo
  i2=i2+1 !start analyzing the next protein
  if (AllHostNN>5) then
    flag=.false.  !cluster is too big for this move
  endif
enddo
if (flag) then !flag, will now be repurposed to test if there is a steric clash
  Move=randiBS(BoxSize)-1
  if (Move(1)==0 .and. Move(2)==0 .and. Move(3)==0) then !Does it move the protein a total distance of zero?  Lets just ignore this move
    Rej(7,2)=Rej(7,2)+1
  endif
  do i2=1,AllHostNN  !run thru all proteins in host
    MPT1=ProteinToType(AllHost(i2))
    do i3=(AllHost(i2)-1)*mVk+1,(AllHost(i2)-1)*mVk+Vk(MPT1)  !move each segment in protein
      x1=PB3(List(i3,:)+Move) !new loc
      if (Lattice(x1(1),x1(2),x1(3))/=0) then !if new loc is occupied
        flag=.false.
        exit
      endif
    enddo
    if (flag.eqv..false.) then
      exit
    endif
  enddo  !end of examining if the move is acceptable
  if (flag) then
    do i2=1,AllHostNN  !run thru all proteins in host and delete them from the lattice
      MPT1=ProteinToType(AllHost(i2))
      do i3=(AllHost(i2)-1)*mVk+1,(AllHost(i2)-1)*mVk+Vk(MPT1)  !move each segment in protein
        Lattice(List(i3,1),List(i3,2),List(i3,3))=0
      enddo
    enddo
    do i2=1,AllHostNN  !run thru all proteins in host and rebuild them on the lattice and list
      MPT1=ProteinToType(AllHost(i2))
      do i3=(AllHost(i2)-1)*mVk+1,(AllHost(i2)-1)*mVk+Vk(MPT1)  !move each segment in protein
        x1=PB3(List(i3,:)+Move)
        List(i3,:)=x1
        Lattice(x1(1),x1(2),x1(3))=i3
      enddo
    enddo
  else
    Rej(7,1)=Rej(7,1)+1
  endif
else
  Rej(7,2)=Rej(7,2)+1
endif
endsubroutine MoveClusterII
!!!!
!!!!
!!!!
!!!!
subroutine MoveGrand

Propose(8)=Propose(8)+1
if (NkN>1) then
    write(*,*) 'Ive only programmed it for one type of protein.  Crashing.'
    stop
endif
if (Vk(1)>1) then
    write(*,*) 'Ive only programmed polymers of length 1.  Crashing.'
    stop
endif

call random_number(rr3)
if (rr3(1)<.5) then !Add a protein
  if (sum(Nk)<mNk) then ! must fit in the list form
    x1=randiBS(BoxSize)
    if (Lattice(x1(1),x1(2),x1(3))==0) then
      MPT=randi(NkN)
      rr1=BoxSize(1)*BoxSize(2)*BoxSize(3)/(sum(Nk)+1)*exp(real(ChemPot(MPT))/1000) ! r
      if (rr3(2)<rr1) then
        Nk(1)=Nk(1)+1
        Lattice(x1(1),x1(2),x1(3))=sum(Nk)
        List(sum(Nk),:)=x1
        BoundTo(sum(Nk))=0
        ProteinToType(sum(Nk))=MPT
      else
        Rej(8,1)=Rej(8,1)+1
      endif
    else 
      Rej(8,2)=Rej(8,2)+1
    endif
  else
    write(*,*) 'warning, sumNk=mNk=',mNk,' and Im trying to add a protein.  t=',t
  endif
else !remove a protein, then shift what ever protein is at the end to fill it's spot
  if (Nk(1)>0) then !can't go below Nk=0
    MP=randi(sum(Nk))
    MPT=ProteinToType(MP)
    flag=.true.
    
    if (all(BoundTo((MP-1)*mVk+1:(MP-1)*mVk+Vk(MPT))==0)) then
      rr1=real(sum(Nk))/BoxSize(1)*BoxSize(2)*BoxSize(3)/exp(real(ChemPot(MPT))/1000) ! 1/r
      if (rr3(2)<rr1) then
        do i1=(MP-1)*mVk+1,(MP-1)*mVk+Vk(MPT)  ! remove the protein
          Lattice(List(i1,1),List(i1,2),List(i1,3))=0
        enddo
        if (MP<sum(Nk)) then !shift the protein at the end to fill what was deleted
          MPT=ProteinToType(sum(Nk))
          ProteinToType(MP)=MPT
          do i1=1,Vk(MPT)  ! rename what's on the lattice
            Lattice(List((sum(Nk)-1)*mVk+i1,1),List((sum(Nk)-1)*mVk+i1,2),List((sum(Nk)-1)*mVk+i1,3))=(MP-1)*mVk+i1
          enddo
          List((MP-1)*mVk+1:(MP-1)*mVk+Vk(MPT),:)=List((sum(Nk)-1)*mVk+1:(sum(Nk)-1)*mVk+Vk(MPT),:)
          BoundTo((MP-1)*mVk+1:(MP-1)*mVk+Vk(MPT))=BoundTo((sum(Nk)-1)*mVk+1:(sum(Nk)-1)*mVk+Vk(MPT))
          do i1=1,Vk(MPT)
          if (BoundTo((MP-1)*mVk+i1)>0) then
            BoundTo(BoundTo((MP-1)*mVk+i1))=(MP-1)*mVk+i1
          endif
          enddo
        endif
        Nk(MPT)=Nk(MPT)-1
      else
        Rej(8,1)=Rej(8,1)+1
      endif
    else
      Rej(8,2)=Rej(8,2)+1
    endif
  else
    Rej(8,2)=Rej(8,2)+1
  endif
endif


endsubroutine MoveGrand
!!!!
!!!!
!!!!
!!!!
subroutine NetworkAnalysis
! constructs how all the proteins are networked
! AllHost gives a list of the proteins in the order they are bound to complexes
! AllHostN gives the list of the start and end position of each cluster
! AllHostNN gives the total list length of AllHostN

IsHost(1:sum(Nk))=.true.  !the name is legacy.  It's actually a logical for if this protein has been seen in the script yet
AllHostNN=0  !number of clusters
flag=.true.  !marks if it is the first cluster on the list because that has to be treated differently
! new code
do i1=1,sum(Nk)
if (IsHost(i1)) then ! if it is a host, collect up what it's bound to, otherwise it's in someone else's cluster
  AllHostNN=AllHostNN+1 !number of clusters, we just found another cluster
  FreeHostN=1  !number of polymers in this cluster
  FreeHost(1)=i1  !list of polymers in this cluster, just added self to the list
  IsHost(i1)=.false.  !prevents it from being picked again in its own cluster
  i2=1  !current number of proteins bound and counted
  do while (i2<=FreeHostN) !indexing list of proteins connected and counted
    MPT=ProteinToType(FreeHost(i2))
    do i3=(FreeHost(i2)-1)*mVk+1,(FreeHost(i2)-1)*mVk+Vk(MPT) !list of beads in the protein
      if (BoundTo(i3)/=0) then
        MP=BeadToProtein(BoundTo(i3))
        if (IsHost(MP)) then !checks if the protein needs to be connected to this cluster or is already counted
          FreeHostN=FreeHostN+1  !add it to the cluster
          FreeHost(FreeHostN)=MP
          IsHost(MP)=.false.   !prevent it from being counted again
        endif
      endif
    enddo
    i2=i2+1
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

if (WriteCoreCluster) then
  IsHalfFull(1:mNk) = .false.
  IsFull(1:mNk) = .false.
  do i1 = 1, sum(Nk)
    bonds = 0
    do i2 = (i1 - 1) * mVk + 1, (i1 - 1) * mVk + Vk(ProteinToType(i1))
      if (BoundTo(i2) /= 0) then
        bonds = bonds + 1
      endif
    end do
    Occupancy(i1) = bonds / Vak(ProteinToType(i1))
    if (Occupancy(i1) >= 0.5) then
      IsHalfFull(i1) = .true.
      if (1-Occupancy(i1) < 1e-5) then
        IsFull(i1) = .true.
      endif
    endif
  end do

  Core = 0 !Core means fully occupied polymers
  CoreHalf = 0!CoreHalf means more than half occupied polymers
  CoreN = 0
  CoreHalfN = 0
  i1 = maxloc(AllHostN(1:AllHostNN, 2) - AllHostN(1:AllHostNN, 1), 1)
  do i2 = AllHostN(i1, 1), AllHostN(i1, 2) ! Runs through all proteins
    if (IsFull(AllHost(i2))) then
      CoreN = CoreN + 1
      Core(CoreN) = AllHost(i2)
    endif
    if (IsHalfFull(AllHost(i2))) then
      CoreHalfN = CoreHalfN + 1
      CoreHalf(CoreHalfN) = AllHost(i2)
    endif
  enddo

endif

endsubroutine NetworkAnalysis
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
  !Calculate Center of Mass of the Cluster
  y1=0
  y2=0
  i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
  do i3=AllHostN(i1,1),AllHostN(i1,2)
    MPT=ProteinToType(AllHost(i3))
    do i4=(AllHost(i3)-1)*mVk+1,(AllHost(i3)-1)*mVk+Vk(MPT)
      y1=y1+sin(pi*2/BoxSize*List(i4,:))
      y2=y2+cos(pi*2/BoxSize*List(i4,:))
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
        MPT=ProteinToType(AllHost(i3))
        if (MPT==i2) then !if right protein type, run through all beads
        do i4=(AllHost(i3)-1)*mVk+1,(AllHost(i3)-1)*mVk+Vk(MPT)
          y1=y1+sin(pi*2/BoxSize*List(i4,:)) ! uses fourier series to find COM of protein type i2
          y2=y2+cos(pi*2/BoxSize*List(i4,:))
        enddo
        endif
      enddo
      COMi(i2,:)=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize)) 
    enddo
  endif
endif

if ( (E_XYZDist_All/=0 .and. mod(t,E_XYZDist_All)==0) .or. &
     (E_RG_All/=0 .and. mod(t,E_RG_All)==0)) then
  ! calc Center of Mass of all polymers
  y1=0
  y2=0
  do i1=1,sum(Nk)
  do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
    y1=y1+sin(pi*2/BoxSize*List(i2,:))
    y2=y2+cos(pi*2/BoxSize*List(i2,:))
  enddo
  enddo
  COM_All=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the system
  !write(*,*) 'Center of Mass is', COM_All
endif

! cluster size of the three largest clusters and the polymers that are in them
if (E_LargestCluster/=0 .and. mod(t,E_LargestCluster)==0)  then
  Analysis_log=.true.
  do i1=1,min(3,AllHostNN) !three largest clusters
    i2=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1,Analysis_log(1:AllHostNN))
    LargestCluster(1+(NkN+1)*(i1-1))=AllHostN(i2,2)-AllHostN(i2,1)+1
    do i3=AllHostN(i2,1),AllHostN(i2,2)
      LargestCluster(1+(NkN+1)*(i1-1)+ProteinToType(AllHost(i3))) = &
          LargestCluster(1+(NkN+1)*(i1-1)+ProteinToType(AllHost(i3)))+1
    enddo
    Analysis_log(i2)=.false.
  enddo
  if(WriteCoreCluster) then
    LargestCluster(1+(NkN+1)*3)=CoreN !Core Cluster detail in second last columns
    do i3=1,CoreN
      LargestCluster(1+(NkN+1)*3+ProteinToType(Core(i3))) = &
          LargestCluster(1+(NkN+1)*3+ProteinToType(Core(i3)))+1
    enddo
        LargestCluster(1+(NkN+1)*4)=CoreHalfN
    do i3=1,CoreHalfN !Core Cluster detail in last columns
      LargestCluster(1+(NkN+1)*4+ProteinToType(CoreHalf(i3))) = &
          LargestCluster(1+(NkN+1)*4+ProteinToType(CoreHalf(i3)))+1
    enddo
  endif
  write(102,*) LargestCluster
  LargestCluster=0
endif
! histogram of the cluster sizes (up to 20, all else gets put into the 21 bin)
if (E_ClusterHist/=0 .and. mod(t,E_ClusterHist)==0)  then
  do i1=1,AllHostNN
    if (AllHostN(i1,2)-AllHostN(i1,1)+1>20) then
      ClusterHist(21)=ClusterHist(21)+1
    else
      ClusterHist(AllHostN(i1,2)-AllHostN(i1,1)+1)= &
          ClusterHist(AllHostN(i1,2)-AllHostN(i1,1)+1)+1
    endif
  enddo
  write(104,*) ClusterHist(:)
  ClusterHist=0
endif


!XYZDistribution for all proteins, recentered on COM
if (E_XYZDist_All/=0 .and. mod(t,E_XYZDist_All)==0)  then
    do i1=1,sum(Nk)
      MPT=ProteinToType(i1)
      do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(MPT)
        !now calculate the histogram
        !x1=PB3(List(i2,:)-floor(COM_Cluster)+BoxSize/2)
        x1=PB3(List(i2,:)-floor(COM_All)+BoxSize/2) !new coordinates re-centered on COM_All
        XDist_All(x1(1),MPT)= XDist_All(x1(1),MPT) + 1
        YDist_All(x1(2),MPT)= YDist_All(x1(2),MPT) + 1
        ZDist_All(x1(3),MPT)= ZDist_All(x1(3),MPT) + 1
      enddo
    enddo
  do i1=1,NkN
      write(Str1,*) i1
      write(200+i1+10,*) XDist_All(:,i1)
      write(200+i1+20,*) YDist_All(:,i1)
      write(200+i1+30,*) ZDist_All(:,i1)
  enddo
  XDist_All=0
  YDist_All=0
  ZDist_All=0
endif

!! !XYZDistribution (spacial organization) for the largest cluster, recentered on COM: 
if (E_XYZDist_Cluster/=0 .and. mod(t,E_XYZDist_Cluster)==0)  then
  i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
  if (AllHostN(i1,2)-AllHostN(i1,1)+1>=50) then ! no histogram if the cluster isn't at least 50 big
     !now calculate the histogram
      do i3=AllHostN(i1,1),AllHostN(i1,2)
        MPT=ProteinToType(AllHost(i3))
        do i4=(AllHost(i3)-1)*mVk+1,(AllHost(i3)-1)*mVk+Vk(MPT)
        x1=PB3(List(i4,:)-floor(COM_Cluster)+BoxSize/2) !new coordinates re-centered on COM
        XDist_Cluster(x1(1),MPT)= XDist_Cluster(x1(1),MPT) + 1
        YDist_Cluster(x1(2),MPT)= YDist_Cluster(x1(2),MPT) + 1
        ZDist_Cluster(x1(3),MPT)= ZDist_Cluster(x1(3),MPT) + 1
        if(WriteIndiDist_Cluster) then
          do i2=1,NkN
            x1=PB3(List(i4,:)-floor(COMi(i2,:))+BoxSize/2) !MPT centered on COM_i2
            XDist_Clusterij(x1(1),MPT,i2)= XDist_Clusterij(x1(1),MPT,i2) + 1
            YDist_Clusterij(x1(2),MPT,i2)= YDist_Clusterij(x1(2),MPT,i2) + 1
            ZDist_Clusterij(x1(3),MPT,i2)= ZDist_Clusterij(x1(3),MPT,i2) + 1
          enddo
        endif
        enddo
      enddo
  endif
  do i1=1,NkN
    write(Str1,*) i1
    write(250+i1+10,*) XDist_Cluster(:,i1)
    write(250+i1+20,*) YDist_Cluster(:,i1)
    write(250+i1+30,*) ZDist_Cluster(:,i1)
  enddo
  XDist_Cluster=0
  YDist_Cluster=0
  ZDist_Cluster=0

  if(WriteIndiDist_Cluster) then
  do i1=1,NkN
    write(Str1,*) i1
    do i2=1,NkN
      write(Str2,*) i2
      write(100*i1+10000*i2+1,*) XDist_Clusterij(:,i1,i2)
      write(100*i1+10000*i2+2,*) YDist_Clusterij(:,i1,i2)
      write(100*i1+10000*i2+3,*) ZDist_Clusterij(:,i1,i2)
    enddo
  enddo
  XDist_Clusterij=0
  YDist_Clusterij=0
  ZDist_Clusterij=0
  endif
endif

!! Radial Histogram: (spacial organization) for the largest cluster
if (E_RadDist_Cluster/=0 .and. mod(t,E_RadDist_Cluster)==0)  then
  i1=maxloc(AllHostN(1:AllHostNN,2)-AllHostN(1:AllHostNN,1),1)
  if (AllHostN(i1,2)-AllHostN(i1,1)+1>=50) then ! no histogram if the cluster isn't at least 50 big
    if(WriteIndiDist_Cluster) then
    do i2=1,NkN !COM of protein i2
      !now calculate the histogram
      do i3=AllHostN(i1,1),AllHostN(i1,2)
        MPT=ProteinToType(AllHost(i3))
        do i4=(AllHost(i3)-1)*mVk+1,(AllHost(i3)-1)*mVk+Vk(MPT)
          y2=real(List(i4,:))
          i5=floor(sqrt(sum(rDistPB(y2,COMi(i2,:))**2))+1)
          RadDistCluster(i5,ProteinToType(AllHost(i3)),i2)= &
            RadDistCluster(i5,ProteinToType(AllHost(i3)),i2) + 1
        enddo
      enddo
    enddo
    endif
    ! COM_Cluster is based on all proteins in the largest cluster
    i2=NkN+1!I did not change the way of assign memory, which is not most effcient in terms of memory usage, but current format is easier in programming.
    do i3=AllHostN(i1,1),AllHostN(i1,2)
      MPT=ProteinToType(AllHost(i3))
      do i4=(AllHost(i3)-1)*mVk+1,(AllHost(i3)-1)*mVk+Vk(MPT)
        y2=real(List(i4,:))
        i5=floor(sqrt(sum(rDistPB(y2,COM_Cluster)**2))+1)
        RadDistCluster(i5,ProteinToType(AllHost(i3)),i2)= &
          RadDistCluster(i5,ProteinToType(AllHost(i3)),i2) + 1
      enddo
    enddo
  endif
  do i1=1,NkN
    write(Str1,*) i1
    if(WriteIndiDist_Cluster) then
    do i2=1,NkN
      write(Str2,*) i2
      write(100*i1+10000*i2,*) RadDistCluster(:,i1,i2)
    enddo
    endif
    write(200+i1,*) RadDistCluster(:,i1,NkN+1)
  enddo
  RadDistCluster=0

  ! Radial Histgram for Core_Cluster, very expensive!
  If(WriteCoreCluster) then
    if(WriteIndiDist_Cluster) then
    do i2=1,NkN !Which specific component is used for the com?
      y1=0
      y2=0
      do i3=1,CoreN !runs through all proteins
        MPT=ProteinToType(Core(i3))
        if (MPT==i2) then !if right protein type, run through all beads
          do i4=(Core(i3)-1)*mVk+1, (Core(i3)-1)*mVk+Vk(MPT)
            y1=y1+sin(pi*2/BoxSize*List(i4,:)) ! uses fourier series to find COM
            y2=y2+cos(pi*2/BoxSize*List(i4,:))
          enddo
        endif
      enddo
      y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the protein type i2 in Core_cluster
      do i3=1,CoreN
        MPT=ProteinToType(Core(i3))
        do i4=(Core(i3)-1)*mVk+1,(Core(i3)-1)*mVk+Vk(MPT)
          y2=real(List(i4,:))
          i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
          RadDistCore(i5,MPT,i2)= &
            RadDistCore(i5,MPT,i2) + 1
        enddo
      enddo
    enddo
    endif

    i2=NkN+1
    y1=0
    y2=0
    do i3=1,CoreN
      MPT=ProteinToType(Core(i3))
      do i4=(Core(i3)-1)*mVk+1,(Core(i3)-1)*mVk+Vk(MPT)
        y1=y1+sin(pi*2/BoxSize*List(i4,:)) ! uses fourier series to find COM
        y2=y2+cos(pi*2/BoxSize*List(i4,:))
      enddo
    enddo
    y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the Core_cluster
    do i3=1,CoreN
      MPT=ProteinToType(Core(i3))
      do i4=(Core(i3)-1)*mVk+1,(Core(i3)-1)*mVk+Vk(MPT)
        y2=real(List(i4,:))
        i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
        RadDistCore(i5,MPT,i2)= &
          RadDistCore(i5,MPT,i2) + 1
      enddo
    enddo
    do i1=1,NkN
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(Str1,*) i1
        write(Str2,*) i2
        write(100*i1+10000*i2+1,*) RadDistCore(:,i1,i2)
      enddo
      endif
      write(200+i1+30,*) RadDistCore(:,i1,NkN+1)
    enddo
    RadDistCore=0

    if(WriteIndiDist_Cluster) then
    do i2=1,NkN !Which specific component is used for the com?
      y1=0
      y2=0
      do i3=1,CoreHalfN !runs through all proteins
        MPT=ProteinToType(CoreHalf(i3))
        if (MPT==i2) then !if right protein type, run through all beads
          do i4=(CoreHalf(i3)-1)*mVk+1, (CoreHalf(i3)-1)*mVk+Vk(MPT)
            y1=y1+sin(pi*2/BoxSize*List(i4,:)) ! uses fourier series to find COM
            y2=y2+cos(pi*2/BoxSize*List(i4,:))
          enddo
        endif
      enddo
      y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the protein type i2 in CoreHalf_cluster
      do i3=1,CoreHalfN
        MPT=ProteinToType(CoreHalf(i3))
        do i4=(CoreHalf(i3)-1)*mVk+1,(CoreHalf(i3)-1)*mVk+Vk(MPT)
          y2=real(List(i4,:))
          i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
          RadDistCoreHalf(i5,MPT,i2)= &
            RadDistCoreHalf(i5,MPT,i2) + 1
        enddo
      enddo
    enddo
    endif

    i2=NkN+1
    y1=0
    y2=0
    do i3=1,CoreHalfN
      MPT=ProteinToType(CoreHalf(i3))
      do i4=(CoreHalf(i3)-1)*mVk+1,(CoreHalf(i3)-1)*mVk+Vk(MPT)
        y1=y1+sin(pi*2/BoxSize*List(i4,:)) ! uses fourier series to find COM
        y2=y2+cos(pi*2/BoxSize*List(i4,:))
      enddo
    enddo
    y3=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the CoreHalf_cluster
    do i3=1,CoreHalfN
      MPT=ProteinToType(CoreHalf(i3))
      do i4=(CoreHalf(i3)-1)*mVk+1,(CoreHalf(i3)-1)*mVk+Vk(MPT)
        y2=real(List(i4,:))
        i5=floor(sqrt(sum(rDistPB(y2,y3)**2))+1)
        RadDistCoreHalf(i5,MPT,i2)= &
          RadDistCoreHalf(i5,MPT,i2) + 1
      enddo
    enddo
    do i1=1,NkN
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        write(Str1,*) i1
        write(Str2,*) i2
        write(100*i1+10000*i2+2,*) RadDistCoreHalf(:,i1,i2)
      enddo
      endif
      write(200+i1+60,*) RadDistCoreHalf(:,i1,NkN+1)
    enddo
    RadDistCoreHalf=0
  endif
endif

if (E_RG_All/=0 .and. mod(t,E_RG_All)==0)  then
  !calc Gyration Tensor
  GT=0
  do i1=1,sum(Nk)
  do i2=(i1-1)*mVk+1,(i1-1)*mVk +Vk(ProteinToType(i1))
    y3=rDistPB(real(List(i2,:)),COM_All)
    do i3=1,3
    do i4=i3,3
      GT(i3,i4)=GT(i3,i4)+y3(i3)*y3(i4)/sum(Nk*Vk)
      GT(i4,i3)=GT(i3,i4)
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
    MPT=ProteinToType(AllHost(i3))
    NB_Cluster=NB_Cluster+Vk(MPT)
  enddo
  !calc Gyration Tensor
  GT=0
  do i3=AllHostN(i1,1),AllHostN(i1,2)
  MPT=ProteinToType(AllHost(i3))
  do i4=(AllHost(i3)-1)*mVk+1,(AllHost(i3)-1)*mVk+Vk(MPT)
    y3=rDistPB(real(List(i4,:)),COM_Cluster)
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
      MPT=ProteinToType(Core(i3))
      do i4=(Core(i3)-1)*mVk+1,(Core(i3)-1)*mVk+Vk(MPT)
        y1=y1+sin(pi*2/BoxSize*List(i4,:))
        y2=y2+cos(pi*2/BoxSize*List(i4,:))
        NB_Cluster=NB_Cluster+1
      enddo
    enddo
    y1=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the cluster
    !write(*,*) 'Center of Mass of Largest Cluster is', y1

    
    !calc Gyration Tensor
    GT=0
    do i3=1,CoreN
    MPT=ProteinToType(Core(i3))
    do i4=(Core(i3)-1)*mVk+1,(Core(i3)-1)*mVk+Vk(MPT)
      y3=rDistPB(real(List(i4,:)),y1)
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
      MPT=ProteinToType(CoreHalf(i3))
      do i4=(CoreHalf(i3)-1)*mVk+1,(CoreHalf(i3)-1)*mVk+Vk(MPT)
        y1=y1+sin(pi*2/BoxSize*List(i4,:))
        y2=y2+cos(pi*2/BoxSize*List(i4,:))
        NB_Cluster=NB_Cluster+1
      enddo
    enddo
    y1=mod(atan2(y1,y2)*BoxSize/2/pi+BoxSize,real(BoxSize))  !center of the cluster
    !write(*,*) 'Center of Mass of Largest Cluster is', y1

    
    !calc Gyration Tensor
    GT=0
    do i3=1,CoreHalfN
    MPT=ProteinToType(CoreHalf(i3))
    do i4=(CoreHalf(i3)-1)*mVk+1,(CoreHalf(i3)-1)*mVk+Vk(MPT)
      y3=rDistPB(real(List(i4,:)),y1)
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
  do i1=1,sum(Nk)
  do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
  if (BoundTo(i2) /=0) then
    BoundSites(BeadToType(i2))=BoundSites(BeadToType(i2))+1
  endif
  enddo
  enddo
  write(105,*) BoundSites
  BoundSites=0
endif

if (E_BondTypes/=0 .and. mod(t,E_BondTypes)==0) then
  BondTypes=0
  do i1=1,sum(Nk)
  do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
  if (BoundTo(i2)/=0) then
    BondTypes(BeadToType(i2),BeadToType(BoundTo(i2)))=BondTypes(BeadToType(i2),BeadToType(BoundTo(i2)))+1
  endif
  enddo
  enddo
  do i1=1,sum(Vk)
    write(106,*) BondTypes(i1,:)
  enddo
endif

if(UseSpringLinkerEq)then
if(E_SpringLinkerVector/=0 .and. mod(t,E_SpringLinkerVector)==0 ) then
  SLVec=0
  do i1=1,sum(Nk)
  do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
  MST1=BeadToType(i2)
  MPT=ProteinToType(BeadToProtein(i2))
  do i4=1,SLConN(MST1)
    if (SLCon(MST1,i4)>MST1) then
    x1=List(i2,:)
    x2=List(i2-MST1+SLCon(MST1,i4),:)
    dx=DistPB(x1,x2)
    SLVec(1:3)=SLVec(1:3)+dx! average values of x y z
    SLVec(4)=SLVec(4)+sqrt(real(sum(dx**2))) ! average length
    SLVec(5)=SLVec(5)+acos(SLVec(1)/SLVec(4)) ! average alpha (angle between x axis)
    SLVec(6)=SLVec(6)+acos(SLVec(2)/SLVec(4))
    SLVec(7)=SLVec(7)+acos(SLVec(3)/SLVec(4))
    endif
  enddo
  enddo
  enddo
  SLVec=SLVec/Nk(MPT)
  write(99,*) SLVec
endif
endif



if(E_SelfLoop/=0 .and. mod(t,E_SelfLoop)==0 ) then
  SelfLoopN=0
  do i1=sum(Nk(1:LinkerType-1))+1, sum(Nk(1:LinkerType-1))+Nk(LinkerType)
    if (BoundTo((i1-1)*mVk+1)/=0 .and. BoundTo((i1-1)*mVk+Vk(LinkerType))/=0) then !both ends bound
      if (BeadToProtein(BoundTo((i1-1)*mVk+1))==BeadToProtein(BoundTo((i1-1)*mVk+Vk(LinkerType)))) then !both ends bound to same protein
        SelfLoopN=SelfLoopN+1
      endif
    endif
  enddo
  write(97,*) SelfLoopN
endif

if(E_Network/=0 .and. mod(t,E_Network)==0 ) then
  LinkerCon=0
  do i1=sum(Nk(1:LinkerType-1))+1, sum(Nk(1:LinkerType-1))+Nk(LinkerType)
    LinkerCon(i1-sum(Nk(1:LinkerType-1)),1)=i1
    i3=1
    do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(LinkerType)
      if (BoundTo(i2)/=0) then
        i3=i3+1
        LinkerCon(i1-sum(Nk(1:LinkerType-1)),i3)=BeadToProtein(BoundTo(i2))
      endif
    enddo
  enddo
  ! write(*,*) 'Linker Type is', LinkerType
  ! write(*,*) 'There are', Nk(LinkerType), 'Linkers'
  do i1=1, Nk(LinkerType)
  write(98,*) LinkerCon(i1,:)
  enddo
endif

if (E_Lattice/=0 .and. mod(t,E_Lattice)==0)  then
  write(1,*) sum(Vk*Nk)
  write(1,*)
  i3=0
  str2='   '
  do i1=1,sum(Nk)
  do i2=1,Vk(ProteinToType(i1))
    i3=(i1-1)*mVk+i2
    write(str1,'(A,I1,A)') '(A,I',floor(log10(real(i2))+1),',A, I5, I5, I5)'
    write(1, str1) Char(64 + ProteinToType(i1)),i2,str2(1:4-floor(log10(real(i2))+1)), List(i3,:)
    write(2,*) BoundTo(i3)
  enddo
  enddo
endif

endsubroutine Analysis
!!!!
!!!!
!!!!
!!!!
subroutine Sanity
call system_clock(tclock3)

! checks that every protein has a type
do i1=1,sum(Nk)
if (ProteinToType(i1)==0) then
  write(*,*) 'Some proteins arent existing!'
  write(*,*) i1
  write(*,*) ProteinToType(i1)
  write(*,*) ProteinToType(max(i1-3,1):min(i1+3,mNk))
  write(*,*) 'positions'
  do i2=(i1-1)*mVk+1,(i1-1)*mVk+1
    write(*,*) i2
    write(*,*) List(i2,:)
    write(*,*) Lattice(List(i2,1),List(i2,2),List(i2,3))
  enddo

  
  write(*,*) MCMove
  write(*,*) iMCMoves
  stop
endif
enddo

! check list->lattice
do i1=1,sum(Nk)
do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
  if (Lattice(List(i2,1),List(i2,2),List(i2,3))/=i2) then
    write(*,*) 'Lattice doesn`t have a segment listed in List'
    write(*,*) MCMove
    write(*,*) iMCMoves
    write(*,*) 'Step Number'
    write(*,*) t
    write(*,*) 'Protein:',BeadToProtein(i2)
    write(*,*) 'Bead:',i2
    write(*,*) 'Listed as at:', List(i2,:)
    write(*,*) 'This Lattice site has:', Lattice(List(i2,1),List(i2,2),List(i2,3))
    write(*,*) BeadToProtein(Lattice(List(i2,1),List(i2,2),List(i2,3)))
    write(*,*) BeadToProtein(i2-1)
    if (Lattice(List(i2,1),List(i2,2),List(i2,3))/=0) then
      write(*,*) List(Lattice(List(i2,1),List(i2,2),List(i2,3)),:)
    endif
    
    do i3=1,BoxSize(1)
    do i4=1,BoxSize(2)
    do i5=1,BoxSize(3)
    if (Lattice(i3,i4,i5)==i2) then
      write(*,*) 'Found segment at a different location:',i2
      write(*,*) i3,i4,i5
      stop
    endif
    enddo
    enddo
    enddo
    write(*,*) 'Segment not found anywhere!',i2
    stop
  endif
enddo
enddo

! check Lattice->List
do i1=1,BoxSize(1)
do i2=1,BoxSize(2)
do i3=1,BoxSize(3)
  if (Lattice(i1,i2,i3)/=0) then
    if (List(Lattice(i1,i2,i3),1)/=i1 .or. List(Lattice(i1,i2,i3),2)/=i2 .or. List(Lattice(i1,i2,i3),3)/=i3) then
      write(*,*) 'List doesn`t match the segment placement on lattice'
      write(*,*) 'Step Number',t
      write(*,*) Lattice(i1,i2,i3)
      write(*,*) List(Lattice(i1,i2,i3),:)
      write(*,*) [i1,i2,i3]
      write(*,*) 'here1'
      write(*,*) sum(Nk)
      write(*,*) rr3(1)
      write(*,*) 
      write(*,*) MCMove
      write(*,*) iMCMoves
      !stop
    endif
  endif
enddo
enddo
enddo

! BoundTo->BoundTo(BoundTo)
do i1=1,sum(Nk)
do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
  if (BoundTo(i2)/=0) then
    if (BoundTo(BoundTo(i2))/=i2) then
      write(*,*) 'Something isn`t bound to what is bound to it'
      write(*,*) 'Step Number:',t
      write(*,*) i2
      write(*,*) BoundTo(i2)
      write(*,*) BoundTo(BoundTo(i2))
      write(*,*) 'here2'
      write(*,*) MCMove
      write(*,*) iMCMoves
      stop
    endif
  endif
enddo
enddo

! Binding Length = 1
do i1=1,sum(Nk)
do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
if (BoundTo(i2)/=0) then
  x1=List(i2,:)
  x2=List(BoundTo(i2),:)
  if (maxval(abs(DistPB(x1,x2)))>1) then
    write(*,*) 'Bound Modules are too far apart to actually bind'
    write(*,*) 'failed on timestep:', t
    write(*,*) i2
    write(*,*) x1
    write(*,*) BoundTo(i2)
    write(*,*) x2
    write(*,*) 'Should be bound to something bound to itself:'
    write(*,*) BoundTo(BoundTo(i2))
    write(*,*) 'Scaled to closest neighbors'
    write(*,*) DistPB(x1,x2)
    write(*,*) 'Last move type:'
    write(*,*) MCMove
    write(*,*) iMCMoves
    write(*,*) 'Bead ',i2,' is on protein ',BeadToProtein(i2)
    write(*,*) 'Bead ',BoundTo(i2),' is on protein ',BeadToProtein(BoundTo(i2))
    stop
  endif
endif
enddo
enddo

!Linker Length
do i1=1,sum(Nk)
do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
x1=List(i2,:)
MST1=BeadToType(i2)
do i4=1,ConnectionN(MST1)
if (Connection(MST1,i4).gt.MST1) then
  x2=List(i2-MST1+Connection(MST1,i4),:)
  x3=DistPB(x1,x2)
  if (maxval(abs(x3)).gt.LL_temp(MST1,i4)) then
    write(*,*) 'Linkers are too far apart, bead and its tethered bead:'
    write(*,*) i2 !bead number
    write(*,*) i2-MST1+Connection(MST1,i4)
    write(*,*) x1 !connected bead
    write(*,*) x2
    write(*,*) 'Move'
    write(*,*) MCmove
    write(*,*) iMCMoves
    write(*,*) 'time'
    write(*,*) t
    write(*,*) 'Distances'
    write(*,*) x1
    write(*,*) x2
    write(*,*) x3
    write(*,*) LL_temp(MST1,i4)
    write(*,*)
    write(*,*) i1,i2,i4
    stop
  endif
endif
enddo
enddo
enddo



! check the total energy
EM=0
do i1=1,sum(Nk)
do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
  MST1=BeadToType(i2)
  if (BoundTo(i2)/=0 .and. BoundTo(i2).gt.i2) then
    i=BeadToType(i2)
    j=BeadToType(BoundTo(i2))
    EM=EM+Energy(i,j)
    if (BeadToProtein(i2)==BeadToProtein(BoundTo(i2))) then
      EM=EM+floor(ESelf*exp(real(-(abs(i-j)-1)/LSteric)))
    endif
    if (Energy(MST1,BeadToType(BoundTo(i2)))==0) then
      write(*,*) 'An interaction has formed between domains that arent allowed to interact'
      write(*,*) t
      write(*,*) iMCMoves
      write(*,*) MCMove
      write(*,*) i1
      write(*,*) i2
      write(*,*) MST1
      write(*,*) BoundTo(i2)
      write(*,*) BeadToType(BoundTo(i2))
      stop
    endif
  endif
  if (UseSpringLinkerEq) then
    do i4=1,SLConN(MST1)
      if (SLCon(MST1,i4)>MST1) then
        x1=List(i2,:)
        x2=List(i2-MST1+SLCon(MST1,i4),:)
        dx=DistPB(x1,x2)
        dx_max=floor(sqrt(real(sum(dx**2)))-SLEq_temp(MST1,i4))!3d distance differnce to Eq Leng for MS1
        EM=EM + SLPot_temp(MST1,i4)*dx_max**2
      endif
    enddo
  endif
enddo
enddo
if (ETC/=EM) then
  write(*,*) 'Energy is wrong.  This is a symptom of a greater problem.'
  write(*,*) iMCMoves
  write(*,*) MCMove
  write(*,*) 'Current Step: ',t
  write(*,*) 'Current Energy EM: ',EM
  write(*,*) 'Energy Saved From Simulation ETC: ',ETC
  stop
endif

call system_clock(tclock4)
tclock_sanity=tclock_sanity+tclock4-tclock3
endsubroutine Sanity
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
endif

if (E_RG_All/=0)  then
close(103)
endif

if (E_RG_Cluster/=0)  then
close(107)
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
if (E_XYZDist_All/=0)  then
  do i1=1,NkN
    do i3=1,3
      close(200*i1+10*i3)
    enddo
  enddo
endif

if (E_XYZDist_Cluster/=0)  then
  do i1=1,NkN
    do i3=1,3
      close(250*i1+10*i3)
    enddo
  enddo
  if(WriteIndiDist_Cluster) then
  do i1=1,NkN
    do i2=1,NkN
      do i3=1,3
        close(100*i1+10000*i2+i3)
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
      close(100*i1+10000*i2)
    enddo
    endif
    close(200+i1)
  enddo
  if(WriteCoreCluster) then
    !Open Core Radial Hist
    do i1=1,NkN
      if(WriteIndiDist_Cluster) then
      do i2=1,NkN
        close(100*i1+10000*i2+1)
        close(100*i1+10000*i2+2)
      enddo
      endif
      close(200+i1+30)
      close(200+i1+60)
    enddo
  endif
endif



if (E_Lattice/=0) then
  close(1)
  close(2)
endif

if (WriteLatticeEnd) then
  open(unit=1,file='EndParms.txt')
  write(1,*) BoxSize
  write(1,*) NkN
  write(1,*) Nk
  write(1,*) Vk
  do i1=1,sum(Vk)
    write(1,*) ConnectionN(i1)
    write(1,*) Connection(i1,1:ConnectionN(i1))
    write(1,*) LL_temp(i1,1:ConnectionN(i1))
  enddo
  close(1)
  
  open(unit=1,file='EndListForm.txt')
  open(unit=2,file='EndBoundTo.txt')
  open(unit=3,file='EndProteinToType.txt')
  do i1=1,sum(Nk)
    do i2=(i1-1)*mVk+1,(i1-1)*mVk+Vk(ProteinToType(i1))
      write(1,*) List(i2,:)
      write(2,*) BoundTo(i2)
    enddo
    write(3,*) ProteinToType(i1)
  enddo
  close(1)
  close(2)
  close(3)
endif



call system_clock(tclock2)
write(*,*)   'Time Spent During Simulation :'
write(*,'(I17,A,I10,A)') tclock2-tclock1, ' milliseconds, ~', (tclock2-tclock1)/3600000, ' hours'
write(*,*)   'Time Spent During Sanity :'
write(*,'(I17,A,I10,A)') tclock_sanity, ' milliseconds, ~', tclock_sanity/3600000, ' hours'
write(*,*)
write(*,*) '        No Moves, Rejected, Proposed, Accepted'
write(*,'(A, I10,I10,I10,I10)') 'TransI ', Rej(1,2), Rej(1,1), Propose(1), Propose(1)-sum(Rej(1,:))
write(*,'(A, I10,I10,I10,I10)') 'TransII', Rej(2,2), Rej(2,1), Propose(2), Propose(2)-sum(Rej(2,:))
write(*,'(A, I10,I10,I10,I10)') 'Rot    ', Rej(3,2), Rej(3,1), Propose(3), Propose(3)-sum(Rej(3,:))
write(*,'(A, I10,I10,I10,I10)') 'SltI   ', Rej(4,2), Rej(4,1), Propose(4), Propose(4)-sum(Rej(4,:))
write(*,'(A, I10,I10,I10,I10)') 'SltII  ', Rej(5,2), Rej(5,1), Propose(5), Propose(5)-sum(Rej(5,:))
write(*,'(A, I10,I10,I10,I10)') 'ClstI  ', Rej(6,2), Rej(6,1), Propose(6), Propose(6)-sum(Rej(6,:))
write(*,'(A, I10,I10,I10,I10)') 'ClstII ', Rej(7,2), Rej(7,1), Propose(7), Propose(7)-sum(Rej(7,:))
write(*,'(A, I10,I10,I10,I10)') 'Grand  ', Rej(8,2), Rej(8,1), Propose(8), Propose(8)-sum(Rej(8,:))
endsubroutine EndofProgram

END PROGRAM First
