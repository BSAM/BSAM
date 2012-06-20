!==========================================================================
! BEARCLAW Boundary Embedded, Adaptive Refinement, Conservation LAW package
!==========================================================================
!
! (c) Copyright Sorin Mitran, 2000
! Department of Applied Mathematics
! University of Washington
! mitran@amath.washington.edu
!
! (c) Copyright Sorin Mitran, 2002
! Department of Mathematics
! University of North Carolina at Chapel Hill
! mitran@amath.unc.edu
!
! Portions of the code
!
! (c) Copyright Steven M. Wise, 2006
! Department of Mathematics
! University of California at Irvine
! swise@math.uci.edu
!
! -----------------------------------------------------------------------
! This software is made available for research and instructional use only.
! You may copy and use this software without charge for these
! non-commercial purposes, provided that the copyright notice and
! associated text is reproduced on all copies. For all other uses,
! including distribution of modified versions, please contact the authors.
!
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             treeops.f90
! Purpose:          Tree operations module.
! Contains:
! Revision History: Ver. 1.0  Oct. 2000 Sorin Mitran
!                   Additions Oct. 2006 Steven Wise
! -----------------------------------------------------------------------
module TreeOps
  !
  ! User definitions of data structures associated with each node
  !
  use NodeInfoDef
  implicit none
  private
  !
  public :: InitForest,KillForest,AddRootLevelNode
  public :: CreateChild,KillNode,DeleteMarkedNode
  public :: ApplyOnLeaves,ApplyOnLevels,ApplyOnForest
  public :: ApplyOnLevel,ApplyOnLevelPairs,ApplyOnChildren
  public :: GetTreeOpsErrCode,GetNodeInfo,SetNodeInfo
  public :: GetCurrentNode,GetNodeNo,GetChildInfo
  public :: GetParentInfo,GetParent,GetSibling,GetChild
  public :: GetLevel,GetSiblingIndex,GetNrOfChildren
  public :: GetRootInfo,SetRootInfo,ExistLevel
  public :: PushForest,PopForest,CreateBelowSeedLevels
  public :: CurrentNodeToYoungest
  !
  ! Error handling
  ! The most recent error code
  integer :: ErrCode
  !
  ! Error codes
  integer, parameter         :: ErrPrint          = 1000
  integer, parameter, public :: err_OK            =    0
  integer, parameter, public :: err_UndefinedNode =  100
  integer, parameter, public :: err_NoParent      =  201
  integer, parameter, public :: err_NoSibling     =  202
  integer, parameter, public :: err_NoChild       =  203
  integer, parameter         :: err_BadLevel      =  301
  integer, parameter         :: IntTrue=0
  integer, parameter         :: IntFalse=-1
  !
  type, public :: Node
    private
    type(NodeInfo), pointer :: Info
    type(Node),     pointer :: Parent
    type(Node),     pointer :: Sibling  ! Next elder sibling
    type(Node),     pointer :: Child
    type(Node),     pointer :: Neighbor ! Nodes not of this parent but
    !
    ! on the same level
    integer :: level         ! This node's level
    integer :: ChildNo       ! Child's ordinal. Eldest=1, Next=2, ...
    integer :: NrOfChildren
    integer :: LeafDist      ! Distance from leaf
    !
    ! This is *not* maintained dynamically as nodes are created/destroyed
    integer :: NodeNo        ! Node number. Useful as ID in debugging
  end type Node
  !
  type :: NodePointer
    type(Node), pointer :: NodePtr
  end type NodePointer
  !
  type, public :: InfoPointer
    type(NodeInfo), pointer :: InfoPtr
  end type InfoPointer
  !
  integer, parameter, public :: FOREST_SEED         = -999
  integer, parameter         :: NO_CHILDREN         = 0
  integer, parameter         :: FIRST_CHILD         = 1
  integer, parameter         :: NOT_A_CHILD         = 0
  integer, parameter         :: BAD_CHILD_NO        = 0
  integer, parameter         :: ZERO_LEAF_NODE_DIST = 0
  !
  logical, parameter :: NoInfoInit   = .false.
  logical, parameter :: PreEvalNext  = .true.
  logical, parameter :: PostEvalNext = .false.
  !
  ! Global variables showing current state of tree traversal
  type(Node), pointer, save                        :: ForestSeed
  type(Node), pointer, save                        :: Root
  type(Node), pointer, save                        :: CurrentNode
  type(NodePointer), dimension(-MaxDepth:MaxDepth) :: Stack
  type(NodePointer), dimension(-MaxDepth:MaxDepth) :: YoungestOnLevel
  type(NodePointer), dimension(-MaxDepth:MaxDepth) :: EldestOnLevel
  integer, dimension(-MaxDepth:MaxDepth)           :: NrOfNodes
  integer                                          :: CurrentLevel
  integer                                          :: GlobalNodeNo
  !
  ! Stack mechanism to allow nested tree traversals
  integer, parameter :: MaxForestStacks=32
  integer, save      :: StackLevel
  type(NodePointer), save, dimension(0:MaxForestStacks) :: ForestStack
  !
  ! Global I/O variables
  integer :: InputUnit,OutputUnit
  !
  ! Global variables showing memory usage (UNUSED)
  !integer, parameter   :: iPrec = selected_int_kind(9)
  !integer (kind=iPrec) :: MemFree, MemAllocated
  !
contains
  ! ---------------------------------------------------------------------------
  ! Interface routines are defined first
  ! ---------------------------------------------------------------------------
  subroutine InitForest(InfoInit)
    implicit none
    logical, intent(in), optional :: InfoInit
    integer                       :: iError
    integer                       :: L
    !
    allocate (Root,ForestSeed,STAT=iError)
    if (iError /= 0) then
      print *,"Error allocating tree/forest roots in InitForest."
      stop
    end if
    StackLevel=0
    ! Seed of all trees
    nullify(ForestSeed%Parent)
    nullify(ForestSeed%Sibling)
    nullify(ForestSeed%Neighbor)
    ForestSeed%Child=>Root
    ForestSeed%level=FOREST_SEED
    ForestSeed%ChildNo=NOT_A_CHILD
    ForestSeed%NrOfChildren=1
    ForestSeed%LeafDist=ZERO_LEAF_NODE_DIST
    ForestSeed%NodeNo=FOREST_SEED
    ! First tree root
    nullify(Root%Sibling)
    nullify(Root%Child)
    nullify(Root%Neighbor)
    Root%Parent=>ForestSeed
    Root%level=rootlevel
    Root%ChildNo=FIRST_CHILD
    Root%NrOfChildren=NO_CHILDREN
    Root%LeafDist=ZERO_LEAF_NODE_DIST
    GlobalNodeNo=1
    Root%NodeNo=GlobalNodeNo
    do L=-MaxDepth,MaxDepth
      nullify(YoungestOnLevel(L)%NodePtr)
      nullify(EldestOnLevel(L)%NodePtr)
      NrOfNodes(L)=0
    end do
    YoungestOnLevel(rootlevel)%NodePtr=>Root
    EldestOnLevel(rootlevel)%NodePtr=>Root
    NrOfNodes(rootlevel)=1
    allocate(Root%Info,STAT=iError)  ! Allocate space for NodeInfo
    if (iError /= 0) then
      print *,"Error allocating tree root Info in InitForest."
      stop
    end if
    ! Set current node pointer
    CurrentNode => Root
    CurrentLevel=rootlevel
    !
    if (.not. present(InfoInit)) return
    if (.not. InfoInit) return
  end subroutine InitForest
  ! ---------------------------------------------------------------------------
  subroutine CreateBelowSeedLevels(MinLevel,InfoInit)
    implicit none
    integer, intent(in)           :: MinLevel
    logical, intent(in), optional :: InfoInit
    integer                       :: iError
    !
    if (MinLevel<-MaxDepth) then
      print *,"Error MinLevel less than MaxDepth in CreateBelowSeedLevels."
      stop
    end if
    !
    allocate(ForestSeed%Info,STAT=iError)  ! Allocate space for NodeInfo
    if (iError /= 0) then
      print *,"Error allocating forest seed Info in CreateBelowSeedLevels."
      stop
    end if
    !
    ForestSeed%Level = -1
    !
    YoungestOnLevel(-1)%NodePtr => ForestSeed
    EldestOnLevel(-1)%NodePtr => ForestSeed
    NrOfNodes(-1) = 1

    if (-1>MinLevel) then
      ForestSeed%ChildNo = 1
      CurrentNode => ForestSeed
      CurrentLevel=-1
      call CreateBelowSeedLevelNode(MinLevel,InfoInit)
      ! Set current node pointer back to root.
      CurrentNode => Root
      CurrentLevel=rootlevel
    end if
    !
    if (.not. present(InfoInit)) return
    if (.not. InfoInit) return
  end subroutine CreateBelowSeedLevels
  ! ---------------------------------------------------------------------------
  recursive subroutine CreateBelowSeedLevelNode(MinLevel,InfoInit)
    implicit none
    integer, intent(in)           :: MinLevel
    logical, intent(in), optional :: InfoInit
    type(Node), pointer           :: BelowSeedNode
    integer                       :: iError, BelowSeedLevel
    !
    ! There are only one of these grids at each level below the root.
    ! Current node is assumed to be root node or lower.
    !
    allocate(BelowSeedNode,STAT=iError)
    if (iError /= err_OK) then
      print *, "Error allocating BelowSeedNode in treeops" // &
               "::CreateBelowSeedLevelNode."
      stop
    end if
    !
    nullify(BelowSeedNode%Parent)
    nullify(BelowSeedNode%Sibling)
    nullify(BelowSeedNode%Child)
    nullify(BelowSeedNode%Neighbor)
    !
    BelowSeedNode%NodeNo = FOREST_SEED
    BelowSeedLevel = CurrentNode%level-1 ! BelowSeedNode is one level up.
    !
    NrOfNodes(BelowSeedLevel) = 1
    BelowSeedNode%Child => CurrentNode   ! Parent is created for CurrentNode
    BelowSeedNode%NrOfChildren = 1       ! BelowSeedNodes have only one child.
    CurrentNode%Parent => BelowSeedNode
    !
    BelowSeedNode%ChildNo = NOT_A_CHILD
    !
    !
    ! Start of level stack
    YoungestOnLevel(BelowSeedLevel)%NodePtr => BelowSeedNode
    ! End of level stack
    EldestOnLevel(BelowSeedLevel)%NodePtr => BelowSeedNode
    !
    ! Set the below-root level counter
    BelowSeedNode%level = BelowSeedLevel
    BelowSeedNode%LeafDist = FOREST_SEED
    !
    ! Finished with internal tree structure maintenance
    ! Allocate space for NodeInfo
    allocate(BelowSeedNode%Info,STAT=iError)
    if (iError /= err_OK) then
      print *,"Error allocating BelowSeedNode Info."
      stop
    end if
    !
    if (BelowSeedLevel>MinLevel) then
      BelowSeedNode%ChildNo = 1
      CurrentNode => BelowSeedNode
      CurrentLevel = BelowSeedLevel
      call CreateBelowSeedLevelNode(MinLevel,InfoInit)
    end if
    !
    if (.not. present(InfoInit)) return
    if (.not. InfoInit) return
  end subroutine CreateBelowSeedLevelNode
  ! ---------------------------------------------------------------------------
  subroutine KillForest
    implicit none
    !
    call KillNode(Root)
  end subroutine KillForest
  ! ---------------------------------------------------------------------------
  subroutine AddRootLevelNode(InfoInit)
    implicit none
    logical, intent(in), optional :: InfoInit
    type(Node), pointer           :: NewRootLevelNode
    integer                       :: iError
    !
    allocate (NewRootLevelNode,STAT=iError)
    if (iError /= 0) then
      print *,"Error allocating node in AddRootLevelNode."
      stop
    end if
    NewRootLevelNode%level=rootlevel
    NrOfNodes(rootlevel)=NrOfNodes(rootlevel)+1
    GlobalNodeNo=GlobalNodeNo+1
    NewRootLevelNode%NodeNo=GlobalNodeNo
    nullify(NewRootLevelNode%Parent)
    nullify(NewRootLevelNode%Child)
    NewRootLevelNode%NrOfChildren=NO_CHILDREN
    NewRootLevelNode%ChildNo=Root%ChildNo+1
    NewRootLevelNode%LeafDist=ZERO_LEAF_NODE_DIST
    ! The newly created node becomes the first node on this level and
    ! the start of the tree
    !  1. Set the Sibling and Neighbor of the newly created node to
    !     point to the old Root node
    NewRootLevelNode%Sibling => Root
    NewRootLevelNode%Neighbor => Root
    !  2. Set the global Root to point to the newly created root level node
    Root => NewRootLevelNode
    ForestSeed%Child=>Root
    CurrentNode=>Root
    !  3. The newly created node starts this level's neighbor list
    YoungestOnLevel(rootlevel)%NodePtr=>NewRootLevelNode
    ! Allocate space for NodeInfo
    allocate(NewRootLevelNode%Info,STAT=iError)
    if (iError /= 0) then
      print *,"Error allocating Info in AddRootLevelNode."
      stop
    end if
    if (.not. present(InfoInit)) return
    if (.not. InfoInit) return
  end subroutine AddRootLevelNode
  ! ---------------------------------------------------------------------------
  subroutine CreateChild(InfoInit,ReadMode)
    implicit none
    logical, intent(in), optional :: InfoInit
    logical, intent(in), optional :: ReadMode
    logical                       :: UpdateYoungest
    type(Node), pointer           :: NewNode
    integer                       :: iError,ThisLevel
    !
    allocate(NewNode,STAT=iError)
    if (iError /= err_OK) then
      print *,"Error allocating NewNode in treeops ::CreateChild."
      stop
    end if
    nullify(NewNode%Parent)
    nullify(NewNode%Sibling)
    nullify(NewNode%Child)
    nullify(NewNode%Neighbor)
    GlobalNodeNo=GlobalNodeNo+1
    NewNode%NodeNo=GlobalNodeNo
    ! Child is one level down
    ThisLevel = CurrentNode%level+1
    if (CurrentLevel /= CurrentNode%level) then
      write(1,*) 'Internal inconsistency in level counters in treeops' // &
                 '::CreateChild'
      write(1,1000) CurrentLevel,CurrentNode%level
      1000 format('Global level=',i2,' CurrentNode%level=',i2)
      stop
    end if
    if (ThisLevel > MaxDepth) then
      print *,'Error in treeops ::CreateChild. Maximum tree depth exceeded'
      stop
    end if
    NrOfNodes(ThisLevel)=NrOfNodes(ThisLevel)+1
    ! Child is born of Current node
    NewNode%Parent => CurrentNode
    CurrentNode%NrOfChildren=CurrentNode%NrOfChildren+1
    ! Point to next elder sibling
    NewNode%Sibling => CurrentNode%Child
    ! First node of this parent?
    if (.not. associated(NewNode%Sibling)) then
      NewNode%ChildNo=FIRST_CHILD
    else
      NewNode%ChildNo=NewNode%Sibling%ChildNo+1
    end if
    ! Determine if we're reading nodes from a file, in which case
    ! we'll be adding elements to the end of the level list
    if (.not. present(ReadMode)) then
      UpdateYoungest = .true.
    else
      if (ReadMode) then
        UpdateYoungest = .false.
      else
        UpdateYoungest = .true.
      end if
    end if
    ! Is this the first child created on this level?
    if (NrOfNodes(ThisLevel)==1) then
      ! Yes. Start stack spanning this level
      YoungestOnLevel(ThisLevel)%NodePtr => NewNode  ! Start of level stack
      EldestOnLevel(ThisLevel)%NodePtr => NewNode    ! End of level stack
      nullify(NewNode%Neighbor)   ! ApplyOnLevel will end on this node
    else
      ! No. Update the Neighbor list spanning this level
      if (UpdateYoungest) then
        ! We're creating a new node during program execution
        ! Set the new node's Neighbor pointer to the previous YoungestOnLevel
        NewNode%Neighbor => YoungestOnLevel(thisLevel)%NodePtr
      else
        ! We're creating a new node while reading from file
        EldestOnLevel(ThisLevel)%NodePtr%Neighbor => NewNode
        nullify(NewNode%Neighbor)
        EldestOnLevel(ThisLevel)%NodePtr => NewNode
      end if
    end if
    CurrentNode%Child => NewNode           ! Parent points to youngest child
    nullify(NewNode%Child)                 ! Just born doesn't have children
    NewNode%NrOfChildren=0
    NewNode%level = ThisLevel              ! Set the child's level counter
    NewNode%LeafDist = ZERO_LEAF_NODE_DIST ! New node is a leaf
    ! Save the most recently created child on this level. Used in maintaining
    ! the Neighbor list spanning a level
    if (UpdateYoungest) YoungestOnLevel(ThisLevel)%NodePtr => NewNode
    ! Finished with internal tree structure maintenance
    ! Allocate space for NodeInfo
    allocate(NewNode%Info,STAT=iError)
    if (iError /= err_OK) then
      print *,"Error allocating NewNode Info."
      stop
    end if
    if (.not. present(InfoInit)) return
    if (.not. InfoInit) return
  end subroutine CreateChild
  ! ---------------------------------------------------------------------------
  ! Kill aNode and all of its children
  ! ---------------------------------------------------------------------------
  recursive subroutine KillNode(aNode)
    implicit none
    type(Node), pointer :: aNode
    type(Node), pointer :: Child,Sibling,Prev
    logical             :: NeighborUpdated
    integer             :: iError,ThisLevel
    !
    ! If aNode has children kill those first
    Child => aNode%Child
    do
      if (.not. associated(Child)) exit
      Sibling => Child%Sibling
      Call KillNode(Child)
      Child => Sibling
    end do
    ! Remove leaf node
    if (associated(aNode%Parent)) then
      aNode%Parent%NrOfChildren = aNode%Parent%NrOfChildren - 1
      if (associated(aNode%Parent%Child,aNode)) then
        ! Update parent child pointer to next eldest
        aNode%Parent%Child => aNode%Sibling
      else
        ! Search for position of this node withing sibling list
        Sibling => aNode%Parent%Child
        do
          if (associated(Sibling%Sibling,aNode)) exit
          Sibling => Sibling%Sibling
        end do
        ! Remove this node from the sibling list
        Sibling%Sibling => aNode%Sibling
      end if
    end if
    ThisLevel = aNode%level
    NrOfNodes(ThisLevel)=NrOfNodes(ThisLevel)-1
    ! Update the Neighbor list spanning this level
    NeighborUpdated=.false.
    ! Was this node the last node on this level?
    if (NrOfNodes(ThisLevel)==0) then
      nullify(YoungestOnLevel(ThisLevel)%NodePtr)
      nullify(EldestOnLevel(ThisLevel)%NodePtr)
      NeighborUpdated=.true.
    end if
    ! Was it the start of the Neighbor list?
    if ((.not. NeighborUpdated) .and. &
      (associated(aNode,YoungestOnLevel(ThisLevel)%NodePtr))) then
      ! Set the start of the list to the next node on this level
      YoungestOnLevel(ThisLevel)%NodePtr => &
                      YoungestOnLevel(ThisLevel)%NodePtr%Neighbor
      NeighborUpdated=.true.
    end if
    if (.not. NeighborUpdated) then
      ! aNode definitely has a previous element. Find it.
      Prev => YoungestOnLevel(ThisLevel)%NodePtr
      do while (.not. associated(Prev%Neighbor,aNode))
        Prev => Prev%Neighbor
      end do
    end if
    ! Was it the end of the Neighbor stack?
    if ((.not. NeighborUpdated) .and. &
      (associated(aNode,EldestOnLevel(ThisLevel)%NodePtr))) then
      nullify(Prev%Neighbor)
      EldestOnLevel(ThisLevel)%NodePtr => Prev
      NeighborUpdated = .true.
    end if
    if (.not. NeighborUpdated) then
      ! If we're here, aNode is neither the last nor the first in Neighbor list
      Prev%Neighbor => aNode%Neighbor ! Skip aNode in Neighbor list
    end if
    deallocate(aNode%Info,STAT=iError)
    if (iError /= 0) then
      print *,"Error deallocating aNode Info."
      stop
    end if
    deallocate(aNode,STAT=iError)
    if (iError /= 0) then
      print *,"Error deallocating aNode."
      stop
    end if
  end subroutine KillNode
  ! ---------------------------------------------------------------------------
  integer function DeleteMarkedNode(info,dummy)
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    DeleteMarkedNode = err_ok
    !
    if (.not. info%tobedeleted) return
    !
    call KillNode(currentnode)
    !
  end function DeleteMarkedNode
  ! ---------------------------------------------------------------------------
  subroutine CurrentNodeToYoungest(level)
    implicit none
    integer, intent(in) :: level
    !
    currentnode => youngestonlevel(level)%nodeptr
    currentlevel = level
    !
  end subroutine CurrentNodeToYoungest
  ! ---------------------------------------------------------------------------
  subroutine ApplyOnLevel(level,f,fparam)
    implicit none
    interface
      integer function f(info,param)
        use NodeInfoDef
        implicit none
        type(nodeinfo)  :: info
        type(funcparam) :: param
      end function f
    end interface
    integer,         intent(in) :: level
    type(funcparam), intent(in) :: fparam
    type(node), pointer         :: nextnode
    integer                     :: ferrcode
    !
    if (abs(level)>maxdepth) then
      print *,'Error in ApplyOnLevel. Maximum tree depth exceeded.'
      stop
    end if
    !
    currentnode => youngestonlevel(level)%nodeptr
    currentlevel = level
    !
    do
      if (.not. associated(currentnode)) exit
      nextnode => currentnode%neighbor
      ferrcode = f(currentnode%info,fparam)
      if (ferrcode/=err_ok) then
        errcode = ferrcode
        if (ferrcode<errprint) print *, 'Error in ApplyOnLevel. ferrcode=', &
          ferrcode
        return
      end if
      currentnode => nextnode
    end do
    !
  end subroutine ApplyOnLevel
  ! ---------------------------------------------------------------------------
  subroutine ApplyOnLevelPairs(L,f,fparam)
    implicit none
    interface
      integer function f(Info1,Info2,Param)
        use NodeInfoDef
        implicit none
        type(NodeInfo)  :: Info1,Info2
        type(FuncParam) :: Param
      end function f
    end interface
    integer,         intent(in) :: L
    type(FuncParam), intent(in) :: fparam
    type(Node), pointer         :: Node1,Node2
    integer                     :: fErrCode
    !
    if (abs(L) > MaxDepth) then
      print *,'Error in ApplyOnLevelPairs. Maximum tree depth exceeded.'
      stop
    end if
    Node1 => YoungestOnLevel(L)%NodePtr
    do
      if (.not. associated(Node1)) exit
      Node2 => Node1%Neighbor
      do
        if (.not. associated(Node2)) exit
        fErrCode=f(Node1%Info,Node2%Info,fparam)
        if (fErrCode /= err_OK) then
          ErrCode = fErrCode
          if (fErrCode<ErrPrint ) &
            print *, 'Error in ApplyOnLevelPairs. fErrCode=',fErrCode
          return
        end if
        Node2=>Node2%Neighbor
      end do
      Node1=>Node1%Neighbor
    end do
  end subroutine ApplyOnLevelPairs
  ! ---------------------------------------------------------------------------
  subroutine ApplyOnChildren(f,fparam)
    implicit none
    interface
      integer function f(Info,Param)
        use NodeInfoDef
        ! Interface declarations
        type(NodeInfo) :: Info
        type(FuncParam) :: Param
      end function f
    end interface
    type(FuncParam), intent(in) :: fparam
    integer                     :: fErrCode
    type(Node), pointer         :: SaveCurrentNode
    !
    ! Apply f on all children of the current node
    SaveCurrentNode=>CurrentNode
    CurrentNode=>CurrentNode%Child
    do
      if (.not. associated(CurrentNode)) exit
      fErrCode=f(CurrentNode%Info,fparam)
      if (fErrCode /= err_OK) then
        ErrCode = fErrCode
        if (fErrCode<ErrPrint) &
          print *, 'Error in ApplyOnChildren. fErrCode=',fErrCode
        return
      end if
      CurrentNode=>CurrentNode%Sibling
    end do
    CurrentNode=>SaveCurrentNode
  end subroutine ApplyOnChildren
  ! ---------------------------------------------------------------------------
  subroutine ApplyOnForest(f,fparam,PreEval)
    !
    ! Apply function f on all nodes within tree
    !
    implicit none
    interface
      integer function f(Info,Param)
        use NodeInfoDef
        ! Interface declarations
        type(NodeInfo) :: Info
        type(FuncParam) :: Param
      end function f
    end interface
    type(FuncParam), intent(in)           :: fparam
    logical,         intent(in), optional :: PreEval
    logical                               :: EvalNextBefore_f
    !
    ! First executable statement
    CurrentNode => Root
    CurrentLevel=rootlevel
    if (present(PreEval)) then
      EvalNextBefore_f=PreEval
    else
      EvalNextBefore_f=.true.
    end if
    call ForestTraversal(Root,f,TrueCond,fparam,EvalNextBefore_f)
  end subroutine ApplyOnForest
  ! ---------------------------------------------------------------------------
  subroutine ApplyOnLevels(f,fparam,PreEval)
    !
    ! Level  0: Root
    ! Level  1: 1st generation children
    ! ....
    ! Level  n: n-th generation children
    ! Level -1: Parents of leaves
    ! Level -2: Grandparents of leaves
    !
    implicit none
    interface
      integer function f(Info,Param)
        use NodeInfoDef
        implicit none
        type(NodeInfo)  :: Info
        type(FuncParam) :: Param
      end function f
    end interface
    type(FuncParam), intent(in)           :: fparam
    logical,         intent(in), optional :: PreEval
    logical                               :: EvalNextBefore_f
    !
    CurrentNode => Root
    CurrentLevel=rootlevel
    if (present(PreEval)) then
      EvalNextBefore_f=PreEval
    else
      EvalNextBefore_f=.true.
    end if
    call ForestTraversal(Root,f,LevelCond,fparam,EvalNextBefore_f)
  end subroutine ApplyOnLevels
  ! ---------------------------------------------------------------------------
  subroutine ApplyOnLeaves(f,fparam,PreEval)
    !
    ! Leaves = youngest generation in existence
    !
    implicit none
    interface
      integer function f(Info,Param)
        use NodeInfoDef
        ! Interface declarations
        type(NodeInfo) :: Info
        type(FuncParam) :: Param
      end function f
    end interface
    type(FuncParam), intent(in)           :: fparam
    logical,         intent(in), optional :: PreEval
    logical                               :: EvalNextBefore_f
    !
    ! First executable statement
    CurrentNode => Root
    CurrentLevel=rootlevel
    if (present(PreEval)) then
      EvalNextBefore_f=PreEval
    else
      EvalNextBefore_f=.true.
    end if
    call ForestTraversal(Root,f,LeafCond,fparam,EvalNextBefore_f)
  end subroutine ApplyOnLeaves
  ! ---------------------------------------------------------------------------
  ! Routines interior to the module
  ! ---------------------------------------------------------------------------
  subroutine ForestTraversal(aNode,f,cond,fparam,EvalNextBefore_f)
    !
    ! Traverse the root level nodes
    !
    implicit none
    interface
      integer function f(Info,Param)
        use NodeInfoDef
        implicit none
        type(NodeInfo)  :: Info
        type(FuncParam) :: Param
      end function f
      logical function cond()
      end function cond
    end interface
    type(Node), pointer, intent(in) :: aNode
    type(FuncParam),     intent(in) :: fparam
    logical,             intent(in) :: EvalNextBefore_f
    type(Node), pointer             :: Next
    !
    Next => aNode
    do
      if (.not. associated(Next)) exit
      call TreeTraversal(Next,f,cond,fparam,EvalNextBefore_f)
      Next=>Next%Sibling
    end do
  end subroutine ForestTraversal
  ! ---------------------------------------------------------------------------
  recursive subroutine TreeTraversal(aNode,f,cond,fparam,EvalNextBefore_f)
    !
    ! The core tree function in terms of which all others are defined:
    !   1) Traverse tree
    !   2) If condition *Cond* is satisfied, apply function *Func* to node
    !      *aNode*
    !   3) Update *Level* counter
    !
    implicit none
    interface
      integer function f(Info,Param)
        use NodeInfoDef
        implicit none
        type(NodeInfo)  :: Info
        type(FuncParam) :: Param
      end function f
      logical function cond()
      end function cond
    end interface
    type(Node), target, intent(in) :: aNode
    type(FuncParam),    intent(in) :: fparam
    logical,            intent(in) :: EvalNextBefore_f
    type(Node), pointer            :: Next
    integer                        :: fErrCode,level
    !
    ! First executable statement
    if (ErrCode /= err_OK) return   ! Go up recursion upon error
    CurrentNode => aNode            ! Set current node and maintain stack
    ! so other functions know what to work on
    Stack(CurrentLevel)%NodePtr => CurrentNode
    if (EvalNextBefore_f) Next => CurrentNode%Child
    ! Do work on this node if condition cond() is satisfied
    if (cond()) then
      fErrCode=f(CurrentNode%Info,fparam)
      if (fErrCode /= err_OK) then
        ErrCode = fErrCode
        if (fErrCode<ErrPrint) &
          print *, 'Error on applying f. fErrCode=', fErrCode
        return
      end if
    end if
    if (.not. EvalNextBefore_f) Next => CurrentNode%Child
    ! Find next node to work on
    do
      if (.not. associated(Next)) then
        ! Reached a leaf
        CurrentLevel=max(CurrentLevel-1,rootlevel)
        level=CurrentLevel
        ! Go up tree, exit clause from do loop
        exit
      else
        CurrentLevel=CurrentLevel+1
        level=CurrentLevel  ! Node has children
        go down
        call TreeTraversal(Next,f,cond,fparam,EvalNextBefore_f)
        ! Reset global context from local context when returning from recursion
        CurrentNode => aNode
        Next => Next%Sibling       ! After children have been exhausted
        ! work on next sibling
      end if
    end do
  end subroutine TreeTraversal
  ! ---------------------------------------------------------------------------
  logical function TrueCond()
    implicit none
    TrueCond = .true.
  end function TrueCond
  ! ---------------------------------------------------------------------------
  logical function LevelCond()
    implicit none
    LevelCond = .true.
  end function LevelCond
  ! ---------------------------------------------------------------------------
  logical function LeafCond()
    implicit none
    if (.not. associated(CurrentNode%Child)) then
      LeafCond = .true.
    else
      LeafCond = .false.
    end if
  end function LeafCond
  ! ---------------------------------------------------------------------------
  integer function SetChildNo(Info,Param)
    implicit none
    type(NodeInfo)      :: Info
    type(FuncParam)     :: Param
    type(Node), pointer :: Youngest
    integer             :: ListNo, NrOfChildren
    !
    if (CurrentNode%Level==rootlevel) then
      CurrentNode%ChildNo=FIRST_CHILD
      SetChildNo=err_OK
      return
    end if
    if (CurrentNode%ChildNo==BAD_CHILD_NO) then
      ! Number all children within this family
      Youngest => CurrentNode%Parent%Child
      ListNo=1
      NrOfChildren=CurrentNode%Parent%NrOfChildren
      Youngest%ChildNo=NrOfChildren
      do
        if (.not. associated(Youngest%Sibling)) exit
        Youngest => Youngest%Sibling
        ListNo=ListNo+1
        Youngest%ChildNo=NrOfChildren-ListNo+1
      end do
    end if
    SetChildNo=err_OK
  end function SetChildNo
  ! ---------------------------------------------------------------------------
  integer function GetTreeOpsErrCode()
    implicit none
    GetTreeOpsErrCode = ErrCode
  end function GetTreeOpsErrCode
  ! ---------------------------------------------------------------------------
  integer function GetNodeInfo(aNode,aNodeInfo)
    implicit none
    type(Node), pointer :: aNode
    type(NodeInfo)      :: aNodeInfo
    !
    if (associated(aNode)) then
      aNodeInfo = aNode%Info
      GetNodeInfo = err_OK
    else
      GetNodeInfo = err_UndefinedNode
    end if
  end function GetNodeInfo
  ! ---------------------------------------------------------------------------
  integer function SetNodeInfo(aNode,aNodeInfo)
    implicit none
    type(Node), pointer :: aNode
    type(NodeInfo)      :: aNodeInfo
    !
    if (associated(aNode)) then
      aNode%Info = aNodeInfo
      SetNodeInfo = err_OK
    else
      SetNodeInfo = err_UndefinedNode
    end if
  end function SetNodeInfo
  ! ---------------------------------------------------------------------------
  integer function GetCurrentNodeInfo(aNodeInfo)
    implicit none
    type(NodeInfo) :: aNodeInfo
    if (associated(CurrentNode)) then
      aNodeInfo = CurrentNode%Info
      GetCurrentNodeInfo = err_OK
    else
      GetCurrentNodeInfo = err_UndefinedNode
    end if
  end function GetCurrentNodeInfo
  ! ---------------------------------------------------------------------------
  integer function GetNodeNo()
    implicit none
    GetNodeNo=CurrentNode%NodeNo
  end function GetNodeNo
  ! ---------------------------------------------------------------------------
  integer function SetCurrentNodeInfo(aNodeInfo)
    implicit none
    type(NodeInfo) :: aNodeInfo
    !
    if (associated(CurrentNode)) then
      CurrentNode%Info = aNodeInfo
      SetCurrentNodeInfo = err_OK
    else
      SetCurrentNodeInfo = err_UndefinedNode
    end if
  end function SetCurrentNodeInfo
  ! ---------------------------------------------------------------------------
  integer function GetRootInfo(aNodeInfo)
    implicit none
    type(NodeInfo), pointer :: aNodeInfo
    !
    if (associated(Root)) then
      aNodeInfo => Root%Info
      GetRootInfo = err_OK
    else
      GetRootInfo = err_UndefinedNode
    end if
  end function GetRootInfo
  ! ---------------------------------------------------------------------------
  integer function SetRootInfo(aNodeInfo)
    implicit none
    type(NodeInfo) :: aNodeInfo
    !
    if (associated(Root)) then
      Root%Info = aNodeInfo
      SetRootInfo = err_OK
    else
      SetRootInfo = err_UndefinedNode
    end if
  end function SetRootInfo
  ! ---------------------------------------------------------------------------
  integer function GetParentInfo(aNodeInfo)
    implicit none
    type(NodeInfo), pointer :: aNodeInfo
    !
    if (associated(CurrentNode) .and. &
      associated(CurrentNode%Parent) .and.  &
      CurrentNode%Parent%level > FOREST_SEED) then
      aNodeInfo => CurrentNode%Parent%Info
      GetParentInfo = err_OK
    else
      GetParentInfo = err_NoParent
    end if
  end function GetParentInfo
  ! ---------------------------------------------------------------------------
  integer function GetCurrentNode(aNode)
    implicit none
    type(Node), pointer :: aNode
    !
    GetCurrentNode = err_OK
    aNode => CurrentNode
  end function GetCurrentNode
  ! ---------------------------------------------------------------------------
  integer function GetParent(aNode,Parent)
    implicit none
    type(Node), pointer :: aNode
    type(Node), pointer :: Parent
    !
    if (associated(aNode)) then
      if (associated(aNode%Parent)) then
        Parent => aNode%Parent
        GetParent = err_OK
      else
        nullify(Parent)
        GetParent = err_NoParent
      end if
    else
      GetParent = err_UndefinedNode
    end if
  end function GetParent
  ! ---------------------------------------------------------------------------
  integer function GetSibling(aNode,Sibling)
    implicit none
    type(Node), pointer :: aNode
    type(Node), pointer :: Sibling
    !
    if (associated(aNode)) then
      if (associated(aNode%Sibling)) then
        Sibling => aNode%Sibling
        GetSibling = err_OK
      else
        nullify(Sibling)
        GetSibling = err_NoSibling
      end if
    else
      GetSibling = err_UndefinedNode
    end if
  end function GetSibling
  ! ---------------------------------------------------------------------------
  integer function GetChild(aNode,Child)
    implicit none
    type(Node), pointer :: aNode
    type(Node), pointer :: Child
    !
    if (associated(aNode)) then
      if (associated(aNode%Child)) then
        Child => aNode%Child
        GetChild = err_OK
      else
        nullify(Child)
        GetChild = err_NoChild
      end if
    else
      GetChild = err_UndefinedNode
    end if
  end function GetChild
  ! ---------------------------------------------------------------------------
  integer function GetChildInfo(aNodeInfo)
    implicit none
    type(NodeInfo), pointer :: aNodeInfo
    !
    if (associated(CurrentNode) .and. associated(CurrentNode%Child)) then
      aNodeInfo => CurrentNode%Child%Info
      GetChildInfo = err_OK
    else
      GetChildInfo = err_NoChild
    end if
  end function GetChildInfo
  ! ---------------------------------------------------------------------------
  integer function GetLevel(ThisLevel)
    implicit none
    integer :: ThisLevel
    !
    ThisLevel = CurrentLevel
    if (ThisLevel /= CurrentNode%level) then
      print *,'Internal inconsistency in level counter'
      stop
    end if
    GetLevel = err_OK
  end function GetLevel
  ! ---------------------------------------------------------------------------
  integer function GetSiblingIndex(SiblingIndex)
    implicit none
    integer             :: SiblingIndex
    integer             :: NrOfSiblings
    type(Node), pointer :: NextSibling,ThisSibling
    !
    ThisSibling => CurrentNode
    SiblingIndex=1
    NrOfSiblings=1
    if (.not. associated(CurrentNode%Parent)) then
      ! We are on the forest seed level
      GetSiblingIndex = err_OK
      return
    else
      ! We are below the forest seed level
      NextSibling => CurrentNode%Parent%Child
      do
        if (associated(NextSibling,ThisSibling)) exit
        SiblingIndex=SiblingIndex+1
        NextSibling => NextSibling%Sibling
      end do
      GetSiblingIndex = err_OK
      NextSibling => CurrentNode%Parent%Child
      do
        if (.not. associated(nextSibling%Sibling)) exit
        NrOfSiblings=NrOfSiblings+1
        NextSibling => NextSibling%Sibling
      end do
    end if
    if ((NrOfSiblings+1-SiblingIndex) /= CurrentNode%ChildNo) then
      print *,'Internal inconsistency in ChildNo counter'
      print *,'Level=',CurrentNode%level
      print *,'ChildNo=',CurrentNode%ChildNo
      print *,'SiblingIndex=',SiblingIndex
      print *,'NrOfSiblings=',NrOfSiblings
      stop
    end if
  end function GetSiblingIndex
  ! ---------------------------------------------------------------------------
  integer function GetNrOfChildren(NrOfChildren)
    implicit none
    integer             :: NrOfChildren
    type(Node), pointer :: NextChild
    !
    NrOfChildren=0
    NextChild => CurrentNode%Child
    do
      if (.not. associated(NextChild)) exit
      NrOfChildren=NrOfChildren+1
      NextChild => NextChild%Sibling
    end do
    if (NrOfChildren /= CurrentNode%NrOfChildren) then
      write(1,1001)
      write(1,1002) NrOfChildren, CurrentNode%NrOfChildren
    end if
    GetNrOfChildren = err_OK
    !
    ! formats
    1001 format('GetNrOfChildren: Internal inconsistency in NrOfChildren')
    1002 format('Computed from links:',i3,' stored in Node:',i3)
  end function GetNrOfChildren
  ! ---------------------------------------------------------------------------
  logical function ExistLevel(level)
    implicit none
    integer :: level
    !
    if (abs(level)>MaxDepth) then
      ErrCode=err_BadLevel
      print *,'Bad level in call to ExistLevel'
      stop
    end if
    if (associated(EldestOnLevel(level)%NodePtr)) then
      ExistLevel=.true.
    else
      ExistLevel=.false.
    end if
  end function ExistLevel
  ! ---------------------------------------------------------------------------
  subroutine PushForest
    implicit none
    !
    ! Save current tree traversal state to allow a subsidiary tree traversal
    ! to take place
    !
    ForestStack(StackLevel)%NodePtr => CurrentNode
    StackLevel=StackLevel+1
    if (StackLevel>MaxForestStacks) then
      print *,'Too many forest traversal recursions. ' // &
              'Increase MaxForestStacks in treeops.f90'
      stop
    end if
    CurrentNode=>Root
  end subroutine PushForest
  ! ---------------------------------------------------------------------------
  subroutine PopForest
    implicit none
    !
    ! Restore tree traversal state on return from a subsidiary tree traversal
    !
    StackLevel=StackLevel-1
    if (StackLevel<0) then
      print *,'PopForest requested without prior PushForest'
      stop
    end if
    CurrentNode => ForestStack(StackLevel)%NodePtr
  end subroutine PopForest
  ! ---------------------------------------------------------------------------
end module TreeOps
