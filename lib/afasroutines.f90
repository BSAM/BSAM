!==========================================================================
! BSAM: Block-Structured Adaptive Multigrid Solver
!==========================================================================
!
! (c) Copyright Steven M. Wise, 2006
! Department of Mathematics
! University of California at Irvine
! swise@math.uci.edu
!
! (c) Copyright Steven M. Wise, 2007
! Department of Mathematics
! University of Tennessee
! swise@math.utk.edu
!
! -----------------------------------------------------------------------
! This software is made available for research and instructional use only.
! You may copy and use this software without charge for these
! non-commercial purposes, provided that the copyright notice and
! associated text is reproduced on all copies. For all other uses,
! including distribution of modified versions, please contact the authors.
!
! Commercial use is strictly forbidden without permission.
!
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             afasroutines.f90
! Purpose:          Adaptive FAS Multigrid module
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
module AFASRoutines
contains
  ! ---------------------------------------------------------------------------
  subroutine MultigridIterations
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel
    use problem, only: MultiGridUserBreak
    implicit none
    type(funcparam) :: dummy
    integer         :: itmg, level
    real            :: residual,oldresidual=1.0e30
    !
    call FillDown(solutionfield)
    !
    ! If the auxiliary fields are updated in the vcycle, perform here:
    do level = finestlevel, minlevel, -1
      call ApplyOnLevel(level,UpdateAuxInVcycle,dummy)
    end do
    !
    ! Maybe implement this later, but it needs more work:
    !call FillDown(auxiliaryfield)
    !
    do level = finestlevel, rootlevel, -1
      call ApplyOnLevel(level,GetSourceFunction,dummy)
    end do
    !
    ! If the source function is updated after each vcycle, perform here first:
    do level = finestlevel, rootlevel, -1
      call ApplyOnLevel(level,UpdateSourceFunction,dummy)
    end do
    !
    call FillDown(sourcefield)
    !
    ! Perform Adaptive Full Aproximation Scheme Vcycles on the grid hierarchy:
    vcycleloop: do itmg = 1, maxvcycles
      !
      call AFASVcycle(finestlevel)
      !
      call FillDown(solutionfield)
      !
      ! If the auxiliary fields are updated in the vcycle, perform here:
      if (modulo(itmg-1, updateauxfreq)==0) then
        do level = finestlevel, minlevel, -1
          call ApplyOnLevel(level, UpdateAuxInVcycle, dummy)
        end do
      end if
      !
      ! Maybe implement this later, but it needs more work:
      !  call FillDown(auxiliaryfield)
      !
      ! If the source function is updated after each vcycle, perform here:
      do level = finestlevel, rootlevel, -1
        call ApplyOnLevel(level, UpdateSourceFunction, dummy)
      end do
      !
      call FillDown(sourcefield)
      !
      ! Note that it is necessary to have an up-to-date source function for
      ! calculating the residual:
      residual = ErrorAFAS(finestlevel)
      print *, 'it =', itmg, ' residual =', residual
      !
      ! Check for convergence
      if (residual < qerrortol) then
        lconverged=.true.
        exit vcycleloop
      endif
      if (MultiGridUserBreak(eqtag, residual, oldresidual)) exit vcycleloop
      oldresidual = residual
    end do vcycleloop
    !
  end subroutine MultigridIterations
  ! ---------------------------------------------------------------------------
  function ErrorAFAS(level) result(errorresult)
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel, GetRootInfo
    implicit none
    integer, intent(in)          :: level
    integer                      :: i, ierror, nrvars
    type(nodeinfo), pointer      :: rootinfo
    type(funcparam)              :: dummy
    real                         :: errorresult
    real, dimension(1:maxnrvars) :: componenterror
    !
    select case(errortype)
    case(1)
      !
      ! Blended error calculation:
      integralresult(1:2) = 0.0
      !
      call ApplyOnLevel(level,L2Error,dummy)
      !
      errorresult = sqrt(integralresult(1)/integralresult(2))
    case(2)
      !
      ! Component errors:
      ierror = GetRootInfo(rootinfo)
      nrvars = rootinfo%nrvars
      integralresult(2) = 0.0
      componentintegral(1:nrvars) = 0.0
      !
      call ApplyOnLevel(level,L2ComponentErrors,dummy)
      !
      componenterror(1:nrvars) = sqrt(componentintegral(1:nrvars) &
                                 / integralresult(2))
      !
      do i = 1, nrvars
        print *, i, componenterror(i)
      end do
      !
      errorresult = maxval(componenterror(1:nrvars))
    case default
      print *, 'ErrorAFAS: Only errortype = 1,2 are supported.'
      stop
      !
    end select
    !
  end function ErrorAFAS
  ! ---------------------------------------------------------------------------
  integer function L2Error(info,dummy)
    !
    ! Adds the L2 error on this grid to the total.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer                       :: nrvars
    integer, dimension(1:maxdims) :: mx
    !
    L2Error = err_ok
    !
    nrvars = info%nrvars
    mx = 1
    mx(1:ndims) = info%mx(1:ndims)
    !
    call Residual(info)
    !
    integralresult(1) = SquareSum(info%rf(1:mx(1),1:mx(2),1:mx(3),1:nrvars)) &
                        + integralresult(1)
    !
    integralresult(2) = real(nrvars*product(mx(1:ndims))) + integralresult(2)
    !
  end function L2Error
  ! ---------------------------------------------------------------------------
  integer function L2ComponentErrors(info,dummy)
    !
    ! Adds the L2 component error on this grid to the component total.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer                       :: i, nrvars
    integer, dimension(1:maxdims) :: mx
    !
    L2ComponentErrors = err_ok
    !
    nrvars = info%nrvars
    mx = 1
    mx(1:ndims) = info%mx(1:ndims)
    !
    call Residual(info)
    !
    do i = 1, nrvars
      componentintegral(i) = SquareSum(info%rf(1:mx(1),1:mx(2),1:mx(3),i:i)) &
                             + componentintegral(i)
    end do
    !
    integralresult(2) = real(product(mx(1:ndims))) + integralresult(2)
    !
  end function L2ComponentErrors
  ! ---------------------------------------------------------------------------
  function SquareSum(rf) result(squaresumresult)
    !
    ! Dimensionally invariant square integral.
    !
    use NodeInfoDef
    implicit none
    real, dimension(:,:,:,:), intent(in) :: rf
    real                                 :: squaresumresult
    !
    squaresumresult = sum(rf*rf)
    !
  end function SquareSum
  ! ---------------------------------------------------------------------------
  subroutine FillDown(field)
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel
    use Boundary, only: SetGhost
    implicit none
    integer, intent(in) :: field
    type(funcparam)     :: dummy
    integer             :: level
    !
    dummy%iswitch = field
    !
    do level = finestlevel, minlevel+1, -1
      call ApplyOnLevel(level,FillDownLevel,dummy)
      if (field==solutionfield) call SetGhost(level-1,0)
    end do
    !
  end subroutine FillDown
  ! ---------------------------------------------------------------------------
  integer function FillDownLevel(info,dummy)
    use NodeInfoDef
    use TreeOps,       only: err_ok, GetParentInfo
    use GridUtilities, only: Restriction2D, Restriction3D
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    type(nodeinfo), pointer           :: parent
    integer                           :: field, ierror, nrvars
    integer, dimension(1:maxdims)     :: cmx, mx
    integer, dimension(1:maxdims,1:2) :: mb
    !
    FillDownLevel = err_ok
    !
    if (info%tobedeleted) return
    !
    field = dummy%iswitch
    !
    nrvars = info%nrvars
    mx     = 1
    cmx    = 1
    mb     = 1
    mx(1:ndims)     = info%mx(1:ndims)
    cmx(1:ndims)    = mx(1:ndims)/2
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    !
    ierror = GetParentInfo(parent)
    !
    select case(field)
    case(sourcefield)
      select case(ndims)
      case(2)
        parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nrvars) &
                 = Restriction2D(info%f(1:mx(1),1:mx(2),1,1:nrvars))
      case(3)
        parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars) &
                 = Restriction3D(info%f(1:mx(1),1:mx(2),1:mx(3),1:nrvars))
      case default
        print *, 'FillDownLevel: Only ndims = 2,3 are supported.'
        stop
      end select
    case(solutionfield)
      select case(ndims)
      case(2)
        parent%q(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nrvars) &
                 = Restriction2D(info%q(1:mx(1),1:mx(2),1,1:nrvars))
      case(3)
        parent%q(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars) &
                 = Restriction3D(info%q(1:mx(1),1:mx(2),1:mx(3),1:nrvars))
      case default
        print *, 'FillDownLevel: Only ndims = 2,3 are supported.'
        stop
      end select
      !  case(auxiliaryfield)
      !    if (maux>0) then
      !      select case(ndims)
      !        case(2)
      !          parent%aux(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:mmaux) &
      !            = Restriction2D(info%aux(1:mx(1),1:mx(2),1,1:mmaux))
      !        case(3)
      !          parent%aux(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:mmaux) &
      !            = Restriction3D(info%aux(1:mx(1),1:mx(2),1:mx(3),1:mmaux))
      !        case default
      !          print *, 'FillDownLevel: Only ndims = 2,3 are supported.'
      !          stop
      !      end select
      !    end if
    end select
    !
  end function FillDownLevel
  ! ---------------------------------------------------------------------------
  integer function GetSourceFunction(info,dummy)
    !
    ! Dimensionally invariant source (or partial source) function routine.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    use Problem, only: Source2D, Source3D
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer                       :: ll, maux, mbc, mmaux, nrvars
    integer, dimension(1:maxdims) :: amx, aul, mx, ul
    real                          :: h
    real, dimension(1:maxdims)    :: xlower
    !
    GetSourceFunction = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx     = 1
    h      = info%dx(1)
    mx(1:ndims)     = info%mx(1:ndims)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amx   = 1
    if (maux>0) amx(1:ndims) = mx(1:ndims)
    !
    mbc = info%mbc
    ll  = 1-mbc
    ul  = 1
    aul = 1
    ul(1:ndims)  = mx(1:ndims)+mbc
    aul(1:ndims) = amx(1:ndims)+mbc
    !
    select case(ndims)
    case(2)
      info%ftmp(1:mx(1),1:mx(2),1,1:nrvars) &
                = Source2D(info%q(   ll: ul(1),ll: ul(2),1,1:nrvars), &
                           info%qold(ll: ul(1),ll: ul(2),1,1:nrvars), &
                           info%aux( ll:aul(1),ll:aul(2),1,1:mmaux ), &
                           ll,mx(1:2),nrvars,maux,mbc,h,xlower(1:2))
    case(3)
      info%ftmp(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
                = Source3D(info%q(   ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                           info%qold(ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                           info%aux( ll:aul(1),ll:aul(2),ll:aul(3),1:mmaux ), &
                           ll,mx(1:3),nrvars,maux,mbc,h,xlower(1:3))
    case default
      print *, 'GetSourceFunction: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end function GetSourceFunction
  ! ---------------------------------------------------------------------------
  integer function UpdateAuxInVcycle(info,dummy)
    !
    ! Used if the auxiliary variables are updated after each Vcycle iteration.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    use Problem, only: UpdateAuxVcycle2D, UpdateAuxVcycle3D
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer                       :: ll, maux, mbc, mmaux, nrvars
    integer, dimension(1:maxdims) :: amx, aul, mx, ul
    real                          :: h
    real, dimension(1:maxdims)    :: xlower
    !
    UpdateAuxInVcycle = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx     = 1
    h      = info%dx(1)
    mx(1:ndims)     = info%mx(1:ndims)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amx   = 1
    if (maux>0) amx(1:ndims) = mx(1:ndims)
    !
    mbc = info%mbc
    ll  = 1-mbc
    ul  = 1
    aul = 1
    ul(1:ndims)  = mx(1:ndims)+mbc
    aul(1:ndims) = amx(1:ndims)+mbc
    !
    if (maux<=0) return
    !
    select case(ndims)
    case(2)
      call UpdateAuxVcycle2D(info%q(   ll: ul(1),ll: ul(2),1,1:nrvars), &
                             info%qold(ll: ul(1),ll: ul(2),1,1:nrvars), &
                             info%aux( ll:aul(1),ll:aul(2),1,1:mmaux ), &
                             ll,mx(1:2),nrvars,maux,mbc,h,xlower(1:2))
    case(3)
      call UpdateAuxVcycle3D( &
                         info%q(   ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                         info%qold(ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                         info%aux( ll:aul(1),ll:aul(2),ll:aul(3),1:mmaux ), &
                         ll,mx(1:3),nrvars,maux,mbc,h,xlower(1:3))
    case default
      print *, 'UpdateAuxInVcycle: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end function UpdateAuxInVcycle
  ! ---------------------------------------------------------------------------
  integer function UpdateSourceFunction(info,dummy)
    !
    ! Dimensionally invariant source update routine.  If no update is needed set
    ! SourceUpdate# = 0.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    use Problem, only: SourceUpdate2D, SourceUpdate3D
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer                       :: ll, maux, mbc, mmaux, nrvars
    integer, dimension(1:maxdims) :: amx, aul, mx, ul
    real                          :: h
    real, dimension(1:maxdims)    :: xlower
    !
    UpdateSourceFunction = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx     = 1
    h      = info%dx(1)
    mx(1:ndims)     = info%mx(1:ndims)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amx   = 1
    if (maux>0) amx(1:ndims) = mx(1:ndims)
    !
    mbc = info%mbc
    ll  = 1-mbc
    ul  = 1
    aul = 1
    ul(1:ndims)  = mx(1:ndims)+mbc
    aul(1:ndims) = amx(1:ndims)+mbc
    !
    select case(ndims)
    case(2)
     !info%f(1:mx(1),1:mx(2),1,1:nrvars)             &
     !       = info%ftmp(1:mx(1),1:mx(2),1,1:nrvars) &
      info%f(1:mx(1),1:mx(2),1,1:nrvars)             &
               = SourceUpdate2D(info%q(   ll: ul(1),ll: ul(2),1,1:nrvars), &
                                info%qold(ll: ul(1),ll: ul(2),1,1:nrvars), &
                                info%aux( ll:aul(1),ll:aul(2),1,1:mmaux ), &
                                ll,mx(1:2),nrvars,maux,mbc,h,xlower(1:2))
    case(3)
      info%f(1:mx(1),1:mx(2),1:mx(3),1:nrvars)             &
             = info%ftmp(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
               + SourceUpdate3D(                           &
                     info%q(   ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                     info%qold(ll: ul(1),ll: ul(2),ll: ul(3),1:nrvars), &
                     info%aux( ll:aul(1),ll:aul(2),ll:aul(3),1:mmaux ), &
                     ll,mx(1:3),nrvars,maux,mbc,h,xlower(1:3))
    case default
      print *, 'UpdateSourceFunction: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end function UpdateSourceFunction
  ! ---------------------------------------------------------------------------
  recursive subroutine AFASVcycle(level)
    !
    ! Dimensionally invariant recursive AFAS vcycle routine.
    !
    use NodeInfoDef
    use TreeOps,  only: ApplyOnLevel
    use Boundary, only: SetGhost, GetCoarseGhostPoints
    implicit none
    integer, intent(in) :: level
    type(funcparam)     :: dummy
    !
    call LevelRelax(level)
    !
    if (level>minlevel) then
      !
      call ApplyOnLevel(level,RestrictSolution,dummy)
      !
      call SetGhost(level-1,0)
      !
      call ApplyOnLevel(level,GetCoarseGhostPoints,dummy)
      !
      call ApplyOnLevel(level,CoarseLoadingFunction,dummy)
      !
      call AFASVcycle(level-1)
      !
      call ApplyOnLevel(level,CorrectFine,dummy)
      !
      call SetGhost(level,0)
      !
      call LevelRelax(level)
      !
    end if
    !
  end subroutine AFASVcycle
  ! ---------------------------------------------------------------------------
  subroutine LevelRelax(level)
    !
    ! Dimensionally invariant relaxation.
    !
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel
    use Boundary, only: SetGhost
    implicit none
    integer, intent(in) :: level
    type(funcparam)     :: dummy
    integer             :: redblack, smoothingpass
    !
    ! Assume Ghost points are set on entry to this routine.
    do smoothingpass = 1, nsmoothingpasses
      do redblack = 1, 2
        !
        dummy%iswitch = redblack
        call ApplyOnLevel(level,RelaxPatch,dummy)
        !
        call SetGhost(level,redblack)
        !
      end do
    end do
    !
  end subroutine LevelRelax
  ! ---------------------------------------------------------------------------
  integer function RelaxPatch(info,dummy)
    !
    ! Dimensionally invariant relaxation of a patch, on RED squares for
    ! dummy%iswitch = 1, and BLACK squares for dummy%iswitch = 2.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    use Problem, only: RelaxGrid2d, RelaxGrid3D
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer                       :: maux, mmaux, nrvars
    integer, dimension(1:maxdims) :: amx, mx
    real                          :: h
    real, dimension(1:maxdims)    :: xlower
    !
    RelaxPatch = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx     = 1
    h      = info%dx(1)
    mx(1:ndims)     = info%mx(1:ndims)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amx   = 1
    if (maux>0) amx(1:ndims) = mx(1:ndims)
    !
    select case(ndims)
    case(2)
      call RelaxGrid2D(info%q(   0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                       info%qold(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                       info%aux( 0:amx(1)+1,0:amx(2)+1,1,1:mmaux ), &
                       info%f(   1: mx(1),  1: mx(2)  ,1,1:nrvars), &
                       mx(1:2),nrvars,maux,h,xlower(1:2),dummy%iswitch)
    case(3)
      call RelaxGrid3D(info%q(   0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                       info%qold(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                       info%aux( 0:amx(1)+1,0:amx(2)+1,0:amx(3)+1,1:mmaux ), &
                       info%f(   1: mx(1),  1: mx(2)  ,1: mx(3)  ,1:nrvars), &
                       mx(1:3),nrvars,maux,h,xlower(1:3),dummy%iswitch)
    case default
      print *, 'RelaxPatch: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end function RelaxPatch
  ! ---------------------------------------------------------------------------
  integer function RestrictSolution(info,dummy)
    !
    ! Dimensionally invariant restriction.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok, GetParentInfo
    use GridUtilities, only: Restriction2D, Restriction3D
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    type(nodeinfo), pointer           :: parent
    integer                           :: ierror, nrvars
    integer, dimension(1:maxdims)     :: cmx, mx
    integer, dimension(1:maxdims,1:2) :: mb
    !
    RestrictSolution = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx     = 1
    cmx    = 1
    mb     = 1
    mx(1:ndims)     = info%mx(1:ndims)
    cmx(1:ndims)    = mx(1:ndims)/2
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    !
    ierror = GetParentInfo(parent)
    !
    select case(ndims)
    case(2)
      info%qc(1:cmx(1),1:cmx(2),1,1:nrvars) &
              = Restriction2D(info%q(1:mx(1),1:mx(2),1,1:nrvars))
    case(3)
      info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
              = Restriction3D(info%q(1:mx(1),1:mx(2),1:mx(3),1:nrvars))
    case default
      print *, 'RestrictSolution: Only ndims = 2,3 are supported.'
      stop
    end select
    !
    parent%q(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars) &
             = info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars)
    !
  end function RestrictSolution
  ! ---------------------------------------------------------------------------
  integer function CoarseLoadingFunction(info,dummy)
    !
    ! Dimensionally invariant coarse loading function routine.
    !
    use NodeInfoDef
    use TreeOps, only: err_ok, GetParentInfo
    use GridUtilities, only: Restriction2D, Restriction3D
    use Problem, only: Operator2D, Operator3D
    implicit none
    type(nodeinfo)                     :: info
    type(funcparam)                    :: dummy
    type(nodeinfo), pointer            :: parent
    integer                            :: ierror, maux, mmaux, nrvars
    integer, dimension(1:maxdims)      :: cmx, mx
    integer, dimension(1:maxdims,1:2)  :: amb, mb
    real                               :: ch
    real, dimension(1:maxdims)         :: xlower
    !
    CoarseLoadingFunction = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx     = 1
    cmx    = 1
    mb     = 1
    mx(1:ndims)     = info%mx(1:ndims)
    cmx(1:ndims)    = mx(1:ndims)/2
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    !
    ierror          = GetParentInfo(parent)
    ch              = parent%dx(1)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amb   = 1
    if (maux>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
    !
    call Residual(info)
    !
    select case(ndims)
    case(2)
      parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nrvars)               &
               = Restriction2D(info%rf(1:mx(1),1:mx(2),1,1:nrvars))      &
                 + Operator2D(info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars), &
                              parent%qold(mb(1,1)-1: mb(1,2)+1,          &
                              mb(2,1)-1: mb(2,2)+1,1,1:nrvars),          &
                              parent%aux(amb(1,1)-1:amb(1,2)+1,          &
                              amb(2,1)-1:amb(2,2)+1,1,1:mmaux ),         &
                              cmx(1:2),nrvars,maux,ch,xlower(1:2))
    case(3)
      parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nrvars)     &
               = Restriction3D(info%rf(1:mx(1),1:mx(2),1:mx(3),1:nrvars))    &
                 + Operator3D(                                               &
                         info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars), &
                         parent%qold(mb(1,1)-1: mb(1,2)+1,                   &
                         mb(2,1)-1: mb(2,2)+1,                               &
                         mb(3,1)-1: mb(3,2)+1,1:nrvars),                     &
                         parent%aux(amb(1,1)-1:amb(1,2)+1,                   &
                         amb(2,1)-1:amb(2,2)+1,                              &
                         amb(3,1)-1:amb(3,2)+1,1:mmaux ),                    &
                         cmx(1:3),nrvars,maux,ch,xlower(1:3))
    case default
      print *, 'CoarseLoadingFunction: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end function CoarseLoadingFunction
  ! ---------------------------------------------------------------------------
  integer function RelativeTruncationError(info,dummy)
    !
    ! Dimensionally invariant relative trunctation error:
    !
    !   L_{2h}(I_h^{2h} q_h)-I_h^{2h}(L_h q_h).
    !
    use NodeInfoDef
    use TreeOps,       only: err_ok, GetParentInfo
    use GridUtilities, only: Restriction2D, Restriction3D
    use Problem,       only: Operator2D, Operator3D
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    type(nodeinfo), pointer           :: parent
    integer                           :: ierror, maux, mmaux, nrvars
    integer, dimension(1:maxdims)     :: amx, cmx, mx
    integer, dimension(1:maxdims,1:2) :: amb, mb
    real                              :: ch, h
    real, dimension(1:maxdims)        :: xlower
    !
    RelativeTruncationError = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx     = 1
    cmx    = 1
    mb     = 1
    h      = info%dx(1)
    mx(1:ndims)     = info%mx(1:ndims)
    cmx(1:ndims)    = mx(1:ndims)/2
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amx   = 1
    if (maux>0) amx(1:ndims) = mx(1:ndims)
    !
    ierror = GetParentInfo(parent)
    ch     = parent%dx(1)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amb   = 1
    if (maux>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
    !
    select case(ndims)
    case(2)
      info%qrte(1:cmx(1),1:cmx(2),1,1:nrvars)                                &
                = Restriction2D(                                             &
                     Operator2D(info%q(   0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                                info%qold(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                                info%aux( 0:amx(1)+1,0:amx(2)+1,1,1:mmaux ), &
                                mx(1:2),nrvars,maux,h,xlower(1:2)))          &
                - Operator2D(info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars),      &
                             parent%qold(mb(1,1)-1: mb(1,2)+1,               &
                             mb(2,1)-1: mb(2,2)+1,1,1:nrvars),               &
                             parent%aux(amb(1,1)-1:amb(1,2)+1,               &
                             amb(2,1)-1:amb(2,2)+1,1,1:mmaux ),              &
                             cmx(1:2),nrvars,maux,ch,xlower(1:2))
    case(3)
      info%qrte(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars)                       &
                = Restriction3D(                                           &
          Operator3D(info%q(   0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                     info%qold(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                     info%aux( 0:amx(1)+1,0:amx(2)+1,0:amx(3)+1,1:mmaux ), &
                     mx(1:3),nrvars,maux,h,xlower(1:3)))                   &
      - Operator3D(info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars),     &
                   parent%qold(mb(1,1)-1: mb(1,2)+1,                       &
                   mb(2,1)-1: mb(2,2)+1,                                   &
                   mb(3,1)-1: mb(3,2)+1,1:nrvars),                         &
                   parent%aux(amb(1,1)-1:amb(1,2)+1,                       &
                   amb(2,1)-1:amb(2,2)+1,                                  &
                   amb(3,1)-1:amb(3,2)+1,1:mmaux ),                        &
                   cmx(1:3),nrvars,maux,ch,xlower(1:3))
    case default
      print *, 'CoarseLoadingFunction: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end function RelativeTruncationError
  ! ---------------------------------------------------------------------------
  subroutine Residual(info)
    use NodeInfoDef
    use Problem, only: Operator2D, Operator3D
    implicit none
    type(nodeinfo)                :: info
    integer                       :: maux, mmaux, nrvars
    integer, dimension(1:maxdims) :: amx, mx
    real                          :: h
    real, dimension(1:maxdims)    :: xlower
    !
    nrvars = info%nrvars
    mx     = 1
    h      = info%dx(1)
    mx(1:ndims)     = info%mx(1:ndims)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    maux  = info%maux
    mmaux = max(maux,1)
    amx   = 1
    if (maux>0) amx(1:ndims) = mx(1:ndims)
    !
    select case(ndims)
    case(2)
      info%rf(1:mx(1),1:mx(2),1,1:nrvars) &
              = info%f(1:mx(1),1:mx(2),1,1:nrvars) &
                - Operator2D(info%q(   0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                             info%qold(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                             info%aux( 0:amx(1)+1,0:amx(2)+1,1,1:mmaux ), &
                             mx(1:2),nrvars,maux,h,xlower(1:2))
    case(3)
      info%rf(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
         = info%f(1:mx(1),1:mx(2),1:mx(3),1:nrvars) &
           - Operator3D(info%q(   0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                        info%qold(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                        info%aux( 0:amx(1)+1,0:amx(2)+1,0:amx(3)+1,1:mmaux ), &
                        mx(1:3),nrvars,maux,h,xlower(1:3))
    case default
      print *, 'Residual: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end subroutine Residual
  ! ---------------------------------------------------------------------------
  integer function CorrectFine(info,dummy)
    !
    ! Dimensionally invariant coarse-grid-correction routine.
    !
    use NodeInfoDef
    use TreeOps,       only: err_ok, GetParentInfo
    use GridUtilities, only: BiLinProlongationP1,  BiLinProlongationP2, &
                             TriLinProlongationP1, TriLinProlongationP2
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    type(nodeinfo), pointer           :: parent
    integer                           :: ierror, mbc, nrvars
    integer, dimension(1:maxdims)     :: cmx, mx
    integer, dimension(1:maxdims,1:2) :: mb
    !
    CorrectFine = err_ok
    !
    if (info%tobedeleted) return
    !
    ierror = GetParentInfo(parent)
    !
    nrvars = info%nrvars
    mbc    = info%mbc
    mx     = 1
    cmx    = 1
    mb     = 1
    mx(1:ndims)     = info%mx(1:ndims)
    cmx(1:ndims)    = mx(1:ndims)/2
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    !
    select case(ndims)
    case(2)
      info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nrvars) &
            = parent%q(mb(1,1)-mbc:mb(1,2)+mbc,             &
                       mb(2,1)-mbc:mb(2,2)+mbc,1,1:nrvars)  &
              -  info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nrvars)
      !
      select case(mbc)
      case(1)
        !
        ! Bilinear prolongation:
        info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars)          &
               = info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars) &
               + BiLinProlongationP1(info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars))
        !
        ! Mass conserving bilinear prolongation:
        !          info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nrvars) &
        !        = info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nrvars) &
        !        + BiLinProlongationP1MC( &
        !          info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nrvars))
        !
        ! Simple injection:
        !          info%q( 1: mx(1),1: mx(2),1,1:nrvars) &
        !        = info%q( 1: mx(1),1: mx(2),1,1:nrvars) &
        !        + Prolongation2D( &
        !          info%qc(1:cmx(1),1:cmx(2),1,1:nrvars))
        !
      case(2)
        !
        ! Bilinear prolongation:
        info%q(-1:mx(1)+2,-1:mx(2)+2,1,1:nrvars)        &
             = info%q(-1:mx(1)+2,-1:mx(2)+2,1,1:nrvars) &
             + BiLinProlongationP2(info%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nrvars))
        !
        ! Mass-conserving bilinear prolongation:
        !          info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
        !        = info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
        !        + BiLinProlongationP2MC( &
        !          info%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nrvars))
      end select
      !
    case(3)
      info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nrvars) &
            = parent%q(mb(1,1)-mbc:mb(1,2)+mbc,                            &
                       mb(2,1)-mbc:mb(2,2)+mbc,                            &
                       mb(3,1)-mbc:mb(3,2)+mbc,1:nrvars)                   &
            - info%qc(1-mbc:cmx(1)+mbc,                                    &
                      1-mbc:cmx(2)+mbc,                                    &
                      1-mbc:cmx(3)+mbc,1:nrvars)
      !
      select case(mbc)
      case(1)
        !
        ! Trilinear prolongation:
        info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars)          &
               = info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars) &
                + TriLinProlongationP1(                         &
                        info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars))
        !
        ! Mass conserving trilinear prolongation:
        !        info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars) &
        !      = info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars) &
        !      + TriLinProlongationP1MC( &
        !        info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars))
        !
        ! Simple injection:
        !        info%q( 1: mx(1),1: mx(2),1: mx(3),1:nrvars) &
        !      = info%q( 1: mx(1),1: mx(2),1: mx(3),1:nrvars) &
        !      + Prolongation3D( &
        !        info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars))
      case(2)
        !
        ! Trilinear prolongation:
        info%q(-1:mx(1)+2,-1:mx(2)+2,-1:mx(3)+2,1:nrvars) &
               = info%q(-1:mx(1)+2,-1:mx(2)+2,-1:mx(3)+2,1:nrvars) &
                + TriLinProlongationP2( &
                        info%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nrvars))
        !
        ! Mass-conserving trilinear prolongation:
        !        info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
        !      = info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
        !      + TriLinProlongationP2( &
        !        info%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nrvars))
      end select
      !
    case default
      print *, 'CorrectFine: Only ndims=2,3 are supported.'
      stop
    end select
    !
  end function CorrectFine
  ! ---------------------------------------------------------------------------
end module AFASRoutines
