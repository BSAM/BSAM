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
! Portions of the code
!
! (c) Copyright Sorin Mitran, 2002
! Department of Mathematics
! University of North Carolina at Chapel Hill
! mitran@amath.unc.edu
!
! -----------------------------------------------------------------------
! This software is made available for research and instructional use only.
! You may copy and use this software without charge for these
! non-commercial purposes, provided that the copyright notice and
! associated text is reproduced on all copies.  For all other uses,
! including distribution of modified versions, please contact the authors.
!
! Commercial use is strictly forbidden without permission.
!
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             bsamstorage.f90
! Purpose:          BSAM memory allocation module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
module BSAMStorage
  !
contains
  ! ---------------------------------------------------------------------------
  subroutine AllocFields(info,parent)
    use NodeInfoDef
    implicit none
    type(nodeinfo)                       :: info
    type(nodeinfo), intent(in), optional :: parent
    integer                              :: ierror,maux,mmaux,mbc,n,nrvars
    integer, dimension(1:maxdims)        :: amx, cmx, mx
    !
    mx = info%mx
    mbc = info%mbc
    nrvars = info%nrvars
    maux = info%maux
    cmx = 1
    do n = 1, ndims
      if (modulo(mx(n),2)==1) then
        cmx(n) = 1
      else
        cmx(n) = mx(n)/2
      end if
    end do
    !
    mmaux = max(maux,1)
    amx = 1
    if (maux>0) amx(1:ndims) = mx(1:ndims)
    !
    ! Allocate space for all pointer components from nodeinfodef:
    select case(ndims)
    case(2)
      allocate( &
 info%errorflags(1    : mx(1)    ,1    : mx(2)    ,1:1         ), &
 info%q(         1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1:1,1:nrvars), &
 info%qold(      1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1:1,1:nrvars), &
 info%qc(        1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1:1,1:nrvars), &
 info%qrte(      1    :cmx(1)    ,1    :cmx(2)    ,1:1,1:nrvars), &
 info%aux(       1-mbc:amx(1)+mbc,1-mbc:amx(2)+mbc,1:1,1:mmaux ), &
 info%f(         1    : mx(1)    ,1    : mx(2)    ,1:1,1:nrvars), &
 info%rf(        1    : mx(1)    ,1    : mx(2)    ,1:1,1:nrvars), &
 info%ftmp(      1    : mx(1)    ,1    : mx(2)    ,1:1,1:nrvars), &
 STAT=ierror)
      !
      if (ierror/=0) then
        print *,'Error in allocation of nodeinfodef components in AllocFields'
        stop
      end if
      !
    case(3)
      allocate( &
 info%errorflags(1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)             ), &
 info%q(         1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1-mbc: mx(3)+mbc,1:nrvars), &
 info%qold(      1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1-mbc: mx(3)+mbc,1:nrvars), &
 info%qc(        1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nrvars), &
 info%qrte(      1    :cmx(1)    ,1    :cmx(2)    ,1    :cmx(3)    ,1:nrvars), &
 info%aux(       1-mbc:amx(1)+mbc,1-mbc:amx(2)+mbc,1-mbc:amx(3)+mbc,1:mmaux ), &
 info%f(         1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1:nrvars), &
 info%rf(        1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1:nrvars), &
 info%ftmp(      1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1:nrvars), &
 STAT=ierror)
      !
      if (ierror/=0) then
        print *,'Error in allocation of nodeinfodef components in AllocFields'
        stop
      end if
      !
    case default
      print *, 'AllocFields: Only ndims = 2,3 are supported.'
      stop
    end select
    !
    ! Initialize solution arrays:
    info%q = 0
    info%qold = 0
    !
    ! Initialize error flags:
    info%errorflags = 0
    !
    info%fieldsallocated = .true.
    !
  end subroutine AllocFields
  ! ---------------------------------------------------------------------------
  subroutine DeAllocFields(info)
    use NodeInfoDef
    implicit none
    type(nodeinfo) :: info
    !
    if (.not. info%fieldsallocated) return
    !
    if (associated(info%errorflags)) deallocate(info%errorflags)
    if (associated(info%q)) deallocate(info%q)
    if (associated(info%qold)) deallocate(info%qold)
    if (associated(info%qc)) deallocate(info%qc)
    if (associated(info%qrte)) deallocate(info%qrte)
    if (associated(info%aux)) deallocate(info%aux)
    if (associated(info%f)) deallocate(info%f)
    if (associated(info%rf)) deallocate(info%rf)
    if (associated(info%ftmp)) deallocate(info%ftmp)
    !
    ! Nullify the pointers:
    nullify(info%errorflags,info%q,info%qold,info%qc,info%qrte, &
            info%aux,info%f,info%rf,info%ftmp)
    !
    info%fieldsallocated = .false.  ! Don't compute on this node:
    info%tobedeleted = .true.       ! Mark for garbage collection:
    !
  end subroutine DeAllocFields
  ! ---------------------------------------------------------------------------
  subroutine AllocPeriodicBCStorage(np,nrvars)
    use NodeInfoDef
    implicit none
    integer, intent(in) :: np
    integer, intent(in) :: nrvars
    integer             :: ierror
    !
    allocate(periodicoffsetindex(1:np),poffset(1:maxdims,1:np), &
             qoffset(1:np,1:nrvars),STAT=ierror)
    !
    if (ierror/=0) then
      print *,'Error in AllocPeriodicBCStorage.'
      stop
    end if
    !
  end subroutine AllocPeriodicBCStorage
  ! ---------------------------------------------------------------------------
  subroutine DeallocPeriodicBCStorage
    use NodeInfoDef
    implicit none
    !
    if (allocated(periodicoffsetindex)) deallocate(periodicoffsetindex)
    if (allocated(poffset)) deallocate(poffset)
    if (allocated(qoffset)) deallocate(qoffset)
    !
  end subroutine DeallocPeriodicBCStorage
  ! ---------------------------------------------------------------------------
  subroutine AllocUniformGrids(mxroot,mbc,nrvars)
    use NodeInfoDef
    implicit none
    integer, dimension(1:maxdims), intent(in) :: mxroot ! root-level grid size.
    integer,                       intent(in) :: mbc
    integer,                       intent(in) :: nrvars
    integer                                   :: ierror, level
    integer, dimension(1:maxdims)             :: high, low, mx
    !
    mx = 1
    low = 1
    high = 1
    mx(1:ndims) = mxroot(1:ndims)
    low(1:ndims) = 1-mbc
    !
    do level = 0, maxlevel
      !
      if (level>0) mx(1:ndims) = mx(1:ndims)*2
      uniformgrid(level)%mx = 1
      uniformgrid(level)%mx(1:ndims) = mx(1:ndims)
      !
      high(1:ndims) = mx(1:ndims)+mbc
      !
      allocate( &
 uniformgrid(level)%q(low(1):high(1),low(2):high(2),low(3):high(3),1:nrvars), &
 STAT=ierror)
      !
      if (ierror/=0) then
        print *, 'AllocUniformGrids: Error allocating uniformgrid on level', &
                 level
        stop
      end if
    end do
    !
  end subroutine AllocUniformGrids
  ! ---------------------------------------------------------------------------
  subroutine DeallocUniformGrids
    use NodeInfoDef
    implicit none
    integer :: level
    !
    do level = rootlevel, maxlevel
      uniformgrid(level)%mx = 0
      if (associated(uniformgrid(level)%q)) deallocate(uniformgrid(level)%q)
    end do
    !
  end subroutine DeallocUniformGrids
  ! ---------------------------------------------------------------------------
end module BSAMStorage
