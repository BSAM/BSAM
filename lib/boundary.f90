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
! File:             boundary.f90
! Purpose:          Ghost-cell interpolation module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
module Boundary
  implicit none
  !
  real, parameter :: a1 = - 1.0/ 5.0
  real, parameter :: a2 =   2.0/ 3.0
  real, parameter :: a3 =   8.0/15.0
  real, parameter :: b1 =   5.0/32.0
  real, parameter :: b2 =  15.0/16.0
  real, parameter :: b3 = - 3.0/32.0
  real, parameter :: c1 =   9.0/16.0
  real, parameter :: c2 =   3.0/16.0
  real, parameter :: c3 =   3.0/16.0
  real, parameter :: c4 =   1.0/16.0
  !
contains
  ! ---------------------------------------------------------------------------
  subroutine SetGhost(level,ipass)
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel, ApplyOnLevelPairs
    implicit none
    integer, intent(in) :: level
    integer, intent(in) :: ipass
    type(funcparam) :: dummy
    !
    dummy%iswitch = ipass
    !
    if (level>rootlevel) then
      ! 0. Interpolate between coarse and fine grids:
      call ApplyOnLevel(level,InterpolateCoarseFine,dummy)
      !
      ! 1. Transfer boundary conditions among grids on same level:
      call ApplyOnLevelPairs(level,TransferBC,dummy)
    endif
    !
    if (periodicboundaryconditions) then
      !
      ! 2. Apply periodic boundary conditions on same grid:
      call ApplyOnLevel(level,PeriodicBC,dummy)
      !
      ! 3. Apply periodic boundary conditions on different grids:
      if (level>rootlevel) &
           call ApplyOnLevelPairs(level,TransferPeriodicBC,dummy)
    end if
    !
    ! 4. Apply physical boundary conditions:
    call ApplyOnLevel(level,SetBC,dummy)
    !
  end subroutine SetGhost
  ! ---------------------------------------------------------------------------
  integer function InterpolateCoarseFine(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok, GetParentInfo
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    type(nodeinfo), pointer           :: parent
    integer                           :: ierror, ipass, nrvars
    integer, dimension(1:maxdims)     :: mx, pmx
    integer, dimension(1:2*maxdims)   :: mthbc
    integer, dimension(1:maxdims,1:2) :: mb
    !
    InterpolateCoarseFine = err_ok
    !
    if (info%tobedeleted) return
    !
    ipass = dummy%iswitch
    ierror = GetParentInfo(parent)
    nrvars = info%nrvars
    mx = 1
    mx(1:ndims) = info%mx(1:ndims)
    pmx = 1
    pmx(1:ndims) = parent%mx(1:ndims)
    mb = 1
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    mthbc = 1
    mthbc(1:2*ndims) = info%mthbc(1:2*ndims)
    !
    select case(ndims)
    case(2)
      call Interpolate2D(  info%q(0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                         parent%q(0:pmx(1)+1,0:pmx(2)+1,1,1:nrvars), &
                         mx(1:2),nrvars,mb(1:2,1:2),mthbc(1:4),ipass)
    case(3)
      call Interpolate3D(  info%q(0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nrvars), &
                         parent%q(0:pmx(1)+1,0:pmx(2)+1,0:pmx(3)+1,1:nrvars), &
                         mx(1:3),nrvars,mb(1:3,1:2),mthbc(1:6),ipass)
    case default
      print *, 'InterpolateCoarseFine: only n=2,3 supported.'
      stop
    end select
    !
  end function InterpolateCoarseFine
  ! ---------------------------------------------------------------------------
  subroutine Interpolate2D(qf,qc,mx,nrvars,mb,mthbc,ipass)
    implicit none
    real, dimension(0:,0:,1:), intent(in out) :: qf
    real, dimension(0:,0:,1:), intent(in)     :: qc
    integer, dimension(1:2), intent(in)       :: mx
    integer, intent(in)                       :: nrvars
    integer, dimension(1:2,1:2), intent(in)   :: mb
    integer, dimension(1:4), intent(in)       :: mthbc
    integer, intent(in)                       :: ipass
    !
    ! Corners:
    !
    ! Red corners:
    if (ipass==1 .or. ipass==0) then
      qf(mx(1)+1,mx(2)+1,1:nrvars) &
         = a1*qf(mx(1)-1,mx(2)-1,1:nrvars) &
         + a2*qf(mx(1)  ,mx(2)  ,1:nrvars) &
         + a3*qc(mb(1,2)+1,mb(2,2)+1,1:nrvars)
      !
      qf(0,0,1:nrvars) &
         = a1*qf(2,2,1:nrvars) &
         + a2*qf(1,1,1:nrvars) &
         + a3*qc(mb(1,1)-1,mb(2,1)-1,1:nrvars)
    end if
    !
    ! Black corners:
    if (ipass==2 .or. ipass==0) then
      qf(mx(1)+1,0,1:nrvars) &
         = a1*qf(mx(1)-1,2,1:nrvars) &
         + a2*qf(mx(1)  ,1,1:nrvars) &
         + a3*qc(mb(1,2)+1,mb(2,1)-1,1:nrvars)
      !
      qf(0,mx(2)+1,1:nrvars) &
         = a1*qf(2,mx(2)-1,1:nrvars) &
         + a2*qf(1,mx(2)  ,1:nrvars) &
         + a3*qc(mb(1,1)-1,mb(2,2)+1,1:nrvars)
    end if
    !
    ! Faces:
    !
    select case(mthbc(1))
    case(2,999) ! Periodic or Internal boundary.
      !
      ! Face #1, Red:
      if (ipass==1 .or. ipass==0) then
        qf(0,2:mx(2)-2:2,1:nrvars) &
           = a1*qf(2,4:mx(2)  :2,1:nrvars) &
             + a2*qf(1,3:mx(2)-1:2,1:nrvars) &
             + a3*qc(mb(1,1)-1,mb(2,1):mb(2,2)-1,1:nrvars)
        !
        qf(0,mx(2),1:nrvars) &
           = a1*qf(2,mx(2),1:nrvars) &
             + a2*qf(1,mx(2),1:nrvars) &
             + a3*(b1*qc(mb(1,1)-1,mb(2,2)+1,1:nrvars) &
             +     b2*qc(mb(1,1)-1,mb(2,2)  ,1:nrvars) &
             +     b3*qc(mb(1,1)-1,mb(2,2)-1,1:nrvars))
      end if
      !
      ! Face #1, Black:
      if (ipass==2 .or. ipass==0) then
        qf(0,3:mx(2)-1:2,1:nrvars) &
           = a1*qf(2,1:mx(2)-3:2,1:nrvars) &
             + a2*qf(1,2:mx(2)-2:2,1:nrvars) &
             + a3*qc(mb(1,1)-1,mb(2,1)+1:mb(2,2),1:nrvars)
        !
        qf(0,1,1:nrvars) &
           =          a1*qf(2,1,1:nrvars) &
             + a2*qf(1,1,1:nrvars) &
             + a3*(b1*qc(mb(1,1)-1,mb(2,1)-1,1:nrvars) &
             +     b2*qc(mb(1,1)-1,mb(2,1)  ,1:nrvars) &
             +     b3*qc(mb(1,1)-1,mb(2,1)+1,1:nrvars))
      end if
      !
    case default
      continue ! Physical boundary.
    end select
    !
    select case(mthbc(2))
    case(2,999) ! Periodic or Internal boundary.
      !
      ! Face #2, Red:
      if (ipass==1 .or. ipass==0) then
        qf(mx(1)+1,3:mx(2)-1:2,1:nrvars) &
           = a1*qf(mx(1)-1,1:mx(2)-3:2,1:nrvars) &
             + a2*qf(mx(1)  ,2:mx(2)-2:2,1:nrvars) &
             + a3*qc(mb(1,2)+1,mb(2,1)+1:mb(2,2):1,1:nrvars)
        !
        qf(mx(1)+1,1,1:nrvars) &
           = a1*qf(mx(1)-1,1,1:nrvars) &
             + a2*qf(mx(1)  ,1,1:nrvars) &
             + a3*(b1*qc(mb(1,2)+1,mb(2,1)-1,1:nrvars) &
             +     b2*qc(mb(1,2)+1,mb(2,1)  ,1:nrvars) &
             +     b3*qc(mb(1,2)+1,mb(2,1)+1,1:nrvars))
      end if
      !
      ! Face #2, Black:
      if (ipass==2 .or. ipass==0) then
        qf(mx(1)+1,2:mx(2)-2:2,1:nrvars) &
           = a1*qf(mx(1)-1,4:mx(2)  :2,1:nrvars) &
             + a2*qf(mx(1)  ,3:mx(2)-1:2,1:nrvars) &
             + a3*qc(mb(1,2)+1,mb(2,1):mb(2,2)-1:1,1:nrvars)
        !
        qf(mx(1)+1,mx(2),1:nrvars) &
           = a1*qf(mx(1)-1,mx(2),1:nrvars) &
             + a2*qf(mx(1)  ,mx(2),1:nrvars) &
             + a3*(b1*qc(mb(1,2)+1,mb(2,2)+1,1:nrvars) &
             +     b2*qc(mb(1,2)+1,mb(2,2)  ,1:nrvars) &
             +     b3*qc(mb(1,2)+1,mb(2,2)-1,1:nrvars))
      end if
      !
    case default
      continue ! Physical boundary.
    end select
    !
    select case(mthbc(3))
    case(2,999) ! Periodic or Internal boundary.
      !
      ! Face #3, Red:
      if (ipass==1 .or. ipass==0) then
        qf(2:mx(1)-2:2,0,1:nrvars) &
           = a1*qf(4:mx(1)  :2,2,1:nrvars) &
             + a2*qf(3:mx(1)-1:2,1,1:nrvars) &
             + a3*qc(mb(1,1):mb(1,2)-1,mb(2,1)-1,1:nrvars)
        !
        qf(mx(1),0,1:nrvars) &
           = a1*qf(mx(1),2,1:nrvars) &
             + a2*qf(mx(1),1,1:nrvars) &
             + a3*(b1*qc(mb(1,2)+1,mb(2,1)-1,1:nrvars) &
             +     b2*qc(mb(1,2)  ,mb(2,1)-1,1:nrvars) &
             +     b3*qc(mb(1,2)-1,mb(2,1)-1,1:nrvars))
      end if
      !
      ! Face #3, Black:
      if (ipass==2 .or. ipass==0) then
        qf(3:mx(1)-1:2,0,1:nrvars) &
           = a1*qf(1:mx(1)-3:2,2,1:nrvars) &
             + a2*qf(2:mx(1)-2:2,1,1:nrvars) &
             + a3*qc(mb(1,1)+1:mb(1,2),mb(2,1)-1,1:nrvars)
        !
        qf(1,0,1:nrvars) &
           = a1*qf(1,2,1:nrvars) &
             + a2*qf(1,1,1:nrvars) &
             + a3*(b1*qc(mb(1,1)-1,mb(2,1)-1,1:nrvars) &
             +     b2*qc(mb(1,1)  ,mb(2,1)-1,1:nrvars) &
             +     b3*qc(mb(1,1)+1,mb(2,1)-1,1:nrvars))
      end if
      !
    case default
      continue
    end select
    !
    select case(mthbc(4))
    case(2,999) ! Periodic or Internal boundary.
      !
      ! Face #4, Red:
      if (ipass==1 .or. ipass==0) then
        qf(3:mx(1)-1:2,mx(2)+1,1:nrvars) &
           = a1*qf(1:mx(1)-3:2,mx(2)-1,1:nrvars) &
             + a2*qf(2:mx(1)-2:2,mx(2)  ,1:nrvars) &
             + a3*qc(mb(1,1)+1:mb(1,2),mb(2,2)+1,1:nrvars)
        !
        qf(1,mx(2)+1,1:nrvars) &
           = a1*qf(1,mx(2)-1,1:nrvars) &
             + a2*qf(1,mx(2)  ,1:nrvars) &
             + a3*(b1*qc(mb(1,1)-1,mb(2,2)+1,1:nrvars) &
             +     b2*qc(mb(1,1)  ,mb(2,2)+1,1:nrvars) &
             +     b3*qc(mb(1,1)+1,mb(2,2)+1,1:nrvars))
      end if
      !
      ! Face #4, Black:
      if (ipass==2 .or. ipass==0) then
        qf(2:mx(1)-2:2,mx(2)+1,1:nrvars) &
           = a1*qf(4:mx(1)  :2,mx(2)-1,1:nrvars) &
             + a2*qf(3:mx(1)-1:2,mx(2)  ,1:nrvars) &
             + a3*qc(mb(1,1):mb(1,2)-1,mb(2,2)+1,1:nrvars)
        !
        qf(mx(1),mx(2)+1,1:nrvars) &
           =          a1*qf(mx(1),mx(2)-1,1:nrvars) &
             + a2*qf(mx(1),mx(2)  ,1:nrvars) &
             + a3*(b1*qc(mb(1,2)+1,mb(2,2)+1,1:nrvars) &
             +     b2*qc(mb(1,2)  ,mb(2,2)+1,1:nrvars) &
             +     b3*qc(mb(1,2)-1,mb(2,2)+1,1:nrvars))
      end if
      !
    case default
      continue ! Physical boundary.
    end select
    !
  end subroutine Interpolate2D
  ! ---------------------------------------------------------------------------
  integer function TransferBC(grid1,grid2,dummy)
    !
    ! Transfer boundary values among grids on same level
    !
    use NodeinfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: grid1, grid2
    type(funcparam) :: dummy
    !
    TransferBC=err_ok
    !
    ! Check for inactive grids awaiting garbage collection
    if (grid1%tobedeleted .or. grid2%tobedeleted) return
    !
    ! Look for overlap of ghost cell regions of grid1 and grid2.  Order doesn't
    ! matter, as transfer goes both ways:
    call GhostOverlap(grid1,grid2)
    !
  end function TransferBC
  ! ---------------------------------------------------------------------------
  subroutine GhostOverlap(grid1,grid2)
    !
    ! Assumes ghost layer of size mbc:
    !
    use NodeinfoDef
    implicit none
    type(nodeinfo)                    :: grid1
    type(nodeinfo)                    :: grid2
    integer                           :: mbc, n, nrvars
    integer, dimension(1:maxdims,1:2) :: mg1, mg2, ml1, ml2, ovg, ovg1, ovg2
    !
    mbc = grid2%mbc
    nrvars = grid1%nrvars
    !
    mg1 = 1
    mg2 = 1
    mg1(1:ndims,1:2) = grid1%mglobal(1:ndims,1:2)
    mg2(1:ndims,1:2) = grid2%mglobal(1:ndims,1:2)
    !
    ! Calculate overlapping region, taking into account the extra ghost layer
    ! with the addition of +/-mbc:
    ovg = 1
    do n = 1, ndims
      ovg(n,1) = max(mg1(n,1)-mbc,mg2(n,1)-mbc)
      ovg(n,2) = min(mg1(n,2)+mbc,mg2(n,2)+mbc)
    end do
    !
    ! Check for the empty set:
    if (any((ovg(1:ndims,2)-ovg(1:ndims,1))<0)) return
    !
    ! Nonempty; begin transfer:
    !
    ! ovg1 is the global overlap of ovg and mg2
    !
    ! ovg2 is the global overlap of ovg and mg1
    !
    ! These must have nonempty overlap:
    ovg1 = 1
    ovg2 = 1
    do n = 1, ndims
      ovg1(n,1) = max(ovg(n,1),mg2(n,1))
      ovg1(n,2) = min(ovg(n,2),mg2(n,2))
      !
      ovg2(n,1) = max(ovg(n,1),mg1(n,1))
      ovg2(n,2) = min(ovg(n,2),mg1(n,2))
    end do
    !
    ! Convert to local coordinates by subtracting off the global offsets:
    ml1 = 1
    ml2 = 1
    ml1(1:ndims,1) = ovg1(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
    ml1(1:ndims,2) = ovg1(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
    ml2(1:ndims,1) = ovg1(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
    ml2(1:ndims,2) = ovg1(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
    !
    ! Transfer:
    grid1%q(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2),1:nrvars) &
      = grid2%q(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2),1:nrvars)
    !
    ! Convert to local coordinates by subtracting off the global offsets:
    ml1(1:ndims,1) = ovg2(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
    ml1(1:ndims,2) = ovg2(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
    ml2(1:ndims,1) = ovg2(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
    ml2(1:ndims,2) = ovg2(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
    !
    ! Transfer:
    grid2%q(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2),1:nrvars) &
      = grid1%q(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2),1:nrvars)
    !
  end subroutine GhostOverlap
  ! ---------------------------------------------------------------------------
  integer function PeriodicBC(info,timestepparam)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: timestepparam
    integer         :: ibc, mx, mbc, n, nrvars
    !
    PeriodicBC = err_OK
    !
    nrvars = info%nrvars
    !
    ! Check for inactive grid awaiting garbage collection
    if (info%tobedeleted) return
    !
    ! Note qpo is the quasi-periodic offset in q along dimension 1.
    !
    select case(ndims)
    case (2)
      if (info%mthbc(1)==2 .and. info%mthbc(2)==2) then
        mbc = info%mbc
        mx = info%mx(1)
        do ibc = 1, mbc
          do n = 1, nrvars
            info%q(mx+ibc,:,1,n) = info%q(ibc,:,1,n)+qpo(n)
            info%q(1-ibc,:,1,n) = info%q(mx-ibc+1,:,1,n)-qpo(n)
          end do
        end do
      end if
      if (info%mthbc(3)==2 .and. info%mthbc(4)==2) then
        mbc = info%mbc
        mx = info%mx(2)
        do ibc = 1, mbc
          info%q(:,mx+ibc,1,:) = info%q(:,ibc,1,:)
          info%q(:,1-ibc,1,:) = info%q(:,mx-ibc+1,1,:)
        end do
      end if
    case (3)
      if (info%mthbc(1)==2 .and. info%mthbc(2)==2) then
        mbc = info%mbc
        mx = info%mx(1)
        do ibc = 1, mbc
          do n = 1, nrvars
            info%q(mx+ibc,:,:,n) = info%q(ibc,:,:,n)+qpo(n)
            info%q(1-ibc,:,:,n) = info%q(mx-ibc+1,:,:,n)-qpo(n)
          end do
        end do
      end if
      if (info%mthbc(3)==2 .and. info%mthbc(4)==2) then
        mbc = info%mbc
        mx = info%mx(2)
        do ibc = 1, mbc
          info%q(:,mx+ibc,:,:) = info%q(:,ibc,:,:)
          info%q(:,1-ibc,:,:) = info%q(:,mx-ibc+1,:,:)
        end do
      end if
      if (info%mthbc(5)==2 .and. info%mthbc(6)==2) then
        mbc = info%mbc
        mx = info%mx(3)
        do ibc = 1, mbc
          info%q(:,:,mx+ibc,:) = info%q(:,:,ibc,:)
          info%q(:,:,1-ibc,:) = info%q(:,:,mx-ibc+1,:)
        end do
      end if
    case default
      print *,'PeriodicBC: Sorry. Code only works up to ndims=3.'
      stop
    end select
    !
  end function PeriodicBC
  ! ---------------------------------------------------------------------------
  integer function TransferPeriodicBC(grid1,grid2,dummy)
    !
    ! Transfer periodic boundary values among grids on same level:
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: grid1, grid2
    type(funcparam) :: dummy
    !
    TransferPeriodicBC = err_ok
    !
    ! Check for inactive grids awaiting garbage collection:
    if (grid1%tobedeleted .or. grid2%tobedeleted) return
    !
    ! Look for periodic overlap of ghost cell regions of grid1 and grid2.  Order
    ! doesn't matter, as transfer goes both ways:
    call PeriodicGhostOverlap(grid1,grid2)
    call PeriodicGhostOverlap(grid2,grid1)
    !
  end function TransferPeriodicBC
  ! ---------------------------------------------------------------------------
  subroutine PeriodicGhostOverlap(grid1,grid2)
    !
    ! Assumes ghost layer of size mbc, and that periodic self-transfer has
    ! occurred on each grid.  The periodic offsets work up to 3D only.
    !
    use NodeinfoDef
    implicit none
    type(nodeinfo)                    :: grid1
    type(nodeinfo)                    :: grid2
    integer                           :: mbc, n, nrvars, offset
    integer, dimension(1:maxdims,1:2) :: mg1, mg2, ml1, ml2, ovg, ovg1, ovg2
    !
    mbc = grid2%mbc
    nrvars = grid1%nrvars
    !
    call GetPeriodicOffsets(grid1)
    !
    mg1 = 1
    mg1(1:ndims,1:2) = grid1%mglobal(1:ndims,1:2)
    !
    offset_loop: do offset = 1, nperiodicoffsets
      !
      mg2 = 1
      do n = 1, ndims
        mg2(n,1:2) = grid2%mglobal(n,1:2)+poffset(n,offset)
      end do
      !
      ! Calculate overlapping region, taking into account the extra ghost layer
      ! with the addition of +/-mbc:
      ovg = 1
      do n = 1, ndims
        ovg(n,1) = max(mg1(n,1)-mbc,mg2(n,1)-mbc)
        ovg(n,2) = min(mg1(n,2)+mbc,mg2(n,2)+mbc)
      end do
      !
      ! Check for the empty set:
      if (any((ovg(1:ndims,2)-ovg(1:ndims,1))<0)) cycle offset_loop
      !
      ! Nonempty; begin transfer:
      !
      !
      ! ovg1 is the global overlap of ovg and mg2
      !
      ! ovg2 is the global overlap of ovg and mg1
      !
      ! These must have nonempty overlap:
      ovg1 = 1
      ovg2 = 1
      do n = 1, ndims
        ovg1(n,1) = max(ovg(n,1),mg2(n,1))
        ovg1(n,2) = min(ovg(n,2),mg2(n,2))
        !
        ovg2(n,1) = max(ovg(n,1),mg1(n,1))
        ovg2(n,2) = min(ovg(n,2),mg1(n,2))
      end do
      !
      ! Convert to local coordinates by subtracting off the global offsets:
      ml1 = 1
      ml2 = 1
      ml1(1:ndims,1) = ovg1(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
      ml1(1:ndims,2) = ovg1(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
      ml2(1:ndims,1) = ovg1(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
      ml2(1:ndims,2) = ovg1(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
      !
      ! Transfer:
      do n = 1, nrvars
        grid1%q(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2),n) &
          = grid2%q(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2),n) &
            + qoffset(offset,n)
      end do
      !
      ! Convert to local coordinates by subtracting off the global offsets:
      ml1(1:ndims,1) = ovg2(1:ndims,1)-mg1(1:ndims,1)+1 ! Grid 1 local
      ml1(1:ndims,2) = ovg2(1:ndims,2)-mg1(1:ndims,1)+1 ! Grid 1 local
      ml2(1:ndims,1) = ovg2(1:ndims,1)-mg2(1:ndims,1)+1 ! Grid 2 local
      ml2(1:ndims,2) = ovg2(1:ndims,2)-mg2(1:ndims,1)+1 ! Grid 2 local
      !
      ! Transfer:
      do n = 1, nrvars
        grid2%q(ml2(1,1):ml2(1,2),ml2(2,1):ml2(2,2),ml2(3,1):ml2(3,2),n) &
          = grid1%q(ml1(1,1):ml1(1,2),ml1(2,1):ml1(2,2),ml1(3,1):ml1(3,2),n) &
            - qoffset(offset,n)
      end do
      !
    end do offset_loop
    !
  end subroutine PeriodicGhostOverlap
  ! ---------------------------------------------------------------------------
  integer function SetBC(info,dummy)
    use NodeinfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                  :: info
    type(funcparam)                 :: dummy
    integer                         :: ll, mbc, nrvars
    integer, dimension(1:maxdims)   :: mx, ul
    real, dimension(1:maxdims)      :: xlower,dx
    integer, dimension(1:2*maxdims) :: mthbc
    integer                         :: ipass
    !
    ! Red-black switch when used:
    ipass = dummy%iswitch
    !
    SetBC = err_ok
    !
    mbc = info%mbc
    ll = 1-mbc
    mx = 1
    mx(1:ndims) = info%mx(1:ndims)
    ul = 1
    ul(1:ndims) = mx(1:ndims)+mbc
    dx = info%dx
    xlower = info%xlower
    nrvars = info%nrvars
    !
    mthbc = 1
    mthbc(1:2*ndims) = info%mthbc(1:2*ndims)
    !
    ! Check for inactive grids awaiting garbage collection
    if (info%tobedeleted) return
    select case(ndims)
    case (2)
      call SetBC2D(info%q(ll:ul(1),ll:ul(2),       1,1:nrvars), &
                   ll,mx(1:2),nrvars,dx(1),xlower,mbc,mthbc(1:4))
    case (3)
      call SetBC3D(info%q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars), &
                   ll,mx(1:3),nrvars,mbc,mthbc(1:6))
    case default
      print *, 'SetBC: Only ndims = 2,3 are supported.'
    end select
    !
  end function SetBC
  ! ---------------------------------------------------------------------------
  subroutine SetBC2D(q,ll,mx,nrvars,h,xlower,mbc,mthbc)
    use NodeinfoDef
    use Problem, only: UserBC2D
    implicit none
    real, dimension(ll:,ll:,1:), intent(in out) :: q
    integer, intent(in)                         :: ll
    integer, dimension(1:2), intent(in)         :: mx
    real, intent(in)                            :: h
    real, dimension(1:2), intent(in)            :: xlower
    integer, intent(in)                         :: nrvars
    integer, intent(in)                         :: mbc
    integer, dimension(1:4), intent(in)         :: mthbc
    integer                                     :: ibc
    integer, dimension(1:2)                     :: ul
    !
    ul(1:2) = mx(1:2)+mbc
    !
    ! Left Boundary along dimension 1:
    select case(mthbc(1))
    case(1)
      do ibc = 1, mbc
        q(1-ibc,ll:ul(2),1:nrvars) = q(ibc,ll:ul(2),1:nrvars)
      end do
    case(3)
      print *, 'SetBC2D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC2D(q(ll:ul(1),ll:ul(2),1:nrvars),ll,mx,nrvars,h,xlower,mbc,1)
    case default
      print 1001, mthbc(1), 1
      stop
    end select
    !
    ! Right Boundary along dimension 1:
    select case(mthbc(2))
    case(1)
      do ibc = 1, mbc
        q(mx(1)+ibc,ll:ul(2),1:nrvars) = q(mx(1)-ibc+1,ll:ul(2),1:nrvars)
      end do
    case(3)
      print *, 'SetBC2D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC2D(q(ll:ul(1),ll:ul(2),1:nrvars),ll,mx,nrvars,h,xlower,mbc,2)
    case default
      print 1001, mthbc(2), 2
      stop
    end select
    !
    ! Left Boundary along dimension 2:
    select case(mthbc(3))
    case(1)
      do ibc = 1, mbc
        q(ll:ul(1),1-ibc,1:nrvars) = q(ll:ul(1),  ibc,1:nrvars)
      end do
    case(3)
      print *, 'SetBC2D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC2D(q(ll:ul(1),ll:ul(2),1:nrvars),ll,mx,nrvars,h,xlower,mbc,3)
    case default
      print 1001, mthbc(3), 3
      stop
    end select
    !
    ! Right Boundary along dimension 2:
    select case(mthbc(4))
    case(1)
      do ibc = 1, mbc
        q(ll:ul(1),mx(2)+ibc  ,1:nrvars) = q(ll:ul(1),mx(2)-ibc+1,1:nrvars)
      end do
    case(3)
      print *, 'SetBC2D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC2D(q(ll:ul(1),ll:ul(2),1:nrvars),ll,mx,nrvars,h,xlower,mbc,4)
    case default
      print 1001, mthbc(4), 4
      stop
    end select
    !
    ! formats
    1001 format('SetBC2D: User boundary condition ', i4, &
            ' on boundary ', i1, //, '         not coded in SetBC2D.')
    !
  end subroutine SetBC2D
  ! ---------------------------------------------------------------------------
  subroutine SetBC3D(q,ll,mx,nrvars,mbc,mthbc)
    use NodeinfoDef
    use Problem, only: UserBC3D
    implicit none
    real, dimension(ll:,ll:,ll:,1:), intent(in out) :: q
    integer, intent(in)                             :: ll
    integer, dimension(1:3), intent(in)             :: mx
    integer, intent(in)                             :: nrvars
    integer, intent(in)                             :: mbc
    integer, dimension(1:6), intent(in)             :: mthbc
    integer                                         :: ibc
    integer, dimension(1:3)                         :: ul
    !
    ul(1:3) = mx(1:3)+mbc
    !
    ! Left Boundary along dimension 1:
    select case(mthbc(1))
    case(1)
      do ibc = 1, mbc
        q(1-ibc,ll:ul(2),ll:ul(3),1:nrvars) &
          = q(  ibc,ll:ul(2),ll:ul(3),1:nrvars)
      end do
    case(3)
      print *, 'SetBC3D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC3D(q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,1)
    case default
      print 1001, mthbc(1), 1
      stop
    end select
    !
    ! Right Boundary along dimension 1:
    select case(mthbc(2))
    case(1)
      do ibc = 1, mbc
        q(mx(1)+ibc  ,ll:ul(2),ll:ul(3),1:nrvars) &
          = q(mx(1)-ibc+1,ll:ul(2),ll:ul(3),1:nrvars)
      end do
    case(3)
      print *, 'SetBC3D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC3D(q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,2)
    case default
      print 1001, mthbc(2), 2
      stop
    end select
    !
    ! Left Boundary along dimension 2:
    select case(mthbc(3))
    case(1)
      do ibc = 1, mbc
        q(ll:ul(1),1-ibc,ll:ul(3),1:nrvars) &
          = q(ll:ul(1),  ibc,ll:ul(3),1:nrvars)
      end do
    case(3)
      print *, 'SetBC3D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC3D(q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,3)
    case default
      print 1001, mthbc(3), 3
      stop
    end select
    !
    ! Right Boundary along dimension 2:
    select case(mthbc(4))
    case(1) ! Zero-order extrapolation:
      do ibc = 1, mbc
        q(ll:ul(1),mx(2)+ibc  ,ll:ul(3),1:nrvars) &
          = q(ll:ul(1),mx(2)-ibc+1,ll:ul(3),1:nrvars)
      end do
    case(3)
      print *, 'SetBC3D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC3D(q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,4)
    case default
      print 1001, mthbc(4), 4
      stop
    end select
    !
    ! Left Boundary along dimension 3:
    select case(mthbc(5))
    case(1)
      do ibc = 1, mbc
        q(ll:ul(1),ll:ul(2),1-ibc,1:nrvars) &
          = q(ll:ul(1),ll:ul(2),  ibc,1:nrvars)
      end do
    case(3)
      print *, 'SetBC3D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC3D(q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,5)
    case default
      print 1001, mthbc(5), 5
      stop
    end select
    !
    ! Right Boundary along dimension 3:
    select case(mthbc(6))
    case(1) ! Zero-order extrapolation:
      do ibc = 1, mbc
        q(ll:ul(1),ll:ul(2),mx(3)+ibc  ,1:nrvars) &
          = q(ll:ul(1),ll:ul(2),mx(3)-ibc+1,1:nrvars)
      end do
    case(3)
      print *, 'SetBC3D: mthbc = 3 not in use.'
      stop
    case(2,999)
      ! Periodic, internal boundary condition. Coded below:
      continue
    case(10:99)
      ! User boundary condition. Call Problem module:
      call UserBC3D(q(ll:ul(1),ll:ul(2),ll:ul(3),1:nrvars),ll,mx,nrvars,mbc,6)
    case default
      print 1001, mthbc(6), 6
      stop
    end select
    !
    ! formats
    1001 format('SetBC3D: User boundary condition ', i4, &
                ' on boundary ', i1, //, '         not coded in SetBC3D.')
    !
  end subroutine SetBC3D
  ! ---------------------------------------------------------------------------
  subroutine PeriodicSetup(rootinfo)
    !
    ! This works for up to three dimensions.  We assume that there are no
    ! inconsistencies in the boundary condition flags mthbc:
    !
    use NodeInfoDef
    use BSAMStorage, only: AllocPeriodicBCStorage
    implicit none
    type(nodeinfo)                  :: rootinfo
    integer                         :: np, nrvars
    integer, dimension(1:2*maxdims) :: mthbc
    !
    mthbc = 1
    mthbc(1:2*ndims) = rootinfo%mthbc(1:2*ndims)
    rootinfo%mthbc = mthbc
    nrvars = rootinfo%nrvars
    !
    periodicboundaryconditions = .true.
    !
    if (mthbc(1)==2) then
      if (mthbc(3)==2) then
        if (mthbc(5)==2) then
          nperiodicoffsets = 13
          np = nperiodicoffsets
          call AllocPeriodicBCStorage(np,nrvars)
          periodicoffsetindex(1:np) = (/1,2,3,4,5,6,7,8,9,10,11,12,13/)
        else
          nperiodicoffsets = 4
          np = nperiodicoffsets
          call AllocPeriodicBCStorage(np,nrvars)
          periodicoffsetindex(1:np) = (/1,2,4,5/)
        end if
      else
        if (mthbc(5)==2) then
          nperiodicoffsets = 4
          np = nperiodicoffsets
          call AllocPeriodicBCStorage(np,nrvars)
          periodicoffsetindex(1:np) = (/1,3,8,9/)
        else
          nperiodicoffsets = 1
          np = nperiodicoffsets
          call AllocPeriodicBCStorage(np,nrvars)
          periodicoffsetindex(1:np) = (/1/)
        end if
      end if
    else
      if (mthbc(3)==2) then
        if (mthbc(5)==2) then
          nperiodicoffsets = 4
          np = nperiodicoffsets
          call AllocPeriodicBCStorage(np,nrvars)
          periodicoffsetindex(1:np) = (/2,3,6,7/)
        else
          nperiodicoffsets = 1
          np = nperiodicoffsets
          call AllocPeriodicBCStorage(np,nrvars)
          periodicoffsetindex(1:np) = (/2/)
        end if
      else
        if (mthbc(5)==2) then
          nperiodicoffsets = 1
          np = nperiodicoffsets
          call AllocPeriodicBCStorage(np,nrvars)
          periodicoffsetindex(1:np) = (/3/)
        else
          periodicboundaryconditions = .false.
        end if
      end if
    end if
    !
  end subroutine PeriodicSetup
  ! ---------------------------------------------------------------------------
  subroutine GetPeriodicOffsets(info)
    use NodeInfoDef
    implicit none
    type(nodeinfo) :: info
    integer        :: level, nrvars, offset
    !
    level = info%level
    nrvars = info%nrvars
    !
    do offset = 1, nperiodicoffsets
      select case(periodicoffsetindex(offset))
      case(1)
        poffset(1,offset) =  mxmax(level,1)
        poffset(2,offset) =  0
        poffset(3,offset) =  0
        qoffset(offset,1:nrvars) =  qpo(1:nrvars)
      case(2)
        poffset(1,offset) =  0
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) =  0
        qoffset(offset,1:nrvars) =  0.0
      case(3)
        poffset(1,offset) =  0
        poffset(2,offset) =  0
        poffset(3,offset) =  mxmax(level,3)
        qoffset(offset,1:nrvars) =  0.0
      case(4)
        poffset(1,offset) =  mxmax(level,1)
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) =  0
        qoffset(offset,1:nrvars) =  qpo(1:nrvars)
      case(5)
        poffset(1,offset) = -mxmax(level,1)
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) =  0
        qoffset(offset,1:nrvars) = -qpo(1:nrvars)
      case(6)
        poffset(1,offset) =  0
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) =  mxmax(level,3)
        qoffset(offset,1:nrvars) =  0.0
      case(7)
        poffset(1,offset) =  0
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) = -mxmax(level,3)
        qoffset(offset,1:nrvars) =  0.0
      case(8)
        poffset(1,offset) =  mxmax(level,1)
        poffset(2,offset) =  0
        poffset(3,offset) =  mxmax(level,3)
        qoffset(offset,1:nrvars) =  qpo(1:nrvars)
      case(9)
        poffset(1,offset) = -mxmax(level,1)
        poffset(2,offset) =  0
        poffset(3,offset) =  mxmax(level,3)
        qoffset(offset,1:nrvars) = -qpo(1:nrvars)
      case(10)
        poffset(1,offset) =  mxmax(level,1)
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) =  mxmax(level,3)
        qoffset(offset,1:nrvars) =  qpo(1:nrvars)
      case(11)
        poffset(1,offset) =  mxmax(level,1)
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) = -mxmax(level,3)
        qoffset(offset,1:nrvars) =  qpo(1:nrvars)
      case(12)
        poffset(1,offset) = -mxmax(level,1)
        poffset(2,offset) =  mxmax(level,2)
        poffset(3,offset) =  mxmax(level,3)
        qoffset(offset,1:nrvars) = -qpo(1:nrvars)
      case(13)
        poffset(1,offset) =  mxmax(level,1)
        poffset(2,offset) = -mxmax(level,2)
        poffset(3,offset) =  mxmax(level,3)
        qoffset(offset,1:nrvars) =  qpo(1:nrvars)
      end select
    end do
    !
  end subroutine GetPeriodicOffsets
  ! ---------------------------------------------------------------------------
  subroutine GetPeriodicTagOffset(coordinate,level,offset,periodicbuffer)
    use NodeInfoDef
    implicit none
    integer, dimension(1:maxdims), intent(in out) :: coordinate
    integer, intent(in)                           :: level
    integer, intent(in)                           :: offset
    logical, intent(out)                          :: periodicbuffer
    integer, dimension(1:3,1:2)                   :: mg
    !
    ! Global coordinates of the buffer patch:
    mg = 1
    mg(1:ndims,1) = coordinate(1:ndims)-ibuffer(level)
    mg(1:ndims,2) = coordinate(1:ndims)+ibuffer(level)
    !
    periodicbuffer = .true.
    !
    select case(offset)
    case(1) ! x1 direction:
      if (mg(1,2)>mxmax(level,1)) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
      else if (mg(1,1)<1) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
      else
        periodicbuffer = .false.
      end if
    case(2) ! x2 direction:
      if (mg(2,2)>mxmax(level,2)) then
        coordinate(2) = coordinate(2)-mxmax(level,2)
      else if (mg(2,1)<1) then
        coordinate(2) = coordinate(2)+mxmax(level,2)
      else
        periodicbuffer = .false.
      end if
    case(3) ! x3 direction:
      if (mg(3,2)>mxmax(level,3)) then
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else if (mg(3,1)<1) then
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(4) ! x1 and x2 directions:
      if (mg(1,2)>mxmax(level,1) .and. mg(2,2)>mxmax(level,2)) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(2) = coordinate(2)-mxmax(level,2)
      else if (mg(1,1)<1 .and. mg(2,1)<1) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(2) = coordinate(2)+mxmax(level,2)
      else
        periodicbuffer = .false.
      end if
    case(5) ! x1 and -x2 directions:
      if (mg(1,2)>mxmax(level,1) .and. mg(2,1)<1) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(2) = coordinate(2)+mxmax(level,2)
      else if (mg(1,1)<1 .and. mg(2,2)>mxmax(level,2)) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(2) = coordinate(2)-mxmax(level,2)
      else
        periodicbuffer = .false.
      end if
    case(6) ! x2 and x3 directions:
      if (mg(2,2)>mxmax(level,2) .and. mg(3,2)>mxmax(level,3)) then
        coordinate(2) = coordinate(2)-mxmax(level,2)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else if (mg(2,1)<1 .and. mg(3,1)<1) then
        coordinate(2) = coordinate(2)+mxmax(level,2)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(7) ! x2 and -x3 directions:
      if (mg(2,2)>mxmax(level,2) .and. mg(3,1)<1) then
        coordinate(2) = coordinate(2)-mxmax(level,2)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else if (mg(2,1)<1 .and. mg(3,2)>mxmax(level,3)) then
        coordinate(2) = coordinate(2)+mxmax(level,2)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(8) ! x1 and x3 directions:
      if (mg(1,2)>mxmax(level,1) .and. mg(3,2)>mxmax(level,3)) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else if (mg(1,1)<1 .and. mg(3,1)<1) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(9) ! -x1 and x3 directions:
      if (mg(1,2)>mxmax(level,1) .and. mg(3,1)<1) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else if (mg(1,1)<1 .and. mg(3,2)>mxmax(level,2)) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(10) ! x1, x2 and x3 directions:
      if (mg(1,2)>mxmax(level,1) .and. mg(2,2)>mxmax(level,2) &
          .and. mg(3,2)>mxmax(level,3)) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(2) = coordinate(2)-mxmax(level,2)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else if (mg(1,1)<1 .and. mg(2,1)<1 .and. mg(3,1)<1) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(2) = coordinate(2)+mxmax(level,2)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(11) ! x1, x2 and -x3 directions:
      if (mg(1,2)>mxmax(level,1) .and. mg(2,2)>mxmax(level,2) &
          .and. mg(3,1)<1) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(2) = coordinate(2)-mxmax(level,2)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else if (mg(1,1)<1 .and. mg(2,1)<1 .and. mg(3,2)>mxmax(level,3)) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(2) = coordinate(2)+mxmax(level,2)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(12) ! -x1, x2 and x3 directions:
      if (mg(1,1)<1 .and. mg(2,2)>mxmax(level,2) &
          .and. mg(3,2)>mxmax(level,3)) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(2) = coordinate(2)-mxmax(level,2)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else if (mg(1,2)>mxmax(level,1) .and. mg(2,1)<1 .and. mg(3,1)<1) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(2) = coordinate(2)+mxmax(level,2)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    case(13) ! x1, -x2 and x3 directions:
      if (mg(1,2)>mxmax(level,1) .and. mg(2,1)<1 &
          .and. mg(3,2)>mxmax(level,3)) then
        coordinate(1) = coordinate(1)-mxmax(level,1)
        coordinate(2) = coordinate(2)+mxmax(level,2)
        coordinate(3) = coordinate(3)-mxmax(level,3)
      else if (mg(1,1)<1 .and. mg(2,2)>mxmax(level,2) .and. mg(3,1)<1) then
        coordinate(1) = coordinate(1)+mxmax(level,1)
        coordinate(2) = coordinate(2)-mxmax(level,2)
        coordinate(3) = coordinate(3)+mxmax(level,3)
      else
        periodicbuffer = .false.
      end if
    end select
    !
  end subroutine GetPeriodicTagOffset
  ! ---------------------------------------------------------------------------
  integer function GetCoarseGhostPoints(info,dummy)
    !
    ! Dimensionally invariant routine to fill coarse grid ghost points.
    !
    use NodeInfoDef
    use TreeOps, only: GetParentInfo, err_ok
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    type(nodeinfo), pointer           :: parent
    integer                           :: ierror, nrvars
    integer, dimension(1:maxdims)     :: cmx, ih, il, jh, jl, mx
    integer, dimension(1:maxdims,1:2) :: mb
    !
    GetCoarseGhostPoints = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    mx = 1
    cmx = 1
    il = 1
    ih = 1
    mb = 1
    jl = 1
    jh = 1
    mx(1:ndims) = info%mx(1:ndims)
    cmx(1:ndims) = mx(1:ndims)/2
    il(1:ndims) = 0
    ih(1:ndims) = cmx(1:ndims)+1
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    jl(1:ndims) = mb(1:ndims,1)-1
    jh(1:ndims) = mb(1:ndims,2)+1
    !
    ierror = GetParentInfo(parent)
    !
    ! Face #1:
    info%qc(il(1),il(2):ih(2),il(3):ih(3),1:nrvars) &
            = parent%q( jl(1),jl(2):jh(2),jl(3):jh(3),1:nrvars)
    !
    ! Face #2:
    info%qc(ih(1),il(2):ih(2),il(3):ih(3),1:nrvars) &
            = parent%q( jh(1),jl(2):jh(2),jl(3):jh(3),1:nrvars)
    !
    if (ndims>=2) then
      !
      ! Face #3:
      info%qc(     1:cmx(1)  ,il(2),il(3):ih(3),1:nrvars) &
              = parent%q(mb(1,1): mb(1,2),jl(2),jl(3):jh(3),1:nrvars)
      !
      ! Face #4:
      info%qc(     1:cmx(1)  ,ih(2),il(3):ih(3),1:nrvars) &
              = parent%q(mb(1,1): mb(1,2),jh(2),jl(3):jh(3),1:nrvars)
      !
      if (ndims>=3) then
        !
        ! Face #5:
        info%qc(     1:cmx(1) ,       1:cmx(2)  ,il(3),1:nrvars) &
                = parent%q(mb(1,1): mb(1,2),mb(2,1): mb(2,2),jl(3),1:nrvars)
        !
        ! Face #6:
        info%qc(     1:cmx(1) ,       1:cmx(2)  ,ih(3),1:nrvars) &
                = parent%q(mb(1,1): mb(1,2),mb(2,1): mb(2,2),jh(3),1:nrvars)
      end if
    end if
    !
  end function GetCoarseGhostPoints
  ! ---------------------------------------------------------------------------
  subroutine GetFaceIndex(mx,mb,dim,parity,fi)
    implicit none
    integer, dimension(1:3), intent(in)     :: mx
    integer, dimension(1:3,1:2), intent(in) :: mb
    integer, intent(in)                     :: dim
    integer, intent(in)                     :: parity
    integer, dimension(1:5), intent(out)    :: fi
    !
    select case(parity)
    case(1)
      fi(1) = 0
      fi(2) = 2
      fi(3) = 1
      fi(4) = mb(dim,parity)-1
      fi(5) = mb(dim,parity)
    case(2)
      fi(1) = mx(dim)+1
      fi(2) = mx(dim)-1
      fi(3) = mx(dim)
      fi(4) = mb(dim,parity)+1
      fi(5) = mb(dim,parity)
    case default
      print *, 'GetFaceIndex: Parity error'
      stop
    end select
    !
  end subroutine GetFaceIndex
  ! ---------------------------------------------------------------------------
  subroutine Interpolate3D(qf,qc,mx,nrvars,mb,mthbc,ipass)
    implicit none
    real, dimension(0:,0:,0:,1:), intent(in out) :: qf
    real, dimension(0:,0:,0:,1:), intent(in)     :: qc
    integer, dimension(1:3), intent(in)          :: mx
    integer, intent(in)                          :: nrvars
    integer, dimension(1:3,1:2), intent(in)      :: mb
    integer, dimension(1:6), intent(in)          :: mthbc
    integer, intent(in)                          :: ipass
    logical, dimension(1:6)                      :: cycle_face
    integer                                      :: fc
    integer                                      :: parity_1,parity_2,parity_3
    integer, dimension(1:3,1:5)                  :: fi
    !
    ! ipass = 1: red (or odd-odd and even-even) squares should be updated,
    ! ipass = 2: black (or even-odd and odd-even) squares should be updated.
    ! ipass = 0: both red and black squares should be updated.
    !
    ! Faces that are on physical boundaries have mthbc(face) equal to a number
    ! other than 2 and 999. If the face is on a physical boundary its ghost
    ! points aren't found by interpolation.
    ! interpol
    !
    ! (1) Interpolate the 8 corners of each box.  We only interpolate on the
    !     physical boundaries if ipass = 0, 2:
    parity_1c_loop: do parity_1 = 1, 2
      if (mthbc(parity_1)/=2 .and. &
          mthbc(parity_1)/=999 .and. ipass==1) cycle parity_1c_loop
      call GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
      parity_2c_loop: do parity_2 = 1, 2
        if (mthbc(parity_2+2)/=2 .and. &
            mthbc(parity_2+2)/=999 .and. ipass==1) cycle parity_2c_loop
        call GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
        parity_3c_loop: do parity_3 = 1, 2
          if (mthbc(parity_3+4)/=2 .and. &
              mthbc(parity_3+4)/=999 .and. ipass==1) cycle parity_3c_loop
          call GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
          qf(fi(1,1),fi(2,1),fi(3,1),1:nrvars) = &
                  a1*qf(fi(1,2),fi(2,2),fi(3,2),1:nrvars) &
                + a2*qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
                + a3*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars)
        end do parity_3c_loop
      end do parity_2c_loop
    end do parity_1c_loop
    !
    ! (2) Interpolate the 12 edges of each box.  We only interpolate on the
    !     physical boundaries if ipass = 0, 2:
    !
    ! Edges 1--4 (coordinates 1 and 3 fixed):
    parity_1e1_loop: do parity_1 = 1, 2
      if (mthbc(parity_1)/=2 .and. &
          mthbc(parity_1)/=999 .and. ipass==1) cycle parity_1e1_loop
      call GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
      parity_3e1_loop: do parity_3 = 1, 2
        if (mthbc(parity_3+4)/=2 .and. &
            mthbc(parity_3+4)/=99 .and. ipass==1) cycle parity_3e1_loop
        call GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
        !
        qf(fi(1,1),      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
               a1*qf(fi(1,2),      1:mx(2)-3:2,fi(3,2),1:nrvars) &
             + a2*qf(fi(1,3),      2:mx(2)-2:2,fi(3,3),1:nrvars) &
             + a3*qc(fi(1,4),mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars)
        !
        qf(fi(1,1),        1,fi(3,1),1:nrvars) = &
               a1*    qf(fi(1,2),        1,fi(3,2),1:nrvars) &
             + a2*    qf(fi(1,3),        1,fi(3,3),1:nrvars) &
             + a3*(b1*qc(fi(1,4),mb(2,1)-1,fi(3,4),1:nrvars) &
             +     b2*qc(fi(1,4),mb(2,1)  ,fi(3,4),1:nrvars) &
             +     b3*qc(fi(1,4),mb(2,1)+1,fi(3,4),1:nrvars))
        !
        qf(fi(1,1),      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
               a1*qf(fi(1,2),      4:mx(2)  :2,fi(3,2),1:nrvars) &
             + a2*qf(fi(1,3),      3:mx(2)-1:2,fi(3,3),1:nrvars) &
             + a3*qc(fi(1,4),mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars)
        !
        qf(fi(1,1),    mx(2),fi(3,1),1:nrvars) = &
               a1*    qf(fi(1,2),    mx(2),fi(3,2),1:nrvars) &
             + a2*    qf(fi(1,3),    mx(2),fi(3,3),1:nrvars) &
             + a3*(b1*qc(fi(1,4),mb(2,2)+1,fi(3,4),1:nrvars) &
             +     b2*qc(fi(1,4),mb(2,2)  ,fi(3,4),1:nrvars) &
             +     b3*qc(fi(1,4),mb(2,2)-1,fi(3,4),1:nrvars))
        !
      end do parity_3e1_loop
    end do parity_1e1_loop
    !
    ! Edges 5--8 (coordinates 1 and 2 fixed):
    parity_1e2_loop: do parity_1 = 1, 2
      if (mthbc(parity_1)/=2 .and. &
          mthbc(parity_1)/= 999 .and. ipass==1) cycle parity_1e2_loop
      call GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
      parity_2e2_loop: do parity_2 = 1, 2
        if (mthbc(parity_2+2)/=2 .and. &
            mthbc(parity_2+2)/=999 .and. ipass==1) cycle parity_2e2_loop
        call GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
        !
        qf(fi(1,1),fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
               a1*qf(fi(1,2),fi(2,2),      1:mx(3)-3:2,1:nrvars) &
             + a2*qf(fi(1,3),fi(2,3),      2:mx(3)-2:2,1:nrvars) &
             + a3*qc(fi(1,4),fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars)
        !
        qf(fi(1,1),fi(2,1),        1,1:nrvars) = &
               a1*    qf(fi(1,2),fi(2,2),        1,1:nrvars) &
             + a2*    qf(fi(1,3),fi(2,3),        1,1:nrvars) &
             + a3*(b1*qc(fi(1,4),fi(2,4),mb(3,1)-1,1:nrvars) &
             +     b2*qc(fi(1,4),fi(2,4),mb(3,1)  ,1:nrvars) &
             +     b3*qc(fi(1,4),fi(2,4),mb(3,1)+1,1:nrvars))
        !
        qf(fi(1,1),fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
               a1*qf(fi(1,2),fi(2,2),      4:mx(3)  :2,1:nrvars) &
             + a2*qf(fi(1,3),fi(2,3),      3:mx(3)-1:2,1:nrvars) &
             + a3*qc(fi(1,4),fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars)
        !
        qf(fi(1,1),fi(2,1),    mx(3),1:nrvars) = &
               a1*    qf(fi(1,2),fi(2,2),    mx(3),1:nrvars) &
             + a2*    qf(fi(1,3),fi(2,3),    mx(3),1:nrvars) &
             + a3*(b1*qc(fi(1,4),fi(2,4),mb(3,2)+1,1:nrvars) &
             +     b2*qc(fi(1,4),fi(2,4),mb(3,2)  ,1:nrvars) &
             +     b3*qc(fi(1,4),fi(2,4),mb(3,2)-1,1:nrvars))
        !
      end do parity_2e2_loop
    end do parity_1e2_loop
    !
    ! Edges 9--12 (coordinates 2 and 3 fixed):
    parity_2e3_loop: do parity_2 = 1, 2
      if (mthbc(parity_2+2)/=2 .and. &
          mthbc(parity_2+2)/=999 .and. ipass==1) cycle parity_2e3_loop
      call GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
      parity_3e3_loop: do parity_3 = 1, 2
        if (mthbc(parity_3+4)/=2 .and. &
            mthbc(parity_3+4)/=999 .and. ipass==1) cycle parity_3e3_loop
        call GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
        !
        qf(      3:mx(1)-1:2,fi(2,1),fi(3,1),1:nrvars) = &
               a1*qf(      1:mx(1)-3:2,fi(2,2),fi(3,2),1:nrvars) &
             + a2*qf(      2:mx(1)-2:2,fi(2,3),fi(3,3),1:nrvars) &
             + a3*qc(mb(1,1)+1:mb(1,2),fi(2,4),fi(3,4),1:nrvars)
        !
        qf(        1,fi(2,1),fi(3,1),1:nrvars) = &
               a1*    qf(        1,fi(2,2),fi(3,2),1:nrvars) &
             + a2*    qf(        1,fi(2,3),fi(3,3),1:nrvars) &
             + a3*(b1*qc(mb(1,1)-1,fi(2,4),fi(3,4),1:nrvars) &
             +     b2*qc(mb(1,1)  ,fi(2,4),fi(3,4),1:nrvars) &
             +     b3*qc(mb(1,1)+1,fi(2,4),fi(3,4),1:nrvars))
        !
        qf(      2:mx(1)-2:2,fi(2,1),fi(3,1),1:nrvars) = &
               a1*qf(      4:mx(1)  :2,fi(2,2),fi(3,2),1:nrvars) &
             + a2*qf(      3:mx(1)-1:2,fi(2,3),fi(3,3),1:nrvars) &
             + a3*qc(mb(1,1):mb(1,2)-1,fi(2,4),fi(3,4),1:nrvars)
        !
        qf(    mx(1),fi(2,1),fi(3,1),1:nrvars) = &
               a1*    qf(    mx(1),fi(2,2),fi(3,2),1:nrvars) &
             + a2*    qf(    mx(1),fi(2,3),fi(3,3),1:nrvars) &
             + a3*(b1*qc(mb(1,2)+1,fi(2,4),fi(3,4),1:nrvars) &
             +     b2*qc(mb(1,2)  ,fi(2,4),fi(3,4),1:nrvars) &
             +     b3*qc(mb(1,2)-1,fi(2,4),fi(3,4),1:nrvars))
        !
      end do parity_3e3_loop
    end do parity_2e3_loop
    !
    ! (3) Interpolate the 6 faces:
    !
    ! Only interpolate periodic or internal interfaces:
    do fc = 1, 6
      select case(mthbc(fc))
      case(2,999)
        cycle_face(fc) = .false.
      case default
        cycle_face(fc) = .true.
      end select
    end do
    !
    ! ipass = 0: Red and Black
    ! ipass = 1: Red
    ! ipass = 2: Black
    ! parity = 1: Yellow, Blue = Red
    !             Green, Red   = Black
    ! parity = 2: Green, Red   = Red
    !             Yellow, Blue = Black
    !
    ! Faces 1 and 2 (coordinate 1 fixed):
    parity_1_loop: do parity_1 = 1, 2
      if (cycle_face(parity_1)) cycle parity_1_loop
      call GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
      !
      if (ipass==0 .or. modulo(parity_1+ipass,2)==0) then
        !
        ! Yellow:
        qf(fi(1,1),      2:mx(2)-2:2,      2:mx(3)-2:2,1:nrvars) = &
              a1*qf(fi(1,2),      4:mx(2)  :2,      4:mx(3)  :2,1:nrvars) &
            + a2*qf(fi(1,3),      3:mx(2)-1:2,      3:mx(3)-1:2,1:nrvars) &
            + a3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1):mb(3,2)-1,1:nrvars)
        !
        qf(fi(1,1),    mx(2),      2:mx(3)-2:2,1:nrvars) = &
              a1*    qf(fi(1,2),    mx(2),      4:mx(3)  :2,1:nrvars) &
            + a2*    qf(fi(1,3),    mx(2),      3:mx(3)-1:2,1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,2)+1,mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,2)  ,mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,2)-1,mb(3,1):mb(3,2)-1,1:nrvars))
        !
        qf(fi(1,1),      2:mx(2)-2:2,    mx(3),1:nrvars) = &
              a1*    qf(fi(1,2),      4:mx(2)  :2,    mx(3),1:nrvars) &
            + a2*    qf(fi(1,3),      3:mx(2)-1:2,    mx(3),1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,2)+1,1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,2)  ,1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,2)-1,1:nrvars))
        !
        ! Blue:
        qf(fi(1,1),      3:mx(2)-1:2,      3:mx(3)-1:2,1:nrvars) = &
              a1*qf(fi(1,2),      1:mx(2)-3:2,      1:mx(3)-3:2,1:nrvars) &
            + a2*qf(fi(1,3),      2:mx(2)-2:2,      2:mx(3)-2:2,1:nrvars) &
            + a3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)+1:mb(3,2),1:nrvars)
        !
        qf(fi(1,1),      3:mx(2)-1:2,        1,1:nrvars) = &
              a1*    qf(fi(1,2),      1:mx(2)-3:2,        1,1:nrvars) &
            + a2*    qf(fi(1,3),      2:mx(2)-2:2,        1,1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)-1,1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)  ,1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1)+1,1:nrvars))
        !
        qf(fi(1,1),        1,      3:mx(3)-1:2,1:nrvars) = &
              a1*    qf(fi(1,2),        1,      1:mx(3)-3:2,1:nrvars) &
            + a2*    qf(fi(1,3),        1,      2:mx(3)-2:2,1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,1)-1,mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,1)  ,mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,1)+1,mb(3,1)+1:mb(3,2),1:nrvars))
      end if
      !
      if (ipass==0 .or. modulo(parity_1+ipass,2)==1) then
        !
        ! Green:
        qf(fi(1,1),      3:mx(2)-1:2,      2:mx(3)-2:2,1:nrvars) = &
              a1*qf(fi(1,2),      1:mx(2)-3:2,      4:mx(3)  :2,1:nrvars) &
            + a2*qf(fi(1,3),      2:mx(2)-2:2,      3:mx(3)-1:2,1:nrvars) &
            + a3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,1):mb(3,2)-1,1:nrvars)
        !
        qf(fi(1,1),        1,      2:mx(3)-2:2,1:nrvars) = &
              a1*    qf(fi(1,2),        1,      4:mx(3)  :2,1:nrvars) &
            + a2*    qf(fi(1,3),        1,      3:mx(3)-1:2,1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,1)-1,mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,1)  ,mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,1)+1,mb(3,1):mb(3,2)-1,1:nrvars))
        !
        qf(fi(1,1),      3:mx(2)-1:2,    mx(3),1:nrvars) = &
              a1*    qf(fi(1,2),      1:mx(2)-3:2,    mx(3),1:nrvars) &
            + a2*    qf(fi(1,3),      2:mx(2)-2:2,    mx(3),1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,2)+1,1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,2)  ,1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,1)+1:mb(2,2),mb(3,2)-1,1:nrvars))
        !
        ! Red:
        qf(fi(1,1),      2:mx(2)-2:2,      3:mx(3)-1:2,1:nrvars) = &
              a1*qf(fi(1,2),      4:mx(2)  :2,      1:mx(3)-3:2,1:nrvars) &
            + a2*qf(fi(1,3),      3:mx(2)-1:2,      2:mx(3)-2:2,1:nrvars) &
            + a3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)+1:mb(3,2),1:nrvars)
        !
        qf(fi(1,1),      2:mx(2)-2:2,        1,1:nrvars) = &
              a1*    qf(fi(1,2),      4:mx(2)  :2,        1,1:nrvars) &
            + a2*    qf(fi(1,3),      3:mx(2)-1:2,        1,1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)-1,1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)  ,1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,1):mb(2,2)-1,mb(3,1)+1,1:nrvars))
        !
        qf(fi(1,1),    mx(2),      3:mx(3)-1:2,1:nrvars) = &
              a1*    qf(fi(1,2),    mx(2),      1:mx(3)-3:2,1:nrvars) &
            + a2*    qf(fi(1,3),    mx(2),      2:mx(3)-2:2,1:nrvars) &
            + a3*(b1*qc(fi(1,4),mb(2,2)+1,mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b2*qc(fi(1,4),mb(2,2)  ,mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b3*qc(fi(1,4),mb(2,2)-1,mb(3,1)+1:mb(3,2),1:nrvars))
      end if
      !
      ! Corner grey cells:
      do parity_2 = 1, 2
        call GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
        do parity_3 = 1, 2
          call GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
          qf(fi(1,1),fi(2,3),fi(3,3),1:nrvars) = &
                a1*    qf(fi(1,2),fi(2,3),fi(3,3),1:nrvars) &
              + a2*    qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
              + a3*(c1*qc(fi(1,4),fi(2,5),fi(3,5),1:nrvars) &
              +     c2*qc(fi(1,4),fi(2,4),fi(3,5),1:nrvars) &
              +     c3*qc(fi(1,4),fi(2,5),fi(3,4),1:nrvars) &
              +     c4*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars))
        end do
      end do
    end do parity_1_loop
    !
    ! Faces 3 and 4 (coordinate 2 fixed):
    parity_2_loop: do parity_2 = 1, 2
      if (cycle_face(2+parity_2)) cycle parity_2_loop
      call GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
      !
      if (ipass==0 .or. modulo(parity_2+ipass,2)==0) then
        !
        ! Yellow:
        qf(      2:mx(1)-2:2,fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
              a1*qf(      4:mx(1)  :2,fi(2,2),      4:mx(3)  :2,1:nrvars) &
            + a2*qf(      3:mx(1)-1:2,fi(2,3),      3:mx(3)-1:2,1:nrvars) &
            + a3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars)
        !
        qf(      2:mx(1)-2:2,fi(2,1),    mx(3),1:nrvars) = &
              a1*    qf(      4:mx(1)  :2,fi(2,2),    mx(3),1:nrvars) &
            + a2*    qf(      3:mx(1)-1:2,fi(2,3),    mx(3),1:nrvars) &
            + a3*(b1*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,2)+1,1:nrvars) &
            +     b2*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,2)  ,1:nrvars) &
            +     b3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,2)-1,1:nrvars))
        !
        qf(    mx(1),fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
              a1*    qf(    mx(1),fi(2,2),      4:mx(3)  :2,1:nrvars) &
            + a2*    qf(    mx(1),fi(2,3),      3:mx(3)-1:2,1:nrvars) &
            + a3*(b1*qc(mb(1,2)+1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b2*qc(mb(1,2)  ,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b3*qc(mb(1,2)-1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars))
        !
        ! Blue:
        qf(      3:mx(1)-1:2,fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
              a1*qf(      1:mx(1)-3:2,fi(2,2),      1:mx(3)-3:2,1:nrvars) &
            + a2*qf(      2:mx(1)-2:2,fi(2,3),      2:mx(3)-2:2,1:nrvars) &
            + a3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars)
        !
        qf(        1,fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
              a1*    qf(        1,fi(2,2),      1:mx(3)-3:2,1:nrvars) &
            + a2*    qf(        1,fi(2,3),      2:mx(3)-2:2,1:nrvars) &
            + a3*(b1*qc(mb(1,1)-1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b2*qc(mb(1,1)  ,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b3*qc(mb(1,1)+1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars))
        !
        qf(      3:mx(1)-1:2,fi(2,1),        1,1:nrvars) = &
              a1*    qf(      1:mx(1)-3:2,fi(2,2),        1,1:nrvars) &
            + a2*    qf(      2:mx(1)-2:2,fi(2,3),        1,1:nrvars) &
            + a3*(b1*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)-1,1:nrvars) &
            +     b2*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)  ,1:nrvars) &
            +     b3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1)+1,1:nrvars))
      end if
      !
      if (ipass==0 .or. modulo(parity_2+ipass,2)==1) then
        !
        ! Green:
        qf(      3:mx(1)-1:2,fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
              a1*qf(      1:mx(1)-3:2,fi(2,2),      4:mx(3)  :2,1:nrvars) &
            + a2*qf(      2:mx(1)-2:2,fi(2,3),      3:mx(3)-1:2,1:nrvars) &
            + a3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars)
        !
        qf(      2:mx(1)-2:2,fi(2,1),        1,1:nrvars) = &
              a1*    qf(      4:mx(1)  :2,fi(2,2),        1,1:nrvars) &
            + a2*    qf(      3:mx(1)-1:2,fi(2,3),        1,1:nrvars) &
            + a3*(b1*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)-1,1:nrvars) &
            +     b2*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)  ,1:nrvars) &
            +     b3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)+1,1:nrvars))
        !
        qf(    mx(1),fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
              a1*    qf(    mx(1),fi(2,2),      1:mx(3)-3:2,1:nrvars) &
            + a2*    qf(    mx(1),fi(2,3),      2:mx(3)-2:2,1:nrvars) &
            + a3*(b1*qc(mb(1,2)+1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b2*qc(mb(1,2)  ,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars) &
            +     b3*qc(mb(1,2)-1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars))
        !
        ! Red:
        qf(      2:mx(1)-2:2,fi(2,1),      3:mx(3)-1:2,1:nrvars) = &
              a1*qf(      4:mx(1)  :2,fi(2,2),      1:mx(3)-3:2,1:nrvars) &
            + a2*qf(      3:mx(1)-1:2,fi(2,3),      2:mx(3)-2:2,1:nrvars) &
            + a3*qc(mb(1,1):mb(1,2)-1,fi(2,4),mb(3,1)+1:mb(3,2),1:nrvars)
        !
        qf(        1,fi(2,1),      2:mx(3)-2:2,1:nrvars) = &
              a1*    qf(        1,fi(2,2),      4:mx(3)  :2,1:nrvars) &
            + a2*    qf(        1,fi(2,3),      3:mx(3)-1:2,1:nrvars) &
            + a3*(b1*qc(mb(1,1)-1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b2*qc(mb(1,1)  ,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars) &
            +     b3*qc(mb(1,1)+1,fi(2,4),mb(3,1):mb(3,2)-1,1:nrvars))
        !
        qf(      3:mx(1)-1:2,fi(2,1),    mx(3),1:nrvars) = &
              a1*    qf(      1:mx(1)-3:2,fi(2,2),    mx(3),1:nrvars) &
            + a2*    qf(      2:mx(1)-2:2,fi(2,3),    mx(3),1:nrvars) &
            + a3*(b1*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,2)+1,1:nrvars) &
            +     b2*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,2)  ,1:nrvars) &
            +     b3*qc(mb(1,1)+1:mb(1,2),fi(2,4),mb(3,2)-1,1:nrvars))
      end if
      !
      ! Corner grey cells:
      do parity_1 = 1, 2
        call GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
        do parity_3 = 1, 2
          call GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
          qf(fi(1,3),fi(2,1),fi(3,3),1:nrvars) = &
                a1*    qf(fi(1,3),fi(2,2),fi(3,3),1:nrvars) &
              + a2*    qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
              + a3*(c1*qc(fi(1,5),fi(2,4),fi(3,5),1:nrvars) &
              +     c2*qc(fi(1,4),fi(2,4),fi(3,5),1:nrvars) &
              +     c3*qc(fi(1,5),fi(2,4),fi(3,4),1:nrvars) &
              +     c4*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars))
        end do
      end do
    end do parity_2_loop
    !
    ! Faces 5 and 6 (coordinate 3 fixed):
    parity_3_loop: do parity_3 = 1, 2
      if (cycle_face(4+parity_3)) cycle parity_3_loop
      call GetFaceIndex(mx,mb,3,parity_3,fi(3,1:5))
      !
      if (ipass==0 .or. modulo(parity_3+ipass,2)==0) then
        !
        ! Yellow:
        qf(      2:mx(1)-2:2,      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
              a1*qf(      4:mx(1)  :2,      4:mx(2)  :2,fi(3,2),1:nrvars) &
            + a2*qf(      3:mx(1)-1:2,      3:mx(2)-1:2,fi(3,3),1:nrvars) &
            + a3*qc(mb(1,1):mb(1,2)-1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars)
        !
        qf(    mx(1),      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
              a1*    qf(    mx(1),      4:mx(2)  :2,fi(3,2),1:nrvars) &
            + a2*    qf(    mx(1),      3:mx(2)-1:2,fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,2)+1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,2)  ,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,2)-1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars))
        !
        qf(      2:mx(1)-2:2,    mx(2),fi(3,1),1:nrvars) = &
              a1*    qf(      4:mx(1)  :2,    mx(2),fi(3,2),1:nrvars) &
            + a2*    qf(      3:mx(1)-1:2,    mx(2),fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,1):mb(1,2)-1,mb(2,2)+1,fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,1):mb(1,2)-1,mb(2,2)  ,fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,1):mb(1,2)-1,mb(2,2)-1,fi(3,4),1:nrvars))
        !
        ! Blue:
        qf(      3:mx(1)-1:2,      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
              a1*qf(      1:mx(1)-3:2,      1:mx(2)-3:2,fi(3,2),1:nrvars) &
            + a2*qf(      2:mx(1)-2:2,      2:mx(2)-2:2,fi(3,3),1:nrvars) &
            + a3*qc(mb(1,1)+1:mb(1,2),mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars)
        !
        qf(      3:mx(1)-1:2,        1,fi(3,1),1:nrvars) = &
              a1*    qf(      1:mx(1)-3:2,        1,fi(3,2),1:nrvars) &
            + a2*    qf(      2:mx(1)-2:2,        1,fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,1)+1:mb(1,2),mb(2,1)-1,fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,1)+1:mb(1,2),mb(2,1)  ,fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,1)+1:mb(1,2),mb(2,1)+1,fi(3,4),1:nrvars))
        !
        qf(        1,      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
              a1*    qf(        1,      1:mx(2)-3:2,fi(3,2),1:nrvars) &
            + a2*    qf(        1,      2:mx(2)-2:2,fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,1)-1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,1)  ,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,1)+1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars))
      end if
      !
      if (ipass==0 .or. modulo(parity_3+ipass,2)==1) then
        !
        ! Green:
        qf(      3:mx(1)-1:2,      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
              a1*qf(      1:mx(1)-3:2,      4:mx(2)  :2,fi(3,2),1:nrvars) &
            + a2*qf(      2:mx(1)-2:2,      3:mx(2)-1:2,fi(3,3),1:nrvars) &
            + a3*qc(mb(1,1)+1:mb(1,2),mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars)
        !
        qf(        1,      2:mx(2)-2:2,fi(3,1),1:nrvars) = &
              a1*    qf(        1,      4:mx(2)  :2,fi(3,2),1:nrvars) &
            + a2*    qf(        1,      3:mx(2)-1:2,fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,1)-1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,1)  ,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,1)+1,mb(2,1):mb(2,2)-1,fi(3,4),1:nrvars))
        !
        qf(      3:mx(1)-1:2,    mx(2),fi(3,1),1:nrvars) = &
              a1*    qf(      1:mx(1)-3:2,    mx(2),fi(3,2),1:nrvars) &
            + a2*    qf(      2:mx(1)-2:2,    mx(2),fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,1)+1:mb(1,2),mb(2,2)+1,fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,1)+1:mb(1,2),mb(2,2)  ,fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,1)+1:mb(1,2),mb(2,2)-1,fi(3,4),1:nrvars))
        !
        ! Red:
        qf(      2:mx(1)-2:2,      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
              a1*qf(      4:mx(1)  :2,      1:mx(2)-3:2,fi(3,2),1:nrvars) &
            + a2*qf(      3:mx(1)-1:2,      2:mx(2)-2:2,fi(3,3),1:nrvars) &
            + a3*qc(mb(1,1):mb(1,2)-1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars)
        !
        qf(      2:mx(1)-2:2,        1,fi(3,1),1:nrvars) = &
              a1*    qf(      4:mx(1)  :2,        1,fi(3,2),1:nrvars) &
            + a2*    qf(      3:mx(1)-1:2,        1,fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,1):mb(1,2)-1,mb(2,1)-1,fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,1):mb(1,2)-1,mb(2,1)  ,fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,1):mb(1,2)-1,mb(2,1)+1,fi(3,4),1:nrvars))
        !
        qf(    mx(1),      3:mx(2)-1:2,fi(3,1),1:nrvars) = &
              a1*    qf(    mx(1),      1:mx(2)-3:2,fi(3,2),1:nrvars) &
            + a2*    qf(    mx(1),      2:mx(2)-2:2,fi(3,3),1:nrvars) &
            + a3*(b1*qc(mb(1,2)+1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
            +     b2*qc(mb(1,2)  ,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars) &
            +     b3*qc(mb(1,2)-1,mb(2,1)+1:mb(2,2),fi(3,4),1:nrvars))
      end if
      !
      ! Corner grey cells:
      do parity_1 = 1, 2
        call GetFaceIndex(mx,mb,1,parity_1,fi(1,1:5))
        do parity_2 = 1, 2
          call GetFaceIndex(mx,mb,2,parity_2,fi(2,1:5))
          qf(fi(1,3),fi(2,3),fi(3,1),1:nrvars) = &
                a1*    qf(fi(1,3),fi(2,3),fi(3,2),1:nrvars) &
              + a2*    qf(fi(1,3),fi(2,3),fi(3,3),1:nrvars) &
              + a3*(c1*qc(fi(1,5),fi(2,5),fi(3,4),1:nrvars) &
              +     c2*qc(fi(1,4),fi(2,5),fi(3,4),1:nrvars) &
              +     c3*qc(fi(1,5),fi(2,4),fi(3,4),1:nrvars) &
              +     c4*qc(fi(1,4),fi(2,4),fi(3,4),1:nrvars))
        end do
      end do
    end do parity_3_loop
    !
  end subroutine Interpolate3D
  ! ---------------------------------------------------------------------------
end module Boundary
