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
! File:             bsaminputoutput.f90
! Purpose:          BSAM I/O module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
module BSAMInputOutput
  use NodeInfoDef
  implicit none
  private
  save
  !
  public WriteQ, ReadQ, WriteUniformMeshQ
contains
  ! ---------------------------------------------------------------------------
  subroutine WriteQ(nframe,time)
    use NodeinfoDef
    use TreeOps, only: ApplyOnForest, ApplyOnLevel
    implicit none
    integer, intent(in) :: nframe
    real,    intent(in) :: time
    type(funcparam)     :: dummy
    character(len=16)   :: filename
    integer             :: level
    !
    if (.not. ldata) return
    !
    ! Matlab format
    write(filename,'("./out/m",i5.5,".dat")') nframe
    open(unit=54,file=filename, &
         action='readwrite',status='replace',form='formatted')
    write(54,'(f25.12)') time
    write(54,'(i3)') finestlevel
    do level = rootlevel, finestlevel
      call ApplyOnLevel(level,OutputQ,dummy)
    end do
    close(54)
    !
    ! If wanted: compress file
    !call system('gzip -f '//filename)
  end subroutine WriteQ
  ! ---------------------------------------------------------------------------
  integer function OutputQ(info,dummy)
    use NodeinfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer                       :: nrvars
    integer, dimension(1:maxdims) :: mx
    !
    Outputq = err_ok
    !
    if (.not. info%fieldsallocated) return
    !
    mx = info%mx
    nrvars = info%nrvars
    !
    write(54,*) ' '
    write(54,'(I3,3(1x,I3))') info%level, ndims, 2, nrvars
    !
    select case(ndims)
    case(2)
      write(54,2001) info%dx(1), info%dx(2)
      write(54,2001) info%xlower(1), info%xlower(2)
      write(54,2001) info%xupper(1), info%xupper(2)
      write(54,2011) info%mx(1), info%mx(2)
      write(54,2012) info%mglobal(1,1), info%mglobal(1,2), &
                     info%mglobal(2,1), info%mglobal(2,2)
      call WriteQ2D(info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars),mx(1:2), &
                    nrvars)
    case(3)
      write(54,2002) info%dx(1), info%dx(2), info%dx(3)
      write(54,2002) info%xlower(1), info%xlower(2), info%xlower(3)
      write(54,2002) info%xupper(1), info%xupper(2), info%xupper(3)
      write(54,2021) info%mx(1), info%mx(2), info%mx(3)
      write(54,2022) info%mglobal(1,1), info%mglobal(1,2), &
                     info%mglobal(2,1), info%mglobal(2,2), &
                     info%mglobal(3,1), info%mglobal(3,2)
      call WriteQ3D(info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars),mx(1:3), &
                    nrvars)
    end select
    !
    ! Formats
    2001 format(f25.12,1x,f25.12)
    2011 format(i8,1x,i8)
    2012 format(i8,3(1x,i8))
    2002 format(f25.12,2(1x,f25.12))
    2021 format(i8,2(1x,i8))
    2022 format(i8,5(1x,i8))
  end function OutputQ
  ! ---------------------------------------------------------------------------
  subroutine WriteQ2D(q,mx,nrvars)
    use NodeInfoDef
    implicit none
    real, dimension(0:,0:,1:),  intent(in) :: q
    integer, dimension(1:2),    intent(in) :: mx
    integer,                    intent(in) :: nrvars
    integer                                :: i, j, k
    real                                   :: rdummy
    !
    do j=0, mx(2)+1
      do i=0, mx(1)+1
        write(54,3001) (q(i,j,k),k=1,nrvars)
        backspace(54)
        read(54,3002) rdummy
      end do
    end do
    !
    ! Formats
    3001 format(10(es20.12,2x))
    3002 format(es20.12)
  end subroutine WriteQ2D
  ! ---------------------------------------------------------------------------
  subroutine WriteQ3D(q,mx,nrvars)
    use NodeInfoDef
    implicit none
    real, dimension(0:,0:,0:,1:),  intent(in) :: q
    integer, dimension(1:3),       intent(in) :: mx
    integer,                       intent(in) :: nrvars
    integer                                   :: i, j, k, l
    real                                      :: rdummy
    !
    do k=0, mx(3)+1
      do j=0, mx(2)+1
        do i=0, mx(1)+1
          write(54,3001) (q(i,j,k,l),l=1,nrvars)
          backspace(54)
          read(54,3002) rdummy
        end do
      end do
    end do
    !
    ! Formats
    3001 format(10(f25.12,2x))
    3002 format(f25.12)
  end subroutine WriteQ3D
  ! ---------------------------------------------------------------------------
  subroutine ReadQ
    use NodeInfoDef
    use BSAMStorage, only: AllocFields
    use TreeOps,     only: CreateChild, CurrentNodeToYoungest, GetChildInfo, &
                           GetRootInfo
    implicit none
    type(nodeinfo), pointer           :: rootinfo, info
    logical                           :: fileexist
    character(len=1)                  :: string
    character(len=16)                 :: filename
    integer                           :: ierror, ioerror, level, lvl
    integer                           :: mbc, restartndims
    integer                           :: npatch, nrvars, restartfinestlevel, r
    integer, dimension(1:maxdims)     :: mx
    integer, dimension(1:maxdims,1:2) :: mg
    real, dimension(1:maxdims)        :: dx, xlower, xupper
    real, parameter                   :: small = 1.0e-08
    !
    ! On read-in the usual parent-child relationship is broken.  For simplicity,
    ! the youngest grid on level=l-1 is the parent of all level=l grids.  In
    ! particular, coarse-fine boundary conditions can not be properly enforced.
    ! The first ghost layer values are intact on read-in.  However, even these
    ! are not needed at restart.
    !
    ! Initialized values for the root grid:
    ierror = GetRootInfo(rootinfo)
    info => rootinfo
    !
    mx = 1
    mx(1:ndims) = rootinfo%mx(1:ndims)
    mbc = rootinfo%mbc
    !
    write(filename,'("./out/m",i5.5,".dat")') restartframe
    !
    ! If compressed, file must be decompressed
    !call system('gunzip -f '//filename//'.gz')
    !
    inquire(file=filename,exist=fileexist)
    if (.not. fileexist) then
      print *, 'Readq: Error input file does not exist: ', filename
      stop
    end if
    open(unit=54,file=filename, &
         status='old',form='formatted',action='read',iostat=ioerror)
    if (ioerror/=0) then
      print *,'Readq: Error opening restart file ', filename
      stop
    end if
    read(54,'(F25.12)') restarttime
    !
    info%time = restarttime
    !
    read(54,'(I3)') restartfinestlevel
    !
    finestlevel = restartfinestlevel
    !
    if (restartfinestlevel>maxlevel) then
      print *, 'Readq: Error restartfinestlevel>maxlevel.'
      stop
    end if
    !
    npatch = 0
    !
    level_loop: do lvl = rootlevel, restartfinestlevel
      !
      if (lvl>rootlevel) call CurrentNodeToYoungest(lvl-1)
      !
      patch_loop: do
        !
        read(54,'(A1)',iostat=ioerror) string
        if (ioerror<0) exit level_loop
        read(54,'(I3,3(1x,I3))') level, restartndims, r, nrvars
        !
        ! Check to see whether we've read all the patches on this level:
        if (level/=lvl) then
          level = lvl
          backspace(54)
          backspace(54)
          exit patch_loop
        end if
        !
        if (level>rootlevel) then
          call CreateChild
          nullify(info)
          ierror = GetChildInfo(info)
        end if
        !
        if (restartndims/=ndims) then
          print *, 'Readq: Error reading data file ', filename
          print *, 'Spatial dimension (ndims) inconsistancy on restart.'
          stop
        end if
        !
        if (r/=2) then
          print *, 'Readq: Error reading data file ', filename
          print *, 'Refinement ratio (r) inconsistancy on restart.'
          print *, 'Currently only r = 2 is supported for restarting.'
          stop
        end if
        !
        select case(ndims)
        case(2)
          read(54,2001) dx(1), dx(2)
          read(54,2001) xlower(1), xlower(2)
          read(54,2001) xupper(1), xupper(2)
          2001 format(F25.12,1x,F25.12)
          read(54,'(I8,1X,I8)') mx(1), mx(2)
          read(54,'(I8,3(1X,I8))') mg(1,1), mg(1,2), mg(2,1), mg(2,2)
        case(3)
          read(54,2002) dx(1), dx(2), dx(3)
          read(54,2002) xlower(1), xlower(2), xlower(3)
          read(54,2002) xupper(1), xupper(2), xupper(3)
          2002 format(F25.12,2(1x,F25.12))
          read(54,'(I8,2(1X,I8))') mx(1), mx(2), mx(3)
          read(54,'(I8,5(1X,I8))') mg(1,1), mg(1,2), mg(2,1), mg(2,2), &
               mg(3,1), mg(3,2)
        end select
        !
        ! Root-level grid constructs have already been set.  Check for errors:
        if (level==rootlevel) then
          !
          if (nrvars/=rootinfo%nrvars) then
            print *, 'Readq: Error reading data file ', filename
            print *, 'Number of variables (nrvars) inconsistancy on restart.'
            stop
          end if
          !
          if (any(dx(1:ndims)-rootinfo%dx(1:ndims)>small)) then
            print *, 'Readq: Error reading data file ', filename
            print *, 'Root-level grid spacing (dx) inconsistancy on restart.'
            stop
          end if
          !
          if (any(abs(xlower(1:ndims)-rootinfo%xlower(1:ndims))>small)) then
            print *, 'Readq: Error reading data file ', filename
            print *, 'Root-level grid location (xlower) inconsistancy on restart.'
            stop
          end if
          !
          if (any(abs(xupper(1:ndims)-rootinfo%xupper(1:ndims))>small)) then
            print *, 'Readq: Error reading data file ', filename
            print *, 'Root-level grid location (xupper) inconsistancy on restart.'
            stop
          end if
          !
          if (any(mx(1:ndims)/=rootinfo%mx(1:ndims))) then
            print *, 'Readq: Error reading data file ', filename
            print *, 'Root-level grid size (mx) inconsistancy on restart.'
            stop
          end if
          !
          if (any(mg(1:ndims,1)/=1) .or. &
              any(mg(1:ndims,2)/=rootinfo%mx(1:ndims))) then
            print *, 'Readq: Error reading data file ', filename
            print *, 'Root-level grid size (mglobal) inconsistancy on restart.'
            stop
          end if
        end if
        !
        info%tobedeleted = .false.
        info%newgrid     = .false.
        info%initialgrid = .false.
        !
        ! These rootlevel constructs should not be changed:
        if (level>rootlevel) then
          info%fieldsallocated = .false.
          info%maxlevel = maxlevel
          info%nsubgrids = 0
          info%level = level
          !
          info%nrvars = nrvars
          !
          info%mbc = rootinfo%mbc
          !
          info%nout = rootinfo%nout
          info%nframe = rootinfo%nframe
          info%outstyle = rootinfo%outstyle
          info%nroutvars = rootinfo%nroutvars
          !
          info%mx = 1
          info%mx(1:ndims) = mx(1:ndims)
          info%maux = rootinfo%maux
          info%mglobal = 1
          info%mglobal(1:ndims,1) = mg(1:ndims,1)
          info%mglobal(1:ndims,2) = mg(1:ndims,2)
          !
          ! These two items are broken on restart, but are not needed:
          info%mthbc = internalbc
          info%mbounds = 1
          !
          npatch = npatch+1
          !      info%ngrid = npatch
          info%ngrid = -13
          !
          info%time = restarttime
          !
          info%xlower = 0.0
          info%xupper = 0.0
          info%xlower(1:ndims) = xlower(1:ndims)
          info%xupper(1:ndims) = xupper(1:ndims)
          info%dx = 0.0
          info%dx(1:ndims) = dx(1:ndims)
          !
          call AllocFields(info)
        end if
        !
        ! Read-in field data:
        select case(ndims)
        case(2)
          call ReadQ2D(info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars),mx(1:2), &
                       nrvars)
        case(3)
          call ReadQ3D(info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars),mx(1:3), &
                       nrvars)
        end select
        !
      end do patch_loop
    end do level_loop
    !
    close(54)
    !
    ! If desired, file can be compressed
    !call system('gzip -f '//filename)
    !
  end subroutine ReadQ
  ! ---------------------------------------------------------------------------
  subroutine ReadQ2D(q,mx,nrvars)
    use NodeInfoDef
    implicit none
    real,    dimension(0:,0:,1:),  intent(out) :: q
    integer, dimension(1:2),       intent(in)  :: mx
    integer,                       intent(in)  :: nrvars
    integer                                    :: i, j, k
    !
    do j=0, mx(2)+1
      do i=0, mx(1)+1
        read(54,3001) (q(i,j,k),k=1,nrvars)
      end do
    end do
    !
    ! Format
    3001 format(10(f25.12,2x))
  end subroutine ReadQ2D
  ! ---------------------------------------------------------------------------
  subroutine ReadQ3D(q,mx,nrvars)
    use NodeInfoDef
    implicit none
    real,    dimension(0:,0:,0:,1:),  intent(out) :: q
    integer, dimension(1:3),          intent(in)  :: mx
    integer,                          intent(in)  :: nrvars
    integer                                       :: i, j, k, l
    !
    do k=0, mx(3)+1
      do j=0, mx(2)+1
        do i=0, mx(1)+1
          read(54,3001) (q(i,j,k,l),l=1,nrvars)
        end do
      end do
    end do
    !
    ! Format
    3001 format(10(F25.12,2X))
  end subroutine ReadQ3D
  ! ---------------------------------------------------------------------------
  subroutine WriteUniformMeshQ(nframe,time)
    use NodeInfoDef
    use TreeOps,       only: ApplyOnLevel, GetRootInfo
    use BSAMStorage,   only: AllocUniformGrids, DeallocUniformGrids
    use GridUtilities, only: BiLinProlongationP1MC, TriLinProlongationP1MC
    implicit none
    integer, intent(in)           :: nframe
    real,    intent(in)           :: time
    type(nodeinfo), pointer       :: rootinfo
    type(funcparam)               :: dummy
    character(len=16)             :: filename
    integer, dimension(1:maxdims) :: high, low, mx, mxuc, mxuf
    real, dimension(1:maxdims)    :: dx, xlower, xupper
    integer                       :: ierror, level, mbc, nrvars
    !
    ierror          = GetRootInfo(rootinfo)
    nrvars          = rootinfo%nrvars
    mbc             = rootinfo%mbc
    dx              = 0.0
    xlower          = 0.0
    xupper          = 0.0
    mx              = 1
    low             = 1
    high            = 1
    dx(1:ndims)     = rootinfo%dx(1:ndims)
    xlower(1:ndims) = rootinfo%xlower(1:ndims)
    xupper(1:ndims) = rootinfo%xupper(1:ndims)
    mx(1:ndims)     = rootinfo%mx(1:ndims)
    low(1:ndims)    = 0
    high(1:ndims)   = mx(1:ndims)+1
    !
    call AllocUniformGrids(mx,mbc,nrvars)
    !
    uniformgrid(rootlevel)%q(low(1):high(1), &
                             low(2):high(2), &
                             low(3):high(3),1:nrvars)     &
                             = rootinfo%q(low(1):high(1), &
                                          low(2):high(2), &
                                          low(3):high(3),1:nrvars)
    !
    mxuf = mx
    !
    do level = 1, finestlevel
      !
      mxuc = mxuf
      mxuf = 1
      mxuf(1:ndims) = 2*mxuc(1:ndims)
      dx(1:ndims) = dx(1:ndims)/2
      !
      select case(ndims)
      case(2)
        uniformgrid(level  )%q(0:mxuf(1)+1,0:mxuf(2)+1,1,1:nrvars) &
              = BiLinProlongationP1MC(uniformgrid(level-1)%q(0:mxuc(1)+1, &
                                      0:mxuc(2)+1,1,1:nrvars))
      case(3)
        uniformgrid(level  )%q(0:mxuf(1)+1,0:mxuf(2)+1,0:mxuf(3)+1,1:nrvars) &
              = TriLinProlongationP1MC(uniformgrid(level-1)%q(0:mxuc(1)+1, &
                                       0:mxuc(2)+1,0:mxuc(3)+1,1:nrvars))
      end select
      !
      call ApplyOnLevel(level,CopyPatchToUniformGrid,dummy)
    end do
    !
    write(filename,'("./out/u",i5.5,".dat")') restartframe
    open(unit=54,file=filename,status='replace',form='formatted')
    !
    write(54,'(f25.12)') time
    write(54,'(i3)') 0
    write(54,*) ' '
    write(54,'(i3,3(1x,i3))') 0, ndims, 2, nrvars
    !
    select case(ndims)
    case(2)
      write(54,2001) dx(1), dx(2)
      write(54,2001) xlower(1), xlower(2)
      write(54,2001) xupper(1), xupper(2)
      write(54,2010) mxuf(1), mxuf(2)
      write(54,2011) 1, mxuf(1), 1, mxuf(2)
      call WriteQ2D(uniformgrid(finestlevel)%q(0:mxuf(1)+1, &
                    0:mxuf(2)+1, &
                    1          , &
                    1:nrvars    ),mxuf(1:2),nrvars)
    case(3)
      write(54,2002) dx(1), dx(2), dx(3)
      write(54,2002) xlower(1), xlower(2), xlower(3)
      write(54,2002) xupper(1), xupper(2), xupper(3)
      write(54,2020) mxuf(1), mxuf(2), mxuf(3)
      write(54,2021) 1, mxuf(1), 1, mxuf(2), 1, mxuf(3)
      call WriteQ3D(uniformgrid(finestlevel)%q(0:mxuf(1)+1, &
                    0:mxuf(2)+1, &
                    0:mxuf(3)+1, &
                    1:nrvars    ),mxuf(1:3),nrvars)
    end select
    !
    call DeallocUniformGrids
    !
    ! Formats
    2001 format(f25.12,1x,f25.12)
    2002 format(f25.12,2(1x,f25.12))
    2010 format(i8,1x,i8)
    2011 format(i8,3(1x,i8))
    2020 format(i8,2(1x,i8))
    2021 format(i8,5(1x,i8))
  end subroutine WriteUniformMeshQ
  ! ---------------------------------------------------------------------------
  integer function CopyPatchToUniformGrid(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    integer, dimension(1:maxdims)     :: h1, h2, l1, l2, mx
    integer, dimension(1:maxdims,1:2) :: mg
    integer                           :: level, nrvars
    !
    CopyPatchToUniformGrid = err_ok
    !
    if (info%tobedeleted) return
    !
    nrvars = info%nrvars
    level = info%level
    !
    mg = 1
    mg(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
    mx = 1
    mx(1:ndims) = info%mx(1:ndims)
    !
    l1 = 1
    l1(1:ndims) = mg(1:ndims,1)-1
    h1 = 1
    h1(1:ndims) = mg(1:ndims,2)+1
    l2 = 1
    l2(1:ndims) = 0
    h2 = 1
    h2(1:ndims) = mx(1:ndims)+1
    !
    uniformgrid(level)%q(l1(1):h1(1),l1(2):h1(2),l1(3):h1(3),1:nrvars) &
                = info%q(l2(1):h2(1),l2(2):h2(2),l2(3):h2(3),1:nrvars)
    !
  end function CopyPatchToUniformGrid
  ! ---------------------------------------------------------------------------
end module BSAMInputOutput
