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
! File:             bsamroutines.f90
! Purpose:          BSAM control module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
module BSAMRoutines
  implicit none
  private
  save
  !
  public :: BSAMSolver
  !
contains
  ! ---------------------------------------------------------------------------
  subroutine BSAMSolver
    use NodeInfoDef
    use TreeOps,         only: AddRootLevelNode, ApplyOnForest, ApplyOnLevel, &
                               CreateBelowSeedLevels, DeleteMarkedNode, &
                               InitForest, KillForest
    use Problem,         only: SetProb, AfterRun
    use Boundary,        only: SetGhost
    use BSAMStorage,     only: DeallocPeriodicBCStorage
    use BSAMInputOutput, only: ReadQ
    implicit none
    type(funcparam)         :: dummy
    character(len=5)        :: zone
    character(len=8)        :: date
    character(len=10)       :: time
    integer                 :: ierror, level, ilevel
    integer, dimension(1:8) :: values
    !
    namelist /rundata/ dt, errortype, getafterstepstats, maxvcycles,          &
                       nsmoothingpasses, omega, outframes, outputuniformmesh, &
                       qerrortol, restart, restartframe, syncelliptic,        &
                       timeiterations, updateauxfreq, eqfn, ldata
    !
    ! Initializations
    errortype         = 1
    dt                = 0.0
    getafterstepstats = .false.
    maxvcycles        = 20
    nsmoothingpasses  = 2
    omega             = 1.0
    outframes         = 1
    outputuniformmesh = .false.
    qerrortol         = 1.0e-06
    restart           = .false.
    restartframe      = 0
    syncelliptic      = .false.
    timeiterations    = 1
    updateauxfreq     = 1
    eqfn              = 6
    ldata             = .true.
    !
    ! Read general input data from rundata.dat
    open(unit=75, file='rundata.dat', &
         status='old', action='read', iostat=ierror)
    if (ierror /= 0) then
      print *, 'Error opening input file rundata.dat. Program stop.'
      stop
    end if
    read(unit=75, nml=rundata)
    close(75)
    !
    ! Write new run date to output.dat
    call date_and_time(date,time,zone,values)
    open(unit=76, file='output.dat', &
         status='unknown', action='write', form='formatted', position='append')
    write(76,1001) date, time
    1001 format(' '/'New run at date ',A8,' and time ',A10/' ')
    write(76, nml=rundata)
    close(76)
    !
    ! By default only one rootlevel grid.  This may change in the future:
    nrootgrids = 1
    !
    ! Initialize a forest of trees.  One root level grid generated in this call:
    call InitForest
    !
    ! Read data file to initialize root level grids:
    call ApplyOnLevel(rootlevel,RootInit,dummy)
    !
    ! Create the levels below the seed needed for multigrid:
    call CreateBelowSeedLevels(minlevel)
    !
    ! Initialize the levels below the forest seed:
    do ilevel = rootlevel-1, minlevel, -1
      call ApplyOnLevel(ilevel,InitSeed,dummy)
    end do
    !
    ! Set user problem parameters:
    call SetProb
    !
    ! Initialize the rootlevel data fields.  In the case of a restart, we read
    ! the fields at all above root levels, saving the data to
    ! uniformgrid(level)%q':
    if (restart) then
      print *, 'Restart of computation from plot frame ', restartframe
      call ReadQ
      outputinitialdata = .false.
      syncelliptic      = .false.
    else
      call ApplyOnLevel(rootlevel,Initialq,dummy)
      outputinitialdata = .true.
    end if
    !
    call SetGhost(rootlevel,0)
    call ApplyOnLevel(rootlevel,SetAuxFields,dummy)
    call ApplyOnLevel(rootlevel,SetSrcFields,dummy)
    call ApplyOnLevel(rootlevel,CopyQToQold,dummy)
    !
    ! Initialization complete. Start run:
    print *, 'BSAM application is running ... '
    print *, ' '
    !
    ! Carry out the time steps:
    call TakeTimeSteps
    !
    ! User-specified actions before program ends:
    call AfterRun
    !
    ! Delete the below-seed-level grids:
    do level = minlevel, maxlevel
      call ApplyOnLevel(level,MarkNodeOld,dummy)
      call ApplyOnLevel(level,ReleaseOldFields,dummy)
      print *,' Level ',level,' fields have been released.'
    end do
    !
    ! Delete the forest of trees:
    call KillForest
    !
    ! Delete forest seed and below-seed levels:
    call ApplyOnLevel(minlevel,MarkNodeToBeDeleted,dummy)
    call ApplyOnLevel(minlevel,DeleteMarkedNode,dummy)
    !
    ! Deallocate storage for periodic boundary conditions
    call DeallocPeriodicBCStorage
    !
  end subroutine BSAMSolver
  ! ---------------------------------------------------------------------------
  integer function MarkNodeToBeDeleted(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    MarkNodeToBeDeleted = err_ok
    !
    info%tobedeleted = .true.
    !
  end function MarkNodeToBeDeleted
  ! ---------------------------------------------------------------------------
  integer function MarkNodeOld(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    MarkNodeOld = err_ok
    !
    info%newgrid = .false.
    info%initialgrid = .false.
    !
  end function MarkNodeOld
  ! ---------------------------------------------------------------------------
  integer function RootInit(rootinfo,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    use Boundary, only: PeriodicSetup
    use BSAMStorage, only: AllocFields
    implicit none
    type(nodeinfo)                    :: rootinfo
    type(funcparam)                   :: dummy
    integer                           :: i, ierror, maux, mbc, nrvars
    integer, dimension(1:maxdims)     :: mx,mxtmp
    integer, dimension(1:2*maxdims)   :: mthbc
    integer, dimension(1:maxdims,1:2) :: mglobal
    real, dimension(1:maxdims)        :: dx, xlower, xupper
    !
    namelist /griddata/ desiredfillratios, errflagopt, ibuffer, maxlevel, &
                        maux, mbc, mglobal, minimumgridpoints, minlevel, &
                        mthbc, mx, ndims, nrvars, qpo, qtolerance, &
                        xlower, xupper
    !
    print *, 'Reading grid data for root-level grid.'
    !
    ! Set default values
    !
    desiredfillratios = 8.0E-01
    errflagopt = errflagdefault
    ibuffer = 0
    maxlevel = 0
    maux = 0
    mbc = 1
    mglobal = 1
    minimumgridpoints = 2
    minlevel = 0
    mthbc = 1
    mx = 1
    ndims = 2
    nrvars = 1
    qpo = 0.0
    qtolerance = 1.0E-06
    xlower = 0.0
    xupper = 1.0
    !
    ! Read from namelist
    !
    open(unit=75,file='griddata.dat',status='old',action='read',iostat=ierror)
    if (ierror/=0) then
      print *,'Error opening input file griddata.dat. Program stop.'
      stop
    end if
    read(75,nml=griddata)
    close(75)
    open(unit=76,file='output.dat', &
         status='old',action='write',form='formatted',position='append')
    write(76,nml=griddata)
    close(76)
    !
    ! Check input variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if (any(desiredfillratios<0.0) .or. any(desiredfillratios>1.0)) then
      print *, 'Error: Code only supports desiredfilratios between 0 and 1.'
      stop
    end if
    !
    if (any(ibuffer<0)) then
      print *, 'Error: Code only supports ibuffers>=0.'
      stop
    end if
    !
    if (maxlevel<0) then
      print *, 'Error: Code only supports maxlevel>=0.'
      stop
    end if
    !
    if (maux<0) then
      print *, 'Error: Code only supports maux>=0.'
      stop
    end if
    !
    if (mbc<1 .or. mbc>2) then
      print *, 'Error: Code only supports mbc=1,2.'
      stop
    end if
    !
    if (minlevel>1) then
      print *, 'Error: Code only supports minlevel<=1.'
      stop
    elseif (minlevel == 1) then
      minlevel=0
      mxtmp=mx
      do
        if (.not. all(mod(mxtmp(1:ndims),2) == 0)) exit
        if (minlevel == 10) exit
        minlevel = minlevel+1
        mxtmp(1:ndims)=mxtmp(1:ndims)/2
      enddo
      minlevel=-minlevel
      print *, 'minlevel chosen automagically: ', minlevel
    end if
    !
    if (ndims<2 .or. ndims>3) then
      print *, 'Error: Code only supports ndims=2,3.'
      stop
    end if
    !
    do i = 1, 2*ndims-1, 2
      if (nrootgrids==1 .and. mthbc(i)==2 .and. mthbc(i+1)/=2) then
        print *, 'Error: Incorrect periodic boundary conditions.'
        stop
      end if
      if (nrootgrids==1 .and. mthbc(i+1)==2 .and. mthbc(i)/=2) then
        print *, 'Error: Incorrect periodic boundary conditions.'
        stop
      end if
    end do
    !
    do i = 1, ndims
      if (mx(i)<=0) then
        print *, 'Error: mx<=0 along dim', i, '.'
        stop
      end if
    end do
    !
    do i = 1, ndims
      if (modulo(mx(i),2)/=0) then
        print *, 'Error: Initial grid must have dimensions which are a multiple'
        print *, '       of the refinement ratio 2.'
        stop
      end if
    end do
    !
    do i = 1, ndims
      if (mglobal(i,2)-mglobal(i,1)+1/=mx(i)) then
        print *, 'Error: mglobal(i,2)-mglobal(i,1)+1/=mx(i), i=', i, '.'
        stop
      end if
    end do
    !
    if (any(minimumgridpoints<1)) then
      print *, 'Error: Code only supports minimumgridpoints>=1.'
      stop
    end if
    !
    if (nrvars>maxnrvars .or. nrvars<1) then
      print *, 'Error: Code only supports nrvars<=maxnrvars and nrvars>=1.'
      stop
    end if
    !
    if (any(qtolerance<=0.0)) then
      print *, 'Error: Code only supports qtolerance>0.0.'
      stop
    end if
    !
    do i = 1, ndims
      if (xupper(i)<=xlower(i)) then
        print *, 'Error: Code only supports xupper > xlower.'
        stop
      end if
    end do
    !
    dx = 0.0
    dx(1:ndims) = (xupper(1:ndims)-xlower(1:ndims))/real(mx(1:ndims))
    !
    do i = 2, ndims
      if (abs(dx(1)-dx(i))>1.0E-10) then
        print *, 'Error: dx(1)\=dx(i) along dim', i, '.'
        stop
      end if
    end do
    !
    mxmax = 1
    mxmax(0,1:ndims) = mx(1:ndims)
    do i = rootlevel+1, maxlevel
      mxmax(i,1:ndims) = 2*mxmax(i-1,1:ndims)
    end do
    !
    ! Copy data to global storage !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    rootinfo%ngrid = 0
    rootinfo%level = rootlevel
    rootinfo%maxlevel = maxlevel
    !
    rootinfo%tobedeleted = .false.
    rootinfo%initialgrid = .true.
    !
    rootinfo%mx = 1
    rootinfo%mx(1:ndims) = mx(1:ndims)
    rootinfo%mglobal = 1
    rootinfo%mglobal(1:ndims,1:2) = mglobal(1:ndims,1:2)
    !
    rootinfo%time = 0.0
    !
    rootinfo%xlower = 0.0
    rootinfo%xlower(1:ndims) = xlower(1:ndims)
    rootinfo%xupper = 0.0
    rootinfo%xupper(1:ndims) = xupper(1:ndims)
    rootinfo%dx = 0.0
    rootinfo%dx(1:ndims) = dx(1:ndims)
    rootinfo%mbc = mbc
    rootinfo%mthbc = 1
    rootinfo%mthbc(1:2*ndims) = mthbc(1:2*ndims)
    !
    rootinfo%nrvars = nrvars
    rootinfo%maux = maux
    rootinfo%nroutvars = rootinfo%nrvars
    !
    rootinfo%mbounds = 1
    rootinfo%mbounds(1:ndims,2) = rootinfo%mx(1:ndims)
    !
    ! Allocate storage !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Set up the periodic offsets used in transferring periodic bc's:
    call PeriodicSetup(rootinfo)
    !
    call AllocFields(rootinfo)
    !
    ! Finished initialization of root node info structure
    RootInit = err_ok
    !
    print *, 'Finished reading root-level grid data.'
    !
  end function RootInit
  ! ---------------------------------------------------------------------------
  integer function InitSeed(seedinfo,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok, GetChildInfo
    use BSAMStorage, only: AllocFields
    implicit none
    type(nodeinfo)          :: seedinfo
    type(funcparam)         :: dummy
    type(nodeinfo), pointer :: child
    integer                 :: ierror, level, mbc, nrvars
    !
    ! Only data needed for multigrid need be copied from the child grid.  Parent
    ! (seedinfo) is born from child:
    !
    ierror = GetChildInfo(child)
    !
    nrvars = child%nrvars
    mbc = child%mbc
    level = child%level-1
    !
    seedinfo%nrvars = nrvars
    seedinfo%mbc = mbc
    seedinfo%maxlevel = maxlevel
    seedinfo%level = level
    !
    print *, 'level=', level, 'ndims=', ndims
    !
    if (any(modulo(child%mx(1:ndims),2)/=0)) then
      print *, 'Error in InitSeed: grid on level', child%level, &
               'will not coarsen.'
      stop
    end if
    !
    seedinfo%mx(1:ndims) = child%mx(1:ndims)/2
    seedinfo%dx(1:ndims) = child%dx(1:ndims)*2.0
    seedinfo%mglobal(1:ndims,1) = 1
    seedinfo%mglobal(1:ndims,2) = child%mglobal(1:ndims,2)/2
    child%mbounds(1:ndims,1) = 1
    child%mbounds(1:ndims,2) = seedinfo%mx(1:ndims)
    !
    seedinfo%mbounds(1:ndims,1) = 1
    seedinfo%mbounds(1:ndims,2) = seedinfo%mx(1:ndims)
    !
    if (any(seedinfo%mx(1:ndims)/=seedinfo%mglobal(1:ndims,2))) then
      print *, 'Error in InitSeed: on level-1', seedinfo%level
      print *, 'mx(:)/=mglobal(:,2)'
      stop
    end if
    !
    seedinfo%ngrid = 0
    !
    seedinfo%tobedeleted = .false.
    seedinfo%level = child%level-1
    seedinfo%initialgrid = .false.
    !
    seedinfo%time = child%time
    !
    seedinfo%xlower = child%xlower
    seedinfo%xupper = child%xupper
    seedinfo%mthbc = child%mthbc
    !
    seedinfo%maux = child%maux
    seedinfo%nroutvars = seedinfo%nrvars
    !
    call AllocFields(seedinfo)
    !
    ! Finished initialization of seed node info structure:
    InitSeed = err_ok
    !
  end function InitSeed
  ! ---------------------------------------------------------------------------
  integer function SetAuxFields(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    use Problem, only: SetAux
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    SetAuxFields = err_ok
    if (info%tobedeleted) return
    !
    if (info%maux>0) call SetAux(info)
    !
  end function SetAuxFields
  ! ---------------------------------------------------------------------------
  integer function SetSrcFields(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    use Problem, only: SetSrc
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    SetSrcFields = err_ok
    if (info%tobedeleted) return
    !
    call SetSrc(info)
    !
  end function SetSrcFields
  ! ---------------------------------------------------------------------------
  integer function Initialq(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    Initialq = err_ok
    if (info%tobedeleted) return
    !
    call QInit(info)
    !
  end function Initialq
  ! ---------------------------------------------------------------------------
  subroutine QInit(info)
    use NodeInfoDef
    use Problem, only: QInit2D, QInit3D
    implicit none
    type(nodeinfo)                :: info
    integer                       :: nrvars
    integer, dimension(1:maxdims) :: mx
    real                          :: h
    real, dimension(1:maxdims)    :: xlower
    !
    nrvars = info%nrvars
    mx(1:ndims) = info%mx(1:ndims)
    h = info%dx(1)
    xlower(1:ndims) = info%xlower(1:ndims)
    !
    select case(ndims)
    case(2)
      !call QInit2D(info%q(1:mx(1),1:mx(2),1        ,1:nrvars), &
                   !mx(1:2),nrvars,h,xlower(1:2))
      call QInit2D(info%q(0:mx(1)+1,0:mx(2)+1,1        ,1:nrvars), &
                   mx(1:2),nrvars,h,xlower(1:2))
    case(3)
      call QInit3D(info%q(1:mx(1),1:mx(2),1:mx(3)+1,1:nrvars), &
                   mx(1:3),nrvars,h,xlower(1:3))
    end select
    !
  end subroutine QInit
  ! ---------------------------------------------------------------------------
  integer function CopyQToQold(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    CopyQToQold = err_ok
    if (info%tobedeleted) return
    !
    info%qold = info%q
    !
  end function CopyQToQold
  ! ---------------------------------------------------------------------------
  subroutine TakeTimeSteps
    !
    ! currenttime is the time of the rootlevel grid.
    ! restarttime is the time at restart.
    ! finaltime is the calculated final time.
    !
    ! With time-subcycling (not currently implemented) a grid whose level is
    ! greater than the rootlevel may exist at a different time than the
    ! rootlevel grid.
    !
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel
    use BSAMInputOutput, only: WriteQ, WriteUniformMeshQ
    use AFASRoutines, only: FillDown, MultigridIterations
    implicit none
    type(funcparam) :: dummy
    logical         :: firstamr
    integer         :: firstframe, it, itperprint, lastframe, level, n
    real            :: starttime, dtsave
    !
    dtsave = dt
    !
    if (restart) then
      if (restartframe>=outframes) then
        print *, 'Error in rundata.dat: outframes=', outframes, &
                 ' restartframe=', restartframe
        print *, 'outframes must be greater than restartframe.'
        stop
      end if
      !
      currenttime = restarttime
      firstframe = restartframe+1
      lastframe = outframes
      print 122, restartframe, outframes
    else
      currenttime = 0.0
      firstframe = 1
      lastframe = outframes
      if (syncelliptic) then
        dt = 0.0
      end if
    end if
    !
    itperprint = timeiterations/outframes
    !
    starttime = currenttime
    finaltime = starttime + real((lastframe-firstframe+1)*itperprint)*dtsave
    !
    firstamr = .true.
    !
    printloop: do n = firstframe, lastframe
      !
      it = 1
      timesteploop: do
        !
        ! For testing purposes only:
        !
        ! Print out grid information for memory checking:
        !    do level = rootlevel+1, maxlevel
        !      print *, 'Grid information, level =', level
        !      gridnumber = 0
        !      call ApplyOnLevel(level,PrintGridInfo,dummy)
        !      read(*,*)
        !    end do
        !
        ! Adapt/create mesh:
        finestlevel = rootlevel
        call AMR(rootlevel)
        !
        ! This ensures that the below rootlevel values of qold are filled after
        ! a clean start or restart:
        if (firstamr) then
          call FillDown(solutionfield)
          firstamr = .false.
        end if
        !
        ! Save a copy of data at the last time step.
        do level = minlevel, finestlevel
          call ApplyOnLevel(level,CopyQToQold,dummy)
        end do
        !
        ! If initial grid, write data:
        if (outputinitialdata .and. (.not. syncelliptic)) then
          !
          outputinitialdata = .false.
          print *, "Writing initial data"
          if (getafterstepstats .and. .false.) then
            integralresult = 0.0
            totalmeshsize = 0
            do level = finestlevel, rootlevel, -1
              call ApplyOnLevel(level,AfterStepStatistics,dummy)
              call ApplyOnLevel(level,GetMeshSize,dummy)
            end do
            !
            ! Write to stats file
            open(unit=65,file='out/stats.dat', &
                 status='unknown',action='write', &
                 form='formatted',position='append')
            write(65,111) currenttime, integralresult(1:2), totalmeshsize
            close(65)
          end if
          totalmeshsize = 0
          do level = finestlevel, rootlevel, -1
            call ApplyOnLevel(level,GetMeshSize,dummy)
          end do
          print *, 'Initial composite mesh size =', totalmeshsize
          print *, 'Finest level =', finestlevel
          call WriteQ(0,currenttime)
          if (outputuniformmesh) call WriteUniformMeshQ(0,currenttime)
          !
        end if
        !
        if (syncelliptic) then
          print *, ''
          print *, ''
          print *, 'Synchronizing elliptic fields at start.'
          print *, ''
        else
          print *, ''
          print *, ''
          print "(' Advancing to time = ', g10.3, ' out of ', g10.3)", &
                currenttime+dt, finaltime
          print *, ''
        end if
        !
        ! For testing purposes only:
        !
        ! Print out grid information for memory checking:
        !    do level = rootlevel+1, maxlevel
        !      print *, 'Grid information, level =', level
        !      gridnumber = 0
        !      call ApplyOnLevel(level,PrintGridInfo,dummy)
        !      read(*,*)
        !    end do
        !
        ! Perform Multigrid on the multilevel mesh:
        do eqtag = 1,eqfn
          write(*,101)
          print 102, eqtag, sigma0
          write(*,101)
          call MultigridIterations
        enddo
        !
        currenttime = starttime+real((n-firstframe)*itperprint+it)*dt
        !
        ! Mark all above-root-level grids old, set current time:
        do level = finestlevel, rootlevel, -1
          call ApplyOnLevel(level,MarkNodeOld,dummy)
          call ApplyOnLevel(level,SetCurrentTime,dummy)
        end do
        !
        ! Calculate after step statistics:
        if (getafterstepstats .and. (.not. syncelliptic)) then
          !
          integralresult = 0.0
          totalmeshsize = 0
          do level = finestlevel, rootlevel, -1
            call ApplyOnLevel(level,AfterStepStatistics,dummy)
            call ApplyOnLevel(level,GetMeshSize,dummy)
          end do
          open(unit=65,file='out/stats.dat',status='unknown',action='write', &
               form='formatted',position='append')
          write(65,111) currenttime, integralresult(1:2), totalmeshsize
          close(65)
        end if
        !
        if (syncelliptic) then
          syncelliptic = .false.
          dt = dtsave
          totalmeshsize = 0
          do level = finestlevel, rootlevel, -1
            call ApplyOnLevel(level,GetMeshSize,dummy)
          end do
          print *, 'Synchronization composite mesh size =', totalmeshsize
        else
          it = it+1
        end if
        !
        if (it>itperprint) exit timesteploop
        !
      end do timesteploop
      !
      totalmeshsize = 0
      do level = finestlevel, rootlevel, -1
        call ApplyOnLevel(level,GetMeshSize,dummy)
      end do
      print *, 'Composite mesh size =', totalmeshsize
      !
      print 121, n, outframes, currenttime
      !
      call WriteQ(n,currenttime)
      if (outputuniformmesh) call WriteUniformMeshQ(n,currenttime)
      !
    end do printloop
    !
    ! Formats for output
    101 format(1x,78("-"))
    102 format(1x,'eqtag = ',i0.2,', sigma0 = ',en12.3)
    111 format(3(f25.12),1x,i8)
    121 format(1x,78("=") /,                     &
               ' Writing frame ',i0.4,' out of ',i0.4,  &
               ' requested frames at t = ',e11.4 /,     &
                                         1x,78("="))
    122 format(1x,78("=") /,                    &
               ' Restart frame ',i0.4,' out of ',i0.4, &
               ' requested frames' /,                  &
               1x,78("="))
  end subroutine TakeTimeSteps
  ! ---------------------------------------------------------------------------
  integer function PrintGridInfo(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    PrintGridInfo = err_ok
    !
    gridnumber = gridnumber+1
    print *, ' '
    print *, 'Grid number =', gridnumber, 'level=', info%level
    print *, 'mx =', info%mx
    print *, 'mb =', info%mbounds
    print *, 'ngrid', info%ngrid
    print *, 'newgrid', info%newgrid
    print *, 'initialgrid', info%initialgrid
    print *, 'tobedeleted', info%tobedeleted
    print *, 'fieldsallocated', info%fieldsallocated
    print *, ' '
    !
  end function PrintGridInfo
  ! ---------------------------------------------------------------------------
  integer function GetMeshSize(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                :: info
    type(funcparam)               :: dummy
    integer, dimension(1:maxdims) :: mx, cmx
    !
    GetMeshSize = err_ok
    if (info%tobedeleted) return
    !
    mx(1:ndims) = info%mx(1:ndims)
    cmx(1:ndims) = mx(1:ndims)/2
    !
    if (info%level==0) then
      totalmeshsize = totalmeshsize+product(mx(1:ndims))
    else
      totalmeshsize = totalmeshsize+product(mx(1:ndims))-product(cmx(1:ndims))
    end if
    !
  end function GetMeshSize
  ! ---------------------------------------------------------------------------
  subroutine ResetGrids
    !
    ! All grids on level>rootlevel are destroyed.
    !
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel, DeleteMarkedNode
    implicit none
    integer         :: level
    type(funcparam) :: dummy
    !
    do level = finestlevel, rootlevel+1, -1
      call ApplyOnLevel(level,MarkNode,dummy)
      call ApplyOnLevel(level,DeleteMarkedNode,dummy)
    end do
    !
  end subroutine ReSetGrids
  ! ---------------------------------------------------------------------------
  integer function MarkNode(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    MarkNode = err_ok
    !
    info%tobedeleted = .true.
    !
  end function MarkNode
  ! ---------------------------------------------------------------------------
  integer function ReleaseOldFields(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    use BSAMStorage, only: DeAllocFields
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    ReleaseOldFields = err_ok
    !
    if (.not. info%newgrid) call DeAllocFields(info)
    !
  end function ReleaseOldFields
  ! ---------------------------------------------------------------------------
  integer function SetCurrentTime(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    !
    SetCurrentTime = err_ok
    !
    info%time = currenttime
    !
  end function SetCurrentTime
  ! ---------------------------------------------------------------------------
  integer function AfterStepStatistics(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok, GetParentInfo
    use Problem, only: AfterStep2D, AfterStep3D
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    type(nodeinfo), pointer           :: parent
    integer                           :: ierror, level, nrvars
    integer, dimension(1:maxdims)     :: mx
    integer, dimension(1:maxdims,1:2) :: mb
    real                              :: h
    real, dimension(1:maxdims)        :: xlower
    !
    AfterStepStatistics = err_ok
    if (info%tobedeleted) return
    !
    mx(1:ndims) = info%mx(1:ndims)
    nrvars = info%nrvars
    h = info%dx(1)
    xlower(1:ndims) = info%xlower(1:ndims)
    level = info%level
    mb = 1
    mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
    !
    ierror = GetParentInfo(parent)
    !
    select case(ndims)
    case(2)
      call AfterStep2D(info%q(0:mx(1)+1,0:mx(2)+1,1,1:nrvars), &
                       parent%q(mb(1,1)-1:mb(1,2)+1, &
                       mb(2,1)-1:mb(2,2)+1,1,1:nrvars), &
                       mx(1:2),nrvars,h,xlower(1:2),level)
    case(3)
      call AfterStep3D(info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars), &
                       parent%q(mb(1,1)-1:mb(1,2)+1, &
                       mb(2,1)-1:mb(2,2)+1, &
                       mb(3,1)-1:mb(3,2)+1,1:nrvars), &
                       mx(1:3),nrvars,h,xlower(1:3),level)
    case default
      print *, 'AfterStepStatistics: Only ndims = 2,3 are supported.'
      stop
    end select
    !
  end function AfterStepStatistics
  ! ---------------------------------------------------------------------------
  recursive subroutine AMR(level)
    use NodeInfoDef
    use TreeOps, only: ExistLevel
    use Boundary, only: SetGhost
    implicit none
    integer, intent(in) :: level
    !
    ! Fill in the ghost points:
    call SetGhost(level,0)
    if (level<maxlevel) then
      !
      ! Tag cells for refinement:
      call EstimateLevelErrors(level)
      !
      ! Create new subgrids:
      call GridAdapt(level)
      !
      ! Recurse to next finer level:
      call AMR(level+1)
      !
    end if
    !
  end subroutine AMR
  ! ---------------------------------------------------------------------------
  subroutine EstimateLevelErrors(level)
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel
    use AFASRoutines, only: RestrictSolution, RelativeTruncationError
    use Boundary, only: SetGhost, GetCoarseGhostPoints
    implicit none
    integer, intent(in) :: level
    type(funcparam)     :: dummy
    !
    ! 1) Calculate the relative truncation error:
    call ApplyOnLevel(level,RestrictSolution,dummy)
    !
    call SetGhost(level-1,0)
    !
    call ApplyOnLevel(level,GetCoarseGhostPoints,dummy)
    !
    call ApplyOnLevel(level,CopyQToQold,dummy)
    !
    call ApplyOnLevel(level-1,CopyQToQold,dummy)
    !
    call ApplyOnLevel(level,RelativeTruncationError,dummy)
    !
    ! 2) Refine based on the size of the relative truncation error.  Make
    !    a linked list of tagged cells:
    allocate(zerothtaggedcell)
    nullify(zerothtaggedcell%prevcell)
    lasttaggedcell => zerothtaggedcell
    ntaggedcells = 0
    call ApplyOnLevel(level,EstimateError,dummy)
    !
    ! 3) Inflate tagged region and destroy the list:
    if (ntaggedcells>0) then
      call InflateEdgeTags(level)
      call DeleteTaggedCellsList
    end if
    nullify(lasttaggedcell)
    nullify(currenttaggedcell)
    deallocate(zerothtaggedcell)
    !
  end subroutine EstimateLevelErrors
  ! ---------------------------------------------------------------------------
  integer function EstimateError(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: info
    type(funcparam) :: dummy
    type(nodeinfo)  :: coarseinfo
    !
    EstimateError = err_ok
    if (info%tobedeleted) return
    !
    call ErrFlag(info,coarseinfo)
    !
  end function EstimateError
  ! ---------------------------------------------------------------------------
  subroutine ErrFlag(info,coarseinfo)
    use NodeInfoDef
    use Problem, only: SetErrFlagsUser2D, SetErrFlagsUser3D
    implicit none
    type(nodeinfo)                    :: info
    type(nodeinfo)                    :: coarseinfo
    integer                           :: level, nrvars
    integer, dimension(1:maxdims)     :: mx, cmx
    integer, dimension(1:maxdims,1:2) :: mglobal
    real                              :: h
    real, dimension(1:maxdims)        :: xlower
    !
    nrvars = info%nrvars
    level = info%level
    h = info%dx(1)
    xlower(1:ndims) = info%xlower(1:ndims)
    mx = 1
    mx(1:ndims) = info%mx(1:ndims)
    cmx = 1
    cmx(1:ndims) = mx(1:ndims)/2
    mglobal = 1
    mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
    !
    select case(errflagopt(level))
    case(errflagdefault)
      !
      ! Default error tagging based on the relative truncation error:
      select case(ndims)
      case(2)
        call SetErrFlags2D(info%qrte(1:cmx(1),1:cmx(2),1,1:nrvars), &
                           info%errorflags(1:mx(1),1:mx(2),1), &
                           mx(1:2),cmx(1:2),nrvars,h,level)
      case(3)
        call SetErrFlags3D(info%qrte(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars), &
                           info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                           mx(1:3),cmx(1:3),nrvars,h,level)
        !
      case default
        print *, 'ErrFlag: only ndims=2d,3d supported'
        stop
      end select
      !
    case(errflaguser)
      !
      ! User chooses the error tagging proceedure:
      select case(ndims)
      case(2)
        call SetErrFlagsUser2D(info%qrte(1:cmx(1)  ,1:cmx(2)  ,1,1:nrvars), &
                               info%q   (0: mx(1)+1,0: mx(2)+1,1,1:nrvars), &
                               info%errorflags(1:mx(1),1:mx(2),1), &
                               mx(1:2),cmx(1:2),nrvars,h,xlower(1:2),level)
      case(3)
        call SetErrFlagsUser3D(info%qrte(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars), &
                               info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nrvars), &
                               info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                               mx(1:3),cmx(1:3),nrvars,h,xlower(1:3),level)
        !
      case default
        print *, 'ErrFlag: only ndims=2d,3d supported'
        stop
      end select
      !
    case default
      print *, 'ErrFlag: No error flagging algorithm selected'
      stop
    end select
    !
    ! Add the tagged cells to a linked list:
    if (ibuffer(level)>0) &
        call BufferAndList(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                           mglobal(1:3,1:2),mx(1:3),level)
    !
  end subroutine ErrFlag
  ! ---------------------------------------------------------------------------
  subroutine SetErrFlags2D(qrte,errorflags,mx,cmx,nrvars,h,level)
    use NodeInfoDef
    implicit none
    real, dimension(1:,1:,1:), intent(in)  :: qrte
    integer, dimension(1:,1:), intent(out) :: errorflags
    integer, dimension(1:2),   intent(in)  :: mx
    integer, dimension(1:2),   intent(in)  :: cmx
    integer,                   intent(in)  :: nrvars
    real,                      intent(in)  :: h
    integer,                   intent(in)  :: level
    integer :: i, j
    real    :: tol
    !
    tol = qtolerance(level)/h/h
    !
    errorflags(1:mx(1),1:mx(2)) = 0
    !
    do i = 1, cmx(1)
      do j = 1, cmx(2)
        if (maxval(abs(qrte(i,j,1:nrvars)))>tol) then
          errorflags(2*i  ,2*j  ) = 1
          errorflags(2*i-1,2*j  ) = 1
          errorflags(2*i  ,2*j-1) = 1
          errorflags(2*i-1,2*j-1) = 1
        end if
      end do
    end do
    !
  end subroutine SetErrFlags2D
  ! ---------------------------------------------------------------------------
  subroutine SetErrFlags3D(qrte,errorflags,mx,cmx,nrvars,h,level)
    use NodeInfoDef
    implicit none
    real, dimension(1:,1:,1:,1:), intent(in)  :: qrte
    integer, dimension(1:,1:,1:), intent(out) :: errorflags
    integer, dimension(1:3),      intent(in)  :: mx
    integer, dimension(1:3),      intent(in)  :: cmx
    integer,                      intent(in)  :: nrvars
    real,                         intent(in)  :: h
    integer,                      intent(in)  :: level
    integer :: i, j, k
    real    :: tol
    !
    tol = qtolerance(level)/h/h/h
    !
    errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0
    !
    do i = 1, cmx(1)
      do j = 1, cmx(2)
        do k = 1, cmx(3)
          if (maxval(abs(qrte(i,j,k,1:nrvars)))>tol) then
            errorflags(2*i  ,2*j  ,2*k  ) = 1
            errorflags(2*i-1,2*j  ,2*k  ) = 1
            errorflags(2*i  ,2*j-1,2*k  ) = 1
            errorflags(2*i-1,2*j-1,2*k  ) = 1
            errorflags(2*i  ,2*j  ,2*k-1) = 1
            errorflags(2*i-1,2*j  ,2*k-1) = 1
            errorflags(2*i  ,2*j-1,2*k-1) = 1
            errorflags(2*i-1,2*j-1,2*k-1) = 1
          end if
        end do
      end do
    end do
    !
  end subroutine SetErrFlags3D
  ! ---------------------------------------------------------------------------
  subroutine BufferAndList(errorflags,mglobal,mx,level)
    use NodeInfoDef
    implicit none
    integer, dimension(1:,1:,1:), intent(in out) :: errorflags
    integer, dimension(1:3,1:2),  intent(in)     :: mglobal
    integer, dimension(1:3),      intent(in)     :: mx
    integer,                      intent(in)     :: level
    integer                                      :: i, j, k
    integer, dimension(1:3)                      :: index
    integer, dimension(1:3,1:2)                  :: mtg
    integer, dimension(1:mx(1),1:mx(2),1:mx(3))  :: errorflagstmp
    !
    errorflagstmp = errorflags
    mtg = 1
    !
    do k = 1, mx(3)
      index(3) = k
      do j = 1, mx(2)
        index(2) = j
        do i = 1, mx(1)
          index(1) = i
          if (errorflagstmp(i,j,k)==1) then
            mtg(1:ndims,1) = index(1:ndims)-ibuffer(level)
            mtg(1:ndims,2) = index(1:ndims)+ibuffer(level)
            if (any(mtg(1:3,1)<1) .or. any(mtg(1:3,2)>mx(1:3))) then
              !
              ! If the buffer area overlaps with any edge of the patch, then
              ! record the global coordinates of the cell:
              ntaggedcells = ntaggedcells+1
              allocate(currenttaggedcell)
              currenttaggedcell%id = ntaggedcells
              currenttaggedcell%coordinate(1:3) = 1
              currenttaggedcell%coordinate(1:ndims) = index(1:ndims) &
                                                      + mglobal(1:ndims,1)-1
              currenttaggedcell%prevcell => lasttaggedcell
              lasttaggedcell => currenttaggedcell
            else
              !
              ! If the buffer area lies within the patch, apply the buffer:
              errorflags(mtg(1,1):mtg(1,2),mtg(2,1):mtg(2,2), &
                         mtg(3,1):mtg(3,2)) = 1
            end if
          end if
        end do
      end do
    end do
    !
  end subroutine BufferAndList
  ! ---------------------------------------------------------------------------
  subroutine InflateEdgeTags(level)
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel
    use Boundary, only: GetPeriodicTagOffset
    implicit none
    integer, intent(in)           :: level
    type(funcparam)               :: dummy
    logical                       :: periodicbuffer
    integer                       :: i
    integer, dimension(1:maxdims) :: coordinate
    !
    currenttaggedcell => lasttaggedcell
    searchloop: do
      if (.not. associated(currenttaggedcell%prevcell)) exit searchloop
      !
      ! Ordinary buffering of an edge tag:
      call ApplyOnLevel(level,BufferTaggedCells,dummy)
      !
      ! Buffering of periodic edge tags. Check to see if the buffer area cuts
      ! across a periodic boundary.  If so, add offset and apply buffer:
      coordinate = currenttaggedcell%coordinate
      if (periodicboundaryconditions) then
        do i = 1, nperiodicoffsets
          currenttaggedcell%coordinate = coordinate
          call GetPeriodicTagOffset(currenttaggedcell%coordinate,level, &
                                    periodicoffsetindex(i),periodicbuffer)
          if (periodicbuffer) call ApplyOnLevel(level,BufferTaggedCells,dummy)
        end do
      end if
      !
      currenttaggedcell => currenttaggedcell%prevcell
    end do searchloop
    !
  end subroutine InflateEdgeTags
  ! ---------------------------------------------------------------------------
  integer function BufferTaggedCells(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    integer                           :: ibuff, level, n
    integer, dimension(1:maxdims)     :: mx
    integer, dimension(1:maxdims,1:2) :: mglobal, mglobaltag, mlocal, moverlap
    !
    BufferTaggedCells = err_ok
    if (info%tobedeleted) return
    !
    level = info%level
    mx = 1
    mx(1:ndims) = info%mx(1:ndims)
    ibuff = ibuffer(level)
    mglobal = 1
    mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
    mglobaltag = 1
    mglobaltag(1:ndims,1) = currenttaggedcell%coordinate(1:ndims)-ibuff
    mglobaltag(1:ndims,2) = currenttaggedcell%coordinate(1:ndims)+ibuff
    !
    ! 1. Find overlap region in global index space:
    moverlap = 1
    moverlap(1:ndims,1) = max(mglobaltag(1:ndims,1),mglobal(1:ndims,1))
    moverlap(1:ndims,2) = min(mglobaltag(1:ndims,2),mglobal(1:ndims,2))
    !
    ! 2. Check for nonempty intersection:
    if (any(moverlap(:,2)-moverlap(:,1)<0)) return
    !
    ! 3. Transform common index space to grid index spaces:
    mlocal = 1
    do n = 1, ndims
      mlocal(n,:) = moverlap(n,:)-mglobal(n,1)+1
    end do
    !
    info%errorflags(mlocal(1,1):mlocal(1,2), &
                    mlocal(2,1):mlocal(2,2), &
                    mlocal(3,1):mlocal(3,2)) = 1
    !
  end function BufferTaggedCells
  ! ---------------------------------------------------------------------------
  subroutine DeleteTaggedCellsList
    use NodeInfoDef
    implicit none
    !
    searchloop: do
      currenttaggedcell => lasttaggedcell%prevcell
      deallocate(lasttaggedcell)
      lasttaggedcell => currenttaggedcell
      if (.not. associated(lasttaggedcell%prevcell)) exit searchloop
    end do searchloop
    !
  end subroutine DeleteTaggedCellsList
  ! ---------------------------------------------------------------------------
  subroutine GridAdapt(level)
    !
    ! Generate new subgrids of level
    !
    use NodeInfoDef
    use TreeOps, only: ApplyOnLevel, ApplyOnLevelPairs, DeleteMarkedNode
    implicit none
    integer, intent(in) :: level
    type(funcparam)     :: dummy
    integer             :: ilevel
    !
    ! Generate new grids in accordance with error flags:
    call ApplyOnLevel(level,RefineGrid,dummy)
    !
    ! Transfer field values from previous grids on this level to the newly
    ! created grids:
    call ApplyOnLevelPairs(level+1,TransferValues,dummy)
    !
    ! Field variables from previous grids on this level are no longer needed.
    ! Deallocate the fields to release memory space and mark nodes for garbage
    ! collection when finest allowed level is reached:
    call ApplyOnLevel(level+1,ReleaseOldFields,dummy)
    !
    ! Check to see if we can invoke garbage collection:
    if (level+1==maxlevel) then
      !
      do ilevel = maxlevel, rootlevel+1, -1
        !
        ! Kill nodes that are marked by tobedeleted in their info structure:
        call ApplyOnLevel(ilevel,DeleteMarkedNode,dummy)
      end do
    end if
    !
  end subroutine GridAdapt
  ! ---------------------------------------------------------------------------
  integer function RefineGrid(info,dummy)
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)                    :: info
    type(funcparam)                   :: dummy
    integer, dimension(1:maxdims,1:2) :: mbounds
    !
    RefineGrid = err_ok
    if (info%tobedeleted) return
    !
    ! Generate new subgrids of the current grid:
    mbounds(1:ndims,1) = 1
    mbounds(1:ndims,2) = info%mx(1:ndims)
    call NewSubGrids(info,mbounds)
    !
  end function RefineGrid
  ! ---------------------------------------------------------------------------
  subroutine NewSubGrids(info,mbounds)
    !
    ! Implementation of Berger-Rigoutsos algorithm
    ! Ref: (IEEE Trans. Systems, Man & Cyber., 21(5):1278-1286, 1991)
    !
    use NodeInfoDef
    implicit none
    type(nodeinfo),                    intent(inout) :: info
    integer, dimension(1:maxdims,1:2), intent(inout) :: mbounds
    logical                                          :: havesplit
    logical, dimension(1:maxsubgrids)                :: cansplit
    integer, parameter                               :: maxsplitpasses = 15
    integer, dimension(1:maxdims)                    :: isplit, mx
    integer, dimension(1:maxdims,1:2,1:maxsubgrids)  :: msubbounds
    integer, dimension(:,:), allocatable             :: signature ,ddsignature
    integer :: del, dist0, dist1, i, i1, i2, ierror, igrid, inflect
    integer :: level, maxm, mgp, minm, n, ngrid, npass
    real    :: fillratio, desfillratio
    !
    mx = info%mx
    level = info%level
    info%nsubgrids = 0
    desfillratio = desiredfillratios(level)
    !
    mgp = minimumgridpoints(level)
    !
    ! Compute fill ratio for this grid:
    fillratio = GridFlagRatio(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                mbounds)
    !
    ! Don't generate a new grid if we don't have flagged points:
    if (fillratio<1.0e-10) return
    !
    ! Allocate space for signatures:
    maxm = maxval(mx(1:ndims))
    allocate(signature(1:maxm,1:ndims),ddsignature(1:maxm,1:ndims),STAT=ierror)
    if (ierror /= 0) then
      print *,'Error allocating signatures arrays in NewSubGrids'
      stop
    end if
    signature=0
    ddsignature=0
    !
    ! Initialize list of subgrids:
    ngrid = 1
    cansplit(:) = .true.
    msubbounds(1:ndims,1:2,ngrid) = mbounds(1:ndims,1:2)
    !
    ! Loop until no better grid splitting can be found:
    igrid = 1
    do while (ngrid<maxsubgrids .and. igrid<=ngrid)
      npass = 0
      do while (cansplit(igrid) .and. npass<maxsplitpasses)
        npass = npass+1
        signature = GetSignatures(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                    msubbounds(:,:,igrid),maxm)
        !
        ! Trim unflagged points on the edges of this grid:
        do n = 1, ndims
          i1 = msubbounds(n,1,igrid)
          i2=msubbounds(n,2,igrid)
          do while(signature(i1,n)==0 .and. i1<msubbounds(n,2,igrid) .and. &
                     i2-i1+1>mgp)
            i1 = i1+1
          end do
          do while(signature(i2,n)==0 .and. i2>msubbounds(n,1,igrid) .and. &
                     i2-i1+1>mgp)
            i2 = i2-1
          end do
          msubbounds(n,1,igrid) = i1
          msubbounds(n,2,igrid) = i2
        end do
        fillratio = GridFlagRatio(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                    msubbounds(:,:,igrid))
        minm = minval(msubbounds(1:ndims,2,igrid)-msubbounds(1:ndims,1,igrid))+1
        if (fillratio<desfillratio) then
          signature = GetSignatures(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                      msubbounds(:,:,igrid),maxm)
          !
          ! Look for holes along which to split grid:
          isplit = 0
          havesplit = .false.
          do n = 1, ndims
            i1 = msubbounds(n,1,igrid)
            i2 = msubbounds(n,2,igrid)
            do i = i1+mgp, i2-mgp
              if (signature(i,n)==0 .and. min(i-i1+1,i2-i)>=mgp) then
                isplit(n) = i
                havesplit = .true.
                exit
              end if
            end do
            if (havesplit) exit
          end do
          if (.not. havesplit) then
            !
            ! No split along a hole. Try split along inflection point:
            do n = 1, ndims
              i1 = msubbounds(n,1,igrid)
              i2 = msubbounds(n,2,igrid)
              do i = i1+1, i2-1
                ddsignature(i,n) = signature(i-1,n)-2*signature(i,n) &
                                   + signature(i+1,n)
              end do
            end do
            inflect = 0
            dist0 = 0
            do n = 1, ndims
              i1 = msubbounds(n,1,igrid)
              i2 = msubbounds(n,2,igrid)
              do i = i1+max(2,mgp), i2-max(1,mgp-1)
                del = abs(ddsignature(i,n)-ddsignature(i-1,n))
                if (del>inflect) then
                  inflect = del
                  isplit = 0
                  isplit(n) = i-1
                  havesplit = .true.
                  dist0 = min(i-i1,i2-i+1)
                else if (del==inflect .and. inflect>0) then
                  dist1 = min(i-i1,i2-i+1)
                  if (dist1>dist0) then
                    isplit = 0
                    isplit(n) = i-1
                    havesplit = .true.
                    dist0 = dist1
                  end if
                end if
              end do
            end do
          end if
          if (havesplit) then
            !
            ! Split the grid along a determined line:
            do n = 1, ndims
              if (isplit(n)>0 .and. min(msubbounds(n,2,igrid)-isplit(n), &
                  isplit(n)-msubbounds(n,1,igrid)+1)>=mgp) then
                !
                ! Add a new subgrid to the end of the grid list:
                ngrid = ngrid+1
                cansplit(ngrid) = .true.
                msubbounds(1:ndims,1:2,ngrid) = msubbounds(1:ndims,1:2,igrid)
                msubbounds(n,1,ngrid) = isplit(n)+1
                !
                ! Replace current grid with a subgrid:
                msubbounds(n,2,igrid) = isplit(n)
                exit
              end if
            end do
          else
            !
            ! Mark grid if no split is possible:
            cansplit(igrid) = .false.
          end if
        else
          cansplit(igrid) = .false.
        end if
      end do
      igrid = igrid+1
    end do
    !
    ! Generate the newly determined subgrids:
    info%nsubgrids = ngrid
    do i = 1, ngrid
      call MakeNewGrid(info,msubbounds(:,:,i))
      if (minval(msubbounds(1:ndims,2,i)-msubbounds(1:ndims,1,i)+1)<mgp) then
        print *, 'Error in NewSubGrids, grid smaller than minimumgridpoints'
        stop
      end if
    end do
    !
    deallocate(signature,ddsignature,STAT=ierror)
    if (ierror/=0) then
      print *,'Error deallocating signatures arrays in NewSubGrids'
      stop
    end if
    !
  end subroutine NewSubGrids
  ! ---------------------------------------------------------------------------
  function GetSignatures(errorflags,msubbounds,maxm) result(gsresult)
    use NodeInfoDef
    implicit none
    integer, dimension(1:,1:,1:),      intent(in) :: errorflags
    integer, dimension(1:maxdims,1:2), intent(in) :: msubbounds
    integer,                           intent(in) :: maxm
    integer, dimension(1:maxm,1:ndims)            :: gsresult
    integer                                       :: i, n
    integer, dimension(1:maxdims)                 :: i1, i2
    !
    i1 = 1
    i2 = 1
    do n = 1, ndims
      i1(1:ndims) = msubbounds(1:ndims,1)
      i2(1:ndims) = msubbounds(1:ndims,2)
      do i = msubbounds(n,1), msubbounds(n,2)
        i1(n) = i
        i2(n) = i
        !
        gsresult(i,n) = sum(errorflags(i1(1):i2(1),i1(2):i2(2),i1(3):i2(3)))
        !
      end do
    end do
    !
  end function GetSignatures
  ! ---------------------------------------------------------------------------
  function GridFlagRatio(errorflags,mbounds) result(gfrresult)
    use NodeInfoDef
    implicit none
    integer, dimension(1:,1:,1:),      intent(in)    :: errorflags
    integer, dimension(1:maxdims,1:2), intent(inout) :: mbounds
    real :: gfrresult
    real :: flagged, total
    !
    total = real(product(mbounds(1:ndims,2)-mbounds(1:ndims,1)+1))
    !
    mbounds(ndims+1:maxdims,1:2) = 1
    !
    flagged = real(sum(errorflags(mbounds(1,1):mbounds(1,2), &
              mbounds(2,1):mbounds(2,2), &
              mbounds(3,1):mbounds(3,2))))
    !
    if (flagged<-1.0e-08) then
      print *, 'Error in GridFlagRatio: flagged < 0.'
      stop
    end if
    !
    gfrresult = flagged/total
    !
  end function GridFlagRatio
  ! ---------------------------------------------------------------------------
  subroutine MakeNewGrid(parent,mbounds)
    !
    ! Generate a new, finer grid within mbounds of the grid info
    !
    use NodeInfoDef
    use TreeOps, only: CreateChild, GetChildInfo
    use BSAMStorage, only: AllocFields
    implicit none
    type(nodeinfo)                                :: parent
    integer, dimension(1:maxdims,1:2), intent(in) :: mbounds
    type(nodeinfo), pointer                       :: child
    integer                                       :: ierror, n, nb
    integer, dimension(1:maxdims,1:2)             :: mglobalbounds
    real                                          :: rand
    !
    ! Create a child of the currentnode:
    call CreateChild
    !
    ! Initialize the grid information for the new node:
    ! Get pointer to youngest child's info:
    ierror = GetChildInfo(child)
    !
    ! Start filling in the fields of child's info:
    child%tobedeleted = .false.
    child%newgrid = .true.
    child%fieldsallocated = .false.
    child%initialgrid = parent%initialgrid
    !
    child%maxlevel = parent%maxlevel
    child%nsubgrids = 0
    child%level = parent%level+1
    if (child%level>finestlevel) finestlevel = child%level
    !
    child%nrvars = parent%nrvars
    !
    child%mbc = parent%mbc
    !
    child%nout = parent%nout
    child%nframe = parent%nframe
    child%outstyle = parent%outstyle
    child%nroutvars = parent%nroutvars
    !
    child%mx = 1
    child%mx(1:ndims) = (mbounds(1:ndims,2)-mbounds(1:ndims,1)+1)*2
    child%maux = parent%maux
    child%mbounds = 1
    child%mbounds(1:ndims,:) = mbounds(1:ndims,:)
    child%mglobal = 1
    mglobalbounds(1:ndims,1) = parent%mglobal(1:ndims,1)+mbounds(1:ndims,1)-1
    mglobalbounds(1:ndims,2) = mglobalbounds(1:ndims,1)+mbounds(1:ndims,2) &
                               - mbounds(1:ndims,1)
    child%mglobal(1:ndims,1) = (mglobalbounds(1:ndims,1)-1)*2+1
    child%mglobal(1:ndims,2) = mglobalbounds(1:ndims,2)*2
    !
    ! First assume all boundaries are internal:
    child%mthbc = internalbc
    !
    ! Now check if we have any physical boundaries:
    do n = 1, ndims
      nb = 2*n-1
      !
      ! If parent boundary condition is physical and child left is same as
      ! parent's left:
      if ((parent%mthbc(nb)<internalbc) .and. &
          (child%mbounds(n,1)==1)) then
        !
        ! Then child left boundary is physical and inherited from parent:
        child%mthbc(nb)=parent%mthbc(nb)
      end if
      nb = nb+1
      !
      ! If parent boundary condition is physical and child right is same as
      ! parent's right:
      if ((parent%mthbc(nb)<internalbc) .and. &
          (child%mbounds(n,2)==parent%mx(n))) then
        !
        ! Then child right boundary is physical and inherited from parent
        child%mthbc(nb) = parent%mthbc(nb)
      end if
    end do
    !
    ! ID the grids randomly:
    call random_number(rand)
    child%ngrid = nint(1000000*rand)
    !
    child%time = parent%time
    !
    child%xlower = 0.0
    child%xupper = 0.0
    child%xlower(1:ndims) = real(mbounds(1:ndims,1)-1) &
                            * parent%dx(1:ndims) + parent%xlower(1:ndims)
    child%xupper(1:ndims) = real(mbounds(1:ndims,2)) &
                            * parent%dx(1:ndims) + parent%xlower(1:ndims)
    child%dx = 0.0
    child%dx(1:ndims) = parent%dx(1:ndims)/2.0
    !
    ! Allocate dynamic space:
    call AllocFields(child,parent)
    !
    ! Initialize field variables with values from parent:
    call InitFields(parent,child)
    !
  end subroutine MakeNewGrid
  ! ---------------------------------------------------------------------------
  subroutine InitFields(parent,child)
    !
    ! Redone to support only bilinear interpolation and 1st layer updating.
    !
    use NodeInfoDef
    use Problem, only: SetAux, SetSrc
    use GridUtilities, only: BiLinProlongationP1MC,  BiLinProlongationP2MC, &
                             TriLinProlongationP1MC, TriLinProlongationP2MC
    implicit none
    type(nodeinfo)                    :: parent, child
    integer                           :: mbc, nrvars
    integer, dimension(1:maxdims)     :: cmx, mx
    integer, dimension(1:maxdims,1:2) :: mb
    !
    nrvars = parent%nrvars
    mbc = parent%mbc
    !
    ! Check whether we can apply user-provided QInit for initialization.
    if (child%initialgrid) then
      call QInit(child)
      child%qold = child%q
      if (child%maux>0) call SetAux(child)
      call SetSrc(child)
      return
    end if
    !
    mx(1:2) = child%mx(1:2)
    cmx(1:2) = mx(1:2)/2
    mb(1:2,1:2) = child%mbounds(1:2,1:2)
    !
    mx = 1
    mx(1:ndims) = child%mx(1:ndims)
    cmx = 1
    cmx(1:ndims) = mx(1:ndims)/2
    mb = 1
    mb(1:ndims,1:2) = child%mbounds(1:ndims,1:2)
    !
    select case(ndims)
    case(2)
      child%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nrvars) &
               = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
                          mb(2,1)-mbc:mb(2,2)+mbc,1,1:nrvars)
      !
      select case(mbc)
      case(1)
        child%q(  0: mx(1)+1, 0: mx(2)+1,1,1:nrvars) &
                = BiLinProlongationP1MC( &
                                  child%qc( 0:cmx(1)+1, 0:cmx(2)+1,1,1:nrvars))
      case(2)
        child%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nrvars) &
                = BiLinProlongationP2MC( &
                                  child%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nrvars))
  case default
        print *, 'InitFields: only supports mbc=1,2.'
      end select
      !
    case(3)
      child%qc(1-mbc:cmx(1)+mbc, &
               1-mbc:cmx(2)+mbc, &
               1-mbc:cmx(3)+mbc,1:nrvars) &
               = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
                          mb(2,1)-mbc:mb(2,2)+mbc, &
                          mb(3,1)-mbc:mb(3,2)+mbc,1:nrvars)
      !
      select case(mbc)
      case(1)
        child%q(  0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nrvars) &
                = TriLinProlongationP1MC( &
                       child%qc( 0:cmx(1)+1, 0:cmx(2)+1, 0:cmx(3)+1,1:nrvars))
      case(2)
        child%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nrvars) &
                = TriLinProlongationP2MC( &
                       child%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nrvars))
      case default
        print *, 'InitFields: only supports mbc=1,2.'
      end select
      !
    case default
      print *, 'InitFields: only supports ndims=2,3.'
    end select
    !
    ! Call user routine to initialize auxilliary array values
    if (child%maux>0) call SetAux(child)
    call SetSrc(child)
    !
    child%qold = child%q
    !
  end subroutine InitFields
  ! ---------------------------------------------------------------------------
  integer function TransferValues(grid1,grid2,dummy)
    !
    ! Transfer values from previous grids on this level to newly created grids
    !
    use NodeInfoDef
    use TreeOps, only: err_ok
    implicit none
    type(nodeinfo)  :: grid1, grid2
    type(funcparam) :: dummy
    !
    TransferValues = err_ok
    !
    if (grid1%newgrid .and. (.not. grid2%newgrid)) then
      !
      ! Look for overlap and transfer grid values from Grid2 to Grid1
      call Transferq(grid2,grid1)
    end if
    !
    if (grid2%newgrid .and. (.not. grid1%newgrid)) then
      !
      ! Look for overlap and transfer grid values from Grid1 to Grid2
      call Transferq(grid1,grid2)
    end if
    !
    ! We either have two new grids or two old grids, nothing to be done:
    return
  end function TransferValues
  ! ---------------------------------------------------------------------------
  subroutine Transferq(sourceinfo,targetinfo)
    !
    ! Look for overlap and transfer grid values from source to target
    !
    use NodeInfoDef
    implicit none
    type(nodeinfo)                    :: sourceinfo
    type(nodeinfo)                    :: targetinfo
    integer, dimension(1:maxdims,1:2) :: mtarget, msource
    !
    if (.not. sourceinfo%fieldsallocated) then
      print *, 'Warning: Internal conistency error, trying to transfer values'
      print *, '         from deallocated sourceinfo in Transferq'
      return
    end if
    !
    if (.not. targetinfo%fieldsallocated) then
      print *, 'Warning: Internal conistency error, trying to transfer values'
      print *, '         to unallocated targetinfo in Transferq.'
      return
    end if
    !
    ! Look for overlap (Later, maybe copy over ghost cells as well):
    msource = sourceinfo%mglobal
    mtarget = targetinfo%mglobal
    !
    call TransferOverlap(msource,mtarget,sourceinfo,targetinfo)
    !
  end subroutine Transferq
  ! ---------------------------------------------------------------------------
  subroutine TransferOverlap(msource,mtarget,sourceinfo,targetinfo)
    use NodeInfoDef
    implicit none
    integer, dimension(1:maxdims,1:2), intent(in) :: msource
    integer, dimension(1:maxdims,1:2), intent(in) :: mtarget
    type(nodeinfo)                                :: sourceinfo
    type(nodeinfo)                                :: targetinfo
    integer                                       :: n
    integer, dimension(1:maxdims,1:2)             :: moverlap, ms, mt
    !
    ! Transfer values from source to target in overlap region:
    ! 1. Find overlap region in global index space:
    moverlap = 1
    moverlap(1:ndims,1) = max(msource(1:ndims,1),mtarget(1:ndims,1))
    moverlap(1:ndims,2) = min(msource(1:ndims,2),mtarget(1:ndims,2))
    !
    ! 2. Check for nonempty intersection:
    if (any(moverlap(:,2)-moverlap(:,1)<0)) return
    !
    ! 3. Transform common index space to grid index spaces:
    ms = 1
    mt = 1
    do n = 1, ndims
      ms(n,:) = moverlap(n,:)-sourceinfo%mglobal(n,1)+1
      mt(n,:) = moverlap(n,:)-targetinfo%mglobal(n,1)+1
    end do
    !
    ! 4. Carry out the transfer:
    targetinfo%q(mt(1,1):mt(1,2), &
                 mt(2,1):mt(2,2), &
                 mt(3,1):mt(3,2),:) &
                 = sourceinfo%q(ms(1,1):ms(1,2), &
                                ms(2,1):ms(2,2), &
                                ms(3,1):ms(3,2),:)
    !
  end subroutine TransferOverlap
  ! ---------------------------------------------------------------------------
end module BSAMRoutines
