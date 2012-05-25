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
! File:             nodeinfodef.f90
! Purpose:          BSAM node data structures and global allocation module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
module NodeInfoDef
  implicit none
  public
  save
  !
  logical :: getafterstepstats
  logical :: outputinitialdata
  logical :: outputuniformmesh
  logical :: periodicboundaryconditions
  logical :: restart
  logical :: syncelliptic
  !
  ! Double (r8) and extended (rext) precision if available
  integer, parameter :: r8   = selected_real_kind(15,307)
  integer, parameter :: rext = selected_real_kind(15,307)
  !
  integer, parameter :: errflagdefault = 1
  integer, parameter :: errflaguser    = 10
  integer, parameter :: internalbc     = 999
  integer, parameter :: maxsubgrids    = 1024
  integer, parameter :: maxdims        = 3
  integer, parameter :: maxdepth       = 10
  integer, parameter :: maxnrvars      = 25
  integer, parameter :: rootlevel      = 0
  integer, parameter :: sourcefield    = 1
  integer, parameter :: solutionfield  = 2
  integer, parameter :: auxiliaryfield = 3
  integer :: errortype
  integer :: finestlevel
  integer :: gridnumber
  integer :: maxvcycles
  integer :: maxlevel
  integer :: minlevel
  integer :: ndims
  integer :: nperiodicoffsets
  integer :: nrootgrids
  integer :: nsmoothingpasses
  integer :: ntaggedcells
  integer :: outframes
  integer :: restartframe
  integer :: timeiterations
  integer :: totalmeshsize
  integer :: updateauxfreq
  integer :: eqtag
  integer :: eqfn
  real(kind=r8) :: sigma0
  integer, dimension(:),   allocatable     :: periodicoffsetindex
  integer, dimension(:,:), allocatable     :: poffset
  integer, dimension(0:maxdepth)           :: errflagopt
  integer, dimension(0:maxdepth)           :: ibuffer
  integer, dimension(0:maxdepth)           :: minimumgridpoints
  integer, dimension(0:maxdepth,1:maxdims) :: mxmax
  !
  real(kind=r8) :: currenttime
  real(kind=r8) :: dt
  real(kind=r8) :: finaltime
  real(kind=r8) :: omega
  real(kind=r8) :: qerrortol
  real(kind=r8) :: restarttime
  real(kind=r8), dimension(1:2)         :: integralresult
  real(kind=r8), dimension(1:maxnrvars) :: componentintegral
  real(kind=r8), dimension(1:maxnrvars) :: qpo
  real(kind=r8), dimension(0:maxdepth)  :: desiredfillratios
  real(kind=r8), dimension(0:maxdepth)  :: qtolerance
  real(kind=r8), dimension(:,:), allocatable :: qoffset
  !
  ! Yao-li, 2007-06: Add 2 glob variables for GetQAt function in gridutilities.
  real(kind=r8), dimension(1:3),         public :: position_to_get_q
  real(kind=r8), dimension(1:maxnrvars), public :: q_at_position
  !
  type :: taggedcell
    integer :: id
    integer, dimension(1:maxdims) :: coordinate
    type(taggedcell), pointer :: prevcell
  end type taggedcell
  type(taggedcell), pointer :: zerothtaggedcell
  type(taggedcell), pointer :: currenttaggedcell
  type(taggedcell), pointer :: lasttaggedcell
  !
  ! Uniform grids used for restarting and output:
  type :: uniformgridtype
    integer, dimension(1:maxdims) :: mx
    real(kind=r8), dimension(:,:,:,:), pointer :: q
  end type uniformgridtype
  type(uniformgridtype), dimension(0:maxdepth) :: uniformgrid
  !
  type :: nodeinfo
    !
    ! This must be the first component to ensure proper parallel communication
    integer :: nodeinfostart
    !
    ! A necessary component. Mark for garbage collection
    logical :: tobedeleted
    !
    ! Flag to show whether grid may accept values from elder siblings on this
    ! level
    logical :: newgrid
    !
    ! Flag to show whether grid has been created during a restart from
    ! checkpoint file
    logical :: restartgrid
    !
    ! Flag showing allocation status of fields within this node
    logical :: fieldsallocated
    !
    ! Flag to show whether this is an initial grid, i.e. created during start-up
    logical :: initialgrid
    !
    integer :: maxlevel   ! Maximum level to which this grid may be refined
    integer :: ngrid      ! Number of this grid
    integer :: nsubgrids  ! Number of child grids
    integer :: level      ! Level on which this node lives
    integer :: nrvars     ! Number of problem field variables
    integer :: mbc        ! Number of ghost cells
    integer :: nout       ! Number of output frames
    integer :: nframe     ! Current frame being output
    integer :: outstyle   ! Style of output frames
    integer :: nroutvars  ! Number of output variables
    integer :: maux       ! Number of auxilliary field variables
    !
    ! Number of grid cells in q along each dimension
    integer, dimension(1:maxdims) :: mx
    !
    ! Index bounds within parent where this child was created
    integer, dimension(1:maxdims,1:2) :: mbounds
    !
    ! Index bounds of this grid in global indexing system
    integer, dimension(1:maxdims,1:2) :: mglobal
    !
    ! Boundary condition codes
    !   1 - back, 2 - front, 3 - left, 4 - right, 5 - bottom, 6 - top
    integer, dimension(1:2*maxdims) :: mthbc
    !
    ! Array of error flags
    integer, dimension(:,:,:), pointer :: errorflags
    !
    ! Lower coordinates for this grid
    real(kind=r8), dimension(1:maxdims) :: xlower
    !
    ! Upper coordinates for this grid
    real(kind=r8), dimension(1:maxdims) :: xupper
    !
    ! The current time at which this grid exists
    real(kind=r8) :: time
    !
    ! Grid spacings
    real(kind=r8), dimension(1:maxdims) :: dx
    !
    ! Pointer to field variable arrays
    real(kind=r8), dimension(:,:,:,:), pointer :: q
    !
    ! Pointer to field variable arrays at previous time
    real(kind=r8), dimension(:,:,:,:), pointer :: qold
    !
    ! Pointer to the coarse-under-fine field variable arrays
    real(kind=r8), dimension(:,:,:,:), pointer :: qc
    !
    ! Pointer to the coarse-under-fine relative truncation error
    real(kind=r8), dimension(:,:,:,:), pointer :: qrte
    !
    ! Pointer to the load function
    real(kind=r8), dimension(:,:,:,:), pointer :: f
    !
    ! Pointer to the residual of solution q
    real(kind=r8), dimension(:,:,:,:), pointer :: rf
    !
    ! Pointer to the temporary function.  Same size as f and rf
    real(kind=r8), dimension(:,:,:,:), pointer :: ftmp
    !
    ! Pointer to auxilliary arrays
    real(kind=r8), dimension(:,:,:,:), pointer :: aux
    !
    ! This must be the last component to ensure proper parallel communication
    integer :: nodeinfoend
  end type nodeinfo
  !
  type :: funcparam
    integer :: iswitch
    type(nodeinfo), pointer :: info
  end type funcparam
  !
end module NodeInfoDef
