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
  logical :: ldata
  logical :: lconverged
  !
  integer, parameter                       :: errflagdefault = 1
  integer, parameter                       :: errflaguser    = 10
  integer, parameter                       :: internalbc     = 999
  integer, parameter                       :: maxsubgrids    = 1024
  integer, parameter                       :: maxdims        = 3
  integer, parameter                       :: maxdepth       = 10
  integer, parameter                       :: maxnrvars      = 25
  integer, parameter                       :: rootlevel      = 0
  integer, parameter                       :: sourcefield    = 1
  integer, parameter                       :: solutionfield  = 2
  integer, parameter                       :: auxiliaryfield = 3
  integer                                  :: errortype
  integer                                  :: finestlevel
  integer                                  :: gridnumber
  integer                                  :: maxvcycles
  integer                                  :: maxlevel
  integer                                  :: minlevel
  integer                                  :: ndims
  integer                                  :: nperiodicoffsets
  integer                                  :: nrootgrids
  integer                                  :: nsmoothingpasses
  integer                                  :: ntaggedcells
  integer                                  :: outframes
  integer                                  :: restartframe
  integer                                  :: timeiterations
  integer                                  :: totalmeshsize
  integer                                  :: updateauxfreq
  integer                                  :: eqtag
  integer                                  :: eqfn
  integer, dimension(:),   allocatable     :: periodicoffsetindex
  integer, dimension(:,:), allocatable     :: poffset
  integer, dimension(0:maxdepth)           :: errflagopt
  integer, dimension(0:maxdepth)           :: ibuffer
  integer, dimension(0:maxdepth)           :: minimumgridpoints
  integer, dimension(0:maxdepth,1:maxdims) :: mxmax
  !
  real                              :: sigma0
  real                              :: currenttime
  real                              :: dt
  real                              :: finaltime
  real                              :: omega
  real                              :: qerrortol
  real                              :: restarttime
  real, dimension(1:2)              :: integralresult
  real, dimension(1:maxnrvars)      :: componentintegral
  real, dimension(1:maxnrvars)      :: qpo
  real, dimension(0:maxdepth)       :: desiredfillratios
  real, dimension(0:maxdepth)       :: qtolerance
  real, dimension(:,:), allocatable :: qoffset
  !
  type :: taggedcell
    integer                       :: id
    integer, dimension(1:maxdims) :: coordinate
    type(taggedcell), pointer     :: prevcell
  end type taggedcell
  !
  type(taggedcell), pointer :: zerothtaggedcell
  type(taggedcell), pointer :: currenttaggedcell
  type(taggedcell), pointer :: lasttaggedcell
  !
  ! Uniform grids used for restarting and output:
  type :: uniformgridtype
    integer, dimension(1:maxdims)     :: mx
    real, dimension(:,:,:,:), pointer :: q
  end type uniformgridtype
  !
  type(uniformgridtype), dimension(0:maxdepth) :: uniformgrid
  !
  type :: nodeinfo
    !
    ! This must be the first component to ensure proper parallel communication
    integer :: nodeinfostart
    !
    ! Logical variables
    ! ------------------------------------------------------------------------
    ! tobedeleted       A necessary component. Mark for garbage collection.
    ! newgrid           Flag to show if grid may accept values from elder
    !                   siblings on this level.
    ! restartgrid       Flag to show if grid has been created during a restart
    !                   from checkpoint file.
    ! fieldsallocated   Flag showing allocation status of fields within this
    !                   node.
    ! initialgrid       Flag to show whether this is an initial grid, i.e.
    !                   created during start-up.
    !
    logical :: tobedeleted
    logical :: newgrid
    logical :: restartgrid
    logical :: fieldsallocated
    logical :: initialgrid
    !
    ! Integers
    ! ------------------------------------------------------------------------
    ! maxlevel      Maximum level to which this grid may be refined
    ! ngrid         Number of this grid
    ! nsubgrids     Number of child grids
    ! level         Level on which this node lives
    ! nrvars        Number of problem field variables
    ! mbc           Number of ghost cells
    ! nout          Number of output frames
    ! nframe        Current frame being output
    ! outstyle      Style of output frames
    ! nroutvars     Number of output variables
    ! maux          Number of auxilliary field variables
    ! mx            Number of grid cells in q along each dimension
    ! mbounds       Index bounds within parent where this child was created
    ! mglobal       Index bounds of this grid in global indexing system
    ! mthbc         Boundary condition codes
    !                 1 - back
    !                 2 - front
    !                 3 - left
    !                 4 - right
    !                 5 - bottom
    !                 6 - top
    ! errorflags    Array of error flags
    !
    integer                            :: maxlevel
    integer                            :: ngrid
    integer                            :: nsubgrids
    integer                            :: level
    integer                            :: nrvars
    integer                            :: mbc
    integer                            :: nout
    integer                            :: nframe
    integer                            :: outstyle
    integer                            :: nroutvars
    integer                            :: maux
    integer, dimension(1:maxdims)      :: mx
    integer, dimension(1:maxdims,1:2)  :: mbounds
    integer, dimension(1:maxdims,1:2)  :: mglobal
    integer, dimension(1:2*maxdims)    :: mthbc
    integer, dimension(:,:,:), pointer :: errorflags
    !
    ! Reals
    ! ------------------------------------------------------------------------
    ! time          The current time at which this grid exists
    ! xlower        Lower coordinates for this grid
    ! xupper        Upper coordinates for this grid
    ! dx            Grid spacings
    ! q             Pointer to field variable arrays
    ! qold          Pointer to field variable arrays at previous time
    ! qc            Pointer to the coarse-under-fine field variable arrays
    ! qrte          Pointer to the coarse-under-fine relative truncation error
    ! f             Pointer to the load function
    ! rf            Pointer to the residual of solution q
    ! ftmp          Pointer to the temporary function.  Same size as f and rf
    ! aux           Pointer to auxilliary arrays
    !
    real                              :: time
    real, dimension(1:maxdims)        :: xlower
    real, dimension(1:maxdims)        :: xupper
    real, dimension(1:maxdims)        :: dx
    real, dimension(:,:,:,:), pointer :: q
    real, dimension(:,:,:,:), pointer :: qold
    real, dimension(:,:,:,:), pointer :: qc
    real, dimension(:,:,:,:), pointer :: qrte
    real, dimension(:,:,:,:), pointer :: f
    real, dimension(:,:,:,:), pointer :: rf
    real, dimension(:,:,:,:), pointer :: ftmp
    real, dimension(:,:,:,:), pointer :: aux
    !
    ! This must be the last component to ensure proper parallel communication
    integer :: nodeinfoend
  end type nodeinfo
  !
  type :: funcparam
    integer                 :: iswitch
    type(nodeinfo), pointer :: info
  end type funcparam
  !
end module NodeInfoDef
