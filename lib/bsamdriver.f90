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
! File:             bsamdriver.f90
! Purpose:          BSAM driver module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
program BSAMDriver
  use BSAMRoutines, only: BSAMSolver
  implicit none
  real :: t0, t1
  !
  call cpu_time(t0)
  call BSAMSolver
  call cpu_time(t1)
  !
  print '("Total execution time=", f10.2, " seconds")', t1-t0
  !
end program BSAMDriver
