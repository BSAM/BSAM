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
! File:             gridutilities.f90
! Purpose:          BSAM grid utilities module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! -----------------------------------------------------------------------
module GridUtilities
  !
  ! Contains the commonly-used 2 and 3D grid utilities.  These only operate on
  ! uniform, rectangular grid patches.
  !
contains
  ! ---------------------------------------------------------------------------
  function ULap2D(a) result(ulapresult)
    !
    ! 2D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1). The result is stored in
    ! ulapresult(1:mx(1),1:mx(2)).
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,0:), intent(in) :: a
    real(kind=r8), dimension(1:size(a,1)-2,1:size(a,2)-2) :: ulapresult
    integer, dimension(1:2) :: mx
    !
    mx(1) = size(a,1)-2
    mx(2) = size(a,2)-2
    !
    ulapresult(1:mx(1),1:mx(2)) =          a(2:mx(1)+1,1:mx(2)  ) &
                                  +        a(0:mx(1)-1,1:mx(2)  ) &
                                  +        a(1:mx(1)  ,2:mx(2)+1) &
                                  +        a(1:mx(1)  ,0:mx(2)-1) &
                                  - 4.0_r8*a(1:mx(1)  ,1:mx(2)  )
    !
  end function ULap2D
  ! ---------------------------------------------------------------------------
  function UDiv2D(f1,f2) result(udivresult)
    !
    ! UNDIVIDED divergence of the 2D flux function
    !
    !  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).
    !
    ! The result is stored in udivresult(1:mx(1),1:mx(2)).
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,1:), intent(in) :: f1
    real(kind=r8), dimension(1:,0:), intent(in) :: f2
    real(kind=r8), dimension(1:size(f2,1),1:size(f1,2)) :: udivresult
    integer, dimension(1:2) :: mx
    !
    mx(1) = size(f2,1)
    mx(2) = size(f1,2)
    !
    udivresult(1:mx(1),1:mx(2)) =   f1(1:mx(1)  ,1:mx(2)  ) &
                                  - f1(0:mx(1)-1,1:mx(2)  ) &
                                  + f2(1:mx(1)  ,1:mx(2)  ) &
                                  - f2(1:mx(1)  ,0:mx(2)-1)
    !
  end function UDiv2D
  ! ---------------------------------------------------------------------------
  function Restriction2D(a) result(restrictionresult)
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(:,:,:), intent(in) :: a
    real(kind=r8), dimension(1:size(a,1)/2, &
                             1:size(a,2)/2, &
                             1:size(a,3)) :: restrictionresult
    integer :: nrvars
    integer, dimension(1:2) :: cmx, mx
    !
    mx(1)  = size(a,1)
    mx(2)  = size(a,2)
    nrvars = size(a,3)
    cmx(1:2) = mx(1:2)/2
    !
    restrictionresult(1:cmx(1),1:cmx(2),1:nrvars) &
                      =   0.25_r8*(a(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
                        +          a(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
                        +          a(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
                        +          a(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars))
    !
  end function Restriction2D
  ! ---------------------------------------------------------------------------
  function Prolongation2D(a) result(presult)
    !
    ! This prolongation algorithm assumes no ghost layers are present.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(:,:,:), intent(in) :: a
    real(kind=r8), dimension(1:size(a,1)*2, &
                             1:size(a,2)*2, &
                             1:size(a,3)) :: presult
    integer :: nrvars
    integer, dimension(1:2) :: cmx, mx
    !
    cmx(1) = size(a,1)
    cmx(2) = size(a,2)
    nrvars = size(a,3)
    mx(1:2) = cmx(1:2)*2
    !
    presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
    presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars) = a(1:cmx(1),1:cmx(2),1:nrvars)
    !
  end function Prolongation2D
  ! ---------------------------------------------------------------------------
  function BiLinProlongationP1(a) result(presult)
    !
    ! This bilinear prolongation algorithm assumes exactly one ghost layer,
    ! i.e., mbc = 1.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,0:,1:), intent(in) :: a
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0:(size(a,2)-2)*2+1, &
                             1: size(a,3)       ) :: presult
    integer :: nrvars
    integer, dimension(1:2) :: cmx, mx
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0: size(a,2)-2   +1, &
                             1: size(a,3)       ) :: b
    !
    cmx(1) = size(a,1)-2
    cmx(2) = size(a,2)-2
    nrvars = size(a,3)
    mx(1:2) = cmx(1:2)*2
    !
    ! Linear interpolation in the x-direction first:
    b(0: mx(1)  :2,0:cmx(2)+1  ,1:nrvars) &
      =    3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
        +         a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
    b(1: mx(1)+1:2,0:cmx(2)+1  ,1:nrvars) &
      =           a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
        +  3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    presult(0: mx(1)+1  ,0: mx(2)  :2,1:nrvars) &
            =   (3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
              +         b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
    presult(0: mx(1)+1  ,1: mx(2)+1:2,1:nrvars) &
            =   (       b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
              +  3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
    !
  end function BiLinProlongationP1
  ! ---------------------------------------------------------------------------
  function BiLinProlongationP2(a) result(presult)
    !
    ! This bilinear prolongation algorithm assumes exactly two ghost layers,
    ! i.e., mbc = 2.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(-1:,-1:,1:), intent(in) :: a
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1:(size(a,2)-4)*2+2, &
                              1: size(a,3)       ) :: presult
    integer :: nrvars
    integer, dimension(1:2) :: cmx, mx
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1: size(a,2)-4   +2, &
                              1: size(a,3)       ) :: b
    !
    cmx(1) = size(a,1)-4
    cmx(2) = size(a,2)-4
    nrvars = size(a,3)
    mx(1:2) = cmx(1:2)*2
    !
    ! Linear interpolation in the x-direction first:
    b(-1: mx(1)+1:2,-1:cmx(2)+2  ,1:nrvars) &
      =    3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
        +         a(-1:cmx(1)    ,-1:cmx(2)+2  ,1:nrvars)
    !
    b( 0: mx(1)+2:2,-1:cmx(2)+2  ,1:nrvars) &
      =    3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
        +         a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    presult(-1: mx(1)+2  ,-1: mx(2)+1:2,1:nrvars) &
            =   (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
              +         b(-1: mx(1)+2  ,-1:cmx(2)    ,1:nrvars))/16.0_r8
    !
    presult(-1: mx(1)+2  , 0: mx(2)+2:2,1:nrvars) &
            =   (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
              +         b(-1: mx(1)+2  , 1:cmx(2)+2  ,1:nrvars))/16.0_r8
    !
  end function BiLinProlongationP2
  ! ---------------------------------------------------------------------------
  function BiLinProlongationP1MC(a) result(presult)
    !
    ! This mass-corrected bilinear prolongation algorithm assumes exactly one
    ! ghost layer, i.e., mbc = 1.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,0:,1:), intent(in) :: a
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0:(size(a,2)-2)*2+1, &
                             1: size(a,3)       ) :: presult
    integer :: nrvars
    integer, dimension(1:2) :: cmx, mx
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0: size(a,2)-2   +1, &
                             1: size(a,3)       ) :: b
    real(kind=r8), dimension(1: size(a,1)-2     , &
                             1: size(a,2)-2     , &
                             1: size(a,3)       ) :: cor
    !
    cmx(1) = size(a,1)-2
    cmx(2) = size(a,2)-2
    nrvars = size(a,3)
    mx(1:2) = cmx(1:2)*2
    !
    ! Linear interpolation in the x-direction first:
    b(0: mx(1)  :2,0:cmx(2)+1  ,1:nrvars) &
      =    3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
        +         a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
    b(1: mx(1)+1:2,0:cmx(2)+1  ,1:nrvars) &
      =           a(0:cmx(1)    ,0:cmx(2)+1  ,1:nrvars) &
        +  3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    presult(0: mx(1)+1  ,0: mx(2)  :2,1:nrvars) &
            =   (3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
              +         b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
    presult(0: mx(1)+1  ,1: mx(2)+1:2,1:nrvars) &
            =   (       b(0: mx(1)+1  ,0:cmx(2)    ,1:nrvars) &
              +  3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,1:nrvars))/16.0_r8
    !
    ! The mass correction is computed only for regular cells, not for ghost
    ! cells.  Compute mass correction:
    cor(1:cmx(1),1:cmx(2),1:nrvars) &
        = a(1:cmx(1),1:cmx(2),1:nrvars) &
            - 0.25_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
            +          presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
            +          presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
            +          presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars))
    !
    ! Add the mass correction:
    presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
  = presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
  = presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
  = presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars) &
  = presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    !
  end function BiLinProlongationP1MC
  ! ---------------------------------------------------------------------------
  function BiLinProlongationP2MC(a) result(presult)
    !
    ! This mass-corrected bilinear prolongation algorithm assumes exactly two
    ! ghost layers, i.e., mbc = 2.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(-1:,-1:,1:), intent(in) :: a
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1:(size(a,2)-4)*2+2, &
                              1: size(a,3)       ) :: presult
    integer :: nrvars
    integer, dimension(1:2) :: cmx, mx
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1: size(a,2)-4   +2, &
                              1: size(a,3)       ) :: b
    real(kind=r8), dimension( 1: size(a,1)-4     , &
                              1: size(a,2)-4     , &
                              1: size(a,3)       ) :: cor
    !
    cmx(1) = size(a,1)-4
    cmx(2) = size(a,2)-4
    nrvars = size(a,3)
    mx(1:2) = cmx(1:2)*2
    !
    ! Linear interpolation in the x-direction first:
    b(-1: mx(1)+1:2,-1:cmx(2)+2  ,1:nrvars) &
      =    3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
        +         a(-1:cmx(1)    ,-1:cmx(2)+2  ,1:nrvars)
    !
    b( 0: mx(1)+2:2,-1:cmx(2)+2  ,1:nrvars) &
      =    3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:nrvars) &
        +         a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    presult(-1: mx(1)+2  ,-1: mx(2)+1:2,1:nrvars) &
            =   (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
              +         b(-1: mx(1)+2  ,-1:cmx(2)    ,1:nrvars))/16.0_r8
    !
    presult(-1: mx(1)+2  , 0: mx(2)+2:2,1:nrvars) &
            =   (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:nrvars) &
              +         b(-1: mx(1)+2  , 1:cmx(2)+2  ,1:nrvars))/16.0_r8
    !
    ! The mass correction is computed only for regular cells, not for ghost
    ! cells.  Compute mass correction:
    cor(1:cmx(1),1:cmx(2),1:nrvars) &
        = a(1:cmx(1),1:cmx(2),1:nrvars) &
            - 0.25_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
            +          presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
            +          presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
            +          presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars))
    !
    ! Add the mass correction:
    presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars) &
  = presult(2:mx(1)  :2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars) &
  = presult(1:mx(1)-1:2,2:mx(2)  :2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars) &
  = presult(2:mx(1)  :2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars) &
  = presult(1:mx(1)-1:2,1:mx(2)-1:2,1:nrvars)+cor(1:cmx(1),1:cmx(2),1:nrvars)
    !
  end function BiLinProlongationP2MC
  ! ---------------------------------------------------------------------------
  function ULap3D(a) result(ulapresult)
    !
    ! 3D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1). The result is
    ! stored in ulapresult(1:mx(1),1:mx(2)).
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,0:,0:), intent(in) :: a
    real(kind=r8), dimension(1:size(a,1)-2, &
                             1:size(a,2)-2, &
                             1:size(a,3)-2) :: ulapresult
    integer, dimension(1:3) :: mx
    !
    mx(1) = size(a,1)-2
    mx(2) = size(a,2)-2
    mx(3) = size(a,3)-2
    !
    ulapresult(1:mx(1),1:mx(2),1:mx(3)) &
               =          a(2:mx(1)+1,1:mx(2)  ,1:mx(3)  ) &
                 +        a(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ) &
                 +        a(1:mx(1)  ,2:mx(2)+1,1:mx(3)  ) &
                 +        a(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ) &
                 +        a(1:mx(1)  ,1:mx(2)  ,2:mx(3)+1) &
                 +        a(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1) &
                 - 6.0_r8*a(1:mx(1)  ,1:mx(2)  ,1:mx(3)  )
    !
  end function ULap3D
  ! ---------------------------------------------------------------------------
  function UDiv3D(f1,f2,f3) result(udivresult)
    !
    ! UNDIVIDED divergence of the 3D flux function
    !
    !  (f1(0:mx(1),1:mx(2),1:mx(3)),
    !   f2(1:mx(1),0:mx(2),1:mx(3)),
    !   f3(1:mx(1),1:mx(2),0:mx(3))).
    !
    ! The result is stored in udivresult(1:mx(1),1:mx(2),1:mx(3)).
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,1:,1:), intent(in) :: f1
    real(kind=r8), dimension(1:,0:,1:), intent(in) :: f2
    real(kind=r8), dimension(1:,1:,0:), intent(in) :: f3
    real(kind=r8), dimension(1:size(f2,1), &
                             1:size(f3,2), &
                             1:size(f1,3)) :: udivresult
    integer, dimension(1:3) :: mx
    !
    mx(1) = size(f2,1)
    mx(2) = size(f3,2)
    mx(3) = size(f1,3)
    !
    udivresult(1:mx(1),1:mx(2),:mx(3)) =   f1(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
                                         - f1(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ) &
                                         + f2(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
                                         - f2(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ) &
                                         + f3(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
                                         - f3(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1)
    !
  end function UDiv3D
  ! ---------------------------------------------------------------------------
  function Restriction3D(a) result(restrictionresult)
    use NodeInfoDef
    implicit none
    !
    real(kind=r8), dimension(:,:,:,:), intent(in) :: a
    real(kind=r8), dimension(1:size(a,1)/2, &
                             1:size(a,2)/2, &
                             1:size(a,3)/2, &
                             1:size(a,4)) :: restrictionresult
    integer :: nrvars
    integer, dimension(1:3) :: cmx, mx
    !
    mx(1)  = size(a,1)
    mx(2)  = size(a,2)
    mx(3)  = size(a,3)
    nrvars = size(a,4)
    cmx(1:3) = mx(1:3)/2
    !
    restrictionresult(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
       =   0.125_r8*(a(2:mx(1)  :2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
         +           a(1:mx(1)-1:2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
         +           a(2:mx(1)  :2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
         +           a(1:mx(1)-1:2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
         +           a(2:mx(1)  :2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
         +           a(1:mx(1)-1:2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
         +           a(2:mx(1)  :2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars) &
         +           a(1:mx(1)-1:2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars))
    !
  end function Restriction3D
  ! ---------------------------------------------------------------------------
  function Prolongation3D(a) result(presult)
    !
    ! This prolongation algorithm assumes no ghost layers are present.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(:,:,:,:), intent(in) :: a
    real(kind=r8), dimension(1:size(a,1)*2, &
                             1:size(a,2)*2, &
                             1:size(a,3)*2, &
                             1:size(a,4)) :: presult
    integer :: nrvars
    integer, dimension(1:3) :: cmx, mx
    !
    cmx(1) = size(a,1)
    cmx(2) = size(a,2)
    cmx(3) = size(a,3)
    nrvars = size(a,4)
    mx(1:3) = cmx(1:3)*2
    !
    presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
            =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    !
  end function Prolongation3D
  ! ---------------------------------------------------------------------------
  function TriLinProlongationP1(a) result(presult)
    !
    ! This trilinear prolongation algorithm assumes exactly one ghost layer,
    ! i.e., mbc = 1.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,0:,0:,1:), intent(in) :: a
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0:(size(a,2)-2)*2+1, &
                             0:(size(a,3)-2)*2+1, &
                             1: size(a,4)) :: presult
    integer :: nrvars
    integer, dimension(1:3) :: cmx, mx
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0: size(a,2)-2   +1, &
                             0: size(a,3)-2   +1, &
                             1: size(a,4)       ) :: b
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0:(size(a,2)-2)*2+1, &
                             0: size(a,3)-2   +1, &
                             1: size(a,4)       ) :: c
    !
    cmx(1) = size(a,1)-2
    cmx(2) = size(a,2)-2
    cmx(3) = size(a,3)-2
    nrvars = size(a,4)
    mx(1:3) = cmx(1:3)*2
    !
    ! Later, use presult in place of b to save storage space:
    !
    ! Linear interpolation in the x-direction first:
    b(0:mx(1)  :2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
      =   3.0_r8*a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
        +        a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    b(1:mx(1)+1:2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
      =          a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
        + 3.0_r8*a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    c(0:mx(1)+1,0:mx(2)  :2,0:cmx(3)+1,1:nrvars) &
      =   3.0_r8*b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
        +        b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    c(0:mx(1)+1,1:mx(2)+1:2,0:cmx(3)+1,1:nrvars) &
      =          b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
        + 3.0_r8*b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    !
    ! Linear interpolation in the z-direction third:
    presult(0:mx(1)+1,0:mx(2)+1,0:mx(3)  :2,1:nrvars) &
            =   (3.0_r8*c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
              +         c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
    presult(0:mx(1)+1,0:mx(2)+1,1:mx(3)+1:2,1:nrvars) &
            =   (       c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
              +  3.0_r8*c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
    !
  end function TriLinProlongationP1
  ! ---------------------------------------------------------------------------
  function TriLinProlongationP2(a) result(presult)
    !
    ! This trilinear prolongation algorithm assumes exactly two ghost layers,
    ! i.e., mbc = 2.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(-1:,-1:,-1:,1:), intent(in) :: a
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1:(size(a,2)-4)*2+2, &
                             -1:(size(a,3)-4)*2+2, &
                              1: size(a,4)       ) :: presult
    integer :: nrvars
    integer, dimension(1:3) :: cmx, mx
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1: size(a,2)-4   +2, &
                             -1: size(a,3)-4   +2, &
                              1: size(a,4)       ) :: b
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1:(size(a,2)-4)*2+2, &
                             -1: size(a,3)-4   +2, &
                              1: size(a,4)       ) :: c
    !
    cmx(1) = size(a,1)-4
    cmx(2) = size(a,2)-4
    cmx(3) = size(a,3)-4
    nrvars = size(a,4)
    mx(1:3) = cmx(1:3)*2
    !
    ! Linear interpolation in the x-direction first:
    b(-1: mx(1)+1:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
        +        a(-1:cmx(1)    ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
    !
    b( 0: mx(1)+2:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
        +        a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    c(-1: mx(1)+2  ,-1: mx(2)+1:2,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
        +        b(-1: mx(1)+2  ,-1:cmx(2)    ,-1:cmx(3)+2,1:nrvars)
    !
    c(-1: mx(1)+2  , 0: mx(2)+2:2,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
        +        b(-1: mx(1)+2  , 1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
    !
    ! Linear interpolation in the z-direction second:
    presult(-1:mx(1)+2,-1:mx(2)+2,-1: mx(3)+1:2,1:nrvars) &
            =   (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
              +         c(-1:mx(1)+2,-1:mx(2)+2,-1:cmx(3)    ,1:nrvars))/64.0_r8
    presult(-1:mx(1)+2,-1:mx(2)+2, 0: mx(3)+2:2,1:nrvars) &
            =   (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
              +         c(-1:mx(1)+2,-1:mx(2)+2, 1:cmx(3)+2  ,1:nrvars))/64.0_r8
    !
  end function TriLinProlongationP2
  ! ---------------------------------------------------------------------------
  function TriLinProlongationP1MC(a) result(presult)
    !
    ! This trilinear prolongation algorithm assumes exactly one ghost layer,
    ! i.e., mbc = 1.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(0:,0:,0:,1:), intent(in) :: a
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0:(size(a,2)-2)*2+1, &
                             0:(size(a,3)-2)*2+1, &
                             1: size(a,4)) :: presult
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0: size(a,2)-2   +1, &
                             0: size(a,3)-2   +1, &
                             1: size(a,4)       ) :: b
    real(kind=r8), dimension(0:(size(a,1)-2)*2+1, &
                             0:(size(a,2)-2)*2+1, &
                             0: size(a,3)-2   +1, &
                             1: size(a,4)       ) :: c
    real(kind=r8), dimension(1: size(a,1)-2     , &
                             1: size(a,2)-2     , &
                             1: size(a,3)-2     , &
                             1: size(a,4)       ) :: cor
    integer, dimension(1:3) :: cmx, mx
    integer :: nrvars
    !
    cmx(1) = size(a,1)-2
    cmx(2) = size(a,2)-2
    cmx(3) = size(a,3)-2
    nrvars = size(a,4)
    mx(1:3) = cmx(1:3)*2
    !
    ! Later, use presult in place of b to save storage space:
    !
    ! Linear interpolation in the x-direction first:
    b(0:mx(1)  :2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
      =   3.0_r8*a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
        +        a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    b(1:mx(1)+1:2,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
      =          a(0:cmx(1)  ,0:cmx(2)+1,0:cmx(3)+1,1:nrvars) &
        + 3.0_r8*a(1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    c(0:mx(1)+1,0:mx(2)  :2,0:cmx(3)+1,1:nrvars) &
      =   3.0_r8*b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
        +        b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    c(0:mx(1)+1,1:mx(2)+1:2,0:cmx(3)+1,1:nrvars) &
      =          b(0:mx(1)+1,0:cmx(2)  ,0:cmx(3)+1,1:nrvars) &
        + 3.0_r8*b(0:mx(1)+1,1:cmx(2)+1,0:cmx(3)+1,1:nrvars)
    !
    ! Linear interpolation in the z-direction third:
    presult(0:mx(1)+1,0:mx(2)+1,0:mx(3)  :2,1:nrvars) &
            =   (3.0_r8*c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
              +         c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
    presult(0:mx(1)+1,0:mx(2)+1,1:mx(3)+1:2,1:nrvars) &
            =   (       c(0:mx(1)+1,0:mx(2)+1,0:cmx(3)  ,1:nrvars) &
              +  3.0_r8*c(0:mx(1)+1,0:mx(2)+1,1:cmx(3)+1,1:nrvars))/64.0_r8
    !
    ! The mass correction is computed only for regular cells, not for ghost
    ! cells.  Compute mass correction:
    cor(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
        = a(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
            - 0.125_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
            +           presult(1:mx(1)-1:2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
            +           presult(2:mx(1)  :2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
            +           presult(1:mx(1)-1:2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
            +           presult(2:mx(1)  :2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
            +           presult(1:mx(1)-1:2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
            +           presult(2:mx(1)  :2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars) &
            +           presult(1:mx(1)-1:2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars))
    !
    ! Add the mass correction:
    presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
            =   presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
            =   presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
            =   presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
            =   presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
            =   presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
            =   presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
            =   presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
            =   presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    !
  end function TriLinProlongationP1MC
  ! ---------------------------------------------------------------------------
  function TriLinProlongationP2MC(a) result(presult)
    !
    ! This trilinear prolongation algorithm assumes exactly two ghost layers,
    ! i.e., mbc = 2.
    !
    use NodeInfoDef
    implicit none
    real(kind=r8), dimension(-1:,-1:,-1:,1:), intent(in) :: a
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1:(size(a,2)-4)*2+2, &
                             -1:(size(a,3)-4)*2+2, &
                              1: size(a,4)       ) :: presult
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1: size(a,2)-4   +2, &
                             -1: size(a,3)-4   +2, &
                              1: size(a,4)       ) :: b
    real(kind=r8), dimension(-1:(size(a,1)-4)*2+2, &
                             -1:(size(a,2)-4)*2+2, &
                             -1: size(a,3)-4   +2, &
                              1: size(a,4)       ) :: c
    real(kind=r8), dimension( 1: size(a,1)-4     , &
                              1: size(a,2)-4     , &
                              1: size(a,3)-4     , &
                              1: size(a,4)       ) :: cor
    integer, dimension(1:3) :: cmx, mx
    integer :: nrvars
    !
    cmx(1) = size(a,1)-4
    cmx(2) = size(a,2)-4
    cmx(3) = size(a,3)-4
    nrvars = size(a,4)
    mx(1:3) = cmx(1:3)*2
    !
    ! Linear interpolation in the x-direction first:
    b(-1: mx(1)+1:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
        +        a(-1:cmx(1)    ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
    !
    b( 0: mx(1)+2:2,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars) &
        +        a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
    !
    ! Linear interpolation in the y-direction second:
    c(-1: mx(1)+2  ,-1: mx(2)+1:2,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
        +        b(-1: mx(1)+2  ,-1:cmx(2)    ,-1:cmx(3)+2,1:nrvars)
    !
    c(-1: mx(1)+2  , 0: mx(2)+2:2,-1:cmx(3)+2,1:nrvars) &
      =   3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2,1:nrvars) &
        +        b(-1: mx(1)+2  , 1:cmx(2)+2  ,-1:cmx(3)+2,1:nrvars)
    !
    ! Linear interpolation in the z-direction third:
    presult(-1:mx(1)+2,-1:mx(2)+2,-1: mx(3)+1:2,1:nrvars) &
            =   (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
              +         c(-1:mx(1)+2,-1:mx(2)+2,-1:cmx(3)    ,1:nrvars))/64.0_r8
    presult(-1:mx(1)+2,-1:mx(2)+2, 0: mx(3)+2:2,1:nrvars) &
            =   (3.0_r8*c(-1:mx(1)+2,-1:mx(2)+2, 0:cmx(3)+1  ,1:nrvars) &
              +         c(-1:mx(1)+2,-1:mx(2)+2, 1:cmx(3)+2  ,1:nrvars))/64.0_r8
    !
    ! The mass correction is computed only for regular cells, not for ghost
    ! cells.  Compute mass correction:
    cor(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
        = a(1:cmx(1),1:cmx(2),1:cmx(3),1:nrvars) &
            - 0.125_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
            +           presult(1:mx(1)-1:2,2:mx(2)  :2,2:mx(3)  :2,1:nrvars) &
            +           presult(2:mx(1)  :2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
            +           presult(1:mx(1)-1:2,1:mx(2)-1:2,2:mx(3)  :2,1:nrvars) &
            +           presult(2:mx(1)  :2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
            +           presult(1:mx(1)-1:2,2:mx(2)  :2,1:mx(3)-1:2,1:nrvars) &
            +           presult(2:mx(1)  :2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars) &
            +           presult(1:mx(1)-1:2,1:mx(2)-1:2,1:mx(3)-1:2,1:nrvars))
    !
    ! Add the mass correction:
    presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
            =   presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
            =   presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
            =   presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
            =   presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
            =   presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
            =   presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
            =   presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
            =   presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:nrvars) &
              +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:nrvars)
    !
  end function TriLinProlongationP2MC
  ! ---------------------------------------------------------------------------
end module GridUtilities
