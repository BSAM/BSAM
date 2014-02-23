module Problem
  !
  ! This module contains routines for initializing the problem and for driving
  ! the multigrid routines in the module AFASRoutines.  In particular, the
  ! routines in AFASRoutines use RelaxGrid#D, Operator#D, Source#D,
  ! SourceUpdate#D, and UpdateAuxVcycle#D, where # = 2 and 3.
  !
  ! Templates for the Cahn-Hilliard Equation.  We use Eyre's gradient stable
  ! scheme for time stepping.
  !
contains
  subroutine SetProb
    use NodeInfoDef
    use ProblemDef
    implicit none
    integer :: ierror
    namelist /problemdata/ alpha, eps
    !
    ! Open file for reading
    open(unit=75,file='problemdata.dat', &
         status='old',action='read',iostat=ierror)
    if (ierror/=0) then
      print *,'Error opening input file problemdata.dat. Program stop.'
      stop
    end if
    read(75,nml=problemdata)
    close(75)
    !
    ! Write data to output file
    open(unit=76,file='output.dat', &
         status='old',action='write',form='formatted', position='append')
    write(76,nml=problemdata)
    close(76)
    !
    eps2 = eps*eps
    !
  end subroutine SetProb
  subroutine AfterRun
    use NodeInfoDef
    use ProblemDef
    implicit none
  end subroutine AfterRun
  subroutine SetAux(info)
    use NodeInfoDef
    use ProblemDef
    implicit none
    type (nodeinfo) :: info
  end subroutine SetAux
  subroutine SetSrc(info)
    use NodeInfoDef
    use ProblemDef
    implicit none
    type (nodeinfo) :: info
  end subroutine SetSrc
  !
  ! 2D Multigrid Routines
  !
  subroutine QInit2D(q,mx,nrvars,h,xlower)
    use NodeInfoDef
    use GridUtilities, only: ULap2D
    use ProblemDef
    implicit none
    real, dimension(0:,0:,1:), intent(out) :: q
    integer, dimension(2),     intent(in)  :: mx
    integer,                   intent(in)  :: nrvars
    real,                      intent(in)  :: h
    real, dimension(2),        intent(in)  :: xlower
    real    :: r1,r2,r3,r4,tmp,x,y,z
    integer :: i,j
    !
    q = 0.0
    !
    do i = 0, mx(1)+1
      x = xlower(1) + (i - 0.5)*h
      do j = 0, mx(2)+1
        y  = xlower(2) + (j - 0.5)*h
        r1 = sqrt((x-1.4)**2+(y-1.6)**2)
        r2 = sqrt((x-1.8)**2+(y-1.6)**2)
        !
        if     (x< 1.7 .and. x>0.8 .and. y<5.0 .and. y>0.5) then
          q(i,j,1) = 1.0
        elseif (x< 3.7 .and. x>2.8 .and. y<5.0 .and. y>0.5) then
          q(i,j,1) = 1.0
        elseif (x< 5.2 .and. x>4.3 .and. y<4.0 .and. y>1.5) then
          q(i,j,1) = 1.0
        elseif (x< 4.5 .and. x>3.5 .and. y<2.0 .and. y>1.5) then
          q(i,j,1) = 1.0
        elseif (x< 6.7 .and. x>5.8 .and. y<5.0 .and. y>0.5) then
          q(i,j,1) = 1.0
        elseif (x< 7.5 .and. x>4.0 .and. y<5.0 .and. y>4.3) then
          q(i,j,1) = 1.0
        elseif (x<10.7 .and. x>9.8 .and. y<6.0 .and. y>0.5) then
          q(i,j,1) = 1.0
        elseif (x<11.5 .and. x>9.0 .and. y<6.0 .and. y>5.3) then
          q(i,j,1) = 1.0
        elseif (y<0.5                                     ) then
          q(i,j,1) = 1.0
        end if
      end do
    end do
  end subroutine QInit2D
  subroutine AfterStep2D(q,qc,mx,nrvars,h,xlower,level)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(0:,0:,1:), intent(in) :: q
    real, dimension(0:,0:,1:), intent(in) :: qc
    integer, dimension(1:2),   intent(in) :: mx
    integer,                   intent(in) :: nrvars
    real,                      intent(in) :: h
    real, dimension(1:2),      intent(in) :: xlower
    integer,                   intent(in) :: level
    integer, dimension(1:2) :: cmx
    integer :: tmp
    real :: ch,ch2,h2
    !
    !
    cmx(1:2) = mx(1:2)/2
    h2 = h*h; ch = h*2.0; ch2 = ch*ch
    !
    if (level==0) then
      integralresult(1) = integralresult(1)+h2*sum(q(1:mx(1),1:mx(2),1))
      integralresult(2) = integralresult(2) &
                          + Energy2D(q(0:mx(1)+1,0:mx(2)+1,1),h)
    else
      integralresult(1) = integralresult(1) &
                          + h2*sum(q( 1: mx(1),1: mx(2),1)) &
                          - ch2*sum(qc(1:cmx(1),1:cmx(2),1))
      integralresult(2) = integralresult(2) &
                          + Energy2D(q( 0: mx(1)+1,0: mx(2)+1,1), h) &
                          - Energy2D(qc(0:cmx(1)+1,0:cmx(2)+1,1),ch)
    end if
    !
  end subroutine AfterStep2D
  subroutine SetErrFlagsUser2D(qrte,qnew,errorflags,mx,cmx,nrvars,h, &
                               xlower,level)
    use ProblemDef
    use NodeInfoDef
    use GridUtilities, only: ULap2D
    implicit none
    real, dimension(:,:,:),   intent(in)  :: qrte
    real, dimension(0:,0:,:), intent(in)  :: qnew
    integer, dimension(:,:),  intent(out) :: errorflags
    integer, dimension(2),    intent(in)  :: mx
    integer, dimension(2),    intent(in)  :: cmx
    integer,                  intent(in)  :: nrvars
    real,                     intent(in)  :: h
    real, dimension(2),       intent(in)  :: xlower
    integer,                  intent(in)  :: level
    real, dimension(1:mx(1),1:mx(2)) :: agp
    !
    errorflags(1:mx(1),1:mx(2)) = 0
    !
    ! Error tagging based on the undivided gradient of the composition:
    where (abs(sqrt(  (qnew(2:mx(1)+1,1:mx(2),1) &
                    -  qnew(0:mx(1)-1,1:mx(2),1))**2 &
                    + (qnew(1:mx(1),2:mx(2)+1,1) &
                    -  qnew(1:mx(1),0:mx(2)-1,1))**2)) &
                    > qtolerance(level))
      errorflags(1:mx(1),1:mx(2)) = 1
    end where
    !
  end subroutine SetErrFlagsUser2D
  subroutine RelaxGrid2D(qnew,qold,aux,f,mx,nrvars,maux,h,xlower,ipass)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(0:,0:,1:), intent(inout) :: qnew
    real, dimension(0:,0:,1:), intent(in)    :: qold
    real, dimension(0:,0:,1:), intent(in)    :: aux
    real, dimension(1:,1:,1:), intent(in)    :: f
    integer, dimension(1:2),   intent(in)    :: mx
    integer,                   intent(in)    :: nrvars
    integer,                   intent(in)    :: maux
    real,                      intent(in)    :: h
    real, dimension(1:2),      intent(in)    :: xlower
    integer, intent(in) :: ipass
    integer :: i,j
    real    :: det,h2,m1,m2,m3,m4,tmp1,tmp2,tmp3
    real, dimension(1:2,1:2) :: a
    real, dimension(1:2)     :: b
    !
    h2 = h*h
    !
    tmp1 = dt/h2
    tmp2 = eps2/h2
    tmp3 = 4.0*eps2/h2
    !
    a(1,1) = 1.0
    a(2,2) = 1.0
    !
    do j = 1, mx(2)
      do i = 1+modulo(j+ipass,2), mx(1), 2
        m1 = Mob(0.5*(qnew(i-1,j,1)+qnew(i,j,1)))
        m2 = Mob(0.5*(qnew(i+1,j,1)+qnew(i,j,1)))
        m3 = Mob(0.5*(qnew(i,j-1,1)+qnew(i,j,1)))
        m4 = Mob(0.5*(qnew(i,j+1,1)+qnew(i,j,1)))
        !
        a(1,2) = tmp1*(m1+m2+m3+m4)
        a(2,1) = -D2fc(qnew(i,j,1))-tmp3
        !
        b(1) = f(i,j,1)+tmp1*(m1*qnew(i-1,j,2)+m2*qnew(i+1,j,2) &
               +                m3*qnew(i,j-1,2)+m4*qnew(i,j+1,2))
        !
        b(2) = f(i,j,2)+Dfc(qnew(i,j,1))-D2fc(qnew(i,j,1))*qnew(i,j,1) &
               - tmp2*(qnew(i+1,j,1)+qnew(i-1,j,1)+qnew(i,j+1,1)+qnew(i,j-1,1))
        !
        ! Solve using Cramer's rule:
        det = 1.0-a(1,2)*a(2,1)
        !
        qnew(i,j,1) = (b(1)-a(1,2)*b(2))/det
        qnew(i,j,2) = (b(2)-a(2,1)*b(1))/det
      end do
    end do
  end subroutine RelaxGrid2D
  subroutine UpdateAuxVcycle2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(ll:,ll:,1:), intent(in)    :: qnew
    real, dimension(ll:,ll:,1:), intent(in)    :: qold
    real, dimension(ll:,ll:,1:), intent(inout) :: aux
    integer,                     intent(in)    :: ll
    integer, dimension(1:2),     intent(in)    :: mx
    integer,                     intent(in)    :: nrvars
    integer,                     intent(in)    :: maux
    integer,                     intent(in)    :: mbc
    real,                        intent(in)    :: h
    real, dimension(1:2),        intent(in)    :: xlower
  end subroutine UpdateAuxVcycle2D
  function Operator2D(qnew,qold,aux,mx,nrvars,maux,h,xlower) result(res)
    use NodeInfoDef
    use GridUtilities, only: ULap2D
    use ProblemDef
    implicit none
    real, dimension(0:,0:,:), intent(in) :: qnew
    real, dimension(0:,0:,:), intent(in) :: qold
    real, dimension(0:,0:,:), intent(in) :: aux
    integer, dimension(2),    intent(in) :: mx
    integer,                  intent(in) :: nrvars
    integer,                  intent(in) :: maux
    real,                     intent(in) :: h
    real, dimension(2),       intent(in) :: xlower
    real, dimension(mx(1),mx(2),nrvars) :: res
    real :: h2,tmp1,tmp2
    h2 = h*h
    tmp1 = dt/h2
    tmp2 = eps2/h2
    !
    res(1:mx(1),1:mx(2),1) = qnew(1:mx(1),1:mx(2),1) &
                             - tmp1*ULapMob2D(qnew(0:mx(1)+1,0:mx(2)+1,1), &
                                              qnew(0:mx(1)+1,0:mx(2)+1,2))
    !
    res(1:mx(1),1:mx(2),2) = qnew(1:mx(1),1:mx(2),2) &
                             - Dfc(qnew(1:mx(1),1:mx(2),1)) &
                             + tmp2*ULap2D(qnew(0:mx(1)+1,0:mx(2)+1,1))
  end function Operator2D
  function Source2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) result(res)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(ll:,ll:,:), intent(in) :: qnew
    real, dimension(ll:,ll:,:), intent(in) :: qold
    real, dimension(ll:,ll:,:), intent(in) :: aux
    integer,                    intent(in) :: ll
    integer, dimension(2),      intent(in) :: mx
    integer,                    intent(in) :: nrvars
    integer,                    intent(in) :: maux
    integer,                    intent(in) :: mbc
    real,                       intent(in) :: h
    real, dimension(2),         intent(in) :: xlower
    real, dimension(mx(1),mx(2),nrvars) :: res
    real :: h2, tmp
    !
    h2 = h*h
    tmp = dt/h2
    !
    res(1:mx(1),1:mx(2),1) = qold(1:mx(1),1:mx(2),1)
    res(1:mx(1),1:mx(2),2) = -Dfe(qold(1:mx(1),1:mx(2),1))
  end function Source2D
  function SourceUpdate2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
        result(res)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(ll:,ll:,:), intent(in) :: qnew
    real, dimension(ll:,ll:,:), intent(in) :: qold
    real, dimension(ll:,ll:,:), intent(in) :: aux
    integer,                    intent(in) :: ll
    integer, dimension(2),      intent(in) :: mx
    integer,                    intent(in) :: nrvars
    integer,                    intent(in) :: maux
    integer,                    intent(in) :: mbc
    real,                       intent(in) :: h
    real, dimension(2),         intent(in) :: xlower
    real, dimension(mx(1),mx(2),nrvars) :: res
    real, dimension(mx(1),mx(2))        :: agp1
    real    :: x,y,hinv,g
    integer :: i,j
    !
    res(1:mx(1),1:mx(2),1) = 0.0
    res(1:mx(1),1:mx(2),2) = 0.0
  end function SourceUpdate2D
  function ULapMob2D(c,a) result(res)
    !
    ! Level independent, undivided laplacian operator for non-constant mobity.
    !
    use NodeInfoDef
    use GridUtilities, only: UDiv2D
    implicit none
    real,    dimension(0:,0:), intent(in):: c
    real,    dimension(0:,0:), intent(in):: a
    real,    dimension(1:size(a,1)-2,1:size(a,2)-2):: res
    real,    dimension(0:size(a,1)-2,1:size(a,2)-2):: f1
    real,    dimension(1:size(a,1)-2,0:size(a,2)-2):: f2
    integer, dimension(1:2):: mx
    !
    mx(1) = size(a,1)-2; mx(2) = size(a,2)-2
    !
    ! Calculate the undivided 2D flux function:
    f1(0:mx(1),1:mx(2)) = Mob(0.5*(c(1:mx(1)+1,1:mx(2))+c(0:mx(1),1:mx(2)))) &
                          * (a(1:mx(1)+1,1:mx(2))-a(0:mx(1),1:mx(2)))
    !
    f2(1:mx(1),0:mx(2)) = Mob(0.5*(c(1:mx(1),1:mx(2)+1)+c(1:mx(1),0:mx(2)))) &
                          * (a(1:mx(1),1:mx(2)+1)-a(1:mx(1),0:mx(2)))
    !
    ! Calculate the undivided divergence of the flux:
    res(1:mx(1),1:mx(2)) = UDiv2D(f1(0:mx(1),1:mx(2)), f2(1:mx(1),0:mx(2)))
  end function ULapMob2D
  function Energy2D(c,h) result(res)
    use NodeInfoDef
    use ProblemDef
    use GridUtilities, only: ULap2D
    implicit none
    real, dimension(0:,0:), intent(in) :: c
    real,                   intent(in) :: h
    real                               :: res
    real, dimension(1:size(c,1)-2,1:size(c,2)-2) :: lap
    real, dimension(0:size(c,1)-2,1:size(c,2)-2) :: f1
    real, dimension(1:size(c,1)-2,0:size(c,2)-2) :: f2
    real, dimension(1:2) :: dc, nrml
    real :: denominator, h2, tmp, gam
    integer, dimension(1:2) :: mx
    integer :: i, j
    !
    mx(1) = size(c,1)-2; mx(2) = size(c,2)-2
    h2 = h*h
    !
    ! Calculate the gradient contributions to the energy:
    f1(0:mx(1),1:mx(2)) = c(1:mx(1)+1,1:mx(2))-c(0:mx(1),1:mx(2))
    f1(0:mx(1),1:mx(2)) = f1(0:mx(1),1:mx(2))*f1(0:mx(1),1:mx(2))
    !
    f2(1:mx(1),0:mx(2)) = c(1:mx(1),1:mx(2)+1)-c(1:mx(1),0:mx(2))
    f2(1:mx(1),0:mx(2)) = f2(1:mx(1),0:mx(2))*f2(1:mx(1),0:mx(2))
    !
    tmp = eps2/h2/4.0
    !
    res = h2*(      sum(Ff(c(1:mx(1),1:mx(2)))) &
              + tmp*sum(f1(1:mx(1),1:mx(2))+f1(0:mx(1)-1,1:mx(2)  ) &
              +         f2(1:mx(1),1:mx(2))+f2(1:mx(1)  ,0:mx(2)-1)))
  end function Energy2D
  !
  ! 3D Multigrid Routines
  !
  subroutine QInit3D(q,mx,nrvars,h,xlower)
    use NodeInfoDef
    use GridUtilities, only: ULap3D
    use ProblemDef
    implicit none
    real, dimension(0:,0:,0:,1:), intent(out) :: q
    integer, dimension(3),        intent(in) :: mx
    integer,                      intent(in) :: nrvars
    real,                         intent(in) :: h
    real, dimension(3),           intent(in) :: xlower
    real    :: r1,r2,r3,r4,tmp,x,y,z
    integer :: i j k
    !
    q = 0.0
    !
    do i = 0, mx(1)+1
      x = (i-0.5)*h+xlower(1)
      do j = 0, mx(2)+1
        y = (j-0.5)*h+xlower(2)
        do k = 0, mx(3)+1
          z = (k-0.5)*h+xlower(3)
          !
          r1 = sqrt((x-1.6)**2+(y-1.6)**2+(z-1.6)**2)
          r2 = sqrt((x-1.6)**2+(y-1.2)**2+(z-1.2)**2)
          !
          tmp =      1.0-0.5*(1.0-tanh((r1-1.0)/(2.0*sqrt(2.0)*eps)))
          tmp = tmp*(1.0-0.5*(1.0-tanh((r2-0.5)/(2.0*sqrt(2.0)*eps))))
          !
          q(i,j,k,1) = 1.0-tmp
        end do
      end do
    end do
  end subroutine QInit3D
  subroutine AfterStep3D(q,qparent,mx,nrvars,h,xlower,level)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(0:,0:,0:,1:), intent(in) :: q
    real, dimension(0:,0:,0:,1:), intent(in) :: qparent
    integer, dimension(1:3),      intent(in) :: mx
    integer,                      intent(in) :: nrvars
    real,                         intent(in) :: h
    real, dimension(1:3),         intent(in) :: xlower
    integer,                      intent(in) :: level
    integer :: tmp
    integer, dimension(1:3) :: cmx
    real :: ch, ch3, h3
    !
    cmx(1:3) = mx(1:3)/2
    h3 = h*h*h; ch = h*2.0; ch3 = ch*ch*ch
    !
    if (level==0) then
      integralresult(1) = integralresult(1) &
                          + h3*sum(q(1:mx(1),1:mx(2),1:mx(3),1))
    else
      integralresult(1) = integralresult(1) &
                          +  h3*sum(q( 1: mx(1),1: mx(2),1: mx(3),1)) &
                          - ch3*sum(qc(1:cmx(1),1:cmx(2),1:cmx(3),1))
    end if
  end subroutine AfterStep3D
  subroutine SetErrFlagsUser3D(qrte,qnew,errorflags,mx,cmx,nrvars,h, &
                               xlower,level)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(:,:,:,:),     intent(in) :: qrte
    real, dimension(0:,0:,0:,:),  intent(in) :: qnew
    integer, dimension(:,:,:),   intent(out) :: errorflags
    integer, dimension(3),        intent(in) :: mx
    integer, dimension(3),        intent(in) :: cmx
    integer,                      intent(in) :: nrvars
    real,                         intent(in) :: h
    real, dimension(3),           intent(in) :: xlower
    integer,                      intent(in) :: level
    !
    errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0
    !
    ! Error tagging based on the undivided gradient of the composition:
    where (abs(sqrt(  (qnew(2:mx(1)+1,1:mx(2),1:mx(3),1) &
                    -  qnew(0:mx(1)-1,1:mx(2),1:mx(3),1))**2 &
                    + (qnew(1:mx(1),2:mx(2)+1,1:mx(3),1) &
                    -  qnew(1:mx(1),0:mx(2)-1,1:mx(3),1))**2 &
                    + (qnew(1:mx(1),1:mx(2),2:mx(3)+1,1) &
                    -  qnew(1:mx(1),1:mx(2),0:mx(3)-1,1))**2)) &
                    > qtolerance(level))
          errorflags(1:mx(1),1:mx(2),1:mx(3)) = 1
    end where
  end subroutine SetErrFlagsUser3D
  subroutine RelaxGrid3D(qnew,qold,aux,f,mx,nrvars,maux,h,xlower,ipass)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(0:,0:,0:,1:), intent(in out) :: qnew
    real, dimension(0:,0:,0:,1:),     intent(in) :: qold
    real, dimension(0:,0:,0:,1:),     intent(in) :: aux
    real, dimension(1:,1:,1:,1:),     intent(in) :: f
    integer, dimension(1:3),          intent(in) :: mx
    integer,                          intent(in) :: nrvars
    integer,                          intent(in) :: maux
    real,                             intent(in) :: h
    real, dimension(1:3),             intent(in) :: xlower
    integer,                          intent(in) :: ipass
    real, dimension(1:2,1:2) :: a
    real, dimension(1:2)     :: b
    real    :: d1,d2,d3,d4,d5,d6,det,h2,m1,m2,m3,m4,m5,m6,tmp1,tmp2,tmp3,tmp4
    integer :: i,j,k
    !
    h2 = h*h
    !
    tmp1 = dt/h2
    tmp2 = eps2/h2
    tmp3 = 6.0*eps2/h2
    !
    a(1,1) = 1.0
    a(2,2) = 1.0
    !
    ! These loops must always be structured this way, i.e., so that i,j,k = 1
    ! is a black cell. (The first cell relaxed is always i = 2, j,k = 1, since
    ! ipass = 1 is called first.)
    !
    do k = 1, mx(3)
      do j = 1, mx(2)
        do i = 1+modulo(k+j+ipass,2), mx(1), 2
          !
          m1 = Mob(0.5*(qnew(i-1,j,k,1)+qnew(i,j,k,1)))
          m2 = Mob(0.5*(qnew(i+1,j,k,1)+qnew(i,j,k,1)))
          m3 = Mob(0.5*(qnew(i,j-1,k,1)+qnew(i,j,k,1)))
          m4 = Mob(0.5*(qnew(i,j+1,k,1)+qnew(i,j,k,1)))
          m5 = Mob(0.5*(qnew(i,j,k-1,1)+qnew(i,j,k,1)))
          m6 = Mob(0.5*(qnew(i,j,k+1,1)+qnew(i,j,k,1)))
          !
          a(1,2) = tmp1*(m1+m2+m3+m4+m5+m6)
          a(2,1) = -D2fc(qnew(i,j,k,1))-tmp3
          !
          b(1) = f(i,j,k,1)+tmp1*(m1*qnew(i-1,j,k,2)+m2*qnew(i+1,j,k,2) &
                 +                  m3*qnew(i,j-1,k,2)+m4*qnew(i,j+1,k,2) &
                 +                  m5*qnew(i,j,k-1,2)+m6*qnew(i,j,k+1,2))
          !
          b(2) = f(i,j,k,2)+Dfc(qnew(i,j,k,1)) &
                 - D2fc(qnew(i,j,k,1))*qnew(i,j,k,1) &
                 - tmp2*(qnew(i+1,j,k,1)+qnew(i-1,j,k,1) &
                 +       qnew(i,j+1,k,1)+qnew(i,j-1,k,1) &
                 +       qnew(i,j,k+1,1)+qnew(i,j,k-1,1))
          !
          ! Solve for Cahn-Hilliard using Cramer's rule:
          det = 1.0-a(1,2)*a(2,1)
          !
          qnew(i,j,k,1) = (b(1)-a(1,2)*b(2))/det
          qnew(i,j,k,2) = (b(2)-a(2,1)*b(1))/det
          !
        end do
      end do
    end do
  end subroutine RelaxGrid3D
  subroutine UpdateAuxVcycle3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, dimension(ll:,ll:,ll:,1:),    intent(in) :: qnew
    real, dimension(ll:,ll:,ll:,1:),    intent(in) :: qold
    real, dimension(ll:,ll:,ll:,1:), intent(inout) :: aux
    integer,                            intent(in) :: ll
    integer, dimension(1:3),            intent(in) :: mx
    integer,                            intent(in) :: nrvars
    integer,                            intent(in) :: maux
    integer,                            intent(in) :: mbc
    real,                               intent(in) :: h
    real, dimension(1:3),               intent(in) :: xlower
  end subroutine UpdateAuxVcycle3D
  function Operator3D(qnew,qold,aux,mx,nrvars,maux,h,xlower) result(res)
    use NodeInfoDef
    use GridUtilities, only: ULap3D
    use ProblemDef
    implicit none
    real, dimension(0:,0:,0:,1:), intent(in) :: qnew
    real, dimension(0:,0:,0:,1:), intent(in) :: qold
    real, dimension(0:,0:,0:,1:), intent(in) :: aux
    integer, dimension(1:3),      intent(in) :: mx
    integer,                      intent(in) :: nrvars
    integer,                      intent(in) :: maux
    real,                         intent(in) :: h
    real, dimension(1:3),         intent(in) :: xlower
    real, dimension(1:mx(1),1:mx(2),1:mx(3),1:nrvars) :: res
    real :: h2, tmp1, tmp2
    !
    h2 = h*h
    tmp1 = dt/h2
    tmp2 = eps2/h2
    !
    res(1:mx(1),1:mx(2),1:mx(3),1) = qnew(1:mx(1),1:mx(2),1:mx(3),1) &
                    - tmp1*ULapMob3D(qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1), &
                                     qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,2))
    !
    res(1:mx(1),1:mx(2),1:mx(3),2) = qnew(1:mx(1),1:mx(2),1:mx(3),2) &
                    - Dfc(qnew(1:mx(1),1:mx(2),1:mx(3),1)) &
                    + tmp2*ULap3D(qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1))
  end function Operator3D
  function Source3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) result(res)
    use NodeInfoDef
    use GridUtilities, only: ULap3D
    use ProblemDef
    implicit none
    real, dimension(ll:,ll:,ll:,1:), intent(in) :: qnew
    real, dimension(ll:,ll:,ll:,1:), intent(in) :: qold
    real, dimension(ll:,ll:,ll:,1:), intent(in) :: aux
    integer,                         intent(in) :: ll
    integer, dimension(1:3),         intent(in) :: mx
    integer,                         intent(in) :: nrvars
    integer,                         intent(in) :: maux
    integer,                         intent(in) :: mbc
    real,                            intent(in) :: h
    real, dimension(1:3),            intent(in) :: xlower
    real, dimension(1:mx(1),1:mx(2),1:mx(3),1:nrvars) :: res
    real :: h2, tmp
    !
    h2 = h*h
    tmp = dt/h2
    !
    res(1:mx(1),1:mx(2),1:mx(3),1) = qold(1:mx(1),1:mx(2),1:mx(3),1)
    res(1:mx(1),1:mx(2),1:mx(3),2) = - Dfe(qold(1:mx(1),1:mx(2),1:mx(3),1))
  end function Source3D
  function SourceUpdate3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
        result(res)
    use NodeInfoDef
    use GridUtilities, only: UDiv3D, ULap3D
    use ProblemDef
    implicit none
    real, dimension(ll:,ll:,ll:,1:), intent(in) :: qnew
    real, dimension(ll:,ll:,ll:,1:), intent(in) :: qold
    real, dimension(ll:,ll:,ll:,1:), intent(in) :: aux
    integer,                         intent(in) :: ll
    integer, dimension(1:3),         intent(in) :: mx
    integer,                         intent(in) :: nrvars
    integer,                         intent(in) :: maux
    integer,                         intent(in) :: mbc
    real,                            intent(in) :: h
    real, dimension(1:3),            intent(in) :: xlower
    real, dimension(1:mx(1),1:mx(2),1:mx(3),1:nrvars) :: res
    !
    res = 0.0
    !
  end function SourceUpdate3D
  function ULapMob3D(c,a) result (res)
    !
    ! Level independent 3D undivided laplacian operator for non-constant mobity.
    !
    use NodeInfoDef
    use GridUtilities, only: UDiv3D
    implicit none
    real, dimension(0:,0:,0:), intent(in):: c
    real, dimension(0:,0:,0:), intent(in):: a
    real, dimension(1:size(a,1)-2, 1:size(a,2)-2, 1:size(a,3)-2) :: res
    real, dimension(0:size(a,1)-2,1:size(a,2)-2,1:size(a,3)-2):: f1
    real, dimension(1:size(a,1)-2,0:size(a,2)-2,1:size(a,3)-2):: f2
    real, dimension(1:size(a,1)-2,1:size(a,2)-2,0:size(a,3)-2):: f3
    integer, dimension(1:3):: mx
    !
    mx(1) = size(a,1)-2; mx(2) = size(a,2)-2; mx(3) = size(a,3)-2
    !
    ! Calculate the 3D flux function:
    f1(0:mx(1),1:mx(2),1:mx(3)) &
       = Mob(0.5_r8*(c(1:mx(1)+1,1:mx(2),1:mx(3))+c(0:mx(1),1:mx(2),1:mx(3)))) &
         * (a(1:mx(1)+1,1:mx(2),1:mx(3))-a(0:mx(1),1:mx(2),1:mx(3)))
    !
    f2(1:mx(1),0:mx(2),1:mx(3)) &
       = Mob(0.5_r8*(c(1:mx(1),1:mx(2)+1,1:mx(3))+c(1:mx(1),0:mx(2),1:mx(3)))) &
         * (a(1:mx(1),1:mx(2)+1,1:mx(3))-a(1:mx(1),0:mx(2),1:mx(3)))
    !
    f3(1:mx(1),1:mx(2),0:mx(3)) &
       = Mob(0.5_r8*(c(1:mx(1),1:mx(2),1:mx(3)+1)+c(1:mx(1),1:mx(2),0:mx(3)))) &
         * (a(1:mx(1),1:mx(2),1:mx(3)+1)-a(1:mx(1),1:mx(2),0:mx(3)))
    !
    ! Calculate the divergence of the flux:
    ulapmobresult(1:mx(1),1:mx(2),1:mx(3)) &
                  = UDiv3D(f1(0:mx(1),1:mx(2),1:mx(3)), &
                           f2(1:mx(1),0:mx(2),1:mx(3)), &
                           f3(1:mx(1),1:mx(2),0:mx(3)))
  end function ULapMob3D
  !
  ! User supplied physical boundary conditions
  !
  subroutine UserBC2D(q,ll,mx,nrvars,h,xlower,mbc,bcnow)
    use NodeinfoDef
    use ProblemDef
    implicit none
    real, dimension(ll:,ll:,1:), intent(inout) :: q
    integer,                        intent(in) :: ll
    integer, dimension(1:2),        intent(in) :: mx
    integer,                        intent(in) :: nrvars
    real,                           intent(in) :: h
    real, dimension(1:2),           intent(in) :: xlower
    integer,                        intent(in) :: mbc
    integer,                        intent(in) :: bcnow
    integer                                    :: ibc
    integer, dimension(1:2)                    :: ul
    !
    ul(1:2) = mx(1:2)+mbc
    !
    select case(bcnow)
    case(1)
      !
      ! dirichlet on bndry 1:
      do ibc = 1, mbc
        q(1-ibc,ll:ul(2),:) = 2.0-q(ibc,ll:ul(2),:)
      end do
      !
    case(2)
      !
      ! dirichlet on bndry 2:
      do ibc = 1, mbc
        q(mx(1)+ibc,ll:ul(2),:) = 2.0-q(mx(1)-ibc+1,ll:ul(2),:)
      end do
      !
    case(3)
      !
      ! dirichlet on bndry 3, variables 3 and 4:
      do ibc = 1, mbc
        q(ll:ul(1),1-ibc,1) = q(ll:ul(1),ibc,1)
        q(ll:ul(1),1-ibc,2) = q(ll:ul(1),ibc,2)
      end do
      !
    case(4)
      !
      ! dirichlet on bndry 4, variables 3 and 4:
      do ibc = 1, mbc
        q(ll:ul(1),mx(2)+ibc,1) = q(ll:ul(1),mx(2)-ibc+1,1)
        q(ll:ul(1),mx(2)+ibc,2) = q(ll:ul(1),mx(2)-ibc+1,2)
      end do          
    end select
    !
  end subroutine UserBC2D
  subroutine UserBC3D(q,ll,mx,nrvars,mbc,bcnow)
    use NodeinfoDef
    implicit none
    real, dimension(ll:,ll:,ll:,1:), intent(inout) :: q
    integer, intent(in) :: ll
    integer, dimension(1:3), intent(in) :: mx
    integer, intent(in) :: nrvars
    integer, intent(in) :: mbc
    integer, intent(in) :: bcnow
    integer, dimension(1:3) :: ul
    integer                 :: ibc
    !
    ul(1:3) = mx(1:3)+mbc
    !
    select case (bcnow)
    case(1)
      !
      ! dirichlet on bndry 1:
      do ibc = 1, mbc
        q(1-ibc,ll:ul(2),ll:ul(3),:) = 2.0-q(ibc,ll:ul(2),ll:ul(3),:)
      end do
      !
    case(2)
      !
      ! dirichlet on bndry 2:
      do ibc = 1, mbc
        q(mx(1)+ibc,ll:ul(2),ll:ul(3),:) = 2.0 &
                                           - q(mx(1)-ibc+1,ll:ul(2),ll:ul(3),:)
      end do
      !
    case(3)
      !
      ! dirichlet on bndry 3:
      do ibc = 1, mbc
        q(ll:ul(1),1-ibc,ll:ul(3),:) = 2.0-q(ll:ul(1),ibc,ll:ul(3),:)
      end do
      !
    case(4)
      !
      ! dirichlet on bndry 4:
      do ibc = 1, mbc
        q(ll:ul(1),mx(2)+ibc,ll:ul(3),:) = 2.0 &
                                           - q(ll:ul(1),mx(2)-ibc+1,ll:ul(3),:)
      end do
      !
    case(5)
      !
      ! dirichlet on bndry 5:
      do ibc = 1, mbc
        q(ll:ul(1),ll:ul(2),1-ibc,:) = 2.0-q(ll:ul(1),ll:ul(2),ibc,:)
      end do
      !
    case(6)
      !
      ! dirichlet on bndry 6:
      do ibc = 1, mbc
        q(ll:ul(1),ll:ul(2),mx(3)+ibc,:) = 2.0 &
                                           - q(ll:ul(1),ll:ul(2),mx(3)-ibc+1,:)
      end do
    end select
  end subroutine UserBC3D
  !
  ! Dimensionally Invariant Multigrid Routines
  !
  elemental function Mob(c) result (res)
    use nodeinfodef
    use ProblemDef
    implicit none
    real, intent(in) :: c
    real :: res
    !
    res = 1.0
  end function Mob
  elemental function Ff(c) result (res)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, intent(in) :: c
    real :: res
    real :: cm1
    !
    cm1 = c-1.0
    !
    res = c*c*cm1*cm1*alpha/4.0
  end function Ff
  elemental function Dfe(c) result (res)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, intent(in) :: c
    real :: res
    !
    res = (c-0.5)*alpha/4.0
  end function Dfe
  elemental function Dfc(c) result (res)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, intent(in) :: c
    real :: res
    real :: cmh
    !
    cmh = c-0.5
    !
    res = cmh*cmh*cmh*alpha
  end function Dfc
  elemental function D2fc(c) result (res)
    use NodeInfoDef
    use ProblemDef
    implicit none
    real, intent(in) :: c
    real :: res
    real :: cmh
    !
    cmh = c-0.5
    !
    res = 3.0*cmh*cmh*alpha
  end function D2fc
  logical function MultiGridUserBreak(eqtag, residual, oldresidual)
    implicit none
    real,    intent(in) :: residual, oldresidual
    integer, intent(in) :: eqtag
    !
    MultiGridUserBreak = (oldresidual <= residual)
  end function MultiGridUserBreak
end module Problem
