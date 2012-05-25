MODULE Problem
!
! This MODULE CONTAINS routines for initializing the problem and for driving
! the Multigrid routines in the MODULE AFASRoutines.  In particular, the
! routines in AFASRoutines USE RelaxGrid#D, Operator#D, Source#D, and
! SourceCorrection#D where # = 2 and 3.
!
! Templates for the Cahn-Hilliard Equation.  We use Eyre's gradient stable
! scheme for time stepping.
!
CONTAINS
!
SUBROUTINE SetProb
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
INTEGER:: ierror
NAMELIST/problemdata/alpha, eps
!
OPEN(UNIT=75,FILE='problemdata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file problemdata.dat. Program stop.'
  STOP
END IF
READ(75,NML=problemdata)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='OLD',ACTION='WRITE',FORM='FORMATTED', &
     POSITION='APPEND')
WRITE(76,NML=problemdata)
CLOSE(76)
!
eps2 = eps*eps
!
END SUBROUTINE SetProb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AfterRun
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
END SUBROUTINE AfterRun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE QInit2D(q,mx,nrvars,h,xlower)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(OUT):: q
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
INTEGER:: i, j
REAL(KIND=r8):: r1, r2, r3, r4, tmp, x, y, z
!
q = 0.0_r8
!
DO i = 1, mx(1)
  x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
  DO j = 1, mx(2)
    y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
!
    r1 = SQRT((x-1.4_r8)**2+(y-1.6_r8)**2)
    r2 = SQRT((x-1.8_r8)**2+(y-1.6_r8)**2)
!    r3 = SQRT((y-1.6_r8)**2+(x-3.2_r8)**2)
!    r4 = SQRT((y-1.6_r8)**2+(x-0.0_r8)**2)
!
!    tmp =      1.0-0.5*(1.0-TANH((r1-1.0)/(2.0*SQRT(2.0)*eps)))
!    tmp = tmp*(1.0-0.5*(1.0-TANH((r2-1.0)/(2.0*SQRT(2.0)*eps))))
!    tmp = tmp*(1.0-0.5*(1.0-TANH((r3-0.5)/(2.0*SQRT(2.0)*eps))))
!    tmp = tmp*(1.0-0.5*(1.0-TANH((r4-0.5)/(2.0*SQRT(2.0)*eps))))
!
!    q(i,j,1) = 1.0_r8-tmp
!
    IF      (x<1.7 .AND. x>0.8 .AND. y<5.0 .AND. y>0.5) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (x<3.7 .AND. x>2.8 .AND. y<5.0 .AND. y>0.5) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (x<5.2 .AND. x>4.3 .AND. y<4.0 .AND. y>1.5) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (x<4.5 .AND. x>3.5 .AND. y<2.0 .AND. y>1.5) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (x<6.7 .AND. x>5.8 .AND. y<5.0 .AND. y>0.5) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (x<7.5 .AND. x>4.0 .AND. y<5.0 .AND. y>4.3) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (x<10.7 .AND. x>9.8 .AND. y<6.0 .AND. y>0.5) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (x<11.5 .AND. x>9.0 .AND. y<6.0 .AND. y>5.3) THEN
      q(i,j,1) = 1.0_r8
    ELSE IF (y<0.5                                    ) THEN
      q(i,j,1) = 1.0_r8
    END IF
!
!    q(i,j,1) = 0.5*(1.0-TANH((x-1.6_r8)/(2.0*SQRT(2.0)*eps)))
!
!    CALL RANDOM_NUMBER(r1)
!    r1 = 0.05_r8*(2.0_r8*r1-1.0_r8)
!    q(i,j,1) = 0.5_r8+r1
!
  END DO
END DO
!
END SUBROUTINE QInit2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE QInit3D(q,mx,nrvars,h,xlower)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(OUT):: q
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
!
INTEGER:: i, j, k
REAL(KIND=r8):: r1, r2, r3, r4, tmp, x, y, z
!
DO i = 1, mx(1)
  x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
  DO j = 1, mx(2)
    y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
    DO k = 1, mx(3)
      z = (REAL(k,KIND=r8)-0.5_r8)*h+xlower(3)
!
      r1 = SQRT((x-1.6_r8)**2+(y-1.6_r8)**2+(z-1.6_r8)**2)
      r2 = SQRT((x-1.6_r8)**2+(y-1.2_r8)**2+(z-1.2_r8)**2)
!
      tmp =      1.0-0.5*(1.0-TANH((r1-1.0)/(2.0*SQRT(2.0)*eps)))
      tmp = tmp*(1.0-0.5*(1.0-TANH((r2-0.5)/(2.0*SQRT(2.0)*eps))))
!
      q(i,j,k,1) = 1.0_r8-tmp
!
    END DO
  END DO
END DO
!
END SUBROUTINE QInit3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AfterStep2D(q,qc,mx,nrvars,h,xlower,level)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: q
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qc
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
INTEGER:: tmp
INTEGER, DIMENSION(1:2):: cmx
REAL(KIND=r8):: ch, ch2, h2
!
cmx(1:2) = mx(1:2)/2
h2 = h*h; ch = h*2.0_r8; ch2 = ch*ch
!
IF(level==0) THEN
  integralresult(1) &
    = integralresult(1)+h2*SUM(q(1:mx(1),1:mx(2),1))
  integralresult(2) &
    = integralresult(2)+Energy2D(q(0:mx(1)+1,0:mx(2)+1,1),h)
ELSE
  integralresult(1) &
    = integralresult(1)+h2*SUM(q( 1: mx(1),1: mx(2),1)) &
    -                  ch2*SUM(qc(1:cmx(1),1:cmx(2),1))
  integralresult(2) &
    = integralresult(2)+Energy2D(q( 0: mx(1)+1,0: mx(2)+1,1), h) &
    -                   Energy2D(qc(0:cmx(1)+1,0:cmx(2)+1,1),ch)
END IF
!
END SUBROUTINE AfterStep2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AfterStep3D(q,qc,mx,nrvars,h,xlower,level)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: q
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qc
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
INTEGER:: tmp
INTEGER, DIMENSION(1:3):: cmx
REAL(KIND=r8):: ch, ch3, h3
!
cmx(1:3) = mx(1:3)/2
h3 = h*h*h; ch = h*2.0_r8; ch3 = ch*ch*ch
!
IF(level==0) THEN
  integralresult(1) &
    = integralresult(1)+h3*SUM(q(1:mx(1),1:mx(2),1:mx(3),1))
ELSE
  integralresult(1) &
    = integralresult(1)+h3*SUM(q( 1: mx(1),1: mx(2),1: mx(3),1)) &
    -                  ch3*SUM(qc(1:cmx(1),1:cmx(2),1:cmx(3),1))
END IF
!
END SUBROUTINE AfterStep3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
SUBROUTINE SetAux(info)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
INTEGER:: i, j, k
REAL(KIND=r8):: time, x, y, z
!
!time = info%time
!
END SUBROUTINE SetAux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
SUBROUTINE SetSrc(info)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!    
TYPE(nodeinfo):: info
!
INTEGER:: i, j, k
REAL(KIND=r8):: time, x, y, z
!
!time = info%time
!
END SUBROUTINE SetSrc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetErrFlagsUser2D(qrte,qnew,errorflags,mx,cmx,nrvars,h,xlower,level)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: qrte
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qnew
INTEGER, DIMENSION(1:,1:), INTENT(OUT):: errorflags
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, DIMENSION(1:2), INTENT(IN):: cmx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
errorflags(1:mx(1),1:mx(2)) = 0
!
! Error tagging based on the undivided gradient of the composition:
WHERE(ABS(SQRT(  (qnew(2:mx(1)+1,1:mx(2),1) &
               -  qnew(0:mx(1)-1,1:mx(2),1))**2 &
               + (qnew(1:mx(1),2:mx(2)+1,1) &
               -  qnew(1:mx(1),0:mx(2)-1,1))**2)) &
               > qtolerance(level))
      errorflags(1:mx(1),1:mx(2)) = 1
END WHERE
!
END SUBROUTINE SetErrFlagsUser2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetErrFlagsUser3D(qrte,qnew,errorflags,mx,cmx,nrvars,h,xlower,level)
USE NodeInfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: qrte
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qnew
INTEGER, DIMENSION(1:,1:,1:), INTENT(OUT):: errorflags
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, DIMENSION(1:3), INTENT(IN):: cmx
INTEGER, INTENT(IN):: nrvars
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
INTEGER, INTENT(IN):: level
!
errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0
!
! Error tagging based on the undivided gradient of the composition:
WHERE(ABS(SQRT(  (qnew(2:mx(1)+1,1:mx(2),1:mx(3),1) &
               -  qnew(0:mx(1)-1,1:mx(2),1:mx(3),1))**2 &
               + (qnew(1:mx(1),2:mx(2)+1,1:mx(3),1) &
               -  qnew(1:mx(1),0:mx(2)-1,1:mx(3),1))**2 &
               + (qnew(1:mx(1),1:mx(2),2:mx(3)+1,1) &
               -  qnew(1:mx(1),1:mx(2),0:mx(3)-1,1))**2)) &
               > qtolerance(level))
      errorflags(1:mx(1),1:mx(2),1:mx(3)) = 1
END WHERE
!
END SUBROUTINE SetErrFlagsUser3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User-Supplied AFAS Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 2D Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine RelaxGrid2D(qnew,qold,aux,f,mx,nrvars,maux,h,xlower,ipass)
  use NodeInfoDef
  use ProblemDef
  implicit none
  !
  real(kind=r8), dimension(0:,0:,1:), intent(in out):: qnew
  real(kind=r8), dimension(0:,0:,1:), intent(in):: qold
  real(kind=r8), dimension(0:,0:,1:), intent(in):: aux
  real(kind=r8), dimension(1:,1:,1:), intent(in):: f
  integer, dimension(1:2), intent(in):: mx
  integer, intent(in):: nrvars
  integer, intent(in):: maux
  real(kind=r8), intent(in):: h
  real(kind=r8), dimension(1:2), intent(in):: xlower
  integer, intent(in):: ipass
  integer:: i, j
  real(kind=r8):: det, h2, m1, m2, m3, m4, tmp1, tmp2, tmp3
  real(kind=r8), dimension(1:2,1:2):: a
  real(kind=r8), dimension(1:2):: b
  !
  h2 = h*h
  tmp1 = dt/h2
  tmp2 = eps2/h2
  tmp3 = 4.0*eps2/h2
  a(1,1) = 1.0
  a(2,2) = 1.0
  !
  do j = 1, mx(2)
    do i = 1+modulo(j+ipass,2), mx(1), 2
      !
      m1 = mob(0.5*(qnew(i-1,j,1)+qnew(i,j,1)))
      m2 = mob(0.5*(qnew(i+1,j,1)+qnew(i,j,1)))
      m3 = mob(0.5*(qnew(i,j-1,1)+qnew(i,j,1)))
      m4 = mob(0.5*(qnew(i,j+1,1)+qnew(i,j,1)))
      !
      a(1,2) = tmp1*(m1+m2+m3+m4)
      a(2,1) = -d2fc(qnew(i,j,1))-tmp3
      !
      b(1) = f(i,j,1)+tmp1*(m1*qnew(i-1,j,2)+m2*qnew(i+1,j,2) &
           +                m3*qnew(i,j-1,2)+m4*qnew(i,j+1,2))
      !
      b(2) = f(i,j,2)+dfc(qnew(i,j,1))-d2fc(qnew(i,j,1))*qnew(i,j,1) &
          - tmp2*(qnew(i+1,j,1)+qnew(i-1,j,1)+qnew(i,j+1,1)+qnew(i,j-1,1))
      !
      ! solve using cramer's rule:
      det = 1.0-a(1,2)*a(2,1)
      !
      qnew(i,j,1) = (b(1)-a(1,2)*b(2))/det
      qnew(i,j,2) = (b(2)-a(2,1)*b(1))/det
      !
    end do
  end do
  !
end subroutine RelaxGrid2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UpdateAuxVcycle2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN OUT):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
!aux(0:mx(1)+1,0:mx(2)+1,1) = 
!
END SUBROUTINE UpdateAuxVcycle2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator2D(qnew,qold,aux,mx,nrvars,maux,h,xlower) RESULT(opresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: aux
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nrvars):: opresult
!
REAL(KIND=r8):: h2, tmp1, tmp2
!
h2 = h*h
tmp1 = dt/h2
tmp2 = eps2/h2
!
opresult(1:mx(1),1:mx(2),1) = qnew(1:mx(1),1:mx(2),1) &
                            - tmp1*ULapMob2D(qnew(0:mx(1)+1,0:mx(2)+1,1), &
                                               qnew(0:mx(1)+1,0:mx(2)+1,2))
!
opresult(1:mx(1),1:mx(2),2) = qnew(1:mx(1),1:mx(2),2) &
                            - Dfc(qnew(1:mx(1),1:mx(2),1)) &
                            + tmp2*ULap2D(qnew(0:mx(1)+1,0:mx(2)+1,1))
!
END FUNCTION Operator2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Source2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(sourceresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nrvars):: sourceresult
!
REAL(KIND=r8):: h2, tmp
!
h2 = h*h
tmp = dt/h2
!
sourceresult(1:mx(1),1:mx(2),1) = qold(1:mx(1),1:mx(2),1)
!
sourceresult(1:mx(1),1:mx(2),2) = -Dfe(qold(1:mx(1),1:mx(2),1))
!
END FUNCTION Source2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SourceUpdate2D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(updateresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv2D, ULap2D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nrvars):: updateresult
!
updateresult(1:mx(1),1:mx(2),1) = 0.0_r8
!
updateresult(1:mx(1),1:mx(2),2) = 0.0_r8
!
END FUNCTION SourceUpdate2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULapMob2D(c,a) RESULT(ulapmobresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv2D
IMPLICIT NONE
!
! Level independent, 2D UNDIVIDED laplacian operator for non-constant mobity.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: c
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2):: ulapmobresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(a,1)-2,1:SIZE(a,2)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,0:SIZE(a,2)-2):: f2
!
mx(1) = SIZE(a,1)-2; mx(2) = SIZE(a,2)-2
!
! Calculate the UNDIVIDED 2D flux function:
f1(0:mx(1),1:mx(2)) &
  = Mob(0.5_r8*(c(1:mx(1)+1,1:mx(2))+c(0:mx(1),1:mx(2)))) &
  *              (a(1:mx(1)+1,1:mx(2))-a(0:mx(1),1:mx(2)))
!
f2(1:mx(1),0:mx(2)) &
  = Mob(0.5_r8*(c(1:mx(1),1:mx(2)+1)+c(1:mx(1),0:mx(2)))) &
  *              (a(1:mx(1),1:mx(2)+1)-a(1:mx(1),0:mx(2)))
!
! Calculate the UNDIVIDED divergence of the flux:
ulapmobresult(1:mx(1),1:mx(2)) = UDiv2D(f1(0:mx(1),1:mx(2)), &
                                          f2(1:mx(1),0:mx(2)))
!
END FUNCTION ULapMob2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Energy2D(c,h) RESULT(energyresult)
USE NodeInfoDef
USE ProblemDef
USE GridUtilities, ONLY: ULap2D
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: c
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8):: energyresult
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: denominator, h2, tmp, gam
REAL(KIND=r8), DIMENSION(1:2):: dc, nrml
REAL(KIND=r8), DIMENSION(1:SIZE(c,1)-2,1:SIZE(c,2)-2):: lap
REAL(KIND=r8), DIMENSION(0:SIZE(c,1)-2,1:SIZE(c,2)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(c,1)-2,0:SIZE(c,2)-2):: f2
!
mx(1) = SIZE(c,1)-2; mx(2) = SIZE(c,2)-2
h2 = h*h
!
! Calculate the gradient contributions to the energy:
f1(0:mx(1),1:mx(2)) = c(1:mx(1)+1,1:mx(2))-c(0:mx(1),1:mx(2))
f1(0:mx(1),1:mx(2)) = f1(0:mx(1),1:mx(2))*f1(0:mx(1),1:mx(2))
!
f2(1:mx(1),0:mx(2)) = c(1:mx(1),1:mx(2)+1)-c(1:mx(1),0:mx(2))
f2(1:mx(1),0:mx(2)) = f2(1:mx(1),0:mx(2))*f2(1:mx(1),0:mx(2))
!
tmp = eps2/h2/4.0_r8
!
energyresult &
  = h2*(    SUM(Ff(c(1:mx(1),1:mx(2)))) &
  +     tmp*SUM(f1(1:mx(1),1:mx(2))+f1(0:mx(1)-1,1:mx(2)  ) &
  +             f2(1:mx(1),1:mx(2))+f2(1:mx(1)  ,0:mx(2)-1)))
!
END FUNCTION Energy2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 3D Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE RelaxGrid3D(qnew,qold,aux,f,mx,nrvars,maux,h,xlower,ipass)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN OUT):: qnew
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: aux
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: f
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
INTEGER, INTENT(IN):: ipass
!
INTEGER:: i, j, k
REAL(KIND=r8):: d1, d2, d3, d4, d5, d6, det, h2, m1, m2, m3, m4, m5, m6, &
                tmp1, tmp2, tmp3, tmp4
REAL(KIND=r8), DIMENSION(1:2,1:2):: a
REAL(KIND=r8), DIMENSION(1:2):: b
!
h2 = h*h
!
tmp1 = dt/h2
tmp2 = eps2/h2
tmp3 = 6.0_r8*eps2/h2
!
a(1,1) = 1.0_r8
a(2,2) = 1.0_r8
!
! These loops must always be structured this way, i.e., so that i,j,k = 1 is a
! black cell. (The first cell relaxed is always i = 2, j,k = 1, since ipass = 1
! is called first.)
!
DO k = 1, mx(3)
  DO j = 1, mx(2)
    DO i = 1+MODULO(k+j+ipass,2), mx(1), 2
!
      m1 = Mob(0.5_r8*(qnew(i-1,j,k,1)+qnew(i,j,k,1)))
      m2 = Mob(0.5_r8*(qnew(i+1,j,k,1)+qnew(i,j,k,1)))
      m3 = Mob(0.5_r8*(qnew(i,j-1,k,1)+qnew(i,j,k,1)))
      m4 = Mob(0.5_r8*(qnew(i,j+1,k,1)+qnew(i,j,k,1)))
      m5 = Mob(0.5_r8*(qnew(i,j,k-1,1)+qnew(i,j,k,1)))
      m6 = Mob(0.5_r8*(qnew(i,j,k+1,1)+qnew(i,j,k,1)))
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
      det = 1.0_r8-a(1,2)*a(2,1)
!
      qnew(i,j,k,1) = (b(1)-a(1,2)*b(2))/det
      qnew(i,j,k,2) = (b(2)-a(2,1)*b(1))/det
!
    END DO
  END DO
END DO
!
END SUBROUTINE RelaxGrid3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UpdateAuxVcycle3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN OUT):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
!
!aux(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1) = 
!
END SUBROUTINE UpdateAuxVcycle3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator3D(qnew,qold,aux,mx,nrvars,maux,h,xlower) RESULT(opresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: aux
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nrvars):: opresult
!
REAL(KIND=r8):: h2, tmp1, tmp2
!
h2 = h*h
tmp1 = dt/h2
tmp2 = eps2/h2
!
opresult(1:mx(1),1:mx(2),1:mx(3),1) &
  = qnew(1:mx(1),1:mx(2),1:mx(3),1) &
  - tmp1*ULapMob3D(qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1), &
                     qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,2))
!
opresult(1:mx(1),1:mx(2),1:mx(3),2) &
  = qnew(1:mx(1),1:mx(2),1:mx(3),2) &
  - Dfc(qnew(1:mx(1),1:mx(2),1:mx(3),1)) &
  + tmp2*ULap3D(qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1))
!
END FUNCTION Operator3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Source3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(sourceresult)
USE NodeInfoDef
USE GridUtilities, ONLY: ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nrvars):: sourceresult
!
REAL(KIND=r8):: h2, tmp
!
h2 = h*h
tmp = dt/h2
!
sourceresult(1:mx(1),1:mx(2),1:mx(3),1) &
  = qold(1:mx(1),1:mx(2),1:mx(3),1)
!
sourceresult(1:mx(1),1:mx(2),1:mx(3),2) &
  =  &
  - Dfe(qold(1:mx(1),1:mx(2),1:mx(3),1))
!
END FUNCTION Source3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SourceUpdate3D(qnew,qold,aux,ll,mx,nrvars,maux,mbc,h,xlower) &
         RESULT(updateresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv3D, ULap3D
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: maux
INTEGER, INTENT(IN):: mbc
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nrvars):: updateresult
!
updateresult(1:mx(1),1:mx(2),1:mx(3),1) = 0.0_r8
!
updateresult(1:mx(1),1:mx(2),1:mx(3),2) = 0.0_r8
!
END FUNCTION SourceUpdate3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULapMob3D(c,a) RESULT(ulapmobresult)
USE NodeInfoDef
USE GridUtilities, ONLY: UDiv3D
IMPLICIT NONE
!
! Level independent 3D UNDIVIDED laplacian operator for non-constant mobity.
!
REAL(KIND=r8), DIMENSION(0:,0:,0:), INTENT(IN):: c
REAL(KIND=r8), DIMENSION(0:,0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
                         1:SIZE(a,2)-2, &
                         1:SIZE(a,3)-2):: ulapmobresult
!
INTEGER, DIMENSION(1:3):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(a,1)-2,1:SIZE(a,2)-2,1:SIZE(a,3)-2):: f1
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,0:SIZE(a,2)-2,1:SIZE(a,3)-2):: f2
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2,0:SIZE(a,3)-2):: f3
!
mx(1) = SIZE(a,1)-2; mx(2) = SIZE(a,2)-2; mx(3) = SIZE(a,3)-2
!
! Calculate the 3D flux function:
f1(0:mx(1),1:mx(2),1:mx(3)) &
  = Mob(0.5_r8*(c(1:mx(1)+1,1:mx(2),1:mx(3))+c(0:mx(1),1:mx(2),1:mx(3)))) &
  *              (a(1:mx(1)+1,1:mx(2),1:mx(3))-a(0:mx(1),1:mx(2),1:mx(3)))
!
f2(1:mx(1),0:mx(2),1:mx(3)) &
  = Mob(0.5_r8*(c(1:mx(1),1:mx(2)+1,1:mx(3))+c(1:mx(1),0:mx(2),1:mx(3)))) &
  *              (a(1:mx(1),1:mx(2)+1,1:mx(3))-a(1:mx(1),0:mx(2),1:mx(3)))
!
f3(1:mx(1),1:mx(2),0:mx(3)) &
  = Mob(0.5_r8*(c(1:mx(1),1:mx(2),1:mx(3)+1)+c(1:mx(1),1:mx(2),0:mx(3)))) &
  *              (a(1:mx(1),1:mx(2),1:mx(3)+1)-a(1:mx(1),1:mx(2),0:mx(3)))
!
! Calculate the divergence of the flux:
ulapmobresult(1:mx(1),1:mx(2),1:mx(3)) &
    = UDiv3D(f1(0:mx(1),1:mx(2),1:mx(3)), &
             f2(1:mx(1),0:mx(2),1:mx(3)), &
             f3(1:mx(1),1:mx(2),0:mx(3)))
!
END FUNCTION ULapMob3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User supplied physical boundary conditions:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UserBC2D(q,ll,mx,nrvars,mbc,bcnow)
USE NodeinfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN OUT):: q
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: mbc
INTEGER, INTENT(IN):: bcnow
!
INTEGER:: ibc
INTEGER, DIMENSION(1:2):: ul
!
ul(1:2) = mx(1:2)+mbc
!
SELECT CASE(bcnow)
  CASE(1)
!
! Dirichlet on bndry 1:
    DO ibc = 1, mbc
      q(1-ibc,ll:ul(2),:) = 2.0_r8-q(ibc,ll:ul(2),:)
    END DO
!
  CASE(2)
!
! Dirichlet on bndry 2:
    DO ibc = 1, mbc
      q(mx(1)+ibc,ll:ul(2),:) = 2.0_r8-q(mx(1)-ibc+1,ll:ul(2),:)
    END DO
!
  CASE(3)
!
! Dirichlet on bndry 3, variables 3 and 4:
    DO ibc = 1, mbc
      q(ll:ul(1),1-ibc,1) = q(ll:ul(1),ibc,1)
      q(ll:ul(1),1-ibc,2) = q(ll:ul(1),ibc,2)
    END DO
!
  CASE(4)
!
! Dirichlet on bndry 4, variables 3 and 4:
    DO ibc = 1, mbc
      q(ll:ul(1),mx(2)+ibc,1) = q(ll:ul(1),mx(2)-ibc+1,1)
      q(ll:ul(1),mx(2)+ibc,2) = q(ll:ul(1),mx(2)-ibc+1,2)
    END DO          
END SELECT
!
END SUBROUTINE UserBC2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE UserBC3D(q,ll,mx,nrvars,mbc,bcnow)
USE NodeinfoDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN OUT):: q
INTEGER, INTENT(IN):: ll
INTEGER, DIMENSION(1:3), INTENT(IN):: mx
INTEGER, INTENT(IN):: nrvars
INTEGER, INTENT(IN):: mbc
INTEGER, INTENT(IN):: bcnow
!
INTEGER:: ibc
INTEGER, DIMENSION(1:3):: ul
!
ul(1:3) = mx(1:3)+mbc
!
SELECT CASE(bcnow)
  CASE(1)
!
! Dirichlet on bndry 1:
    DO ibc = 1, mbc
      q(1-ibc,ll:ul(2),ll:ul(3),:) = 2.0_r8-q(ibc,ll:ul(2),ll:ul(3),:)
    END DO
!
  CASE(2)
!
! Dirichlet on bndry 2:
    DO ibc = 1, mbc
      q(mx(1)+ibc,ll:ul(2),ll:ul(3),:) = 2.0_r8-q(mx(1)-ibc+1,ll:ul(2),ll:ul(3),:)
    END DO
!
  CASE(3)
!
! Dirichlet on bndry 3:
    DO ibc = 1, mbc
      q(ll:ul(1),1-ibc,ll:ul(3),:) = 2.0_r8-q(ll:ul(1),ibc,ll:ul(3),:)
    END DO
!
  CASE(4)
!
! Dirichlet on bndry 4:
    DO ibc = 1, mbc
      q(ll:ul(1),mx(2)+ibc,ll:ul(3),:) = 2.0_r8-q(ll:ul(1),mx(2)-ibc+1,ll:ul(3),:)
    END DO
!
  CASE(5)
!
! Dirichlet on bndry 5:
    DO ibc = 1, mbc
      q(ll:ul(1),ll:ul(2),1-ibc,:) = 2.0_r8-q(ll:ul(1),ll:ul(2),ibc,:)
    END DO
!
  CASE(6)
!
! Dirichlet on bndry 6:
    DO ibc = 1, mbc
      q(ll:ul(1),ll:ul(2),mx(3)+ibc,:) = 2.0_r8-q(ll:ul(1),ll:ul(2),mx(3)-ibc+1,:)
    END DO
END SELECT
!
END SUBROUTINE UserBC3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Dimensionally Invariant Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Mob(c) RESULT(mobresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: mobresult
!
!IF(c<0.0_r8) THEN
!  mobresult = 1.0E-5_r8
!ELSE IF(c>1.0_r8) THEN
!  mobresult = 1.0E-5_r8
!ELSE
!  mobresult = 16.0_r8*c*c*(1.0_r8-c)*(1.0_r8-c)+1.0E-6_r8
!!  mobresult = 4.0_r8*c*(1.0_r8-c)+1.0E-6_r8
!END IF
!
mobresult = 1.0_r8
!
END FUNCTION Mob
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Ff(c) RESULT(ffresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: ffresult
!
REAL(KIND=r8):: cm1
!
cm1 = c-1.0_r8
!
ffresult = c*c*cm1*cm1*alpha/4.0_r8
!
END FUNCTION Ff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Dfe(c) RESULT(dferesult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: dferesult
!
dferesult = (c-0.5_r8)*alpha/4.0_r8
!
END FUNCTION Dfe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Dfc(c) RESULT(dfcresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: dfcresult
!
REAL(KIND=r8):: cmh
!
cmh = c-0.5_r8
!
dfcresult = cmh*cmh*cmh*alpha
!
END FUNCTION Dfc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION D2fc(c) RESULT(d2fcresult)
USE NodeInfoDef
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: c
REAL(KIND=r8):: d2fcresult
!
REAL(KIND=r8):: cmh
!
cmh = c-0.5_r8
!
d2fcresult = 3.0_r8*cmh*cmh*alpha
!
END FUNCTION D2fc
END MODULE Problem
