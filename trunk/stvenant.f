************************************************************************
*             1-D UNSTEADY NONLINEAR GRADUALLY VARIED FLOW             *
************************************************************************
*       This program solves the non-linear long wave equations         *
*       using the numerical scheme of Koutitas (1988, p. 68)           *
*       as rewritten in FORTRAN by R. Slingerland.                     *
*                                                                      *
************************************************************************
*                The set up constants are:                             *
*   IM     --- Number of cross sections                                *
*   DT,DX  --- Time (s) and space (m) descretisation steps             *
*   C      --- Chezy friction coefficient (m**1/2 s**-1)               *
*   PER,ZO --- Period (s) and amplitude (m) of incoming waves          *
*   NM     --- Number of time steps desired                            *
*   BK     --- Type of downstream boundary: 1 --- free transmission,   *
*                 2 --- full reflection                                *
*   MARK   --- Designates variable to be saved:                        *
*                 1 --- H                                              *
*                 2 --- Q/A                                            *
*                 3 --- Z                                              *
*   NWRITE --- Interval between the timesteps when dependent variables *
*              are written out to a file                               *
*   B(I)   --- Channel width at the Ith cross section (m)              *
*   HO(I)  --- Still water depth (m)                                   *
*                Some internal variables are:                          *
*   Z      --- free surface elevation with respect to the SWL (m)      *
*   H,HN   --- past and present water depths (m)                       *
*   Q,QN   --- past and present discharge values (m**3 s**-1)          *
************************************************************************
************************************************************************
      IMPLICIT NONE
      REAL*8 Q(50),QN(50),H(50),Z(50),R(50),B(50),A(50),HO(50),Z1,L,CO,
     #DT,DX,C,PR,ZO,T,VV,TEMP,SN,PI,G,SUR(100000),SAR(100,1000),Z2
      INTEGER*4 IM,NM,BK,MARK,N,I,II,NC,NUM,NWRITE
      PARAMETER (PI=3.14159,G=9.8)
C-----------------------------------------------------------------------
C    Open the input files
C-----------------------------------------------------------------------
      OPEN(9, FILE='./input.dat')
      READ (9, *) IM, DT, DX, C, PR, ZO, NM, BK, MARK, NWRITE
      READ (9, *) (B(I), I = 1, IM)
      READ (9, *) (HO(I), I = 1, IM)
C-----------------------------------------------------------------------
C  Perform calculations
C-----------------------------------------------------------------------
      II = 0
      DO 1 I = 1, IM
         H(I) = HO(I)
         Q(I) = 0.0
         Z(I) = 0.0
    1 CONTINUE
      N = 0
      T = 0
C-----------------------------------------------------------------------
C  This is the main time loop
C-----------------------------------------------------------------------
    2 CONTINUE
      N = N + 1
      T = T + DT
C-----------------------------------------------------------------------
C  Update water depths and hydraulic geometries
C-----------------------------------------------------------------------
      H(1) = HO(1) + Z(1)
      H(IM) = HO(IM) + Z(IM-1)
      DO 3 I = 2, IM - 1
         H(I) = HO(I) + (Z(I)+Z(I-1))/2.0
    3 CONTINUE
      DO 4 I = 1, IM
         A(I) = B(I)*H(I)
         R(I) = A(I)/(B(I)+2.0*H(I))
    4 CONTINUE
C-----------------------------------------------------------------------
C  Calculate new zeta at node 1, accounting for reflected waves
C-----------------------------------------------------------------------
      CO = SQRT(9.8*H(2))
      L = CO*PR
      IF ((N-1)*DT .GT. DX/CO) THEN
         Z1 = Z(1) - ZO*SIN(2.0*3.14159*(T-DT)/PR)
         Z2 = Z(2) - ZO*SIN(2.0*3.14159*(T-DT)/PR-DX/L)
         Z1 = Z1 + (DT/DX)*CO*(Z2-Z1)
      ELSE
         Z1 = 0
      ENDIF
      Z(1) = ZO*SIN(2.0*3.14159*T/PR) + Z1
C-----------------------------------------------------------------------
C  Solve the continuity equation for the other new zeta(i)
C-----------------------------------------------------------------------
      DO 5 I = 2, IM - 1
         Z(I) = Z(I) - DT/DX*(Q(I+1)-Q(I))/(B(I)+B(I+1))*2.0
    5 CONTINUE
C-----------------------------------------------------------------------
C  Calculate the local energy loss term due to abrupt channel
C  enlargements
C-----------------------------------------------------------------------
      DO 8 I = 2, IM - 1
         VV = 0.0
         IF (B(I+1).GT.B(I-1) .AND. Q(I).GT.0.0) GO TO 6
         IF (B(I+1).LT.B(I-1) .AND. Q(I).LT.0.0) GO TO 6
         GO TO 7
    6    CONTINUE
         VV = (ABS(Q(I+1)/A(I+1))-ABS(Q(I-1)/A(I-1)))**2/4.0/9.8/DX
    7    CONTINUE
         TEMP = Q(I)
         SN = SIGN(1D0,TEMP)
C-----------------------------------------------------------------------
C  Solve the general law of motion for the new discharge, Q(i)
C-----------------------------------------------------------------------
         QN(I) = Q(I) - DT*(Q(I+1)**2/A(I+1)-Q(I-1)**2/A(I-1))/2.0/DX - 
     #      DT*9.8*A(I)*(Z(I)-Z(I-1))/DX - DT*9.8*A(I)*((Q(I)/A(I))**2/C
     #      **2/R(I)+VV)*SN
    8 CONTINUE
C-----------------------------------------------------------------------
C  Define the discharge at section 1 as the discharge at section 2.
C  Define the discharge at section im based on the equation for a moving
C  surge or, if the specified boundary condition is a reflecting wall,
C  set Q at im equal to 0.
C-----------------------------------------------------------------------
      QN(1) = QN(2)
      IF (BK .EQ. 1) THEN
         QN(IM) = Z(IM-1)*SQRT(9.8*B(IM)*A(IM))
      ELSE
         QN(IM) = 0.0
      ENDIF
      DO 9 I = 1, IM
         Q(I) = QN(I)
    9 CONTINUE
C-----------------------------------------------------------------------
C  Store the dependent variable of choice
C-----------------------------------------------------------------------
      NC = NC + 1
      IF (NC .EQ. NWRITE) THEN
         DO 10 I = 1, IM
            II = II + 1
            IF (MARK .EQ. 1) THEN
               SUR(II) = H(I)
            ELSE IF (MARK .EQ. 2) THEN
               SUR(II) = Q(I)/A(I)
            ELSE
               SUR(II) = Z(I)
            ENDIF
   10    CONTINUE
         NC = 1
      ENDIF
      IF (N .LT. NM) GO TO 2
C----------------------------------------------------------------------
C  Write out dependent variable for plotting
C----------------------------------------------------------------------
      NUM = II/IM
      CALL SUB (SUR, SAR, IM, NUM, II)
      STOP 
      END
C----------------------------------------------------------------------
C  Subroutine for the simple purpose of passing the exact dimensions of
C  array VAR
C----------------------------------------------------------------------
      SUBROUTINE SUB(SUR,VAR,IM,NM,II)
      INTEGER*4 IM,NM
      REAL *8 VAR(IM,NM),SUR(II)
      IJ = 0
      DO J = 1, NM
         DO I = 1, IM
            IJ = IJ + 1
            VAR(I,J) = SUR(IJ)
         END DO
      END DO
      OPEN(2, FILE='./out.dat')
      DO J = 1, NM
         WRITE (2, 100) (VAR(I,J), I = 1, IM)
  100 FORMAT(2X,7F11.4)
      END DO
      CLOSE(2)
      RETURN 
      END
