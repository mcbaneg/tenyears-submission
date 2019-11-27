c  this file holds functions stolen from Molscat that are
c  needed for evaluating 4D angular quantum corrections

      FUNCTION YRR(L1,L2,L,CT1,CT2,DP)
C
C     BISPHERICAL HARMONIC ANGULAR FUNCTIONS FOR TWO DIATOMS
C     CT1, CT2 ARE COS(THETA-1) AND COS(THETA-2), AND
C          DP IS DELTA(PHI), I.E., PHI2-PHI1, IN RADIANS
C     CF. GREEN, JCP 62, 2271 (1975) APPENDIX.
C     N.B. P(L,M;X) THERE IS (2*PI)**-1/2 NORMALIZED P(L,M;X)
C          MOLSCAT PLM(L,M,CT) ROUTINE IS NORMALIZED ON CT, AND
C               PLM(L,0,1.D0)=SQRT((2L+1)/2) .
C          THUS, MUST MULT EACH PLM BY (2*PI)**-1/2
C
C     ODD  L1+L2+L  *NOT* ALLOWED; TRAPPED W/MESSAGE AND STOP
C
C     NEEDS ROUTINES THRJ(XJ1,XJ2,XJ3,XM1,XM2,XM3)
C                    PLM(L,M,COSTH)
C                    PARITY(J)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ODD
! added 13 Nov 2019 to avoid collision with new reserved name
      external parity
      DATA PI/3.14159 26535 89793 D0/
      ODD(I)=2*(I/2)-I.NE.0
C
      IF (ODD(L1+L2+L)) GO TO 9999
C
      XL1=L1
      XL2=L2
      XL=L
C     SQRT(4*PI) FROM Y(L,M,THETA=0), 2*PI FOR TWO PLM'S
      DEN=SQRT(4.D0*PI)*2.D0*PI
      FACT=((2.D0*XL+1.D0)/DEN)*PARITY(L1+L2)
      MTOP=MIN(L1,L2)
      M=0
      XM=0.D0
      SUM=THRJ(XL1,XL2,XL,0.D0,0.D0,0.D0)*PLM(L1,0,CT1)*PLM(L2,0,CT2)
 2000 M=M+1
      IF (M.GT.MTOP) GO TO 3000
      XM=XM+1.D0
      SUM=SUM+2.D0*PARITY(M)*THRJ(XL1,XL2,XL,XM,-XM,0.D0)*
     1             PLM(L1,M,CT1)*PLM(L2,M,CT2)*COS(XM*DP)
      GO TO 2000
 3000 YRR=FACT*SUM
      RETURN
 9999 WRITE(6,699) L1,L2,L
  699 FORMAT('0 YRR *** ERROR.  ODD ARGUMENTS NOT ALLOWED',3I5)
      STOP
      END
      FUNCTION THRJ(F1,F2,F3,G1,G2,G3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SMALL CHANGES 31 JUL 95 (SG)
      SAVE MUNG,X,Y
      PARAMETER (MXIX=302)
      DIMENSION X(MXIX),Y(MXIX)
      DATA MUNG/0/
      IF (MUNG.EQ.21) GO TO 69
      MUNG = 21
      X(1) = 0.D0
      DO 100 I = 1, MXIX-1
      A = I
      X(I+1) = LOG(A) +X(I)
      Y(I+1) = LOG(A)
  100 CONTINUE
   69 IF(F1-ABS(G1)) 1,13,13
   13 IF(F2-ABS(G2))1,14,14
   14 IF(F3-ABS(G3))1,15,15
   15 SUM=F1+F2+F3
      NSUM=SUM+.001D0
      IF(SUM-NSUM)2,2,1
    1 THRJ=0.D0
      RETURN
    2 IF(ABS(G1+G2+G3)-1.D-08)3,3,1
    3 IF(F1+F2-F3)1,4,4
    4 IF(F1+F3-F2)1,5,5
    5 IF(F2+F3-F1)1,6,6
    6 J1=2.D0*F3+2.001D0
      J2=F1+F2-F3+1.001D0
      J3=F1-F2+F3+1.001D0
      J4=-F1+F2+F3+1.001D0
      J5=F1+F2+F3+2.001D0
      J6=F1+G1+1.001D0
      J7=F1-G1+1.001D0
      J8=F2+G2+1.001D0
      J9=F2-G2+1.001D0
      J10=F3+G3+1.001D0
      J11=F3-G3+1.001D0
      IF(J5.GT.MXIX) THEN
        WRITE(6,601) J5,MXIX
  601   FORMAT(' *** DIMENSION ERROR IN THRJ - INDEX.GT.MXIX',2I5)
        STOP
      ENDIF
      R=0.5D0*(Y(J1)+X(J2)+X(J3)+X(J4)-X(J5)
     1+X(J6)+X(J7)+X(J8)+X(J9)+X(J10)+X(J11))
      SUM=0.D0
      F=-1
      KZ=-1
    7 KZ=KZ+1
      F=-F
      J1=KZ+1
      J2=F1+F2-F3-KZ+1.001D0
      IF(J2)20,20,8
    8 J3=F1-G1-KZ+1.001D0
      IF(J3)20,20,9
    9 J4=F2+G2-KZ+1.001D0
      IF(J4)20,20,10
   10 J5=F3-F2+G1+KZ+1.001D0
      IF(J5)7,7,11
   11 J6=F3-F1-G2+KZ+1.001D0
      IF(J6)7,7,12
   12 JMAX=MAX(J1,J2,J3,J4,J5,J6)
      IF(JMAX.GT.MXIX) THEN
        WRITE(6,601) JMAX,MXIX
        STOP
      ENDIF
      S=-(X(J1)+X(J2)+X(J3)+X(J4)+X(J5)+X(J6))
      SUM=SUM+F*EXP(R+S)
      GO TO 7
   20 INT=ABS(F1-F2-G3)+0.0001D0
      VAL=((-1.D0)**INT)*SUM/SQRT(2.D0*F3+1.D0)
      IF(ABS(VAL).LE.1.D-6) VAL=0.D0
      THRJ=VAL
      RETURN
      END

      FUNCTION PLM(LIN,MIN,COSTH)
C
C     COMPUTES NORMALIZED ASSOC. LEGENDRE POLYNOMIALS BY RECURSION.
C     THE VALUES RETURNED ARE NORMALIZED FOR INTEGRATION OVER X
C     (I.E. INTEGRATION OVER COS THETA BUT NOT PHI).
C     NOTE THAT THE NORMALIZATION GIVES
C           PLM(L,0,1)=SQRT(L+0.5)
C           PLM(L,0,X)=SQRT(L+0.5) P(L,X)
C     FOR M.NE.0, THE VALUE RETURNED DIFFERS FROM THE USUAL
C           DEFINITION OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C           (E.G. EDMONDS PAGES 23-24)
C           BY A FACTOR OF (-1)**M*SQRT(L+0.5)*SQRT((L-M)!/(L+M)!)
C     THUS THE SPHERICAL HARMONICS ARE
C          CLM = PLM * EXP(I*M*PHI) / SQRT(L+0.5)
C          YLM = PLM * EXP(I*M*PHI) / SQRT(2*PI)
C     THIS ROUTINE ALWAYS RETURNS THE VALUE FOR ABS(MIN); NOTE THAT
C          FOR MIN.LT.0 THIS VALUE SHOULD BE MULTIPLIED BY PARITY(MIN)
C
C     FUNCTION PM1(LIN,MIN,COSTH)
C     This routine appears to be much more stable for large l, m than
C       the routine from Nerf/ modified according to R.T Pack
C     It was obtained:
C     From: Marie-Lise Dubernet <mld@ipp-garching.mpg.de>
C     Date: Mon, 19 Jun 1995 12:48:11 +0200 (MET DST)
C     Some mods 27-28 June 95 by SG for speed and to accord w/ MOLSCAT 
C     Bugs fixed 21 Sept 95 (SG)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
! added here 13 Nov 2019
      external parity
C     CHECK FOR ABS(COSTH).LE.1.D0 ... 
      IF (ABS(COSTH).GT.1.D0) THEN
        WRITE(6,*) ' *** ILLEGAL ARGUMENT TO PLM. X =',COSTH
        STOP
      ENDIF
C     SAVE ARGUMENTS IN LOCAL VARIABLES
      L=LIN
      M=ABS(MIN)
      X=COSTH
C
C  IF M>L PLM=0 !
      IF(M.GT.L) THEN
        PLM=0.D0
        RETURN
      ENDIF
      LMAX=L
C
      IF (M.GT.0) GO TO 5
C  HERE FOR REGULAR LEGENDRE POLYNOMIALS
      PLM=1.D0
      PM2=0.D0
      XL=0.D0
      DO 2 L=1,LMAX
      XL=XL+1.D0
      PP=((2.D0*XL-1.D0)*X*PLM-(XL-1.D0)*PM2)/XL
      PM2=PLM
    2 PLM=PP
      GO TO 9000
C
C  HERE FOR ALEXANDER-LEGENDRE POLYNOMIALS
C
    5 IMAX=2*M
      RAT=1.D0
      AI=0.D0
      DO 6 I=2,IMAX,2
      AI=AI+2.D0
    6 RAT=RAT*((AI-1.D0)/AI)
C     Y=SIN(THETA)
      Y=SQRT(1.D0-X*X)
      PLM=SQRT(RAT)*(Y**M)
      PM2=0.D0
      LOW=M+1
      XL=LOW-1
      DO 10 L=LOW,LMAX
      XL=XL+1.D0
      AL=DBLE((L+M)*(L-M))
      AL=1.D0/AL
      AL2=(DBLE((L+M-1)*(L-M-1)))*AL
      AL=SQRT(AL)
      AL2=SQRT(AL2)
      PP=(2.D0*XL-1.D0)*X*PLM*AL-PM2*AL2
      PM2=PLM
   10 PLM=PP
      PLM=PLM*PARITY(MIN)
C
C     CONVERT TO MOLSCAT'S IDIOSYNCRATIC NORMALIZATION
9000  PLM=PLM*SQRT(XL+0.5D0)
      RETURN
      END

      FUNCTION PARITY(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARITY=1.D0
      IF((I/2)*2-I.NE.0) PARITY=-1.D0
      RETURN
      END

