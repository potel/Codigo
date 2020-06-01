      SUBROUTINE LOGFAC(L)
	use factorials
C  FACT(I)= LN((I-1)!)
C  THE VALUE OF L CORRESPONDS TO THE DIMENSION OF (FACT(
      FACT(1)=0.
      DO 100 J=2,L
  100 FACT(J)=FACT(J-1)+LOG(J-1D0)
      RETURN
      END

      FUNCTION WIG3J(A,B,C,AL,BE,GA)
	use factorials
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8,intent(IN):: A,B,C,AL,BE,GA
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(AL+BE+GA) 11,10,11
11    WIG3J = 0.0
      RETURN
10    IF(FAIL3(C,A,B)) GOTO 11
      IF(A-ABS(AL)) 11,14,14
14    IF(B-ABS(BE)) 11,15,15
15    IF(C-ABS(GA)) 11,13,13
13    IA = C-B+AL
      IB = C-A-BE
      IF(IA) 20,21,21
21    IF(IB) 20,23,23
23    MIN = 0
      GO TO 24
20    MIN = -IA
      IF(MIN+IB) 25,24,24
25    MIN = -IB
24    IC = A-AL
      ID = B+BE
      IE = A+B-C
      NIN = MIN
      T=PHASE(MIN)
      S = T
30    MIN = MIN+1
      IZ = IC-MIN+1
      IF(IZ) 29,29,26
26    IY = ID-MIN+1
      IF(IY) 29,29,27
27    IX = IE-MIN+1
      IF(IX) 29,29,28
28    TA = real(IX)*real(IY)*real(IZ)
      TB = MIN*real(IA+MIN)*real(IB+MIN)
      T = -T*TA/TB
      S = S+T
      GO TO 30
29    IF(S.EQ.0.0) GO TO 11
      I = B-A+GA
      IH = A+AL
      II = B-BE
      IJ = C+GA
      IK = C-GA
      IL = A+C-B
      IM = B+C-A
      IN = A+B+C+1.0
      XDDD = 0.5*(FACT(IH+1)+FACT(IC+1)+FACT(ID+1)
     $           +FACT(II+1)+FACT(IJ+1)+FACT(IK+1)
     $           +FACT(IE+1)+FACT(IM+1)+FACT(IL+1)-FACT(IN+1))
     $      - ( FACT(IC-NIN+1)+FACT(IA+NIN+1)+FACT(ID-NIN+1)+
     $          FACT(IB+NIN+1)+FACT(NIN+1)+FACT(IE-NIN+1))
      WIG3J = PHASE(I)  *  EXP(XDDD) * S
      RETURN
      END
C     FUNCTION WIG6J
C
C     THIS CALCULATES WIGNER 6-J COEFFICIENTS AS DEFINED IN BRINK AND
C     SATCHLER. IT REQUIRES TWO SUBPROGRAMS, FACT AND COMP... SEE
C     LISTING OF 3-J SYMBOL. IT USES THE CLOSED FORM FOR W-COEFICIENTS
C     DUE TO RACAH ALSO IN BRINK AND SATCHLER. THE SAME RULES APPLY FOR
C     THE INPUT PARAMETERS AS STATED IN THE 3-J PROGRAM.
C
      FUNCTION WIG6J(A,B,C,D,E,F)
	use factorials
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,intent(IN):: A,B,C,D,E,F
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(FAIL3(A,B,C)) GO TO 10
      IF(FAIL3(A,E,F)) GO TO 10
      IF(FAIL3(B,D,F)) GO TO 10
      IF(FAIL3(C,D,E)) GO TO 10
      GO TO 14
   10 WIG6J=0.
      RETURN
   14 IC=A+B-C
      ID=E+D-C
      IE=A+E-F
      IG=B+D-F
      IA=C+F-A-D
      IB=C+F-B-E
      IH=A+B+E+D+1.
      M=MIN(IH,IC,ID,IE,IG)
      IF(M)10,17,17
   17 MUP=MIN(IA,IB,0)
      T=PHASE(M)
      N = M
      S=T
   18 M=M-1
      IF(M+MUP)24,16,16
   16 TA=real(IA+M+1)*real(IB+M+1)*real(IH-M)*real(M+1)
      TB=real(IC-M)*real(ID-M)*real(IE-M)*real(IG-M)
      T=-T*TA/TB
      S=S+T
      GO TO 18
   24 IT=A+B+C+1.
      IU=A+E+F+1.
      IV=B+D+F+1.
      IW=C+D+E+1.
      XD  =  .5*(FACT(1+IC)+FACT(1+IE+IB)+FACT(1+IA+IG)
     1+FACT(1+IE)+FACT(1+IB+IC)+FACT(1+IA+ID)+FACT(1+IG)+FACT(1+IC+IA)
     1+FACT(1+ID+IB)
     2+FACT(1+ID)+FACT(1+IA+IE)+FACT(1+IB+IG)-FACT(1+IT)-FACT(1+IU)
     3-FACT(1+IV)-FACT(1+IW))
     4+FACT(1+IH-N)-FACT(1+N)-FACT(1+IA+N)-FACT(1+IB+N)-FACT(1+IC-N)
     5-FACT(1+ID-N)-FACT(1+IE-N)-FACT(1+IG-N)
      WIG6J=PHASE(IH-1)  *S*EXP(XD)
      RETURN
 99    WIG6J = 1.0
      RETURN
      END
C     FUNCTION WIG9J(A,B,C,D,E,F,G,H,Z)
C
C     THIS CALCULATES 9-J SYMBOLS, OR X-COEFICIENTS AS DEFINED IN BRINK
C     AND SATCHLER. IT USES THE FORMULA FOR 9-J SYMBOLS IN TERMS OF 6-J
C     SYMBOLS. IT THEREFORE NEEDS THE 6-J SUBPROGRAM, AND THE CONDITION
C     ON THE INPUT PARAMETERS IS THE SAME AS FOR THE 3-J PROGRAM#      0
C
       FUNCTION WIG9J(A,B,C,D,E,F,G,H,Z)
       IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,intent(IN):: A,B,C,D,E,F,G,H,Z
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      IF(FAIL3(A,B,C)) GO TO 20
      IF(FAIL3(D,E,F)) GO TO 20
      IF(FAIL3(C,F,Z)) GO TO 20
      IF(FAIL3(A,D,G)) GO TO 20
      IF(FAIL3(B,E,H)) GO TO 20
      IF(FAIL3(G,H,Z)) GO TO 20
      GO TO 26
   20 WIG9J=0.
      RETURN
   26 S=0.
      XA=ABS(A-Z)
      XB=ABS(D-H)
      XC=ABS(B-F)
      X=XA
      IF(X-XB)10,11,11
   10 X=XB
   11 IF(X-XC)12,13,13
   12 X=XC
   13 IF(X-A-Z)14,14,15
   14 IF(X-D-H)16,16,15
   16 IF(X-B-F)17,17,15
   17 S=S+(2.*X+1.)*WIG6J(A,Z,X,H,D,G)*WIG6J(B,F,X,D,H,E)*WIG6J(A,Z,X,F,
     1B,C)
      X=X+1.
      GO TO 13
   15 IF(S-0.)18,20,18
   18 K=2.0*(A+B+D+F+H+Z)
      WIG9J=PHASE(K)*S
      RETURN
      END
      FUNCTION CLEB6(A,AL,B,BE,C,M)
	use factorials
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8,intent(IN):: A,AL,B,BE,C,M
      LOGICAL FAIL3,FRAC
      PHASE(I) = (-1)**I
      FRAC(X) = ABS(X-INT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
      GA = -M
C
      IF(AL+BE+GA) 11,10,11
C11    WIG3J = 0.0
11    CLEB6 = 0.0
      RETURN
10    IF(FAIL3(C,A,B)) GO TO 11
      IF(A-ABS(AL)) 11,14,14
14    IF(B-ABS(BE)) 11,15,15
15    IF(C-ABS(GA)) 11,13,13
13    IA = C-B+AL
      IB = C-A-BE
      IF(IA) 20,21,21
21    IF(IB) 20,23,23
23    MIN = 0
      GO TO 24
20    MIN = -IA
      IF(MIN+IB) 25,24,24
25    MIN = -IB
24    IC = A-AL
      ID = B+BE
      IE = A+B-C
      NIN = MIN
      T=PHASE(MIN)
      S = T
30    MIN = MIN+1
      IZ = IC-MIN+1
      IF(IZ) 29,29,26
26    IY = ID-MIN+1
      IF(IY) 29,29,27
27    IX = IE-MIN+1
      IF(IX) 29,29,28
28    TA = real(IX)*real(IY)*real(IZ)
      TB = MIN*real(IA+MIN)*real(IB+MIN)
      T = -T*TA/TB
      S = S+T
      GO TO 30
29    I = B-A+GA
      IF(S.EQ.0.0) GO TO 11
      IH = A+AL
      II = B-BE
      IJ = C+GA
      IK = C-GA
      IL = A+C-B
      IM = B+C-A
      IN = A+B+C+1.0
      XDDD = 0.5*(FACT(IH+1)+FACT(IC+1)+FACT(ID+1)
     $           +FACT(II+1)+FACT(IJ+1)+FACT(IK+1)
     $           +FACT(IE+1)+FACT(IM+1)+FACT(IL+1)-FACT(IN+1))
     $      - ( FACT(IC-NIN+1)+FACT(IA+NIN+1)+FACT(ID-NIN+1)+
     $          FACT(IB+NIN+1)+FACT(NIN+1)+FACT(IE-NIN+1))
C     WIG3J = (-1.0)**I *  EXP(XDDD) * S
      CLEB6 = SQRT(2*C+1.)*EXP(XDDD) * S
      RETURN
      END
      FUNCTION RACAH(A,B,C,D,E,F)
      IMPLICIT REAL*8(A-H,O-Z)
	real*8,intent(IN):: A,B,C,D,E,F
      PHASE(I) = (-1)**I
      Z = ABS(A+B+C+D)
      I = Z + 0.5
      RACAH = PHASE(I) * WIG6J(A,B,E,D,C,F)
      RETURN
      END