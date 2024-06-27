
      SUBROUTINE NEUTR (ST,DE,DEL,REV,ITER)

C****************************************************************
C     NEUTRAL LOADING SURFACES
C     FROM PROF. STOLLE
C     TRYING TO MODEL SAND UNDER UNDRAINED CONDITIONS
C     FOR MORE INFORMATION OF THEORY, REFER TO PROF. STOLLE NOTE 
C****************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ST(4),DE(4),DEL(4,4),REV(2)
      DIMENSION DS(4),DSE(4),SS(4),GV(4),GD(4),DEP(4,4)
      COMMON/PROP/EMOD,ANU,E,BULKW,PORE,RLAM,RC,RF,ACON,G0,B0

      SPHI = 3.D0*RF/(6.D0+RF)
      CALL INVAR (SPHI,ST,V1,V2,G,TH)
      IF(ITER.EQ.1) THEN
      REV(1)=V2/(DSQRT(2.D0)/3.D0*V1*G)
      REV(2)=(V1*V1+4.5D0*(V2/(G*RF))**2)/(2.D0*V1)
      ENDIF

      IF ((V1.LT.0.D0).OR.(V2.LT.0.1D-3)) THEN
      WRITE(6,*) 'V1,V1',V1,V2
      STOP
      ENDIF

      GTH=(3.D0-SPHI)*(DSQRT(3.D0)*DSIN(TH)+SPHI*DCOS(TH))
      GTH=GTH/(DSQRT(3.D0)*DCOS(TH)-SPHI*DSIN(TH))**2

C  COMPUTES ELASTIC STRESS INCREMENT

      DO 10 I=1,4
      DSE(I)=.0D0
      DO 10 J=1,4
  10  DSE(I)=DSE(I)+DEL(I,J)*DE(J)

C  CHECKS FOR REVERSE LOADING

      Q0=DSQRT(2.D0)/3.D0*REV(1)*G*V1
      P0=2.D0*REV(2)/(1.D0+4.5D0*(V2/(V1*G*RF))**2)
      LQ=0
      LV=0
      IF (V2.LT.0.999D0*Q0) LQ=1
      IF (V1.LT.0.999D0*P0) LV=1
      FV1=.0D0
      FV2=1.D0
      FV3=.0D0
      CALL GRAD(V1,V2,FV1,FV2,FV3,ST,GD)
      A=(V1*V1+4.5D0*(V2/(G*RF))**2)/(2.D0*V1)
      FV1=2.D0*(V1-A)
      FV2=9.D0*V2/(G*RF)**2
      FV3=-9.D0*V2*V2/(RF*RF*G*G*G)*GTH
      CALL GRAD (V1,V2,FV1,FV2,FV3,ST,GV)
      FD=.0D0
      DO 20 I=1,4
  20  FD=FD+DSE(I)*GD(I)
      IF ((LQ.EQ.0).AND.(FD.LT.0.D0)) LQ=1
      FV=.0D0
      DO 30 I=1,4
  30  FV=FV+DSE(I)*GV(I)
      IF ((LV.EQ.0).AND.(FV.LT.0.D0)) LV=1

C  CALCULATES HD

      C1=RF*G*V1
      IF(LQ.EQ.1) GOTO 40
      HQ=C1*ACON/(DSQRT(2.D0/3.D0)*C1-DSQRT(3.D0)*V2)**2
      GOTO 50
  40  CONTINUE

c  STRESS REVERSAL
      G1=(3.D0-SPHI)/(2.D0*(DSQRT(3.D0)*DCOS(-TH)-SPHI*DSIN(-TH)))
      QD=DSQRT(2.D0)/3.D0*REV(1)*G1*V1
      D=V2+QD
      D1=Q0+QD
      R=D/D1
      IF(FD.LT.0.D0) R=(1.D0-D/D1)
      IF(R.LT.0.D0) R=0.D0
      HQ1=C1*ACON/(DSQRT(2.D0/3.D0)*C1-DSQRT(3.D0)*Q0)**2
      HQ=HQ1*R**G0
  50  CONTINUE

C  CALCULATES HV

      GVM=FV1/DSQRT(FV1*FV1+FV2*FV2)
      IF (LV.EQ.1) GOTO 60
      HV=RLAM/((1.D0+E)*V1*GVM)*(RC-V2/(DSQRT(2.D0)/3.D0*V1*G))/RC
      GOTO 100
  60  CONTINUE
C  STRESS REVERSAL
      D=V1
      D1=P0
      R=(1.D0-D/D1)
      IF(FV.GT.0.D0) R=D/D1
      IF(R.LT.0.D0) R=0.D0
      HV1=RLAM/((1.D0+E)*P0*GVM)*(RC-V2/(DSQRT(2.D0)/3.D0*V1*G))/RC
      HV=HV1*DEXP(D1/D-1.D0)*R**B0
C  REFLECTION
      DO 70 I=1,4
  70  SS(I)=ST(I)+DSE(I)
      CALL INVAR (SPHI,SS,VV1,VV2,GG,THH)
      FCT=(VV2-V2)-(V2/V1)*(VV1-V1)
      IF (FCT.GT.0.D0) GOTO 100
      FCT=DSQRT(ST(1)*ST(1)+ST(2)*ST(2)+2.D0*ST(3)*ST(3)+ST(4)*ST(4))
      DO 80 I=1,4
  80  SS(I)=ST(I)/FCT
      FCT=SS(1)*GV(1)+SS(2)*GV(2)+2.D0*SS(3)*GV(3)+SS(4)*GV(4)
      DO 90 I=1,4
  90  GV(I)=2.D0*FCT*SS(I)-GV(I)
 100  CONTINUE

C  CALCULATES STRESS RATES

      DO 110 I=1,4
      IF(I.EQ.3) GOTO 110
      SS(I)=ST(I)-(ST(1)+ST(2)+ST(4))/3.D0
 110  CONTINUE
      SS(3)=ST(3)/2.D0
      PO=DEL(1,2)/(2.D0*(DEL(1,2)+DEL(3,3)))
      YO=2.D0*DEL(3,3)*(1.D0+PO)
      DO 120 I=1,4
      DO 120 J=1,4
      F1=1.D0
      F2=1.D0
      IF (I.NE.J) F1=-PO
      IF((I.EQ.3).OR.(J.EQ.3)) F1=.0D0
      IF((I.EQ.3).AND.(J.EQ.3)) F1=1.D0+PO
      IF(I.EQ.3) F2=0.D0
 120  DEP(I,J)=F1/YO-F2*HV*GV(J)/3.D0+HQ/V2*SS(I)*GD(J)
      CALL INVERT(DEP,4)
      DO 130 I=1,4
      DS(I)=.0D0
      DO 130 J=1,4
 130  DS(I)=DS(I)+DEP(I,J)*DE(J)

C  CALCULATES TOTAL STRESS AND UPDATES STRESS REVERSAL SURFACE

      DO 140 I=1,4
 140  SS(I)=ST(I)+DS(I)
      CALL INVAR (SPHI,SS,V1,V2,G,TH)
      FR=V2/(DSQRT(2.D0)/3.D0*V1*G)
      IF(FR.GT.RF) GOTO 160
      DO 150 I=1,4
 150  ST(I)=SS(I)
 160  CONTINUE
C
      DNV=.0D0
      DND=.0D0
      DO 170 I=1,4
      DNV=DNV+GV(I)*DS(I)
 170  DND=DND+GD(I)*DS(I)
      DEV=HV*DNV
      DEQ=DSQRT(2.D0/3.D0)*HQ*DND
      !IF (FR.GT.RF) CALL INVAR (SPHI,V1,V2,G,TH) ! Here, the S being passed in INVAR function, is not given. Added in the next function
      IF (FR.GT.RF) CALL INVAR (SPHI,SS,V1,V2,G,TH)
      IF(LQ.EQ.0) REV(1)=V2/(DSQRT(2.D0)/3.D0*V1*G)
      IF(LV.EQ.0) REV(2)=(V1*V1+4.5D0*(V2/(G*RF))**2)/(2.D0*V1)
      IF(LQ.EQ.1) REV(1)=REV(1)+(RF-REV(1))**2/(ACON*RF)*DEQ
      IF(LV.EQ.1) REV(2)=REV(2)+(1.D0+E)/RLAM*DEV*REV(2)
      RETURN
      END

C.....................................................................

      SUBROUTINE GRAD (V1,V2,FV1,FV2,FV3,ST,GF)
C
C***********************************************
C     CALCULATES GRADIENT TENSOR
C***********************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ST(4),GF(4)
      DIMENSION A1(4),A2(4),A3(4),SD(4)
C
      DO 5 I=1,4
   5  SD(I)=ST(I)-(ST(1)+ST(2)+ST(4))/3.D0
      SD(3)=ST(3)
      V3=SD(4)*(SD(4)*SD(4)-.5D0*V2*V2)
C
      DO 10 I=1,4
      A1(I)=-1.D0/DSQRT(3.D0)
      A2(I)=.0D0
      A3(I)=.0D0
  10  CONTINUE
      A1(3)=.0D0
C
      IF(V2.LT.0.1D-3) GOTO 40
      DO 20 I=1,4
  20  A2(I)=SD(I)/V2
C
      SINT3=-3.D0*DSQRT(6.D0)*V3/(V2**3)
      IF (SINT3.LT.-1.D0) SINT3=-1.D0
      IF (SINT3.GT. 1.D0) SINT3=1.D0
      COST3=DSQRT(1.D0-SINT3*SINT3)
      IF((DABS(SINT3)-1.D0).LT.0.1D-3) GOTO 40
      V2S=V2*V2/6.D0
      A3(1)=SD(2)*SD(4)+V2S
      A3(2)=SD(4)*SD(1)+V2S
      A3(3)=SD(3)*SD(4)
      A3(4)=SD(1)*SD(2)-SD(3)*SD(3)+V2S
      FACT=DSQRT(6.D0)/(V2**3*COST3)
      DO 30 I=1,4
      A3(I)=FACT*(3.D0*V3*A2(I)/V2-A3(I))
  30  CONTINUE
  40  CONTINUE
      DO 50 I=1,4
  50  GF(I)=FV1*A1(I)+FV2*A2(I)+FV3*A3(I)
      GN=DSQRT(GF(1)*GF(1)+GF(2)*GF(2)+2.D0*GF(3)*GF(3)+GF(4)*GF(4))
      DO 60 I=1,4
  60  GF(I)=GF(I)/GN
      RETURN
      END

C.....................................................................

      SUBROUTINE INVAR (SPHI,S,V1,V2,G,TH)

C*****************************************
C     CALCULATES STRESS INVARIANTS
C*****************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(4),SD(4)

      V1=-(S(1)+S(2)+S(4))/DSQRT(3.D0)
      DO 10 I=1,4
      IF(I.EQ.3) GOTO 10
      SD(I)=S(I)-(S(1)+S(2)+S(4))/3.D0
  10  CONTINUE
      V2=(SD(1)*SD(1)+SD(2)*SD(2)+2.D0*S(3)*S(3)+SD(4)*SD(4))
      V3=SD(4)*(SD(4)*SD(4)-0.5D0*V2)
      V2=DSQRT(V2)
      IF(V2.LT.0.1D-3) GOTO 20
      SINT3=-3.D0*DSQRT(6.D0)*V3/(V2**3)
      GOTO 30
  20  SINT3=.0D0
  30  CONTINUE
      IF(SINT3.LT.-1.D0) SINT3=-1.D0
      IF(SINT3.GT.1.D0)  SINT3=1.D0
      TH=DASIN(SINT3)/3.D0
      G=(3.D0-SPHI)/(2.D0*(DSQRT(3.D0)*DCOS(TH)-SPHI*DSIN(TH)))
      RETURN
      END
	  
C.....................................................................

      SUBROUTINE INVERT (A,N)

C*****************************************
C     COMPUTES THE INVERSE MATRIX OF A
C*****************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),K(25)

      DET=1.D0
      DO 5 I=1,N
   5  K(I)=I
      DO 80 II=1,N
      DO 10 I=II,N
      PIV=A(I,II)
      IF (DABS(PIV).GT.1.D-12) GOTO 20
  10  CONTINUE
      WRITE(6,*)'MATRIX IS SINGULAR'
      STOP
  20  DET=DET*PIV
      IF(I.EQ.II) GOTO 40
      I1=K(II)
      K(II)=K(I)
      K(I)=I1
      DO 30 J=1,N
      C=A(I,J)
      A(I,J)=A(II,J)
  30  A(II,J)=C
      DET=-DET
  40  C=1.D0/PIV
      A(II,II)=1.D0
      DO 50 J=1,N
  50  A(II,J)=A(II,J)*C
      DO 70 I=1,N
      IF(I.EQ.II) GOTO 70
      C=A(I,II)
      A(I,II)=.0D0
      DO 60 J=1,N
  60  A(I,J)=A(I,J)-C*A(II,J)
  70  CONTINUE
  80  CONTINUE
      DO 120 J=1,N
      DO 90 J1=J,N
      JJ=K(J1)
      IF(JJ.EQ.J) GOTO 100
  90  CONTINUE
 100  IF (J.EQ.J1) GOTO 120
      K(J1)=K(J)
      DO 110 I=1,N
      C=A(I,J)
      A(I,J)=A(I,J1)
 110  A(I,J1)=C
 120  CONTINUE
      RETURN
      END
