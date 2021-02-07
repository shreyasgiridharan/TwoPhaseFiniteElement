        ! IF(COH(I).LT.1.D-08) COH(I) = 1.D-08
        ! PHI(I)  = SIN(PI*PHI(I)/180.D0)
        ! PSI(I)  = SIN(PI*PSI(I)/180.D0)
        ! COH(I)  = COH(I)*COS(PI*PHI(I)/180.D0)
        !
        !
        !
        !
	       !
        !
        !ST(1:3) = ST(1:3) + DST(1:3)
        !ST(4) = ST(4) + anv(is)*(DST(1)+DST(2))
        !
        !call MOHRC(PHI(IS),PSI(IS),COH(IS),EMOD(IS),ANV(IS),ST,NPLAST,fbar(ip))
!***********************************************************************************
!     MOHR Coulomb Plaxis
!
!***********************************************************************************
      SUBROUTINE MOHRC(SPHI,SPSI,COHS,EMOD,ANV,ST,NPLAST,fbar,plastind)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ST(4)

      LITER = 0
      GMOD = 0.5d0*EMOD/(1.d0+ANV)

      SXE = 0.5 * (ST(1)-ST(2))
      SSTAR = 0.5 * (ST(1)+ST(2))
      TSTAR = DSQRT(SXE**2+ST(3)**2)
      SINALF=0.0
      COSALF=1.0
      IF (TSTAR > 0.0) THEN
        SINALF = ST(3)/TSTAR
        COSALF = SXE/TSTAR
      ENDIF
      SIG1 = SSTAR - TSTAR
      SIG2 = ST(4)
      SIG3 = SSTAR + TSTAR

      SIG1T  = SIG1
      SIG2T  = SIG2
      SIG3T  = SIG3
      ISIGZZ = 2
      IF (SIG2T > SIG3T) THEN
        ISIGZZ = 3
        SIG3 = SIG2T
        SIG2 = SIG3T
      ELSE
        IF (SIG2T < SIG1T) THEN
          ISIGZZ = 1
          SIG1 = SIG2T
          SIG2 = SIG1T
        END IF
      END IF

      F21 = 0.5*(SIG2-SIG1) + 0.5*(SIG2+SIG1)*SPHI - COHS
      F32 = 0.5*(SIG3-SIG2) + 0.5*(SIG3+SIG2)*SPHI - COHS
      F31 = 0.5*(SIG3-SIG1) + 0.5*(SIG3+SIG1)*SPHI - COHS
      fbar =  -0.5*(SIG3-SIG1)/(0.5*(SIG3+SIG1)*SPHI - COHS)

      IF(F31 > 0.D0) THEN
      nplast = nplast + 1
      LITER = 1


C         *** GAUSS INTEGRATION POINTS IN PLASTIC STATE ***
C         *** STRESSES ARE BROUGHT BACK TO THE YIELD SURFACE ***

      VNUI1 = 1.0D0 - 2.0D0*ANV
      VNUI2 = 1.0D0 + 2.0D0*ANV
      VNUIQ = VNUI2 / VNUI1

      DUM = SPSI / VNUI1
      PSIMIN = GMOD*(-1.0D0+DUM)
      PSIMET = GMOD*( 1.0D0+DUM)
      PSINU = 2.0D0*GMOD*ANV*DUM
C
      HA= ( 1.0D0 - SPHI + SPSI - SPHI * SPSI) * SIG1 +
     *    (-2.0D0 - 2.0/VNUI1*SPHI*SPSI) * SIG2 +
     *    ( 1.0D0 - SPHI - SPSI + VNUIQ * SPHI * SPSI) * SIG3 +
     *    2.0*(1.0D0 + SPSI) * COHS

      HB=-( 1.0D0 + SPHI + SPSI + VNUIQ * SPHI * SPSI) * SIG1 -
     *    (-2.0D0 - 2.0/VNUI1*SPHI*SPSI) * SIG2 -
     *    ( 1.0D0 + SPHI - SPSI - SPHI * SPSI) * SIG3 +
     *    2.0*(1.0D0 - SPSI) * COHS

C
      IF (F31 > 0.0) THEN
        IF (HA < 0.0) THEN
          IAREA=3
        ELSE
          IF (HB < 0.0) THEN
            IAREA=1
          ELSE
            IAREA=2
          ENDIF
        ENDIF
      ELSE
        WRITE(*,*) 'WARNING: NO ACTIVE RETURN AREA'
      ENDIF
      DSP1=0.
      DSP2=0.
      DSP3=0.

C         *** EXTENSION POINT ***

      IF (IAREA.EQ.1) THEN
        A11 = GMOD*(1.0D0+SPHI*SPSI/VNUI1)
        A12 = 0.5*GMOD*(1.0D0+SPHI+SPSI+VNUIQ*SPHI*SPSI)
        DETER  = A11*A11-A12*A12
        RLAM31 = (F31*A11-F32*A12) / DETER
        RLAM32 = (F32*A11-F31*A12) / DETER
        RLAM21 = 0.0D0
        DSP1 = RLAM31*PSIMIN + RLAM32*PSINU
        DSP2 = RLAM31*PSINU  + RLAM32*PSIMIN
        DSP3 = RLAM31*PSIMET + RLAM32*PSIMET
        plastind = 2.0
      ENDIF

C         *** REGULAR YIELD SURFACE ***

      IF (IAREA.EQ.2) THEN
        A11 = GMOD*(1.0D0+SPHI*SPSI/VNUI1)
        RLAM31 = F31 / A11
        RLAM21 = 0.0D0
        RLAM32 = 0.0D0
        DSP1 = RLAM31*PSIMIN
        DSP2 = RLAM31*PSINU
        DSP3 = RLAM31*PSIMET
        plastind = 1.0
      END IF

C         *** COMPRESSION POINT ***

      IF (IAREA.EQ.3) THEN
        A11 = GMOD*(1.0D0+SPHI*SPSI/VNUI1)
        A12 = 0.5*GMOD*(1.0D0-SPHI-SPSI+VNUIQ*SPHI*SPSI)
        DETER  = A11*A11-A12*A12
        RLAM31 = (F31*A11-F21*A12) / DETER
        RLAM21 = (F21*A11-F31*A12) / DETER
        RLAM32 = 0.0D0
        DSP1 = RLAM31*PSIMIN + RLAM21*PSIMIN
        DSP2 = RLAM31*PSINU  + RLAM821*PSIMET
        DSP3 = RLAM31*PSIMET + RLAM21*PSINU
        plastind = 1.0
      END IF

C         *** TEST ON APEX POINTS ***

      IF (SIG3 > 0.D0) NTEN = NTEN + 1
      IF (SPHI > 1.0E-6) THEN
        SFAIL = COHS/SPHI
        HC = SIG1+SIG2+SIG3-(DSP1+DSP2+DSP3)-3*SFAIL
        IF (HC.GT.0.0) THEN
          IAPEX=1
          DSP1 = SIG1-SFAIL
          DSP2 = SIG2-SFAIL
          DSP3 = SIG3-SFAIL
          plastind = 2.0
        END IF
      END IF

C         *** COMPUTATION OF CARTESIAN STRESS COMPONENTS ***

      IF (ISIGZZ.EQ.1) THEN
        DTSTAR = 0.5D0 * (DSP3-DSP2)
        DSSTAR = 0.5D0 * (DSP3+DSP2)
        ST(4) = ST(4) - DSP1
        DSPZZ = DSP1
      END IF
      IF (ISIGZZ.EQ.2) THEN
        DTSTAR = 0.5D0 * (DSP3-DSP1)
        DSSTAR = 0.5D0 * (DSP3+DSP1)
        ST(4) = ST(4) - DSP2
        DSPZZ = DSP2
      END IF
      IF (ISIGZZ.EQ.3) THEN
        DTSTAR = 0.5D0 * (DSP2-DSP1)
        DSSTAR = 0.5D0 * (DSP2+DSP1)
        ST(4) = ST(4) - DSP3
        DSPZZ = DSP3
      END IF

      DSPXX = (DSSTAR+DTSTAR*COSALF)
      DSPYY = (DSSTAR-DTSTAR*COSALF)
      DSPXY = DTSTAR*SINALF
      ST(1) = ST(1) - DSPXX
      ST(2) = ST(2) - DSPYY
      ST(3) = ST(3) - DSPXY

C------- Calculate TAUMAX for checking accuracy on plastic point
      TAUMAX=0.5*(SIG3-DSP3-SIG1+DSP1)
      IF (TAUMAX < -1E-6) WRITE(*,*)' TAUMAX IS NEGATIVE!!!'
	IF (TAUMAX < -1E-6) stop
      ENDIF
      return
      END SUBROUTINE 
