        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 21 09:22:29 2017
        MODULE NEUTRAL_LIQUEFACTION__genmod
          INTERFACE 
            SUBROUTINE NEUTRAL_LIQUEFACTION(SIG,DEPS,DEL_MAT,REV_VECT,E,&
     &ANU)
              REAL(KIND=8), INTENT(INOUT) :: SIG(4)
              REAL(KIND=8), INTENT(IN) :: DEPS(4)
              REAL(KIND=8), INTENT(INOUT) :: DEL_MAT(4,4)
              REAL(KIND=8), INTENT(INOUT) :: REV_VECT(2)
              REAL(KIND=8), INTENT(IN) :: E
              REAL(KIND=8), INTENT(IN) :: ANU
            END SUBROUTINE NEUTRAL_LIQUEFACTION
          END INTERFACE 
        END MODULE NEUTRAL_LIQUEFACTION__genmod
