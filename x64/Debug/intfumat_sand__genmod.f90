        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 27 12:53:52 2017
        MODULE INTFUMAT_SAND__genmod
          INTERFACE 
            SUBROUTINE INTFUMAT_SAND(DT,PROPS,DEPS,SIG,STATV)
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: PROPS(26)
              REAL(KIND=8), INTENT(IN) :: DEPS(4)
              REAL(KIND=8), INTENT(INOUT) :: SIG(4)
              REAL(KIND=8), INTENT(INOUT) :: STATV(*)
            END SUBROUTINE INTFUMAT_SAND
          END INTERFACE 
        END MODULE INTFUMAT_SAND__genmod
