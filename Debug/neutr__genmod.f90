        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 21 10:19:30 2017
        MODULE NEUTR__genmod
          INTERFACE 
            SUBROUTINE NEUTR(ST,DE,DEL,REV,ITER)
              REAL(KIND=8) :: ST(4)
              REAL(KIND=8) :: DE(4)
              REAL(KIND=8) :: DEL(4,4)
              REAL(KIND=8) :: REV(2)
              INTEGER(KIND=4) :: ITER
            END SUBROUTINE NEUTR
          END INTERFACE 
        END MODULE NEUTR__genmod
