        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 21 10:22:14 2017
        MODULE MOVE_SIG_H__genmod
          INTERFACE 
            SUBROUTINE MOVE_SIG_H(STRESS,NTENS,PORE,SIG)
              INTEGER(KIND=4) :: NTENS
              REAL(KIND=8) :: STRESS(NTENS)
              REAL(KIND=8) :: PORE
              REAL(KIND=8) :: SIG(6)
            END SUBROUTINE MOVE_SIG_H
          END INTERFACE 
        END MODULE MOVE_SIG_H__genmod
