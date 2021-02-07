        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 21 10:22:14 2017
        MODULE RHS_H__genmod
          INTERFACE 
            SUBROUTINE RHS_H(Y,NY,NASV,PARMS,NPARMS,DEPS,KRK,NFEV,ERROR)
              INTEGER(KIND=4) :: NPARMS
              INTEGER(KIND=4) :: NASV
              INTEGER(KIND=4) :: NY
              REAL(KIND=8) :: Y(NY)
              REAL(KIND=8) :: PARMS(NPARMS)
              REAL(KIND=8) :: DEPS(6)
              REAL(KIND=8) :: KRK(NY)
              INTEGER(KIND=4) :: NFEV
              INTEGER(KIND=4) :: ERROR
            END SUBROUTINE RHS_H
          END INTERFACE 
        END MODULE RHS_H__genmod
