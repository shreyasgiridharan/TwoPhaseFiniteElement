        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 20 15:19:18 2017
        MODULE SOLOUT_H__genmod
          INTERFACE 
            SUBROUTINE SOLOUT_H(STRESS,NTENS,ASV,NASV,DDSDDE,Y,NYDIM,   &
     &PORE,DEPSV_NP1,PARMS,NPARMS,DD)
              INTEGER(KIND=4) :: NPARMS
              INTEGER(KIND=4) :: NYDIM
              INTEGER(KIND=4) :: NASV
              INTEGER(KIND=4) :: NTENS
              REAL(KIND=8) :: STRESS(NTENS)
              REAL(KIND=8) :: ASV(NASV)
              REAL(KIND=8) :: DDSDDE(NTENS,NTENS)
              REAL(KIND=8) :: Y(NYDIM)
              REAL(KIND=8) :: PORE
              REAL(KIND=8) :: DEPSV_NP1
              REAL(KIND=8) :: PARMS(NPARMS)
              REAL(KIND=8) :: DD(6,6)
            END SUBROUTINE SOLOUT_H
          END INTERFACE 
        END MODULE SOLOUT_H__genmod
