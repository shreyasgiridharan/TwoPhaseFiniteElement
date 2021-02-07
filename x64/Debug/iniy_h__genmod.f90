        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 21 10:22:14 2017
        MODULE INIY_H__genmod
          INTERFACE 
            SUBROUTINE INIY_H(Y,NYDIM,NASV,NTENS,SIG,QQ)
              INTEGER(KIND=4) :: NTENS
              INTEGER(KIND=4) :: NASV
              INTEGER(KIND=4) :: NYDIM
              REAL(KIND=8) :: Y(NYDIM)
              REAL(KIND=8) :: SIG(NTENS)
              REAL(KIND=8) :: QQ(NASV)
            END SUBROUTINE INIY_H
          END INTERFACE 
        END MODULE INIY_H__genmod
