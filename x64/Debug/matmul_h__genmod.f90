        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 21 10:22:14 2017
        MODULE MATMUL_H__genmod
          INTERFACE 
            SUBROUTINE MATMUL_H(A,B,C,L,M,N)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: L
              REAL(KIND=8) :: A(L,M)
              REAL(KIND=8) :: B(M,N)
              REAL(KIND=8) :: C(L,N)
            END SUBROUTINE MATMUL_H
          END INTERFACE 
        END MODULE MATMUL_H__genmod
