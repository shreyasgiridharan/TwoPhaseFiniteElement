        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 20 15:19:18 2017
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
