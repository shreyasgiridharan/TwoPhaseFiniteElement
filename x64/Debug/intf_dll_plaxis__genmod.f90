        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 27 12:53:51 2017
        MODULE INTF_DLL_PLAXIS__genmod
          INTERFACE 
            SUBROUTINE INTF_DLL_PLAXIS(DT,MATPROP,DEPSINC,SIGINOUT,STATV&
     &,POREP,TIMEIN)
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: MATPROP(40)
              REAL(KIND=8), INTENT(IN) :: DEPSINC(3)
              REAL(KIND=8), INTENT(INOUT) :: SIGINOUT(4)
              REAL(KIND=8), INTENT(INOUT) :: STATV(*)
              REAL(KIND=8), INTENT(IN) :: POREP(3)
              REAL(KIND=8), INTENT(IN) :: TIMEIN
            END SUBROUTINE INTF_DLL_PLAXIS
          END INTERFACE 
        END MODULE INTF_DLL_PLAXIS__genmod
