    subroutine IntfUMAT_Sand(dt, props, dEps, Sig, Statv)
    !**********************************************************************
    !
    !   Function : 
    !
    !
    !**********************************************************************
     
    implicit none
        
    double precision, intent(in) :: dt, props(26), dEps(4)
    double precision, intent(inout) :: Sig(4), Statv(*)
    ! local variables
    character*80 C80dum
    double precision :: dtime,props16(16),dEps6(6),Sig6(6),Statv14(14), &
        Ddum,Ddum2(2),Ddum3(3),Ddum33(3,3),Ddum6(6),Ddum66(6,6)

    Ddum=0.d0; Ddum2=0.d0; Ddum3=0.d0; Ddum33=0.d0; Ddum6=0.d0
    Ddum66=0.d0

   props16(1:14)=props(11:24)
   props16(15)=props(26)
   props16(16)=props(25)
   
   props16(15) = 2200000.d0
            
   Sig6(1:2)=Sig(1:2); Sig6(4)=Sig(3); Sig6(3)=Sig(4)
   Sig6(5:6)=0.d0

   dEps6(1:2)=dEps(1:2);dEps6(4)=dEps(3);dEps6(3)=dEps(4)
   dEps6(5:6)=0.d0

   Statv14(1:4)=Statv(1:4)      ! delta tensor
   Statv14(5:6)=0.d0           ! delta tensor 3D
   Statv14(7)=Statv(5)          ! void ratio
   Statv14(8)=Statv(10)             ! not used, pore
   Statv14(9)=(Sig(1)+Sig(2)+Sig(3))/3.d0  ! mean stress
   Statv14(13)=Statv(8)                     ! time step

   dtime=dt
   call umat(Sig6,Statv14,Ddum66,Ddum,Ddum,Ddum,Ddum,Ddum6,Ddum6,&
       Ddum,Ddum6,dEps6,Ddum2,dtime,Ddum,Ddum,Ddum,Ddum,C80dum,&
       3,1,4,14,props16,16,Ddum3,Ddum33,1._8,Ddum,Ddum33,Ddum33,&
       1,1,0,0,0,0)

   Sig(1:2)=Sig6(1:2); Sig(4)=Sig6(3); Sig(3)=Sig6(4)
   Statv(1:4)=Statv14(1:4)! delta tensor
   Statv(5)=Statv14(7)    ! void ratio
   Statv(6)=Statv14(11)   ! mobilized friction angle
   Statv(7)=Statv14(12)   ! normalised length rho
   Statv(8)=Statv14(13)   ! time step
   Statv(9)=Statv14(10)   ! number of evaluation
   Statv(10) = Statv14(8)   !pore

    end subroutine IntfUMAT_Sand

      subroutine Intf_DLL_Plaxis(dt, MatProp, dEpsInc, SigInOut, Statv, PoreP, timein)
    
    !**********************************************************************
    !
    !   Function : 
    !   DLL used for evaluating stresses and updating state variables
    !   attach ubc_sand64.lib (64bit) file to the source
    !   compile the program in x64 platform
    !   attaching ubc_sand.lib or compiling in Win32 prevents execution
    !   
    !   Interface file allows reassigning stresses and strain parameters
    !   The state of state parameters yet to be assessed
    !   
    !**********************************************************************
    
    !**********************************************************************
    !   Last changed :
    !   GS. 09.02.17 - Tested for element test against Plaxis
    !   Results similar in both drained and undrained implementations
    !   GS. 13.02.17 - Tested for a shear strain controlled cyclic 
    !   test. Results similar to Plaxis
    !**********************************************************************
    
    implicit none
    
    ! Global InOut variable declarations
    
    double precision, intent(in) :: dt, MatProp(40), dEpsInc(3), PoreP(3), timein
    double precision, intent(inout) :: SigInOut(4), Statv(*)
    
    ! Local Variables
    
    integer :: IDTask, iMod, IsUndr, iStep, iTer, iEl, Int, ipl, nStat, NonSym,&
        iStrsDep, iTimeDep, iAbort, iTang, iPrjDir(255), iPrjLen
    
    double precision :: X, Y, Z, Time0, dTime, Props(15), Sig0(6), Swp0, StVar0(21),&
        dEps(6), D(6,6), BulkW, Sig(6), Swp, StVar(21)
   
    
    ! The !DEC provides an interface for calling the user_mod function
    ! Tells the linker that the symbol user_mod is a function and is to be 
    ! called. Reference is a guesswork. Unsure as to how the user_mod is built
    
    !DEC$ ATTRIBUTES STDCALL,REFERENCE :: USER_MOD
    
    IDtask = 1  ! Identification of the task
    iMod = 1    ! Fails for iMod>1, assumed dll has only one model
    IsUndr = 1  ! 0 = Drained condition, 1 = Undrained condition
                ! IsUndr = 1 does not yield any result that is different, guessing Plaxis specific
    iStep = 0   ! Current calculation step number
    iTer = 0    ! Current iteration number
    iEl = 1     ! Current element number
    Int = 1     ! Current local stress point number
   
    ! iAbort = 0    ! Parameter forcing calculation to stop
    
    X = 0.d0 ; Y = 0.d0 ; Z = 0.d0  ! Global coordinates of current stress point
    Time0 = timein    ! Time at the start of the current step
    Swp0 = 0.d0     ! Previous excess pore pressure at the stress point
    D(6,6) = 0.d0   ! Effective material stiffness matrix at the stress point
    BulkW = 0.d0   ! Bulk modulus of water at current stress point
    Swp = 0.d0      ! Resulting excess pore pressure at current stress point
    
    Props(1:15) = MatProp(26:40)    ! Material properties reassignment
    
    ! Global to local stress components reassignment
    Sig0(1:2)=SigInOut(1:2); Sig0(4)=SigInOut(3); Sig0(3)=SigInOut(4)
    Sig0(5:6)=0.d0
   
    ! Global to local strain components reassignment
    ! Here, shear component multiplied by two to yield *gamma_xy*
    dEps(1:2)=dEpsInc(1:2);dEps(4)=dEpsInc(3)*2.d0;dEps(3)=0.d0
    dEps(5:6)=0.d0
   
    ! Global to local state variable reassignment
    StVar0(1:21) = Statv(1:21)

    dTime = dt  ! Incremental time step
    
    ! Calculation of stresses and state parameters start here
    ! Do not reassign values until end of this block
    
    ! ******** BLOCK START ******* 
    
    ! Sets up the state parameters
    IDTask = 1
    call user_mod(IDTask, iMod, IsUndr, iStep, iTer,& 
        iEl, Int, X, Y, Z, Time0, dTime, Props, Sig0, Swp0, StVar0,& 
        dEps, D, BulkW, Sig, Swp, StVar, ipl, nStat, NonSym, iStrsDep,&
        iTimeDep, iTang, iPrjDir, iPrjLen, iAbort)
    
    ! Calculates the stress per strain increment and state parameters
    IDTask = 2
    call user_mod(IDTask, iMod, IsUndr, iStep, iTer,& 
        iEl, Int, X, Y, Z, Time0, dTime, Props, Sig0, Swp0, StVar0,& 
        dEps, D, BulkW, Sig, Swp, StVar, ipl, nStat, NonSym, iStrsDep,&
        iTimeDep, iTang, iPrjDir, iPrjLen, iAbort)

    ! ******** BLOCK END ******* 
    
    ! Local to global stress reassignment
    SigInOut(1:2)=Sig(1:2); SigInOut(4)=Sig(3); SigInOut(3)=Sig(4)
    
    ! Local to global state variable reassignment
    Statv(1:21) = StVar(1:21)
    
    end subroutine Intf_DLL_Plaxis
    