    module modulefem
!*********************************************************************
!    Function:
!**********************************************************************
    implicit none
    logical, dimension(:), allocatable :: IsFixDof
    integer :: ndtn = 0, ITime = 0
    integer :: NDivX = 0, NDivY = 0
    integer :: NEl = 0, NNod = 0
    integer, dimension(:,:), allocatable :: ICon
    double precision :: gx = 0.D0, gy = 0.D0, rhoD, rhoS, rhoSat, poro
    double precision :: Meshdx = 0.D0, Meshdy = 0.D0, rhoW = 1.d0, kappa
    double precision, dimension(:), allocatable :: GrvF, InrF, MasS, Area, Areai, ExtF, v0, V, dis, DisW
    double precision, dimension(:), allocatable :: MasW, MasbW, InrFW, ExtFW, DragF, GrvFW, ExtP, W, w0
    double precision, dimension(:,:), allocatable :: NodCo, NodCoAvg
    double precision, dimension(:,:,:), allocatable :: Sig0, SigG, SigE, EpsG, F, HS, PoreP, StateVar, StatevarPRINT
    double precision, dimension(:,:,:,:), allocatable :: B, BNEW, TEMPB
    
    double precision, dimension(:,:,:,:), allocatable :: DEL_Mat
    double precision, dimension(:,:,:), allocatable :: Rev_Vect, FBARNEW, PHIMOB
    integer :: Iteration = 0, NPLAST = 0
    double precision :: MatProp(40), dt = 0.D0, FBAR, PLASTIND, RES = 0.D0
    integer, parameter :: DATUnit = 2
    integer, parameter :: LOGUnit = 2
    integer, parameter :: MSHUnit = 3
    integer, parameter :: RESUnit = 4
    integer, parameter :: OUTUnit = 5
    logical :: TERMINATION = .false.
    
    contains

    subroutine solve()
    !**********************************************************************
    !    Function:
    !**********************************************************************

    implicit none
    integer :: IIter, IStep = 0, IPTime
    integer :: iprint, imeshwrite, notoprint, notomesh
    real :: start, finish
    integer :: clock_rate, clock_Max, st, en
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values
    
    notoprint = 1000
    notomesh = 100
    
    iprint = ndtn / notoprint
    imeshwrite = ndtn / notomesh 
    
    open(LOGUnit, file = 'output.log')
    open(OUTUnit, file = 'proglog.log')
    
    call initial()
    call WtMSH(MSHUnit)

    call date_and_time(date,time,zone,values)
    call date_and_time(VALUES=values)


    call CPU_TIME(start)
    call SYSTEM_CLOCK ( st, clock_rate, clock_max )

    do ITime = 1, ndtn ! physical time scale

            call Map2nod()
            call update()
       
    
       !if(itime*dt.gt.10.d0) then     
       if(itime.eq.1 .or. (itime/iprint)*iprint == itime .or. itime == ndtn) then
           
           !write (*,*) ITime, ITime*dt, -SigE(2,1,1) + SigE(1,1,1)
           write (*,*) ITime, ITime*dt, SigE(3,1,1)
           
           
           !write (LOGUnit,*) SigE(2,1,1),',', PoreP(2,1,1), ',', EpsG(2,1,1)
           !write (LOGUnit,*) SigE(1,1,1) - SigE(2,1,1),',', PoreP(2,1,1), ',', EpsG(2,1,1)
           !write (LOGUnit,*)  -1*EpsG(2,1,1),',',SigE(1,1,1) - SigE(2,1,1),',', PoreP(2,1,1)
           !write (LOGUnit,*)  EpsG(3,1,1)*2.D0,',',SigE(3,1,1),',', (SigE(1,1,1) + SigE(2,1,1))/2.d0
           write (LOGUnit,*)  EpsG(3,1,1)*2.D0,',',SigE(3,1,1),',', PoreP(2,4,1)
           !write (LOGUnit,*)  PoreP(2,1,1), ',', PoreP(2,2,1), ',', PoreP(2,3,1) 
           !write (LOGUnit,*)  PoreP(2,4,1)
           !write (LOGUnit,*)  EpsG(3,1,1),',',SigE(3,1,1),',', (SigE(1,1,1) + SigE(2,1,1))/2.d0
           !write (LOGUnit,"(17(E10.4,x))")  EpsG(3,1,1),PoreP(2,1,1), SigE(3,1,1), (SigE(1,1,1) + SigE(2,1,1))/2.d0
           !write (LOGUnit,"(17(E10.4,x))") V(1), V(2), v(3), v(4), v(5), v(6), v(7), v(8)
           !write (LOGUnit,"(17(E10.4,x))") W(1)-V(1), W(2)-V(2), W(3)-V(3), W(4)-V(4), W(5)-V(5), W(6)-V(6), W(7)-V(7), W(8)-V(8)
           !write (LOGUnit,"(17(E10.4,x))") w(1), w(2), w(3), w(4), W(5), w(6), w(7), w(8)
           
           !write (LOGUnit,"(17(E10.4,x))") SigE(1,1,1), SigE(2,1,1), SigE(3,1,1), PoreP(2,1,1), EpsG(2,1,1)
       endif
       !endif
       
      
       if(itime.eq.1 .or. (itime/imeshwrite)*imeshwrite == itime .or. itime == ndtn) then
           IStep = IStep + 1
           call WtRES(IStep)
       endif

       if(itime.eq.ndtn) TERMINATION = .true.          
    enddo

    call CPU_TIME(finish)
    call SYSTEM_CLOCK ( en, clock_rate, clock_max )

    write(*,*) ' '
    write(*,*) 'Run time...'
    write ( *, * ) 'Elapsed real time = ', real ( en - st ) / real ( clock_rate )

    IF(TERMINATION) THEN
        WRITE(OUTUnit,*) 'Program Exectuted in full...'
    ELSE
        WRITE(OUTUnit,*) 'Program terminated'
    END IF 
     WRITE(OUTUnit,*) 'Program start'
    write(OUTUnit,"(8i5)") values
    call date_and_time(date,time,zone,values)
    !call date_and_time(DATE=date,ZONE=zone)
    !call date_and_time(TIME=time)
    call date_and_time(VALUES=values)
    WRITE(OUTUnit,*) 'Program End'
    write(OUTUnit,"(8i5)") values
    WRITE(OUTUnit,*) 'RUN TIME'
    WRITE(OUTUnit,*) 'Time = ',finish-start,' seconds'
    WRITE(OUTUnit,*) 'Time = ',(finish-start)/60.d0,' minutes'
    WRITE(OUTUnit,*) 'Time = ',(finish-start)/3600.d0,' hours'
    write(OUTUnit, * ) 'Elapsed real time = ', real ( en - st ) / real ( clock_rate )
    close(2)
    close(5)

    end subroutine solve

    subroutine readdata()
    !**********************************************************************
    !    Function: Make a file
    !**********************************************************************

    implicit none
    double precision :: r1(2), r2(2), LPos(2), temp
    integer :: I, J, IEl, INod(4)

    open(DATUnit, file = 'input.dat')
    read(DATUnit, *) NDivX, NDivY
    read(DATUnit, *) Meshdx, Meshdy
        NNod = (NDivX + 1)*(NDivY + 1)
        NEl = NDivX * NDivY
    allocate (NodCo(2, NNod)); NodCo = 0.d0
    allocate (NodCoAvg(2, NNod)); NodCoAvg = 0.d0
    allocate (ICon(4, NEl)); icon = -1
    allocate (IsFixDof(2 * NNod))
    allocate(MasS(2 * NNod)); MasS = 0.d0
    allocate(GrvF(2 * NNod)); GrvF = 0.d0
    allocate(InrF(2 * NNod)); InrF = 0.d0
    allocate(ExtF(2 * NNod)); ExtF = 0.d0
    allocate(ExtP(2 * NNod)); ExtP = 0.d0
    allocate(InrFW(2 * NNod)); InrFW = 0.d0
    allocate(GrvFW(2 * NNod)); GrvFW = 0.d0
    allocate(MasW(2 * NNod)); MasW = 0.d0
    allocate(MasbW(2 * NNod)); MasbW = 0.d0
    allocate(DragF(2 * NNod)); DragF = 0.d0

    allocate(v(2 * NNod)); V = 0.d0
    allocate(v0(2 * NNod)); V0 = 0.d0
    allocate(W(2 * NNod)); W = 0.d0
    allocate(W0(2 * NNod)); W0 = 0.d0
    
    allocate(dis(2 * NNod)); dis = 0.d0
    allocate(DisW(2 * NNod)); DisW = 0.d0
    allocate(B(2, 4, NEl, 4)); B = 0.d0
    allocate(BNEW(3, 8, NEl, 4)); BNEW = 0.D0
    allocate(TEMPB(2, 4, NEl, 4)); TEMPB = 0.D0
    allocate(Area(NEl)); Area = 0.d0
    allocate(Areai(NEl)); Areai = 0.d0
    allocate(SigG(4, 4, NEl)); SigG = 0.d0
    allocate(Sig0(4, 4, NEl)); Sig0 = 0.d0
    allocate(SigE(4, 4, NEl)); SigE = 0.d0
    allocate(EpsG(3, 4, NEl)); EpsG = 0.d0
    allocate(PoreP(4, 4, NEl)); PoreP = 0.d0
    
    allocate(StateVar(32,4,NEl)); StateVar = 0.d0   ! State Variable
    allocate(StatevarPRINT(3,4,NEl)); StatevarPRINT = 0.d0   ! State Variable
    
    allocate(DEL_Mat(4, 4, NEl, 4)); DEL_Mat = 0.d0
    allocate(Rev_Vect(2, 4, NEl)); Rev_Vect = 0.d0
    allocate(FBARNEW(1, 4, NEL)); FBARNEW = 0.d0
    allocate(PHIMOB(1, 4, NEL)); PHIMOB= 0.d0

    allocate(F(4, 4, nel)); F(1,:,:) = 1.d0; F(2,:,:) = 1.d0; F(3,:,:) = 0.d0; F(4,:,:) = 0.d0
    allocate(HS(4, NEL, 4)); HS = 0.d0

    do I = 1, NNod
        read(DATUnit, *) temp, NodCo(1, I), NodCo(2, I), temp
    end do
    do I = 1, NEl
        read(DATUnit, *) temp, ICon(1, I), ICon(2, I), ICon(3, I), ICon(4, I)
    end do
    read(DATUnit, *) MatProp(1), MatProp(2), MatProp(3), MatProp(4), MatProp(5)  ! E,nu,rhoS,n,kappa(conductivity)
    read(DATUnit, *) MatProp(6), MatProp(7), MatProp(8), MatProp(9), MatProp(10) ! Dummy values added for debugging 
    read(DATUnit, *) MatProp(11), MatProp(12), MatProp(13), MatProp(14), MatProp(15), MatProp(16), MatProp(17) ! Phic, Pt, hs, n, e_d0, e_c0 e_i0
    read(DATUnit, *) MatProp(18), MatProp(19), MatProp(20), MatProp(21), MatProp(22), MatProp(23), MatProp(24), MatProp(25) ! alpha, beta, m_r, m_t, R_max, beta_r, chi, e_0
    !read(DATUnit, *) MatProp(26), MatProp(27), MatProp(28), MatProp(29), MatProp(30), MatProp(31), MatProp(32) !for clay only
    read(DATUnit, *) MatProp(26), MatProp(27), MatProp(28), MatProp(29), MatProp(30), MatProp(31)
    read(DATUnit, *) MatProp(32), MatProp(33), MatProp(34), MatProp(35), MatProp(36), MatProp(37), MatProp(38), MatProp(39), MatProp(40)
    read(DATUnit, *) gx, gy
    read(DATUnit, *) ndtn, dt
    close(1)
    
    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        Area(IEl) = abs(NodCo(1,INod(3))-NodCo(1,INod(1)))* abs(NodCo(2,INod(3))-NodCo(2,INod(1)))
    end do

    ! Boundary Conditions
    IsFixDof = .false.
    
    !IsFixDof(1) = .true.
    !IsFixDof(2) = .true.
    !IsFixDof(3) = .true.
    !IsFixDof(4) = .true.
    !IsFixDof(5) = .true.
    !IsFixDof(7) = .true.
          
    
    do I = 1, NNod ! vertical bar problem
        if ((NodCo(2, I) .eq. 0.d0)) then
            IsFixDof((I - 1) * 2 + 2) = .true. 
            IsFixDof((I - 1) * 2 + 1) = .true.
        end if
        if ((NodCo(1, I) .eq. 0.0d0)) then
            IsFixDof((I - 1) * 2 + 1) = .true.
        end if
        if ((NodCo(1, I) .eq. 0.1d0)) then
            IsFixDof((I - 1) * 2 + 1) = .true.
        end if ! Base Fixed
    !    
    !    !if ((NodCo(1, I) .eq. 1.d0).and.((NodCo(2, I) .eq. 0.d0))) then
    !    !    IsFixDof((I - 1) * 2 + 1) = .true.
    !    !end if
    !    !if ((NodCo(1, I) .eq. 1.d0).and.((NodCo(2, I) .eq. 5.d0))) then
    !    !    IsFixDof((I - 1) * 2 + 1) = .true.
    !    !end if
    !    
     end do
     
     WRITE(*,*) '***************************************'
     WRITE(*,*) '           TWO-PHASE FEM               '
     WRITE(*,*) '***************************************'
     write(*,*) ' '
     !write(*,*) 'Press Enter to continue...'
     !READ(*,*)

    end subroutine readdata

subroutine Initial()

    !**********************************************************************
    !    Function: Calculates B-Strain Displacement matrix (1 and 4 gauss points)
    !**********************************************************************

    implicit none

    integer :: IEl, I, J, K, INod(4), Id, ig
    double precision :: LPos(2, 1), dNxi(2, 4), Ja(2, 2), JaI(2, 2), A
    double precision :: xi, eta, rp, rm, sp, sm, temp
    INTEGER :: IGRAD = 2
    double precision :: AREAAVG = 0.d0

    !Porosity and Soild density
    rhoS = MatProp(3)
    Poro = MatProp(4)
    kappa = MatProp(5) ! conductivity or Darcy Permeability
    
    !dry and saturated density
    rhoD= (1.d0-poro)*rhoS
    rhoSat= rhoD + poro * rhoW
    
    !********************************************************
    temp = 1.D0/sqrt(3.d0)
    
    !temp=0.d0
    
    !temp = 0.5D0
    !********************************************************

    B = 0.0d0
    BNEW  = 0.0D0

    do IEl = 1, nel
        INod(:) = ICon(:, IEl)

        do ig = 1, 4

            select case (ig)
            case(1)
                xi = -temp
                eta = -temp

            case (2)
                xi = temp
                eta = -temp

            case(3)
                xi = temp
                eta = temp

            case(4)
                xi = -temp
                eta = temp
            end select

            rp = 1.0 + xi
            rm = 1.0 - xi;
            sp = 1.0 + eta;
            sm = 1.0 - eta;

            dNxi(1, 1) = -0.25D0 * sm; dNxi(1, 2) = +0.25D0 * sm; dNxi(1, 3) = +0.25D0 * sp
            dNxi(1, 4) = -0.25D0 * sp; dNxi(2, 1) = -0.25D0 * rm; dNxi(2, 2) = -0.25D0 * rp
            dNxi(2, 3) = +0.25D0 * rp; dNxi(2, 4) = +0.25D0 * rm

            HS(1, iel, ig) = (1.D0 - xi)*(1.D0 - eta)/4.D0
            HS(2, iel, ig) = (1.D0 + xi)*(1.D0 - eta)/4.D0
            HS(3, iel, ig) = (1.D0 + xi)*(1.D0 + eta)/4.D0
            HS(4, iel, ig) = (1.D0 - xi)*(1.D0 + eta)/4.D0

            Area(IEl) = abs(NodCo(1, INod(3)) - NodCo(1, INod(1))) * abs(NodCo(2, INod(3)) - NodCo(2, INod(1)))

            Ja = 0.0D0

            do I = 1, 2
                do J = 1, 2
                    do K = 1, 4
                        Ja(I, J) = Ja(I, J) + dNxi(I, K) * NodCo(J, INod(K))
                    end do
                end do
            end do
            A = Ja(1, 1) * Ja(2, 2) - Ja(1, 2) * Ja(2, 1)

            if (A .gt. 0.D0) then
                JaI(1, 1) = +Ja(2, 2)/A; JaI(1, 2) = -Ja(1, 2)/A
                JaI(2, 1) = -Ja(2, 1)/A; JaI(2, 2) = +Ja(1, 1)/A
            else
                write(LOGUnit, *) 'negative or zero Jacobian !!'; stop
            end if

            do J = 1, 4
                B(1, J, IEl, ig) = dNxi(1, J) * JaI(1, 1) + dNxi(2, J) * JaI(1, 2)
                B(2, J, IEl, ig) = dNxi(1, J) * JaI(2, 1) + dNxi(2, J) * JaI(2, 2)
            end do
            
            AreaAvg = abs(NodCo(1, INod(3)) - NodCo(1, INod(1))) * abs(NodCo(2, INod(4)) - NodCo(2, INod(2))) +&
                        abs(NodCo(1, INod(2)) - NodCo(1, INod(4))) * abs(NodCo(2, INod(3)) - NodCo(2, INod(1)))
                        
            NodCoAvg(1,INod(1)) = (NodCo(2,INod(2)) - NodCo(2,INod(4))) / AreaAvg
            NodCoAvg(1,INod(2)) = (NodCo(2,INod(3)) - NodCo(2,INod(1))) / AreaAvg
            NodCoAvg(1,INod(3)) = (NodCo(2,INod(4)) - NodCo(2,INod(2))) / AreaAvg
            NodCoAvg(1,INod(4)) = (NodCo(2,INod(1)) - NodCo(2,INod(3))) / AreaAvg
            
            NodCoAvg(2,INod(1)) = (NodCo(1,INod(4)) - NodCo(1,INod(2))) / AreaAvg
            NodCoAvg(2,INod(2)) = (NodCo(1,INod(1)) - NodCo(1,INod(3))) / AreaAvg
            NodCoAvg(2,INod(3)) = (NodCo(1,INod(2)) - NodCo(1,INod(4))) / AreaAvg
            NodCoAvg(2,INod(4)) = (NodCo(1,INod(3)) - NodCo(1,INod(1))) / AreaAvg
            
            
            !reassigning B matrix to a new structure to carry out selective reduced integration
            do J = 1 ,4
                BNEW(1, J*2-1, IEL, ig) = B(1, J, IEl, ig)
                BNEW(2, J*2, IEL, ig)   = B(2, J, IEl, ig)
                BNEW(3, J*2-1, IEL, ig) = B(2, J, IEl, ig)
                BNEW(3, J*2, IEL, ig)   = B(1, J, IEl, ig)
            end do
            
            do J = 1 ,4
                TEMPB(1,J,Iel,ig) = B(1, J, IEl, ig) - BNEW(1,J*2-1,IEL, ig)  
                TEMPB(2,J,Iel,ig) = B(2, J, IEl, ig) - BNEW(1,J*2  ,IEL, ig)
                RES = RES + TEMPB(1,J,Iel,ig) + TEMPB(2,J,Iel,ig)
            end do
            
            
            IF(IGRAD.EQ.1) THEN !selective integration for shear
            do J = 1, 4
                BNEW(3, J*2-1, IEL, ig) = NodCoAvg(2,INod(J))
                BNEW(3, J*2, IEL, ig)   = NodCoAvg(1,INod(J))
            end do
            END IF
            
            IF(IGRAD.EQ.2) THEN !selective integration for volume change
            do J = 1, 4
                BNEW(1,J*2-1,IEl,ig) = BNEW(1,J*2-1,IEl,ig) - (BNEW(1,J*2-1,IEl,ig) - NodCoAvg(1,INod(J)))/3.D0
                BNEW(1,J*2,IEl,ig)   =                      - (BNEW(2,J*2  ,IEl,ig) - NodCoAvg(2,INod(J)))/3.D0 !BNEW(1,J*2  ,IEl,ig)
                BNEW(2,J*2-1,IEl,ig) =                      - (BNEW(1,J*2-1,IEl,ig) - NodCoAvg(1,INod(J)))/3.D0
                BNEW(2,J*2,IEl,ig)   = BNEW(2,J*2  ,IEl,ig) - (BNEW(2,J*2  ,IEl,ig) - NodCoAvg(2,INod(J)))/3.D0
            end do
            END IF
                
        enddo
    enddo

end subroutine Initial
    
    !**********************************************************************
    !    Function: Average gradient calculation function
    !**********************************************************************
  
    !
    !subroutine avgrad()
    !
    !implicit none
    !double precision :: area = 0.D0
    !
    !!A2 = (X(3)-X(1))*(Y(4)-Y(2))+(X(2)-X(4))*(Y(3)-Y(1))
    !
    !print(*,*) NodCo, NodCoAvg
    !
    !end subroutine avgrad
subroutine Map2Nod()

    !**********************************************************************
    !    Function: Internal and External Forces, Mass at each nodes
    !**********************************************************************

    implicit none
    integer :: I, IEl, Id, INod(4), ig, J
   double precision :: factor, temp, timer, switch = 0.d0, facmul
   
  

    GrvF = 0.D0; InrF = 0.D0; ExtF = 0.D0; MasS = 0.D0
    GrvFW = 0.D0; InrFW = 0.D0; ExtP = 0.D0; MasW = 0.D0; MasbW = 0.d0; DragF = 0.d0
    
    if(itime*dt.gt.(1.0d0)) then 
       factor = 1.d0
    else
       factor = itime * dt / 1.0d0    
    endif 
    
    !factor = itime * dt
    !if(itime*dt.le.1.0D0) factor = itime * dt / 1.D0
    !if(itime*dt.gt.1.0D0) factor = 1.d0
    
   !if(itime*dt.lt.(0.5d0)) then 
   !     ExtF(2) = -0.3
   !     ExtF(4) = -0.3
   !     ExtP(2) = -0.3
   !     ExtP(4) = -0.3
   !else
   !     factor =  1.d0
   !     ExtF(2) = -0.3
   !     ExtF(4) = -0.3
   !     ExtP(2) = 0
   !     ExtP(4) = 0
   ! endif
        !if(itime*dt.le.10.D0) then
        !ExtF(3) = -50.D0 ; !ExtP(3) = 0.D0
        !ExtF(5) = -50.D0 ; !ExtP(5) = 0.D0
        !ExtF(6) = -50.D0 ; !ExtP(6) = 0.D0
        !ExtF(8) = -50.D0 ; !ExtP(8) = 0.D0
        !ExtF(1) =  50.D0 ;
        !ExtF(7) =  50.D0 ;
        !endif
        !
        !if(itime*dt.gt.10.d0) then
        !  ExtF(3) = -50.D0 ;  !ExtP(3) = PoreP(2,1,1)
        !  ExtF(1) =  50.D0 ;  !ExtP(1) =  PoreP(2,1,1)
        !  ExtF(5) = -50.D0 ;  !ExtP(5) = PoreP(2,1,1)
        !  ExtF(7) =  50.D0 ;  !ExtP(7) = PoreP(2,1,1)
        !  ExtF(8) = -50.D0
        !  ExtF(6) = -50.D0
        !  ExtF(7) =  50.D0
        !  ExtF(5) = -50.D0
        !endif  
                
    
        timer = itime*dt
        
        !if(itime*dt.gt.10.D0) then
        !    ExtF(6) = -50.D0 
        !    ExtF(8) = -50.D0 
            !ExtF(3) = -50.D0
            !ExtF(5) = -50.D0 
            !ExtF(1) =  50.D0
            !ExtF(7) =  50.D0
        !end if 
        
        !facmul = 0.d0
        !
        !if(timer.gt.10.d0) then
        !    
        !    if(timer.ge.10.d0.and.timer.lt.12.5d0) facmul = (timer-10.d0)/2.5d0
        !    if(timer.ge.12.5d0.and.timer.lt.15.d0) facmul = ((timer-12.5d0)/2.5d0 - 1.d0)
        !    if(timer.ge.15.d0.and.timer.lt.17.5d0) facmul = (timer-15.d0)/2.5d0
        !    if(timer.ge.17.5d0.and.timer.lt.20.d0) facmul = ((timer-17.5d0)/2.5d0 - 1.d0)
        !    
        !    if(timer.gt.10.d0.and.timer.le.12.5d0.or.timer.gt.17.5d0) then
        !        switch = +1.D0
        !    else
        !        switch = -1.D0
        !    endif
        !endif 
        !
        !if(timer.ge.20.d0) then
        !    
        !    if(timer.ge.20.d0.and.timer.lt.22.5d0) facmul = (timer-20.d0)/2.5d0
        !    if(timer.ge.22.5d0.and.timer.lt.25.d0) facmul = ((timer-22.5d0)/2.5d0 - 1.d0)
        !    if(timer.ge.25.d0.and.timer.lt.27.5d0) facmul = (timer-25.d0)/2.5d0
        !    if(timer.ge.27.5d0.and.timer.lt.30.d0) facmul = ((timer-27.5d0)/2.5d0 - 1.d0)
        !    
        !    if(timer.ge.20.d0.and.timer.le.22.5d0.or.timer.gt.27.5d0) then
        !        switch = +1.D0
        !    else
        !        switch = -1.D0
        !    endif        
        !end if 
        !
        !
        !if(timer.ge.30.d0) then
        !    
        !    if(timer.ge.30.d0.and.timer.lt.32.5d0) facmul = (timer-30.d0)/2.5d0
        !    if(timer.ge.32.5d0.and.timer.lt.35.d0) facmul = ((timer-32.5d0)/2.5d0 - 1.d0)
        !    if(timer.ge.35.d0.and.timer.lt.37.5d0) facmul = (timer-35.d0)/2.5d0
        !    if(timer.ge.37.5d0.and.timer.lt.40.d0) facmul = ((timer-37.5d0)/2.5d0 - 1.d0)
        !    
        !    if(timer.ge.30.d0.and.timer.le.32.5d0.or.timer.gt.37.5d0) then
        !        switch = +1.D0
        !    else
        !        switch = -1.D0
        !    endif
        !endif 
        !
        !if(timer.ge.40.d0) then
        !    
        !    if(timer.ge.40.d0.and.timer.lt.42.5d0) facmul = (timer-40.d0)/2.5d0
        !    if(timer.ge.42.5d0.and.timer.lt.45.d0) facmul = ((timer-42.5d0)/2.5d0 - 1.d0)
        !    if(timer.ge.45.d0.and.timer.lt.47.5d0) facmul = (timer-45.d0)/2.5d0
        !    if(timer.ge.47.5d0.and.timer.lt.50.d0) facmul = ((timer-47.5d0)/2.5d0 - 1.d0)
        !    
        !    if(timer.ge.40.d0.and.timer.le.42.5d0.or.timer.gt.47.5d0) then
        !        switch = +1.D0
        !    else
        !        switch = -1.D0
        !    endif
        !endif 
        
        
        
    !if(timer.gt.10.d0) then 
    !    if(switch.gt.0.d0) then
    !        ExtF(7) = +5.D0 * facmul; !ExtP(7) = PoreP(2,1,1)
    !        ExtF(5) = +5.D0 * facmul; !ExtP(5) = PoreP(2,1,1)
    !    end if
    !    if(switch.lt.0.d0) then
    !        ExtF(7) = -5.D0 * facmul; !ExtP(7) = -PoreP(2,1,1)
    !        ExtF(5) = -5.D0 * facmul; !ExtP(5) = -PoreP(2,1,1)
    !    end if
    !end if 
        
            
        do IEl = 1, nel
        INod(:) = ICon(:, IEl)
        do ig = 1, 4
            do I = 1, 4
                Id = (INod(I) - 1) * 2
                !! Soil part
                InrF(Id + 1) = InrF(Id + 1)+(Sigg(1,ig,IEl) * B(1,I,iel,ig) + Sigg(3,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                InrF(Id + 2) = InrF(Id + 2)+(Sigg(3,ig,IEl) * B(1,I,iel,ig) + Sigg(2,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                GrvF(Id + 1) = GrvF(Id + 1) + Area(iel) * rhoSat * hs(i, iel, ig) * gx /4.d0 * factor
                GrvF(Id + 2) = GrvF(Id + 2) + Area(iel) * rhoSat * hs(i, iel, ig) * gy /4.d0 * factor
                MasS(Id + 1) = MasS(Id + 1) + Area(iel) * (1 - poro) * rhoS * hs(i, iel, ig) /4.d0
                MasS(Id + 2) = MasS(Id + 2) + Area(iel) * (1 - poro) * rhoS * hs(i, iel, ig) /4.d0
                
                !!! Soil part with new B Matrix
                !InrF(Id + 1) = InrF(Id + 1)+(Sigg(1,ig,IEl) * BNEW(1,I*2-1,iel,ig) + Sigg(3,ig,IEl) * BNEW(3, I*2, iel, ig)) * Area(iel)/4.d0
                !InrF(Id + 2) = InrF(Id + 2)+(Sigg(3,ig,IEl) * BNEW(3,I*2-1,iel,ig) + Sigg(2,ig,IEl) * BNEW(2, I*2, iel, ig)) * Area(iel)/4.d0
                !GrvF(Id + 1) = GrvF(Id + 1) + Area(iel) * rhoSat * hs(i, iel, ig) * gx /4.d0 * factor
                !GrvF(Id + 2) = GrvF(Id + 2) + Area(iel) * rhoSat * hs(i, iel, ig) * gy /4.d0 * factor
                !MasS(Id + 1) = MasS(Id + 1) + Area(iel) * (1 - poro) * rhoS * hs(i, iel, ig) /4.d0
                !MasS(Id + 2) = MasS(Id + 2) + Area(iel) * (1 - poro) * rhoS * hs(i, iel, ig) /4.d0
                
                ! Water part
                InrFW(Id + 1) = InrFW(Id + 1) + (PoreP(1,ig,IEl) * B(1,I,iel,ig) + PoreP(3,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                InrFW(Id + 2) = InrFW(Id + 2) + (PoreP(3,ig,IEl) * B(1,I,iel,ig) + PoreP(2,ig,IEl) * B(2, I, iel, ig)) * Area(iel)/4.d0
                GrvFW(Id + 1) = GrvFW(Id + 1) + Area(iel) * rhoW * hs(i, iel, ig) * gx /4.d0 * factor
                GrvFW(Id + 2) = GrvFW(Id + 2) + Area(iel) * rhoW * hs(i, iel, ig) * gy /4.d0 * factor
                MasW(Id + 1)  = MasW(Id + 1) + Area(iel) *  rhoW * hs(i, iel, ig) /4.d0
                MasW(Id + 2)  = MasW(Id + 2) + Area(iel) *  rhoW * hs(i, iel, ig) /4.d0
                MasbW(Id + 1) = MasbW(Id + 1) + Area(iel) * (poro) * rhoW * hs(i, iel, ig) /4.d0
                MasbW(Id + 2) = MasbW(Id + 2) + Area(iel) * (poro) * rhoW * hs(i, iel, ig) /4.d0
                DragF(Id + 1) = DragF(Id + 1) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) !* (W0(Id+1) - V0(Id+1))
                DragF(Id + 2) = DragF(Id + 2) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) !* (W0(Id+2) - V0(Id+2))
                
                !! Water part
                !InrFW(Id + 1) = InrFW(Id + 1) + (PoreP(1,ig,IEl) * BNEW(1,I*2-1,iel,ig) + PoreP(3,ig,IEl) * BNEW(3, I*2, iel, ig)) * Area(iel)/4.d0
                !InrFW(Id + 2) = InrFW(Id + 2) + (PoreP(3,ig,IEl) * BNEW(3,I*2-1,iel,ig) + PoreP(2,ig,IEl) * BNEW(2, I*2, iel, ig)) * Area(iel)/4.d0
                !GrvFW(Id + 1) = GrvFW(Id + 1) + Area(iel) * rhoW * hs(i, iel, ig) * gx /4.d0 * factor
                !GrvFW(Id + 2) = GrvFW(Id + 2) + Area(iel) * rhoW * hs(i, iel, ig) * gy /4.d0 * factor
                !MasW(Id + 1)  = MasW(Id + 1) + Area(iel) *  rhoW * hs(i, iel, ig) /4.d0
                !MasW(Id + 2)  = MasW(Id + 2) + Area(iel) *  rhoW * hs(i, iel, ig) /4.d0
                !MasbW(Id + 1) = MasbW(Id + 1) + Area(iel) * (poro) * rhoW * hs(i, iel, ig) /4.d0
                !MasbW(Id + 2) = MasbW(Id + 2) + Area(iel) * (poro) * rhoW * hs(i, iel, ig) /4.d0
                !DragF(Id + 1) = DragF(Id + 1) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) * (W0(Id+1) - V0(Id+1))
                !DragF(Id + 2) = DragF(Id + 2) + ( Area(iel) *  rhoW * poro * 9.81d0 * hs(i, iel, ig) / (kappa * 4.d0)) * (W0(Id+2) - V0(Id+2))
            end do
        enddo !gauss / particles
    enddo !nel  
	
	do IEl = 1, nel
	INod(:) = ICon(:, IEl)
	do ig = 1, 4
	do I = 1, 4
                Id = (INod(I) - 1) * 2
				DragF(Id + 1) = DragF(Id + 1)*(W0(Id+1) - V0(Id+1))
				DragF(Id + 2) = DragF(Id + 2)*(W0(Id+2) - V0(Id+2))
	end do
	end do
	end do

end subroutine Map2Nod

subroutine Update()
    !**********************************************************************
    !    Function: Small deformation formulation
    !*************************** *******************************************

    implicit none
    integer :: IEl, INod(4), Id, I, ig, J
    double precision :: delV(4), deps(4), Bulk = 2e6, delW(4), dpore, aW, tem, dampf = 0.D0, temp, timer, switch = 0.d0
    
    
    integer :: dummy1
    double precision :: dummy2, dummy3
    !Initial Velocity
    V = 0.d0
    W = 0.d0
   
    !if(itime*dt.eq.10.d0) EpsG = 0.d0 !Resetting strain to 0 after initial loading
    
    !if(itime*dt.gt.10.d0) dampf = 0.0001D0
    
    
    !if(itime*dt.ge.10.D0) then
    !    IsFixDOF    = .false. 
    !    IsFixDOF(1) = .true.
    !    IsFixDOF(2) = .true.
    !    IsFixDOF(3) = .true.
    !    IsFixDOF(4) = .true.
    !    !IsFixDOF(6) = .true.
    !    !IsFixDOF(8) = .true.
    !endif 
    
    
    do I = 1, NNod
        Id = (I - 1) * 2
        if (.not.IsFixDof(Id + 1).and.(MasS(Id + 1) .gt. 0.D0)) then
            W(Id + 1) = W0(Id + 1) + (GrvFW(Id + 1) + ExtP(Id + 1) - InrFW(Id + 1) - DragF(Id + 1)) / MasW (Id + 1) * dt
            aW = (W(Id + 1) - W0(Id + 1)) / dt
            tem = V0(Id + 1) + (GrvF(Id + 1) + ExtF(Id + 1) - InrF(Id + 1) - (aW * MasbW (Id + 1)) ) / MasS(Id + 1) * dt
            !V(Id + 1) = V0(Id + 1) + (GrvF(Id + 1) + ExtF(Id + 1) - InrF(Id + 1) - (aW * MasbW (Id + 1)) ) / MasS(Id + 1) * dt
            V(Id + 1) = tem  - sign(1.D0,tem) * dampf * abs(tem)
            end if

        if (.not.IsFixDof(Id + 2).and.(MasS(Id + 2) .gt. 0.D0)) then
            W(Id + 2) = W0(Id + 2) + (GrvFW(Id + 2) + ExtP(Id + 2) - InrFW(Id + 2) - DragF(Id + 2)) / MasW (Id + 2) * dt
            aW = (W(Id + 2) - W0(Id + 2)) / dt
            !V(Id + 2) = V0(Id + 2) + (GrvF(Id + 2) + ExtF(Id + 2) - InrF(Id + 2) - (aW * MasbW (Id + 2)) ) / MasS(Id + 2) * dt
            tem = V0(Id + 2) + (GrvF(Id + 2) + ExtF(Id + 2) - InrF(Id + 2) - (aW * MasbW (Id + 2)) ) / MasS(Id + 2) * dt
            V(Id + 2) = tem - sign(1.D0,tem) * dampf * abs(tem)
            end if
    end do
    
    
    !if(itime*dt.GT.10.D0) then
    !
    !    V(6) = -0.008D0; W(6) = -0.008D0 
    !    V(8) = -0.008D0; W(8) = -0.008D0
    !
    !    W(1) = V(1) 
    !    W(3) = V(3) 
    !    W(5) = V(5) 
    !    W(7) = V(7) 
    !
    !endif 
    
    timer = itime*dt
    
    !if(timer.gt.10.d0) then
    !    if(timer.gt.10.d0.and.timer.le.12.5d0.or.timer.gt.17.5d0.and.timer.le.20.d0.or.timer.gt.20.d0.and.timer.le.22.5d0.or.timer.gt.27.5d0) then
    !        switch = +1.D0
    !    else
    !        switch = -1.D0
    !    endif
    !endif 
    !
    !if(timer.ge.30.d0) then
    !    if(timer.ge.30.d0.and.timer.le.32.5d0.or.timer.gt.37.5d0) then
    !        switch = +1.D0
    !    else
    !        switch = -1.D0
    !    endif
    !endif 
    !
    !if(timer.ge.40.d0) then
    !    if(timer.ge.40.d0.and.timer.le.42.5d0.or.timer.gt.47.5d0) then
    !        switch = +1.D0
    !    else
    !        switch = -1.D0
    !    endif
    !endif 
    
    !W(6) = V(6)
    !W(8) = V(8)
    
    !W(5) = V(6)
    !W(7) = V(7)
    
    !if(timer.gt.10.d0) then 
    !    if(switch.gt.0.d0) then
    !        !V(6) = 0.d0   ; W(6) = 0.d0
    !        !V(8) = 0.d0   ; W(8) = 0.d0
    !        !W(6) = V(6)
    !        !W(8) = V(8)
    !        
    !        V(7) = 0.002d0 ; W(7) = 0.002d0
    !        V(5) = 0.002d0 ; W(5) = 0.002d0
    !        !W(5) = V(5)
    !        !W(7) = V(7)
    !        !W(3) = V(3)
    !        !W(1) = V(1)
    !    end if
    !    if(switch.lt.0.d0) then
    !        !V(6) = 0.d0   ; W(6) = 0.d0
    !        !V(8) = 0.d0   ; W(8) = 0.d0
    !        !W(6) = V(6)
    !        !W(8) = V(8)
    !        
    !        V(7) = -0.002d0 ; W(7) = -0.002d0
    !        V(5) = -0.002d0 ; W(5) = -0.002d0
    !        !W(5) = V(5)
    !        !W(7) = V(7)
    !        !W(3) = V(3)
    !        !W(1) = V(1)
    !    end if
    !end if 
    


    do IEl = 1, NEl
        INod(:) = ICon(:, IEl)
        do ig = 1, 4
            delV = 0.0
            delW = 0.0
            do I = 1, 4
                Id = (INod(I) - 1) * 2
                !!! Old implementation
                delV(1) = delV(1) + B(1, i, iel, ig) * V(Id + 1)
                delV(2) = delV(2) + B(2, i, iel, ig) * V(Id + 2)
                delV(3) = delV(3) + B(2, i, iel, ig) * V(Id + 1)
                delV(4) = delV(4) + B(1, i, iel, ig) * V(Id + 2)
                delW(1) = delW(1) + B(1, i, iel, ig) * W(Id + 1)
                delW(2) = delW(2) + B(2, i, iel, ig) * W(Id + 2)
                delW(3) = delW(3) + B(2, i, iel, ig) * W(Id + 1)
                delW(4) = delW(4) + B(1, i, iel, ig) * W(Id + 2)
                !
                !!! New implementation
                !delV(1) = delV(1) + BNEW(1, i*2-1, iel, ig) * V(Id + 1)
                !delV(2) = delV(2) + BNEW(2, i*2, iel, ig) * V(Id + 2)
                !delV(3) = delV(3) + BNEW(3, i*2-1, iel, ig) * V(Id + 1)
                !delV(4) = delV(4) + BNEW(3, i*2, iel, ig) * V(Id + 2)
                !delW(1) = delW(1) + BNEW(1, i*2-1, iel, ig) * W(Id + 1)
                !delW(2) = delW(2) + BNEW(2, i*2, iel, ig) * W(Id + 2)
                !delW(3) = delW(3) + BNEW(3, i*2-1, iel, ig) * W(Id + 1)
                !delW(4) = delW(4) + BNEW(3, i*2, iel, ig) * W(Id + 2)
            end do
            ! Strain increment
            dEps(1:2) = delV(1:2) * dt
            dEps(3) = (delV(3) + delV(4)) * dt        
            
            ! Total strain
            EpsG(1:2, ig, IEl) = EpsG(1:2, ig, IEl) +  dEps(1:2)
            EpsG(3, ig, IEl) = EpsG(3, ig, IEl) +  dEps(3)
            
            ! Pressure increment
            dPore = (((1-poro) * delV(1) + poro * delW(1)) + ((1-poro) * delV(2) + poro * delW(2))) * ((dt * Bulk) / poro )
           
            ! Total pressure
            PoreP(1:2, ig, IEl) = PoreP(1:2, ig, IEl) +  dPore 
            PoreP(3, ig, IEl) = 0.d0
        
        enddo !gauss
 
        do ig = 1, 4
            call Elastic(MatProp(1), MatProp(2), Epsg(:,ig,IEL), SigE(:,ig,IEl))
            
            !call Neutral_Liquefaction(SigE(:,ig,iel), dEps, DEL_Mat(:,:,ig,iel), Rev_Vect(:,ig,iel), MatProp(1), MatProp(2))
            
            !call Elastic(MatProp(6), MatProp(7), Epsg(:,ig,IEL), SigE(:,ig,IEl))
            !call MOHRC(SIN(MatProp(9)*0.0174533D0), SIN(MatProp(10)*0.0174533D0), MatProp(8), MatProp(6), MatProp(7), SigE(:,ig,iel), NPLAST, FBAR, PLASTIND)
            
            !call IntfUMAT_Sand(dt, matprop, dEps, SigE(:,ig,iel), Statevar(:,ig,iel)) 
            
            !call Intf_DLL_Plaxis(dt, matprop, dEps, SigE(:,ig,iel), Statevar(:,ig,iel),PoreP(:,ig,iel), itime*dt)
            Sigg(:,ig,IEL) = SigE(:,ig,IEL) + PoreP(:,ig,IEL)
            
            !StatevarPRINT(1,1:4,1:NEl) = ASIND(Statevar(1,1:4,1:NEl)/0.0174533D0)
            !StatevarPRINT(2,1:4,1:NEl) = ASIND(Statevar(2,1:4,1:NEl)/0.0174533D0)
            !StatevarPRINT(3,1:4,1:NEl) = ASIND(Statevar(5,1:4,1:NEl)/0.0174533D0)
            
            !PHIMOB(1,ig,iel) = ASIND(Statevar(1,ig,iel))
            
            !call smooth(Sigg(:,ig,iel),Area(IEl))
        enddo !gauss
     
    enddo ! elements

    dis = dis + v * dt
    !DisW = DisW + w * dt
    v0 = v
    w0 = w

end subroutine Update


subroutine smooth(SIG,AREA)
    
    implicit none
        
    double precision, intent(in) :: AREA
    double precision, intent(inout) :: SIG(4)
    !Local variables
    double precision :: AR,SIGSMOOTH(4)
    !
    AR = 0.D0
    SIGSMOOTH = 0.D0
    !
    SIGSMOOTH(1) = SIGSMOOTH(1) + SIG(1) * AREA!xx 
    SIGSMOOTH(2) = SIGSMOOTH(2) + SIG(2) * AREA!yy
    SIGSMOOTH(3) = SIGSMOOTH(3) + SIG(3) * AREA!xy
    SIGSMOOTH(4) = SIGSMOOTH(4) + SIG(4) * AREA!zz
    AR = AR + AREA
    
    IF(AR.GT.1.D-12) THEN
        SIGSMOOTH(1) = SIGSMOOTH(1)/AR
        SIGSMOOTH(2) = SIGSMOOTH(2)/AR
        SIGSMOOTH(3) = SIGSMOOTH(3)/AR
        SIGSMOOTH(4) = SIGSMOOTH(4)/AR
    END IF
    
    SIG(1) = SIGSMOOTH(1)/AR
    SIG(2) = SIGSMOOTH(2)/AR
    SIG(3) = SIGSMOOTH(3)/AR
    SIG(4) = SIGSMOOTH(4)/AR
    
        
end subroutine 

!*********************************************************
!       Elastic Hookes model
!**********************************************************
subroutine Elastic(E, nu, eps, Sig)

    implicit none
    
    double precision, intent(in) :: E, nu
    double precision, intent(out) :: Sig(4)
    double precision, intent(in) :: eps(3)
 
    ! local variables
    double precision :: G_mod, K_mod, Eps_tr

    G_mod = E / (2.d0 * (1 + nu))
    K_mod = E / (3.d0 * (1 - 2 * nu))
    
    Eps_tr = eps(1) + eps(2)
    
    Sig(1) = ((K_mod * Eps_tr) + 2 * G_mod * (eps(1) - (Eps_tr/3.d0)))
    Sig(2) = ((K_mod * Eps_tr) + 2 * G_mod * (eps(2) - (Eps_tr/3.d0)))
    Sig(4) = ((K_mod * Eps_tr) + 2 * G_mod * (0.d0 - (Eps_tr/3.d0)))
    Sig(3) = ( G_mod * eps(3))
                
endsubroutine elastic


subroutine MkOpFiles()
    !**********************************************************************
    !    Function: Output files
    !**********************************************************************
    implicit none

    call MkFile(MSHUnit, 'Model.post.msh')
    call MkFile(RESUnit, 'Model.post.res')

end subroutine MkOpFiles


subroutine MkFile(Unit, flNam)
    !**********************************************************************
    !
    !    Function: Make a file
    !
    !**********************************************************************

    implicit none

    integer Unit
    character flNam * (*)

    if (FlExist(flNam)) then
        open(Unit, file = flNam)
        close(Unit, Status = 'Delete')
    endif
    open(Unit, file = flNam)

end subroutine MkFile

logical function FlExist(flNam)
    !**********************************************************************
    !
    !    Function: To check the existence of a file
    !
    !**********************************************************************

    implicit none

    logical lfil
    character flNam * (*)

    lfil = .false.
    inquire(file = flNam, exist = lfil)
    if (lfil) then
        FlExist = .true.
    else
        FlExist = .false.
    endif

end function FlExist


subroutine WtMSH(Unit)
    !**********************************************************************
    !    Function: Writing GiD *.msh file
    !!**********************************************************************

    implicit none
    integer Unit
    ! local variables
    integer :: IPart, INod, IEl, J, K
    double precision, dimension(2) :: Pos, r1, r2, VPos

    write(Unit, *) 'MESH dimension 2 ElemType Quadrilateral Nnode 4'
    write(Unit, *) 'Coordinates'
    do INod = 1, NNod
        write(Unit, 111) INod, NodCo(1, INod), NodCo(2, INod)
    end do
    write(Unit, *) 'End Coordinates'
    write(Unit, *) 'Elements'
    do IEl = 1, NEl
        write(Unit, "(5I7)") IEl, ICon(1, IEl), ICon(2, IEl), ICon(3, IEl), ICon(4, IEl)
    end do
    write(Unit, *) 'End Elements'
    close(Unit)
    111 format(I7, 2E16.6E3)

end subroutine WtMSH

subroutine WtRES(IStep)
    !**********************************************************************
    !
    !    Function:
    !
    !**********************************************************************

    implicit none

    integer IStep
    ! local variables
    integer :: IEl, J, K, INod, Id, Unit, ig

    Unit = ResUnit
    if (IStep .eq. 1) then
        write(Unit, *) 'GiD Post Results File 1.0'
        write(Unit, *) 'GaussPoints "Material_Point" Elemtype Quadrilateral'
        write(Unit, *) 'Number of Gauss Points: 4'
        write(Unit, *) 'Natural Coordinates: Internal'
        write(Unit, *) 'end gausspoints'
        write(Unit, *) 'Result "Boundary" "MPM"', IStep, 'Vector OnNodes'
        write(Unit, *) 'ComponentNames "X-fix", "Y-fix"'
        write(Unit, *) 'values'
        do INod = 1, NNod
            Id = (INod - 1) * 2
            J = 0; K = 0
            if (IsFixDof(Id + 1)) J = 1
            if (IsFixDof(Id + 2)) K = 1
            write(Unit, "(3I7)") INod, J, K
        end do
        write(Unit, *) 'end values'
    end if
    write(Unit, *) 'Result "displacement" "MPM"', IStep, 'Vector OnNodes'
    write(Unit, *) 'ComponentNames "comp. x", "comp. y"'
    write(Unit, *) 'values'
    do INod = 1, NNod
        Id = (INod - 1) * 2
        write(Unit, "(I7, 2E16.6E3)") INod, Dis(Id + 1), Dis(Id + 2)
    end do
    write(Unit, *) 'end values'
    
    !write(Unit, *) 'Result "displacementWater" "MPM"', IStep, 'Vector OnNodes'
    !write(Unit, *) 'ComponentNames "compW. x", "compW. y"'
    !write(Unit, *) 'values'
    !do INod = 1, NNod
    !    Id = (INod - 1) * 2
    !    write(Unit, "(I7, 2E16.6E3)") INod, DisW(Id + 1), DisW(Id + 2)
    !end do
    !write(Unit, *) 'end values'
    
    write(Unit, *) 'Result "Total Stress" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "sigma xx", "sigma yy", "sigma xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 3E16.6E3)") IEl, Sigg(1, 1, IEl), Sigg(2, 1, IEl), Sigg(3, 1, IEl)
        write(Unit, "(3E16.6E3)") Sigg(1, 2, IEl), Sigg(2, 2, IEl), Sigg(3, 2, IEl)
        write(Unit, "(3E16.6E3)") Sigg(1, 3, IEl), Sigg(2, 3, IEl), Sigg(3, 3, IEl)
        write(Unit, "(3E16.6E3)") Sigg(1, 4, IEl), Sigg(2, 4, IEl), Sigg(3, 4, IEl)
    end do
    write(Unit, *) 'end values'
    
    write(Unit, *) 'Result "Eff Stress" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "sigmaE xx", "sigmaE yy", "sigmaE xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 3E16.6E3)") IEl, SigE(1, 1, IEl), SigE(2, 1, IEl), SigE(3, 1, IEl)
        write(Unit, "(3E16.6E3)") SigE(1, 2, IEl), SigE(2, 2, IEl), SigE(3, 2, IEl)
        write(Unit, "(3E16.6E3)") SigE(1, 3, IEl), SigE(2, 3, IEl), SigE(3, 3, IEl)
        write(Unit, "(3E16.6E3)") SigE(1, 4, IEl), SigE(2, 4, IEl), SigE(3, 4, IEl)
    end do
    write(Unit, *) 'end values'
    
    write(Unit, *) 'Result "Strain" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "eps xx", "eps yy", "eps xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 3E16.6E3)") IEl, EpsG(1, 1, IEl), EpsG(2, 1, IEl), EpsG(3, 1, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 2, IEl), EpsG(2, 2, IEl), EpsG(3, 2, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 3, IEl), Epsg(2, 3, IEl), Epsg(3, 3, IEl)
        write(Unit, "( 3E16.6E3)") Epsg(1, 4, IEl), Epsg(2, 4, IEl), Epsg(3, 4, IEl)
    end do
    write(Unit, *) 'end values'
    
    write(Unit, *) 'Result "Pore Pressure" "MPM"', IStep, 'Vector OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "Pore xx", "Pore yy", "Pore xy"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 3E16.6E3)") IEl, PoreP(1, 1, IEl), PoreP(2, 1, IEl), PoreP(3, 1, IEl)
        write(Unit, "(3E16.6E3)") PoreP(1, 2, IEl), PoreP(2, 2, IEl), PoreP(3, 2, IEl)
        write(Unit, "(3E16.6E3)") PoreP(1, 3, IEl), PoreP(2, 3, IEl), PoreP(3, 3, IEl)
        write(Unit, "(3E16.6E3)") PoreP(1, 4, IEl), PoreP(2, 4, IEl), PoreP(3, 4, IEl)
    end do
    write(Unit, *) 'end values'

    write(Unit, *) 'Result "Phi_Mob" "MPM"', IStep, 'Scalar OnGaussPoints "Material_Point"'
    write(Unit, *) 'ComponentNames "Failure Indicator"'
    write(Unit, *) 'values'
    do IEl = 1, NEl
        write(Unit, "(I7, 4E16.6E4)") IEl, PHIMOB(1,1,IEl)
        write(Unit, "(4E16.6E4)") PHIMOB(1,2,IEl)
        write(Unit, "(4E16.6E4)") PHIMOB(1,3,IEl)
        write(Unit, "(4E16.6E4)") PHIMOB(1,4,IEl)
    end do
    write(Unit, *) 'end values'    
end subroutine WtRES

end module modulefem

!************program *****************
program MPM2D

    use modulefem
    use omp_lib
    implicit none

    call MkOpFiles()
    call readdata()

    call solve()
end program
!*************program***************** 