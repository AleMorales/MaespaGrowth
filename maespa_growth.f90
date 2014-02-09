Program maespa_growth
Use Initialize, only: maespa_initialize, maespa_finalize
USE maindeclarations, only: istart, iday, iend, mflag, nstep, RYTABLE1, RXTABLE1, RZTABLE1, FOLTABLE1, ZBCTABLE1, totRespf, totCO2, TAIR, KHRS ! This contains all variables that were defined in the program unit of maespa.
Use growth_module ! This containts the parameter and state variables and the growth module
Implicit None
Double precision :: Assimilation, Pool, Allocation_leaf, Allocation_shoots, Allocation_stem, Allocation_froots, Allocation_croots, Allocation_fruits, Allocation_reserves, PC_fruits, PV_fruits, reallocation, Tf, RmD, senescence,ratio_leaf_shoots, ratio_leaf_stem, LADv
Integer          :: Year, Doy, PhenStage
Logical :: pruning

! This reads all the input files and initializes all arrays. It corresponds to the all the code that appear before the daily loop in maespa
call maespa_initialize
! This reads all the input files for the growth equations
call read_growth_inputs

! The daily loop and its variables (Istart, Iday, Iend) correspond to the original daily loop in Maestra
IDAY = 1
DO WHILE (ISTART + IDAY <= IEND) ! start daily loop
! Determine phenological stage
    DOY = mod(iday, 365) + 1
    Year = iday/365 + 1
!    write(*,*) 'Year = ', Year, 'DOY = ', DOY
    If(DOY == 1) Then
        Age = Age + 1
        biomass_leaf2 = biomass_leaf1
        biomass_leaf2T = biomass_leaf2
        biomass_leaf1 = biomass_leaf0
        biomass_leaf0 = 0.0
    End if
    if(DOY .ge. DOYsenescence1 .AND. DOY .LT. DOYsenescence2) Then
        senescence = 1.0
    else
        senescence = 0.0
    end if

    if(DOY .LT. DOYPhen1) Then
        PhenStage = 1
    else if(DOY .LT. DOYPhen2) Then
        PhenStage = 2
    else if(DOY .LT. DOYPhen3) Then
        PhenStage = 3
    else if(DOY .LT. DOYPhen4) Then
        PhenStage = 4
    else if(DOY == DOYPhen4) Then
        PhenStage = 1
        Biomass_fruits = 0.0 ! Harvest
        ! Give the opportunity of not having pruning
        if(OptPrun == 1) Then
            ! Pruning for high density orchards is always applied but differs when H > Hmax or H <= Hmax
            if(trim(DensOpt) == 'H') Then
            ! Prunning may want to control maximum height
                if(H > Hmax) Then                                                      
                    H = Hmax  
                    Volume = Hmax*Rx*Ry*4.0/3.0/2.0*3.1416
                    Rx = min((0.75*Volume/cp/3.1416)**(1.0/3.0), d_alley/3.0)
                    Ry = min((0.75*Volume/cp/3.1416)**(1.0/3.0), d_row/2.0)                                                   
                    Volume = Hmax*Rx*Ry*4.0/3.0/2.0*3.1416 ! Because Rx and Ry may be limited by tree spacing  
                end if
            ! Calculate LADv and prune LAD if bigger. Use empirical rule from Orgaz et al. (2007)    
                if(Volume/(D_row*D_alley) < 0.5) then
                    LADv = 2.0
                else
                    LADv = 2.0 - 0.8*(Volume/(d_row*d_alley) - 0.5)/1.5
                end if      
                if(LAD > LADv) LAD = LADv              
            ! Update all state variables and canopy dimensions 
                LAI = Volume/(d_alley*d_row)*LAD  
                ratio_leaf_shoots =    Biomass_leaf/Biomass_shoots   
                ratio_leaf_stem   =    Biomass_leaf/Biomass_stem 
                Biomass_leaf = LAI/specific_leaf_area                                                
                Biomass_shoots = Biomass_leaf/ratio_leaf_shoots
                Biomass_stem = Biomass_leaf/ratio_leaf_stem      
            ! Prunning for super-high density orchards only applied when H > Hmax   
            else if (trim(DensOpt) == 'SH' .AND. H > Hmax) Then
            ! Update all state variables and canopy dimensions    
                H = Hmax  
                Volume = Hmax*4.0*Rx*Ry
                Rx = min((0.125*Volume/cp)**(1.0/3.0), d_alley/3.0)
                Ry = min((0.125*Volume/cp)**(1.0/3.0), d_row/2.0)                                                       
                Volume = Hmax*4.0*Rx*Ry ! Because Rx and Ry may be limited by tree spacing                                                
                LAI = Volume*LAD/d_alley/d_row 
                ratio_leaf_shoots =    Biomass_leaf/Biomass_shoots   
                ratio_leaf_stem   =    Biomass_leaf/Biomass_stem                                  
                Biomass_leaf = LAI/specific_leaf_area                                                
                Biomass_shoots = Biomass_leaf/ratio_leaf_shoots
                Biomass_stem = Biomass_leaf/ratio_leaf_stem     
            end if
        end if
    else
        PhenStage = 1
    end if

! Override the arrays with dimensions of all the trees in the orchard
    RYTABLE1(1,:) = Rx ! Radius of green crown in y direction
    RXTABLE1(1,:) = Ry ! Radius of green crown in x direction
    RZTABLE1(1,:) = H ! Diameter of green crown in z direction   
    FOLTABLE1(1,:) = LAD*Volume ! Total leaf area of a crown.
    ZBCTABLE1(1,:) = max(H/2.0,0.5) ! Trunk height (i.e. from groun to first green branch)

! Call the main code of maespa. This corresponds to all the code that was located within the daily loop of maespa
    call run_maespa

! Carbon allocation
    Select Case (PhenStage)
        Case (1)
            Allocation_fruits = 0.0
            Allocation_leaf = 0.0
            Allocation_shoots = 0.0
            Allocation_stem = 0.0
            Allocation_froots = 0.0
            Allocation_croots = 0.0
            Allocation_reserves = 1.0
            reallocation        = 0.0
        Case (2)
            Allocation_fruits = 0.
            Allocation_leaf = PC_leaf
            Allocation_shoots = PC_shoots
            Allocation_stem = PC_stem
            Allocation_froots = PC_froots
            Allocation_croots = PC_croots
            Allocation_reserves = 0.0
            reallocation        = 1.0            
        Case (3)
            if(Age .GE. Adult) Then        
                Allocation_fruits  = PCfr
            else
                Allocation_fruits = 0.0
            End if             
            Allocation_leaf = PC_leaf*(1.0 - Allocation_fruits)
            Allocation_shoots = PC_shoots*(1.0 - Allocation_fruits)
            Allocation_stem = PC_stem*(1.0 - Allocation_fruits)
            Allocation_froots = PC_froots*(1.0 - Allocation_fruits)
            Allocation_croots = PC_croots*(1.0 - Allocation_fruits)
            Allocation_reserves = 0.0
            PV_fruits = PVfr
            reallocation        = 0.0            
        Case (4)
            if(Age > Adult) Then
                Allocation_fruits = PCoil
            else
                Allocation_fruits = 0.0
            End if
            Allocation_leaf = PC_leaf*(1.0 - Allocation_fruits)
            Allocation_shoots = PC_shoots*(1.0 - Allocation_fruits)
            Allocation_stem = PC_stem*(1.0 - Allocation_fruits)
            Allocation_froots = PC_froots*(1.0 - Allocation_fruits)
            Allocation_croots = PC_croots*(1.0 - Allocation_fruits)
            Allocation_reserves = 0.0
            PV_fruits           = PVoil
            reallocation        = 0.0            
    End Select

! Pool of assimilates produce on a given day
    Assimilation = (totCO2(1) + totRespf(1))/d_alley/d_row*12.0
! Maintenance respiration. Temperature factor is the avereage temperature effect over the day
    Tf           = sum(Q10**((TAIR(1:KHRS) - 25.0)/10.0))/KHRS
    RmD          = (Biomass_leaf*RmRef_leaf + Biomass_shoots*RmRef_shoots + Biomass_stem*RmRef_stem + Biomass_froots*RmRef_froots + Biomass_croots*RmRef_croots + Biomass_fruits*RmRef_fruits + Reserves*RmRef_reserves)*Tf*24.0
    If(RmD > Assimilation + reallocation*ReservesT/(DOYPhen2 - DOYPhen1)*CCres) RmD = Assimilation + reallocation*ReservesT/(DOYPhen2 - DOYPhen1)*CCres ! Source limitation of photosynthesis 
    Pool         = Assimilation + reallocation*ReservesT/(DOYPhen2 - DOYPhen1)*CCres - RmD

! Update the state variables 
    Biomass_leaf = Biomass_leaf + Pool*Allocation_leaf*PV_leaf/d_alley/d_row  - senescence*biomass_leaf2T/(DOYsenescence2 - DOYsenescence1)
    biomass_leaf0       = Biomass_leaf + Pool*Allocation_leaf*PV_leaf/d_alley/d_row 
    Volume       = Volume       + Pool*Allocation_leaf*PV_leaf*specific_leaf_area/LAD
    LAI          = Volume/d_alley/d_row*LAD
    If(trim(DensOpt) == trim('SH')) THEN
      Rx = min((0.125*Volume/cp)**(1.0/3.0), d_alley/3.0)
      Ry = min((0.125*Volume/cp)**(1.0/3.0), d_row/2.0)
      H  = Volume/4.0/Rx/Ry
    else if(trim(DensOpt) == 'H') Then      
     Rx  = min((0.75*Volume/cp/3.1416)**(1.0/3.0), d_alley/3.0)
     Ry  = min((0.75*Volume/cp/3.1416)**(1.0/3.0), d_row/2.0)
     H   = 0.75*Volume/3.1416/Rx/Ry*2.0
    End if

    Biomass_shoots = Biomass_shoots + Pool*Allocation_shoots*PV_shoots/d_alley/d_row
    Biomass_stem   = Biomass_stem   + Pool*Allocation_stem*PV_stem/d_alley/d_row
    Biomass_froots = Biomass_froots + Pool*Allocation_froots*PV_froots/d_alley/d_row - Biomass_froots*Kroots
    Biomass_croots = Biomass_croots + Pool*Allocation_croots*PV_croots/d_alley/d_row
    Biomass_fruits = Biomass_fruits + Pool*Allocation_fruits*PV_fruits/d_alley/d_row
    Reserves       = Reserves       + Pool*Allocation_reserves*PVres/d_alley/d_row   - reallocation*ReservesT/(DOYPhen2 - DOYPhen1)
    if(PhenStage == 1) ReservesT = Reserves

! Write something to the output
    Call  write_growth_outputs(Year, DOY)

! Proceed to the next step. Borrowed from the original daily loop of maespa
    IDAY = IDAY + NSTEP
    IF ((ISTART+IDAY) <= IEND) THEN
        CALL SKIPMET(MFLAG,NSTEP)
    END IF
        
END DO ! End daily loop

! This closes files and does other finalization operations. It corresponds to the all the code that appear after the daily loop in maespa
call maespa_finalize
call growth_finalize

End Program maespa_growth

