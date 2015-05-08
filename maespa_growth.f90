Program maespa_growth
Use Initialize, only: maespa_initialize, maespa_finalize
USE maindeclarations, only: istart, iday, iend, mflag, nstep, RYTABLE1, RXTABLE1, RZTABLE1, FOLTABLE1, ZBCTABLE1, totRespf, totCO2, TAIR, KHRS, WSOILMETHOD, ISMAESPA, TDYAB, RADABV, SPERHR ! This contains all variables that were defined in the program unit of maespa.
Use growth_module ! This containts the parameter and state variables and the growth module
Implicit None
Double precision :: Assimilation, Pool, Allocation_leaf, Allocation_stem, Allocation_froots, Allocation_croots, Allocation_fruits, Allocation_reserves, PC_fruits, PV_fruits, reallocation, Tf, RmD, senescence, ratio_leaf_stem, LADv, cohort0, cohort1, cohort2, ChillingHours_day, ratio_leaf_froots, root_loss, residues
Integer          :: Year, Doy, PhenStage, i
Logical :: pruning
double precision, parameter :: pi = 3.14159265359

! This reads all the input files and initializes all arrays. It corresponds to the all the code that appear before the daily loop in maespa
call maespa_initialize
WSOILMETHOD = 1 ! Just to make sure
ISMAESPA = .FALSE. ! Just to make sure
! This reads all the input files for the growth equations
call read_growth_inputs

! The daily loop and its variables (Istart, Iday, Iend) correspond to the original daily loop in Maestra
IDAY = 0
DO WHILE (ISTART + IDAY <= IEND) ! start daily loop

! Determine phenological stage
    DOY = mod(iday, 365) + 1
    Year = iday/365 + 1

    ! Update age and leaf cohorts
    If(DOY == 1) Then
        Age = Age + 1
        biomass_leaf2 = biomass_leaf1
        biomass_leaf2T = biomass_leaf2
        biomass_leaf1 = biomass_leaf0
        biomass_leaf0 = 0.0
    End if

    ! Turn on senescence in the case when senescence concentrates on a period of the year
    if(DOY .ge. DOYsenescence1 .AND. DOY .LT. DOYsenescence2) Then
        senescence = 1.0
    else
        senescence = 0.0
    end if

    ! Restart phenological variables used for flowering date
    if (DOY == DOYPhen4) Then
        ChillingHours = 0.0
        ThermalTime = 0.0
    end if

    ! Phenological stages plus harvest/pruning.
    ! Note that thermaltime and chilling hours are calculated after this code, so they refer to the previous day of simulation
    ! One could think of this code happening at the beginning of each day and temperature-related and carbon allocation happening at the end of the day
    ! It also means that harvesting and pruning is done "at the beginning of the day"
    ! Note that harvesting and vernalization start
    if(DOY .LT. DOYPhen1) Then
        PhenStage = 1
        residues = 0
    else if(DOY .LT. DOYPhen2) Then
    !else if((ChillingHours < ColdRequirement .OR. ThermalTime < ThermalTimeRequirement) .AND. DOY < DOYPhen3) then
        PhenStage = 2
        residues = 0
    else if(DOY .LT. DOYPhen3) Then
    !else if(ChillingHours > ColdRequirement .AND. ThermalTime > ThermalTimeRequirement .AND. DOY < DOYPhen3) then
        PhenStage = 3
        residues = 0
    else if(DOY .LT. DOYPhen4) Then
    !else if(DOY .GE. DOYPhen3 .AND. DOY .LT. DOYPhen4) Then
        PhenStage = 4
        residues = 0
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
                    Volume = Hmax*Rx*Ry*4.0/3.0/2.0*pi
                    Rx = min((0.75*Volume/cp/pi)**(1.0/3.0), d_alley/3.0)
                    Ry = min((0.75*Volume/cp/pi)**(1.0/3.0), d_row/2.0)
                    Volume = Hmax*Rx*Ry*4.0/3.0/2.0*pi ! Because Rx and Ry may be limited by tree spacing
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
                cohort0 = biomass_leaf0/biomass_leaf
                cohort1 = biomass_leaf1/biomass_leaf
                cohort2 = biomass_leaf2/biomass_leaf
                ratio_leaf_stem   =    Biomass_leaf/Biomass_stem
                ratio_leaf_froots   =    Biomass_leaf/Biomass_froots
                residues = Biomass_leaf - LAI/specific_leaf_area
                Biomass_leaf = LAI/specific_leaf_area
                Biomass_leaf0 = cohort0*Biomass_leaf
                Biomass_leaf1 = cohort1*Biomass_leaf
                Biomass_leaf2 = cohort2*Biomass_leaf
                if(1.0/ratio_leaf_stem > ActiveWood) then
                  residues = residues + Biomass_stem - Biomass_leaf*ActiveWood
                  Biomass_stem = Biomass_leaf*ActiveWood
                end if
                residues = residues + Biomass_froots - Biomass_leaf/ratio_leaf_froots
                Biomass_froots = Biomass_leaf/ratio_leaf_froots
            ! Prunning for super-high density orchards only applied when H > Hmax
            else if (trim(DensOpt) == 'SH' .AND. H > Hmax) Then
            ! Update all state variables and canopy dimensions
                H = Hmax
                Volume = Hmax*4.0*Rx*Ry
                Rx = min((0.125*Volume/cp)**(1.0/3.0), d_alley/3.0)
                Ry = min((0.125*Volume/cp)**(1.0/3.0), d_row/2.0)
                Volume = Hmax*4.0*Rx*Ry ! Because Rx and Ry may be limited by tree spacing
                LAI = Volume*LAD/d_alley/d_row
                cohort0 = biomass_leaf0/biomass_leaf
                cohort1 = biomass_leaf1/biomass_leaf
                cohort2 = biomass_leaf2/biomass_leaf
                ratio_leaf_stem   =    Biomass_leaf/Biomass_stem
                ratio_leaf_froots   =    Biomass_leaf/Biomass_froots
                residues = Biomass_leaf - LAI/specific_leaf_area
                Biomass_leaf = LAI/specific_leaf_area
                Biomass_leaf0 = cohort0*Biomass_leaf
                Biomass_leaf1 = cohort1*Biomass_leaf
                Biomass_leaf2 = cohort2*Biomass_leaf
                residues = residues + Biomass_stem - Biomass_leaf/ratio_leaf_stem + Biomass_froots - Biomass_leaf/ratio_leaf_froots
                Biomass_stem = Biomass_leaf/ratio_leaf_stem
                Biomass_froots = Biomass_leaf/ratio_leaf_froots
            end if
        end if
    else if (DOY .GT. DOYPhen4) then
        PhenStage = 1
    end if

! Override the arrays with dimensions of all the trees in the orchard
    RYTABLE1(1,:) = Ry ! Radius of green crown in y direction
    RXTABLE1(1,:) = Rx ! Radius of green crown in x direction
    RZTABLE1(1,:) = H ! Diameter of green crown in z direction
    FOLTABLE1(1,:) = LAD*Volume ! Total leaf area of a crown.
    ZBCTABLE1(1,:) = min(H/2.0,0.5) ! Trunk height (i.e. from ground to first green branch). Not really much effect on simulation (unless we want info on the ground)

! Call the main code of maespa. This corresponds to all the code that was located within the daily loop of maespa
! Note that this code will update the arrays containing weather information, so all calculations dependent on temperature should appear afterwards.
    call run_maespa

    ! Update chilling hours
    if(ChillingHours < ColdRequirement) Then
        ! Code to update chilling hours
        ChillingHours_day = 0
        Do i = 1, KHRS ! KHRS are the number of steps within a day (e.g 48 -> dt = 0.5 h)
            if(Tair(i) > 0.0 .AND. Tair(i) .LE. Phen_T0) Then
                ChillingHours_day = ChillingHours_day + Tair(i)/Phen_T0*24.0/KHRS ! Note 24.0 and not 24 to coerce to floating point.
            else if (Tair(i) .LE. Phen_Tx) Then
                ChillingHours_day = ChillingHours_day + (1.0 - (Tair(i) - Phen_T0)*(1.0 - Phen_a)/(Phen_Tx - Phen_T0))*24.0/KHRS
            else if (Tair(i) > Phen_Tx) Then
                ChillingHours_day = ChillingHours_day + Phen_a*24.0/KHRS
            else
                ChillingHours_day = ChillingHours_day
            end if
        End Do
        ChillingHours = ChillingHours + ChillingHours_day
        IF(ChillingHours < 0.0) ChillingHours = 0.0 ! They cannot be negative according to this model
    else if(ChillingHours > ColdRequirement .AND. ThermalTime < ThermalTimeRequirement) Then
        ! Code to update thermal time. Parameter based on daily average temperature, so use that
        ThermalTime = ThermalTime + max(sum(Tair(1:KHRS))/KHRS - Phen_Tb, 0.0)
    end if
! Carbon allocation
    Select Case (PhenStage)
        Case (1)
            Allocation_fruits = 0.0
            Allocation_leaf = 0.0
            Allocation_stem = 0.0
            Allocation_froots = 0.0
            Allocation_croots = 0.0
            Allocation_reserves = 1.0
            reallocation        = 0.0
        Case (2)
            Allocation_fruits = 0.
            Allocation_leaf = PC_leaf
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
            Allocation_stem = PC_stem*(1.0 - Allocation_fruits)
            Allocation_froots = PC_froots*(1.0 - Allocation_fruits)
            Allocation_croots = PC_croots*(1.0 - Allocation_fruits)
            Allocation_reserves = 0.0
            PV_fruits = PVfr
            reallocation        = 0.0
        Case (4)
            if(Age .GE. Adult) Then
                Allocation_fruits = PCoil
            else
                Allocation_fruits = 0.0
            End if
            Allocation_leaf = PC_leaf*(1.0 - Allocation_fruits)
            Allocation_stem = PC_stem*(1.0 - Allocation_fruits)
            Allocation_froots = PC_froots*(1.0 - Allocation_fruits)
            Allocation_croots = PC_croots*(1.0 - Allocation_fruits)
            Allocation_reserves = 0.0
            PV_fruits           = PVoil
            reallocation        = 0.0
    End Select

! Pool of assimilates produce on a given day (g C (m2 ground)-2)
    Assimilation = (totCO2(1) + totRespf(1))/d_alley/d_row*12.0
! Maintenance respiration. ! Average maintenance respiration on a daily basis (g C (m2 ground)-2)
    RmD          = sum((Biomass_leaf*RmRef_leaf + min(Biomass_stem, Biomass_leaf*ActiveWood)*RmRef_stem + Biomass_froots*RmRef_froots + min(Biomass_croots, Biomass_froots*ActiveWood)*RmRef_croots + Biomass_fruits*RmRef_fruits + Reserves*RmRef_reserves)*Q10**((TAIR(1:KHRS) - 25.0)/10.0))*24.0/KHRS
! Source-limited maintenanced respiration
    If(RmD > Assimilation + reallocation*ReservesT/(DOYPhen2 - DOYPhen1)*CCres) Then
      RmD = Assimilation + reallocation*ReservesT/(DOYPhen2 - DOYPhen1)*CCres
! Pool of assimilates (g C (m ground)-2)
      Pool      = 0.0
    else
      Pool = Assimilation + reallocation*ReservesT/(DOYPhen2 - DOYPhen1)*CCres - RmD
    end if
! Update the state variables
! Leaf biomass (g dm (m ground)-2)
    Biomass_leaf = Biomass_leaf + Pool*Allocation_leaf*PV_leaf - senescence*biomass_leaf2T/(DOYsenescence2 - DOYsenescence1)
    biomass_leaf0  = biomass_leaf0 + Pool*Allocation_leaf*PV_leaf
    biomass_leaf2 = biomass_leaf2 - senescence*biomass_leaf2T/(DOYsenescence2 - DOYsenescence1)
! Volume of an individual crown (m3)
    Volume       = Volume       + (Pool*Allocation_leaf*PV_leaf - senescence*biomass_leaf2T/(DOYsenescence2 - DOYsenescence1))*specific_leaf_area/LAD*d_alley*d_row
! Leaf area index (m2 leaf (m ground)-2)
    LAI          = Volume/d_alley/d_row*LAD
! Radius x and y(m). Diameter z (m)
    If(trim(DensOpt) == trim('SH')) THEN
      Rx = min((0.125*Volume/cp)**(1.0/3.0), d_alley/3.0)
      Ry = min((0.125*Volume/cp)**(1.0/3.0), d_row/2.0)
      H  = Volume/4.0/Rx/Ry
    else if(trim(DensOpt) == 'H') Then
     Rx  = min((0.75*Volume/cp/pi)**(1.0/3.0), d_alley/3.0)
     Ry  = min((0.75*Volume/cp/pi)**(1.0/3.0), d_row/2.0)
     H   = 0.75*Volume/pi/Rx/Ry*2.0
    End if
! Biomass of stem, froots, croots, fruits and reserves (g C (m ground)-2)
    Biomass_stem   = Biomass_stem   + Pool*Allocation_stem*PV_stem
    Biomass_froots = Biomass_froots + Pool*Allocation_froots*PV_froots - Biomass_froots*Kroots
    root_loss      = Biomass_froots*Kroots
    Biomass_croots = Biomass_croots + Pool*Allocation_croots*PV_croots
    Biomass_fruits = Biomass_fruits + Pool*Allocation_fruits*PV_fruits
    Reserves       = Reserves       + Pool*Allocation_reserves*PVres  - reallocation*ReservesT/(DOYPhen2 - DOYPhen1)
    if(PhenStage == 1) ReservesT = Reserves

! Write state variables and fluxes to the output
    Call  write_growth_outputs(Year, DOY, Assimilation, RmD, TDYAB(1,1)/d_alley/d_row, sum(RADABV(1:KHRS, 1))*SPERHR/1e6, residues, root_loss, sum(Tair(1:KHRS))/KHRS)

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
