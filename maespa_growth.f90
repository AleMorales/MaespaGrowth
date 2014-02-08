Program maespa_growth
Use Initialize, only: maespa_initialize, maespa_finalize
USE maindeclarations, only: istart, iday, iend, mflag, nstep, RYTABLE1, RXTABLE1, RZTABLE1, FOLTABLE1, ZBCTABLE1, totRespf, totCO2 ! This contains all variables that were defined in the program unit of maespa.
Use growth_module ! This containts the parameter and state variables and the growth module
Implicit None
Double precision :: Assimilation

! This reads all the input files and initializes all arrays. It corresponds to the all the code that appear before the daily loop in maespa
call maespa_initialize
! This reads all the input files for the growth equations
call read_growth_inputs

! The daily loop and its variables (Istart, Iday, Iend) correspond to the original daily loop in Maestra
IDAY = 0
DO WHILE (ISTART + IDAY <= IEND) ! start daily loop

! Override the arrays with dimensions of all the trees in the orchard
    RYTABLE1(1,:) = Rx ! Radius of green crown in y direction
    RXTABLE1(1,:) = Ry ! Radius of green crown in x direction
    RZTABLE1(1,:) = H ! Diameter of green crown in z direction   
    FOLTABLE1(1,:) = LAD*Volume ! Total leaf area of a crown.
    ZBCTABLE1(1,:) = max(H/2.0,0.5) ! Trunk height (i.e. from groun to first green branch)
! Call the main code of maespa. This corresponds to all the code that was located within the daily loop of maespa
    call run_maespa

! Update biomass
    Assimilation = (totCO2(1) + totRespf(1))/d_alley/d_row*12.0
! Update the state variables    
    Biomass_leaf = Biomass_leaf + Assimilation*PC_leaf*PV_leaf/d_alley/d_row    
    Volume       = Volume       + Assimilation*PC_leaf*PV_leaf*specific_leaf_area/LAD
    LAI          = LAI          + Assimilation*PC_leaf*PV_leaf/d_alley/d_row *specific_leaf_area
    If(trim(DensOpt) == trim('H')) THEN
      Rx = min((0.125*Volume/cp)**(1.0/3.0), d_alley/3.0)
      Ry = min((0.125*Volume/cp)**(1.0/3.0), d_row/2.0)
      H  = Volume/4.0/Rx/Ry
    End if
! Write something to the output
    Call  write_growth_outputs
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

