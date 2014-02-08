Program maespa_growth
Use Initialize, only: maespa_initialize, maespa_finalize
USE maindeclarations ! This contains all variables that were defined in the program unit of maespa.
Use growth_module ! This containts the parameter and state variables and the growth module
Implicit None

! This reads all the input files and initializes all arrays. It corresponds to the all the code that appear before the daily loop in maespa
call maespa_initialize
! This reads all the input files for the growth equations
call read_growth_inputs

! The daily loop and its variables (Istart, Iday, Iend) correspond to the original daily loop in Maestra
IDAY = 0
DO WHILE (ISTART + IDAY <= IEND) ! start daily loop

! Override the arrays with dimensions of all the trees in the orchard
    RYTABLE1(1,:) = 0.37050426492714 ! Radius of green crown in y direction
    RXTABLE1(1,:) = 0.37050426492714 ! Radius of green crown in x direction
    RZTABLE1(1,:) = 0.852159809332423 ! Diameter of green crown in z direction   
    FOLTABLE1(1,:) = 0.49 ! Total leaf area of a crown
    ZBCTABLE1(1,:) = 0.5 ! Trunk height (i.e. from groun to first green branch)

! Call the main code of maespa. This corresponds to all the code that was located within the daily loop of maespa
    call run_maespa

! Just to check: retrieve daily variables and print them
    write(*,*) "Daily gross photos is", totCO2(1) + totrespf(1)
    write(*,*) "Daily respiration is", totrespf(1)
    write(*,*) "Daily absorbed PAR is", tdyab(1,1)

! Proceed to the next step. Borrowed from the original daily loop of maespa
    IDAY = IDAY + NSTEP
    IF ((ISTART+IDAY) <= IEND) THEN
        CALL SKIPMET(MFLAG,NSTEP)
    END IF
        
END DO ! End daily loop

write(*,*) OptPrun
write(*,*) DOYPhen2
write(*,*) RmRef_stem, Biomass_leaf
write(*,*) DensOpt

! This closes files and does other finalization operations. It corresponds to the all the code that appear after the daily loop in maespa
call maespa_finalize
call growth_finalize

End Program maespa_growth

