Program maespa_growth
Use Initialize, only: maespa_initialize, maespa_finalize
USE maindeclarations
Implicit None

call maespa_initialize

IDAY = 0
DO WHILE (ISTART + IDAY <= IEND) ! start daily loop
! Override the dimensions of all the trees in the orchard
    RYTABLE1(1,:) = 0.37050426492714 ! Radius of green crown in y direction
    RXTABLE1(1,:) = 0.37050426492714 ! Radius of green crown in x direction
    RZTABLE1(1,:) = 0.852159809332423 ! Diameter of green crown in z direction   
    FOLTABLE1(1,:) = 0.49 ! Total leaf area of a crown
    ZBCTABLE1(1,:) = 0.5 ! Trunk height (i.e. from groun to first green branch)
! Calculate photosynthesis
    call run_maespa
! Retrieve daily variables
    write(*,*) "Daily gross photos is", totCO2(1) + totrespf(1)
    write(*,*) "Daily respiration is", totrespf(1)
    write(*,*) "Daily absorbed PAR is", tdyab(1,1)

    IDAY = IDAY + NSTEP
    IF ((ISTART+IDAY) <= IEND) THEN
        CALL SKIPMET(MFLAG,NSTEP)
    END IF
        
END DO ! End daily loop

call maespa_finalize

End Program maespa_growth