Module growth_module
Implicit None
INTEGER :: IOERROR
double precision :: RmRef_leaf, RmRef_shoots, RmRef_stem, RmRef_froots, RmRef_croots, RmRef_fruits, PC_leaf, PC_shoots, PC_stem, PC_froots, PC_croots, RF_leaf, RF_shoots, RF_stem, RF_froots, RF_croots, PV_leaf, PV_shoots, PV_stem, PV_froots, PV_croots, CC_leaf, CC_shoots, CC_stem, CC_froots, CC_croots, sla, LAD, Q10, cp, k, RUEg, Hmax, Kroots, CCfr, CCoil, RFfr, RFoil, PGoil, PGfr, Kf, PVfr, PVoil, FractionReproductive, OneLeaf,  PVres, RFres, CCres, PCfr, PCoil, Biomass_leaf, Biomass_shoots, Biomass_stem, Biomass_froots, Biomass_croots, Biomass_fruits, V, LAI, H, Rx, Ry, Reserves, ReservesT, Biomass_leaf0, Biomass_leaf1, Biomass_leaf2, Biomass_leaf2T
integer :: OptPrun, DOYPhen1, DOYPhen2, DOYPhen3, DOYPhen4, Adult, DOYsenescence1, DOYsenescence2, Age
character(LEN = 10) :: DensOpt, OptFruit

Namelist /growth_parameters/ RmRef_leaf, RmRef_shoots, RmRef_stem, RmRef_froots, RmRef_croots, RmRef_fruits, PC_leaf, PC_shoots, PC_stem, PC_froots, PC_croots, RF_leaf, RF_shoots, RF_stem, RF_froots, RF_croots, PV_leaf, PV_shoots, PV_stem, PV_froots, PV_croots, CC_leaf, CC_shoots, CC_stem, CC_froots, CC_croots, sla, LAD, Q10, cp, k, RUEg, Hmax, DensOpt, DOYPhen1, DOYPhen2, DOYPhen3, DOYPhen4, Kroots, CCfr, CCoil, RFfr, RFoil, PGoil, PGfr, Kf, PVfr, PVoil, Adult, FractionReproductive, OneLeaf, DOYsenescence1, DOYsenescence2, PVres, RFres, CCres, PCfr, PCoil, OptFruit, OptPrun

Namelist /growth_initial/ Biomass_leaf, Biomass_shoots, Biomass_stem, Biomass_froots, Biomass_croots, Biomass_fruits, V, LAI, H, Rx, Ry, Reserves, ReservesT, Biomass_leaf0, Biomass_leaf1, Biomass_leaf2, Biomass_leaf2T, Age

contains


subroutine read_growth_inputs
USE maestcom, only: ugrowth, out_growth, iwarn
Implicit none

! Open the input file called 'Growth.dat'   
    OPEN(UGROWTH, FILE= 'Growth.dat', STATUS = 'OLD', IOSTAT = IOERROR)
! Trigger a warning if the file does not exist (model runs with default values)
    IF(IOERROR.NE.0) THEN
        close(ugrowth)
        CALL SUBERROR('ERROR: Growth.dat DOES NOT EXIST', IWARN,IOERROR)
    ENDIF
! Read the parameters
    READ(UGROWTH, nml = growth_parameters, IOSTAT = IOERROR)
! Trigger a warning if there is a problem reading the namelist
    IF(IOERROR.NE.0) THEN
        close(ugrowth)
        CALL SUBERROR('ERROR: namelist growth_parameters was not red properly', IWARN,IOERROR)
    ENDIF       
! Read the initial values for the state and auxilliary variables
	Read(ugrowth, nml = growth_initial, IOSTAT = IOERROR)    
! Trigger a warning if there is a problem reading the namelist
    IF(IOERROR.NE.0) THEN
        close(ugrowth)
        CALL SUBERROR('ERROR: namelist growth_initial was not red properly', IWARN,IOERROR)
    ENDIF    
    close(ugrowth)

! Open the output file called 'out_growth.dat'   
    OPEN(out_growth, FILE= 'out_growth.dat', STATUS = 'REPLACE', IOSTAT = IOERROR)
    IF(IOERROR.NE.0) THEN
    	close(out_growth)
        CALL SUBERROR('ERROR: out_growth.dat does not exist', IWARN,IOERROR)
    ENDIF 
    Write(out_growth, '(14A15)') 'Year','DOY','Leaf','Shoots','Stem','Fine_roots','Coarse_roots','Fruits','Reserves','Volume','LAI','Height','Radius_x','Radius_y'
end subroutine read_growth_inputs

subroutine growth_finalize
Use maestcom, only: out_growth
Implicit None
	close(out_growth)
end subroutine growth_finalize

End module growth_module

