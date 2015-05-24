Module growth_module
Implicit None
INTEGER :: IOERROR
double precision :: RmRef_leaf, RmRef_stem, RmRef_froots, RmRef_croots, RmRef_fruits, PC_leaf, PC_stem, PC_froots, PC_croots, RF_leaf, RF_stem, RF_froots, RF_croots, PV_leaf, PV_stem, PV_froots, PV_croots, CC_leaf, CC_stem, CC_froots, CC_croots, LAD, Q10, cp,Hmax, Kroots, CCfr, CCoil, RFfr, RFoil, PGoil, PGfr, Kf, PVfr, PVoil, FractionReproductive, OneLeaf,  PVres, CCres, PCfr, PCoil, d_alley, d_row, Biomass_leaf, Biomass_stem, Biomass_froots, Biomass_croots, Biomass_fruits, Volume, LAI, H, Rx, Ry, Reserves, ReservesT, Biomass_leaf0, Biomass_leaf1, Biomass_leaf2, Biomass_leaf2T, specific_leaf_area, RmRef_reserves, ColdRequirement, Phen_T0, Phen_Tx, Phen_a, ThermalTimeRequirementFl, ThermalTimeRequirementFr, Phen_TbFl, Phen_TbFr, ChillingHours, ThermalTime, ActiveWood, FruitMax, MaxPhotos
double precision :: TMaxDay, TMinDay, DayLength
integer :: OptPrun, DOYPhen1, DOYPhen2, DOYPhen3, DOYPhen4, Adult, DOYsenescence1, DOYsenescence2, Age, FruitOpt, WinterOpt, PhenSim
character(LEN = 10) :: DensOpt, OptFruit

Namelist /growth_parameters/ RmRef_leaf, RmRef_stem, RmRef_froots, RmRef_croots, RmRef_fruits, PC_leaf, PC_stem, PC_froots, PC_croots, PV_leaf, PV_stem, PV_froots, PV_croots, LAD, Q10, cp, Hmax, DensOpt, DOYPhen1, DOYPhen2, DOYPhen3, DOYPhen4, Kroots, PVfr, PVoil, Adult, DOYsenescence1, DOYsenescence2, PVres, CCres, PCfr, PCoil, OptPrun, d_alley, d_row, specific_leaf_area, RmRef_reserves, ThermalTimeRequirementFl, ThermalTimeRequirementFr, Phen_TbFl, Phen_TbFr, ColdRequirement, Phen_T0, Phen_Tx, Phen_a, ActiveWood, FruitOpt, FruitMax, WinterOpt, MaxPhotos, PhenSim

Namelist /growth_initial/ Biomass_leaf, Biomass_stem, Biomass_froots, Biomass_croots, Biomass_fruits, Volume, LAI, H, Rx, Ry, Reserves, ReservesT, Biomass_leaf0, Biomass_leaf1, Biomass_leaf2, Biomass_leaf2T, Age, ChillingHours, ThermalTime

contains


subroutine read_growth_inputs
USE maestcom, only: ugrowth, out_growth, ifatal
Implicit none

! Open the input file called 'Growth.dat'
    OPEN(UGROWTH, FILE= 'Growth.dat', STATUS = 'OLD', IOSTAT = IOERROR)
! Trigger a warning if the file does not exist (model runs with default values)
    IF(IOERROR.NE.0) THEN
        close(ugrowth)
        CALL SUBERROR('ERROR: Growth.dat DOES NOT EXIST', IFATAL,IOERROR)
    ENDIF
! Read the parameters
    READ(UGROWTH, nml = growth_parameters, IOSTAT = IOERROR)
! Trigger a warning if there is a problem reading the namelist
    IF(IOERROR.NE.0) THEN
        close(ugrowth)
        CALL SUBERROR('ERROR: namelist growth_parameters was not red properly', IFATAL,IOERROR)
    ENDIF
    rewind(UGROWTH)
! Read the initial values for the state and auxilliary variables
	Read(ugrowth, nml = growth_initial, IOSTAT = IOERROR)
! Trigger a warning if there is a problem reading the namelist
    IF(IOERROR.NE.0) THEN
        close(ugrowth)
        CALL SUBERROR('ERROR: namelist growth_initial was not red properly', IFATAL,IOERROR)
    ENDIF
    close(ugrowth)

! Open the output file called 'out_growth.dat'
    OPEN(out_growth, FILE= 'out_growth.dat', STATUS = 'OLD', IOSTAT = IOERROR,ACTION = 'WRITE')
    IF(IOERROR.NE.0) THEN
    	close(out_growth)
        CALL SUBERROR('ERROR: out_growth.dat does not exist', IFATAL,IOERROR)
    ENDIF
    Write(out_growth, '(23(A14,3x))') 'Year','DOY','Leaf','Stem','Fine_roots','Coarse_roots','Fruits', &
    'Reserves','Volume','LAI','Height','Radius_x','Radius_y', 'Assimilation','M.Respiration', 'Leaf_Loss', &
    'Residues','Root_Loss','aPAR','PAR','Chilling','Thermal','Tavg'
end subroutine read_growth_inputs

subroutine write_growth_outputs(Year, DOY, Assimilation, RmD, aPAR, PAR, residues, root_loss, Tavg)
USE maestcom, only: out_growth
Implicit None
integer, intent(in) :: Year, DOY
double precision, intent(in) :: Assimilation, RmD, aPAR, residues, root_loss, Tavg
real, intent(in) :: PAR
	Write(out_growth, '(I2,3x,I3,2x,21(ES16.6E3,3x))') Year, DOY, Biomass_leaf, Biomass_stem, Biomass_froots, Biomass_croots, Biomass_fruits, Reserves, Volume, LAI, H, Rx, Ry, Assimilation, RmD, Biomass_leaf2T, residues, root_loss, aPAR, PAR, ChillingHours, ThermalTime, Tavg
end subroutine

subroutine growth_finalize
Use maestcom, only: out_growth
Implicit None
	close(out_growth)
end subroutine growth_finalize

End module growth_module
