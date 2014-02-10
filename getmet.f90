!**********************************************************************
! GETMET.FOR
!
! This file contains all the subroutines required to read the met data
! file. The main routines (called externally to this file) are:
! OPENMETF - opens the met file & reads the format information
! RESTARTMETF - finds the start of the actual met data in the met file
! GETMET - calls the appropriate subroutine to read the met data, based
!          on the flag MFLAG describing the format
! SKIPMET - skips any unwanted met data
! ALTERMETOTC - alters met data to account for effects of OTC
! ALTERMETCC - alters met data in accordance with climate change scenario
! GETWIND - calculates wind speed through the canopy using exponential formula
!
! Subsidiary routines are:
! 1. Subroutines to read in data
! GETMETDAY - reads the standard daily met data format
! GETMETHR - reads the standard hourly met data format
! READENVIRON - reads in constant met data values
! READLAT - reads latitude and longitude
! READDELTAT - reads in information needed to calculate incident radiation, if needed
! 2. Subroutines to process temperatures & rainfall
! CALCTHRLY - calculates hourly air temps from Tmax & Tmin
! CALCTSOIL - sets hourly soil temp equal to average daily air temp
! ASSIGNRAIN - assigns daily rainfall to hours in the day
! 3. Subroutines to process VPD/RH
! CALCRH - calculates hourly RH from hourly air temp & the minimum temp
! RHTOVPD - converts RH to VPD
! VPDTORH - converts VPD to RH
! TDEWTORH - calculates RH from air temperature & dewpoint temperature
! MFDTORH - converts from water vapour mole fraction to RH
! VPDTOMFD - converts from VPD to water vapour mole fraction
! CALCAH - converts relative humidity to absolute humidity (called by ALTERMETCC)
! AHTORH - converts absolute humidity to relative humidity (called by ALTERMETCC)
! 4. Subroutines to process radiation
! ESTPARIN - estimates daily incident radiation from temperature
! CALPARHHRLY - calculates hourly incident PAR from daily totals of beam & diffuse
! CALCFBMD - calculates diffuse fraction of daily incident radiation
! CALCFBMH - calculates diffuse fraction of hourly incident radiation
! CALCFSUN - calculate sunlit fraction of hour (presently unused)
! CALCNIR - calculates incident NIR from incident PAR
! THERMAL - calculates incident thermal radiation from air temperature
! CALCFBMWN - an alternative method to calculate the diffuse fraction (currently unused)
! 5. Additionally, the following functions are provided:
! SATUR - calculate water vapour saturation pressure at given T MOVED TO UTILS
! TK - calculate T in Kelvin given T in Celsius MOVED TO UTILS
! ETRAD - calculate radiation incident on the atmosphere
!**********************************************************************

!**********************************************************************
SUBROUTINE OPENMETF(ISTART,IEND,CAK,PRESSK,SWMIN,SWMAX,USEMEASET,DIFSKY,ALAT,TTIMD,&
                    DELTAT,MFLAG,METCOLS,NOMETCOLS,MTITLE,MSTART,in_path)
! This subroutine opens the meteorological data file and reads the
! data describing it. It checks that it covers the dates of the simulation
! and returns the format of the met file.
! INPUTS:
! ISTART, IEND - Desired start and end of simulation (in days-since-1950) - from confile. If not defined,
!	program uses all met data. 
! OUTPUTS:
! CAK - Constant CO2 concentration (umol mol-1), if specified; 0 otherwise
! PRESSK - Constant air pressure (Pa), if specified; PATM otherwise
! SWMIN - Minimum soil water content
! SWMAX - Maximum soil water content - in same units as soil water content in met data
! DIFSKY - Parameter indicating light distribution over cloudy sky, if specified; default 0
! ALAT - latitude, in radians (-ve for S hemisphere)
! TTIMD - time difference between longitude of plot & longitude of time zone, in hours
! MFLAG - indicates whether daily (0) or hourly (1) met data
! METCOLS - array containing the number of the column corresponding to each met data element used by program
! NOMETCOLS - number of columns of met data
! MTITLE - title of met data file
! MSTART - start date of met data (in days-since-1950 format)
! DELTAT - monthly mean daily temperature amplitude
!**********************************************************************
    USE maestcom
    USE metcom
    IMPLICIT NONE
    CHARACTER*8 COLUMNS(MAXMET)
    CHARACTER*256 STARTDATE,ENDDATE
    CHARACTER(LEN=*), INTENT(INOUT) :: MTITLE
    CHARACTER(LEN=*), INTENT(IN) :: in_path
    REAL DELTAT(12)
    INTEGER METCOLS(MAXMET),KHRSPERDAY,NOCOLUMNS,IOERROR
    INTEGER, INTENT(IN) :: USEMEASET
    INTEGER, EXTERNAL :: IDATE50
    INTEGER DAYORHR,MSTART,MEND,ISTART,IEND,ICOL,MFLAG,I,NOMETCOLS
    NAMELIST /METFORMAT/ DAYORHR, KHRSPERDAY, NOCOLUMNS, COLUMNS, STARTDATE, ENDDATE
    REAL CAK,PRESSK,DIFSKY,SWMIN,SWMAX,ALAT,TTIMD
     
    990   FORMAT (A60)     ! For reading titles in input files.
    991   FORMAT (A12,A60) ! For writing comments to output files.
    
    ! Open input file with met data
    OPEN (UMET, FILE =trim(in_path)//'met.dat', STATUS='OLD', IOSTAT=IOERROR)
    IF (IOERROR.NE.0) THEN
        CALL SUBERROR('ERROR OPENING MET FILE',IFATAL,IOERROR)
    END IF
    
    READ (UMET, 990, IOSTAT=IOERROR) MTITLE
    IF (IOERROR.NE.0) THEN
        CALL SUBERROR('ERROR READING MET FILE',IFATAL,IOERROR)
    END IF
    
    ! Read in unchanging met data
    CALL READENVIRON(UMET, CAK, PRESSK, DIFSKY, SWMIN, SWMAX)

    ! Read in latitude & longitude
    CALL READLAT(UMET, ALAT, TTIMD)
    
    ! Read format
    KHRSPERDAY = 0  ! Init.
    REWIND(UMET)
    READ (UMET,METFORMAT,IOSTAT = IOERROR)
    IF (IOERROR.NE.0) CALL SUBERROR('ERROR IN MET FILE: MISSING FORMAT INFORMATION', IFATAL,IOERROR)
    
    ! Check dates of simulation are covered
    MSTART = IDATE50(STARTDATE)
    MEND = IDATE50(ENDDATE)
    
    IF ((ISTART.EQ.0).AND.(IEND.EQ.0)) THEN
        ISTART = MSTART
        IEND = MEND
    ELSE
    
        IF ((ISTART.LT.MSTART).OR.(IEND.GT.MEND))CALL SUBERROR('ERROR: MET FILE DOES NOT COVER DATES',IFATAL,IOERROR)
    END IF
      
    ! Process number of hours per day - the old KHRS
    IF (KHRSPERDAY.EQ.0) CALL SUBERROR('ERROR: ADD KHRSPERDAY TO MET.DAT FORMAT',IFATAL,KHRSPERDAY)
    IF (KHRSPERDAY.GT.MAXHRS) CALL SUBERROR('ERROR: KHRSPERDAY MUST NOT EXCEED MAXHRS',IFATAL,KHRSPERDAY)
    
    KHRS = KHRSPERDAY
    HHRS = (KHRS) / 2.0
    SPERHR = 3600 * 24.0 / KHRS

    ! Process columns - assign column numbers to METCOLS array
    IF (NOCOLUMNS.GT.MAXMET) THEN
        CALL SUBERROR('ERROR: TOO MANY MET DATA COLUMNS: DATA LOST',IWARN,IOERROR)
        NOCOLUMNS = MAXMET
    END IF

    ! Initialise array of met columns
    DO ICOL = 1,MAXMET
        METCOLS(ICOL) = MISSING
    END DO

    MFLAG = DAYORHR

    ! If MFLAG = 0, look for daily values; if MFLAG = 1, hourly values
    IF (MFLAG.EQ.0) THEN ! Daily values
        DO I = 1,NOCOLUMNS
            ICOL = MISSING
            IF (COLUMNS(I).EQ.'DATE')     THEN
                ICOL = MDDATE
            ELSEIF (COLUMNS(I).EQ.'WIND') THEN
                ICOL = MDWIND
            ELSEIF (COLUMNS(I).EQ.'TMIN') THEN
                ICOL = MDTMIN
            ELSEIF (COLUMNS(I).EQ.'TMAX') THEN
                ICOL = MDTMAX
            ELSEIF (COLUMNS(I).EQ.'PAR')  THEN
                ICOL = MDPAR
            ELSEIF (COLUMNS(I).EQ.'SI')  THEN
                ICOL = MDSI
            ELSEIF (COLUMNS(I).EQ.'FBEAM')  THEN
                ICOL = MDFBEAM
            ELSEIF (COLUMNS(I).EQ.'CA')  THEN
                ICOL = MDCA
            ELSEIF (COLUMNS(I).EQ.'PRESS')  THEN
                ICOL = MDPRESS
            ELSEIF (COLUMNS(I).EQ.'PPT')  THEN
                ICOL = MDPPT
            ELSEIF (COLUMNS(I).EQ.'SMD')  THEN
                ICOL = MDSMD
            ELSEIF (COLUMNS(I).EQ.'SWP')  THEN
                ICOL = MDSWP
            ELSEIF (COLUMNS(I).EQ.'SWC') THEN
                ICOL = MDSWC
            ELSEIF (COLUMNS(I).EQ.'ET')  THEN
                ICOL = MDET
            ELSEIF (COLUMNS(I).EQ.'EA')  THEN
                ICOL = MDEA                     
            ENDIF
            IF (ICOL.NE.MISSING) METCOLS(ICOL) = I
        END DO
        ! Check to see if any of the essential information is missing.
        IF (METCOLS(MDTMIN).EQ.MISSING) CALL SUBERROR('ERROR: NEED VALUES OF TMIN IN MET FILE',IFATAL,0)
        IF (METCOLS(MDTMAX).EQ.MISSING) CALL SUBERROR('ERROR: NEED VALUES OF TMAX IN MET FILE',IFATAL,0)
        IF ((METCOLS(MDPAR).EQ.MISSING).AND.(METCOLS(MDSI).EQ.MISSING)) CALL READDELTAT(UMET,DELTAT,IOERROR)
        IF ((CAK.EQ.0).AND.(METCOLS(MDCA).EQ.MISSING)) CALL SUBERROR('ERROR: NO VALUE FOR CA IN MET FILE',IFATAL,0)
        IF (USEMEASET.EQ.1.AND.(METCOLS(MDET).EQ.MISSING)) &
            CALL SUBERROR('ERROR: NEED MEASURED ET WHEN USEMEASET=1',IFATAL,0)  ! RAD
    ELSE                 ! Hourly values
        DO I = 1,NOCOLUMNS
            ICOL = MISSING
            IF (COLUMNS(I).EQ.'DATE') THEN
                ICOL = MHDATE
            ELSEIF (COLUMNS(I).EQ.'WIND') THEN
                ICOL = MHWIND
            ELSEIF (COLUMNS(I).EQ.'TAIR') THEN
                ICOL = MHTAIR
            ELSEIF (COLUMNS(I).EQ.'TSOIL') THEN
                ICOL = MHTSOIL
            ELSEIF (COLUMNS(I).EQ.'RH')  THEN
                ICOL = MHRH
            ELSEIF (COLUMNS(I).EQ.'RH%')  THEN
                ICOL = MHRHP
            ELSEIF (COLUMNS(I).EQ.'VPD')  THEN
                ICOL = MHVPD
            ELSEIF (COLUMNS(I).EQ.'PAR')  THEN
                ICOL = MHPAR
            ELSEIF (COLUMNS(I).EQ.'RAD')  THEN
                ICOL = MHRAD
            ELSEIF (COLUMNS(I).EQ.'FBEAM')  THEN
                ICOL = MHFBEAM
            ELSEIF (COLUMNS(I).EQ.'CA')  THEN
                ICOL = MHCA
            ELSEIF (COLUMNS(I).EQ.'MFD')  THEN
                ICOL = MHMFD
            ELSEIF (COLUMNS(I).EQ.'PRESS')  THEN
                ICOL = MHPRESS
            ELSEIF (COLUMNS(I).EQ.'TDEW')  THEN
                ICOL = MHTDEW
            ELSEIF (COLUMNS(I).EQ.'PPT')  THEN
                ICOL = MHPPT
            ELSEIF (COLUMNS(I).EQ.'ET')  THEN
                ICOL = MHET
            ELSEIF (COLUMNS(I).EQ.'SMD')  THEN
                ICOL = MHSMD
            ELSEIF (COLUMNS(I).EQ.'SWC')  THEN
                ICOL = MHSWC
            ELSEIF (COLUMNS(I).EQ.'SWP')  THEN
                ICOL = MHSWP
            ENDIF
            IF (ICOL.NE.MISSING) METCOLS(ICOL) = I
        END DO
        ! Check to see if any of the essential information is missing.
        IF (METCOLS(MHTAIR).EQ.MISSING) CALL SUBERROR('ERROR: NEED VALUES OF TAIR IN MET FILE',IFATAL,0)
        IF ((METCOLS(MHPAR).EQ.MISSING).AND.(METCOLS(MHRAD).EQ.MISSING))CALL READDELTAT(UMET,DELTAT,IOERROR)
        IF ((CAK.EQ.0).AND.(METCOLS(MHCA).EQ.MISSING))CALL SUBERROR('ERROR: NO VALUE FOR CA IN MET FILE',IFATAL,0)
    END IF
    NOMETCOLS = NOCOLUMNS
    RETURN
END SUBROUTINE OPENMETF


!**********************************************************************
SUBROUTINE RESTARTMETF(ISTART,MSTART,MFLAG)
! Position met file at start of met data
! BM Changed 10/00: Start of data must be directly specified
! INPUTS:
! ISTART - Date for which met data is wanted
! MSTART - Date on which met data starts 
! MFLAG - indicates whether daily or hourly data
!**********************************************************************
    USE maestcom
    IMPLICIT NONE
    INTEGER IOERROR,MFLAG,LINESTOSKIP,ISTART,MSTART,I
    CHARACTER*60 TMPTXT, DATASTART
    
    DATASTART = 'DATA STARTS'
    REWIND(UMET)

    ! Read through format information to find start of data
    990   FORMAT (A60)
    30    READ (UMET,990,IOSTAT = IOERROR) TMPTXT
    IF (IOERROR.NE.0) CALL SUBERROR('ERROR: COULD NOT FIND START OF MET DATA',IFATAL,IOERROR)
    IF (TMPTXT.NE.DATASTART) GOTO 30
   
    ! Read data until the start date
    IF (MFLAG.EQ.0) THEN
        LINESTOSKIP = ISTART-MSTART
    ELSE
        LINESTOSKIP = KHRS*(ISTART-MSTART)
    END IF

    DO I = 1,LINESTOSKIP
        READ (UMET,990,IOSTAT = IOERROR) TMPTXT
        IF (IOERROR.NE.0) CALL SUBERROR('READ PAST EOF IN RESTARTMETF', IFATAL,IOERROR)
    END DO
    RETURN
END SUBROUTINE RESTARTMETF


!**********************************************************************
SUBROUTINE READENVIRON(UFILE,CAK,PRESSK,DIFSKYI,SWMINI,SWMAXI)
! Read in environmental conditions which do not change with time.
! These must be at the start of the met.dat file.
! INPUT:
! UFILE - file number of met data file
! OUTPUTS:
! CAK - Constant CO2 concentration (umol mol-1), if specified; 0 otherwise
! PRESSK - Constant air pressure (Pa), if specified; PATM otherwise
! DIFSKYI - Parameter indicating light distribution over cloudy sky, if specified; default 0
! SWMINI - Minimum soil water content
! SWMAXI - Maximum soil water content - in same units as soil water content in met data
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    NAMELIST /ENVIRON/ CA,PRESS,DIFSKY, SWMIN, SWMAX
    INTEGER UFILE,IOERROR
    REAL CA,PRESS,DIFSKY, SWMIN, SWMAX
    REAL CAK,PRESSK,DIFSKYI, SWMINI, SWMAXI
    
    ! Default values
    PRESS = PATM
    CA = 0
    DIFSKY = 0.0
    SWMIN = 0.0
    SWMAX = 1.0

    ! Read namelist
    READ (UFILE, ENVIRON, IOSTAT = IOERROR)
    IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: DEFAULT VALUES: PRESS = 101.25 kPa, DIFSKY = 0', IWARN,IOERROR)
    END IF
    CAK = CA
    PRESSK = PRESS
    DIFSKYI = DIFSKY
    SWMINI = SWMIN
    SWMAXI = SWMAX

    RETURN
END SUBROUTINE READENVIRON

!**********************************************************************
SUBROUTINE READDELTAT(UFILE,DELTATI,IOERROR)
! Read in mean monthly daily temperature amplitudes
! Needed to calculate incident PAR (if not specified) using routine of Bristow & Campbell 1984
! If not present, program will stop (only called if PAR not specified). 
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    NAMELIST /BRISTO/ DELTAT
    REAL DELTAT(12),DELTATI(12)
    INTEGER UFILE,IOERROR,I

    ! Read namelist
    REWIND(UFILE)
    READ (UFILE, BRISTO, IOSTAT = IOERROR)
    IF (IOERROR.NE.0) THEN
        CALL SUBERROR('ERROR: NEED VALUES OF PAR IN MET FILE OR DELTAT', IFATAL,0)
    END IF

    DO I = 1,12
        DELTATI(I) = DELTAT(I)
    END DO

    RETURN
END SUBROUTINE READDELTAT

!**********************************************************************
SUBROUTINE GETMET(IDATE,MFLAG,ZEN,METCOLS,NOMETCOLS,CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,WINDAH,&
                    TSOIL,TAIR,RADABV,FBEAM,RH,VPD,VMFD,CA,PRESS,PPT,SOILMOIST,SOILDATA,TSOILDATA,ETMEAS)

! According to the value of MFLAG, this subroutine calls the appropriate
! subroutine to read in the meteorological data for use in MAESTRO.
! It is important that all handling of meteorological data occurs within
! this function. In particular, the flag MFLAG should not be referred to
! within the main body of the program. This maintains modularity & makes it
! easier to add new met data formats.

! The values passed to this function are:
!   IDATE - date of simulation, in years-since-1950 format
!   MFLAG - flag indicating the format of the input file
!   ZEN - the zenith angle of the sun
!   NOMETCOLS - the number of columns in the met file
!   METCOLS - the format of the columns in the met file
!   CAK - Constant CO2 concentration (umol mol-1), if specified; 0 otherwise
!   PRESSK - Constant air pressure (Pa), if specified; PATM otherwise
!   SWMINI - Minimum soil water content
!   SWMAXI - Maximum soil water content - in same units as soil water content in met data


! The values which must be returned are:
! arrays of KHRSly values of the following met variables:
!   WINDAH - wind speed in m s-1
!   TSOIL - soil surface temperature (degrees !)
!   TAIR - air temperature (degrees !)
!   RADABV(3) - incident radiation in three wavelengths (J m-2 s-1)
!   FBEAM(3) - beam fractions of radiation in 3 wavelengths (fraction)
!   RH - relative humidity of air (fraction)
!   VPD - vapour pressure deficit of air (Pa)
!   VMFD - mole fraction deficit (mmol mol-1)
!   CA - the atmospheric CO2 concentration (umol mol-1)
!   PRESS - the atmospheric pressure (Pa)
!   PPT - rainfall over the time period (mmol H2O m-2)
!   SOILMOIST - soil moisture data - either potential (MPa) or deficit (dimnless)
!   SOILDATA - type of soil moisture data
!   ET - measured evapotranspiration (mm t-1)
!**********************************************************************

    USE maestcom
    IMPLICIT NONE

    REAL ZEN(MAXHRS)
    REAL WINDAH(MAXHRS),TSOIL(MAXHRS),TAIR(MAXHRS)
    REAL RADABV(MAXHRS,3),FBEAM(MAXHRS,3)
    REAL RH(MAXHRS),VPD(MAXHRS),CA(MAXHRS),VMFD(MAXHRS),PRESS(MAXHRS)
    REAL PPT(MAXHRS), SOILMOIST(MAXHRS), ETMEAS(MAXHRS)
    REAL DELTAT(12)
    INTEGER METCOLS(MAXMET),SOILDATA,TSOILDATA,MFLAG,IDATE,NOMETCOLS
    REAL CAK,PRESSK,SWMIN,SWMAX,ALAT,DEC,DAYL

    ! Call appropriate function to read met file.
    IF (MFLAG.EQ.0) THEN
        CALL GETMETDAY(IDATE,ZEN,NOMETCOLS,METCOLS,CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,  &
                        DAYL,WINDAH,TSOIL,TAIR,RADABV,FBEAM,RH,VPD,VMFD,CA,PRESS,PPT,       &
                        SOILMOIST,SOILDATA,TSOILDATA)

    ELSE IF (MFLAG.EQ.1) THEN
        CALL GETMETHR(IDATE,ZEN,NOMETCOLS,METCOLS,CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,   &
                        DAYL,WINDAH,TSOIL,TAIR,RADABV,FBEAM,RH,VPD,VMFD,CA,PRESS,PPT,       &
                        SOILMOIST,SOILDATA,TSOILDATA,ETMEAS)

    !      ELSE
    ! insert other formats if required
    END IF
    RETURN
END SUBROUTINE GETMET


!**********************************************************************
SUBROUTINE GETMETDAY(IDATE,ZEN,NOMETCOLS,METCOLS,CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,   &
                        WINDAH,TSOIL,TAIR,RADABV,FBEAM,RH,VPD,VMFD,CA,PRESS,PPT,SOILMOIST,      &
                        SOILDATA,TSOILDATA)
! Read daily met data: see function GETMET for parameter definitions
!**********************************************************************

    USE maestcom
    USE metcom
    IMPLICIT NONE
    INTEGER METCOLS(MAXMET),SOILDATA,TSOILDATA,IOERROR,I,NOMETCOLS
    INTEGER IHR,IDATE
    REAL ZEN(MAXHRS)
    REAL DATAIN(MAXMET)
    REAL windah(MAXHRS),tsoil(MAXHRS),TAIR(MAXHRS),PPT(MAXHRS)
    REAL SOILMOIST(MAXHRS)
    REAL RADABV(MAXHRS,3),FBEAM(MAXHRS,3),FSUN(MAXHRS)
    REAL RH(MAXHRS),VPD(MAXHRS),VMFD(MAXHRS),CA(MAXHRS),PRESS(MAXHRS)
    REAL DELTAT(12)
    REAL WIND,PRESSK,TMAX,TMIN,DAYL,PRECIP,ALAT,DEC,PAR,FBM
    REAL RADBM,RADDF,SMD,CAK,SWMIN,SWMAX
! Parameters related to the disaggregation of weather and the growth subroutines (added by A. Morales on Feb 2014)
    REAL c, taun, p1, p2, Ik    
! Parameters related to the disaggregation of weather and the growth subroutines (added by A. Morales on Feb 2014)    
    INTEGER OpMie, OpRayleigh, OpEphrat, OpSpitters
! Daily average vapour pressure and function to calculate saturated vapour pressure (added by A. Morales on Feb 2014) 
    REAL EA
    REAL, external :: SATUR

    ! Read in day's data
    READ (UMET,*,IOSTAT = IOERROR) (DATAIN(I), I = 1,NOMETCOLS)
    IF (IOERROR.NE.0) CALL SUBERROR('ERROR READING MET DATA',IFATAL,IOERROR)

! Get input from growth file (contains additional paramters). Added by A. Morales on Feb 2014
    CALL InputDisag(OpMie, OpRayleigh, OpEphrat, OpSpitters, c, taun, p1, p2, Ik)

    ! Set up wind speed array
    IF (METCOLS(MDWIND).EQ.MISSING) THEN
        WIND = DEFWIND  ! Very default value - CHECK!!
    ELSE IF (DATAIN(METCOLS(MDWIND)).LE.0) THEN
        WIND = DEFWIND
    ELSE
        WIND = DATAIN(METCOLS(MDWIND))
    END IF
    DO I = 1,KHRS
      WINDAH(I) = WIND
    END DO

    ! Set up pressure array
    IF (METCOLS(MDPRESS).NE.MISSING) THEN
        DO IHR = 1,KHRS
            PRESS(IHR) = DATAIN(METCOLS(MDPRESS))
        END DO
    ELSE
        DO IHR = 1,KHRS
            PRESS(IHR) = PRESSK
        END DO
    END IF

    ! Calculate hourly temperatures
    TMAX = DATAIN(METCOLS(MDTMAX))
    TMIN = DATAIN(METCOLS(MDTMIN))
! Allow the user to select either Ephrat algorithm or original sinusoidal wave (added by A. Morales on Feb 2014)   
    if(OpEphrat == 1) then
        CALL Ephrat(TMAX, TMIN, taun, p1, p2, Ik, DAYL, TAIR)
    else
        CALL CALCTHRLY(TMAX,TMIN,DAYL,TAIR)
    end if
    ! Calculate soil temperatures
    !      CALL CALCTSOIL(TAIR,TSOIL)
    TSOILDATA = 0  ! Can't read in daily values.

! Added by A. Morales on Feb 2014
! Calculate VPD (Pa) throughout the day using the average daily vapor pressure as input (if available)
! or estimate the vapor pressure from temperature
      if (METCOLS(MDEA).EQ.MISSING) THEN
        ! Calculate relative humidities
        CALL CALCRH(TMIN,TAIR,RH)
        ! Calculate VPDs
        CALL RHTOVPD(RH,TAIR,VPD)
      else
        EA = DATAIN(METCOLS(MDEA))
        VPD = SATUR(TAIR) - EA*1000.0
      end if
! Calculate VMFDs
      CALL VPDTOMFD(VPD,PRESS,VMFD)

    ! Calculate incident PAR
    IF (METCOLS(MDSI).EQ.MISSING) THEN
        IF (METCOLS(MDPAR).EQ.MISSING) THEN
            PRECIP = 0.0
            IF (METCOLS(MDPPT).NE.MISSING) PRECIP = DATAIN(METCOLS(MDPPT))
            CALL BRISTO(IDATE,TMAX,TMIN,PRECIP,DELTAT,ALAT,DEC,DAYL,PAR)
        ELSE
            PAR = DATAIN(METCOLS(MDPAR))
        END IF
    END IF

    ! Calculate diffuse fraction
    IF (METCOLS(MDFBEAM).EQ.MISSING) THEN
        IF (METCOLS(MDSI).EQ.MISSING) THEN
            CALL CALCFBMD(IDATE,ZEN,PAR,FBM)
        ELSE 
            CALL CALCFBMD(IDATE,ZEN,DATAIN(METCOLS(MDSI))*FPAR,FBM)
        END IF
    ELSE
        FBM = DATAIN(METCOLS(MDFBEAM))
    END IF

    ! Calculate PAR
    IF (METCOLS(MDSI).EQ.MISSING) THEN
        RADBM = PAR*FBM
        RADDF = PAR*(1.-FBM)
    ELSE
        RADBM = DATAIN(METCOLS(MDSI))*FBM*FPAR
        RADDF = DATAIN(METCOLS(MDSI))*(1.-FBM)*FPAR
    END IF
! Allow the user to select full Spitters or original Maestra method (added by A. Morales on Feb 2014)  
    If (OpSpitters == 1) then
        CALL Spitters(IDATE,DATAIN(METCOLS(MDSI)),ZEN,RADABV,FBEAM,OpMie,OpRayleigh,DAYL,c,Alat)
    else    
        CALL CALCPARHRLY(RADBM,RADDF,ZEN,RADABV,FBEAM)
    end if

! Calculate NIR. If We are using the full Spitters method this is already calculated. (Added by A. Morales on Feb 2014)
    If (OpSpitters == 0) then
        CALL CALCNIR(RADABV,FBEAM)
    end if

    ! Calculate FSUN
    CALL CALCFSUN(FBEAM,FSUN)

    ! Calculate thermal radiation
    CALL THERMAL(TAIR,VPD,FSUN,RADABV)

    ! Read in value of CA if present - NB could add daily variation if wished?
    IF (METCOLS(MDCA).NE.MISSING) THEN
        DO IHR = 1,KHRS
            CA(IHR) = DATAIN(METCOLS(MDCA))
        END DO
    ELSE
        DO IHR = 1,KHRS
            CA(IHR) = CAK
        END DO
    END IF

    ! Read in rainfall data
    DO IHR = 1,KHRS
        PPT(IHR) = 0.0  ! Initially set all values to zero
    END DO
    IF (METCOLS(MDPPT).NE.MISSING) THEN
        CALL ASSIGNRAIN(DATAIN(METCOLS(MDPPT)),PPT)
    END IF

    ! Three options for soil water content: SMD, SWC, SWP
    ! For SWP, data is soil water potential (MPa)
    ! For SMD, data is soil moisture deficit (0-1)
    ! For SWC, data is soil water content (m3 m-3)
    IF (METCOLS(MDSWP).NE.MISSING) THEN
        SOILDATA = POTENTIAL
        SMD = DATAIN(METCOLS(MDSWP))
        ! Data given is soil water deficit;
        ! Changed from Belinda's version: no conversion done here.
    ELSE IF (METCOLS(MDSMD).NE.MISSING) THEN
        SOILDATA = DEFICIT
        SMD = DATAIN(METCOLS(MDSMD))
        !	  IF (SMD.GT.SWMAX) THEN
        !	    CALL SUBERROR('SW CONTENT EXCEEDS MAX',IWARN,IERROR)
        !	    SMD = 0.0
        !	  ELSE IF (SMD.LT.SWMIN) THEN
        !	    CALL SUBERROR('SW CONTENT IS LESS THAN MIN',IWARN,IERROR)
        !	    SMD = 1.0
        !        ELSE
        !          SMD = (SWMAX - DATAIN(METCOLS(MDSW)))/(SWMAX - SWMIN)
        !	  END IF
    ELSE IF (METCOLS(MDSWC).NE.MISSING) THEN
        SOILDATA = CONTENT
        SMD = DATAIN(METCOLS(MDSWC))
        ! Not 'sw', not 'swp', not 'swc'
    ELSE
        SOILDATA = NONE
        SMD = 0.0
    END IF
	
    DO IHR = 1,KHRS
        SOILMOIST(IHR) = SMD
    END DO

    RETURN
END SUBROUTINE GETMETDAY


!**********************************************************************
SUBROUTINE GETMETHR(IDATE,ZEN,NOMETCOLS,METCOLS,CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,WINDAH, &
                    TSOIL,TAIR,RADABV,FBEAM,RH,VPD,VMFD,CA,PRESS,PPT,SOILMOIST,SOILDATA,            &
                    TSOILDATA,ETMEAS)
! Read hourly met data: see function GETMET for parameter definitions
!**********************************************************************

    USE maestcom
    USE metcom
    IMPLICIT NONE
    INTEGER IHR,IOERROR,I,NOMETCOLS
    INTEGER METCOLS(MAXMET),SOILDATA, TSOILDATA
    INTEGER IDATES(MAXHRS),IDATE
    REAL DATAIN(MAXHRS,MAXMET)
    REAL WINDAH(MAXHRS),TSOIL(MAXHRS),TAIR(MAXHRS)
    REAL RADABV(MAXHRS,3),FBEAM(MAXHRS,3),FSUN(MAXHRS)
    REAL RH(MAXHRS),VPD(MAXHRS),VMFD(MAXHRS),TDEW(MAXHRS)
    REAL ETMEAS(MAXHRS)
    REAL CA(MAXHRS),PRESS(MAXHRS),ZEN(MAXHRS),PPT(MAXHRS)
    REAL SOILMOIST(MAXHRS)
    REAL DELTAT(12),PRESSK,TMIN,TMAX,DAYPPT,ALAT,DEC,DAYL
    REAL PAR,FBM,RADBM,RADDF,CALCFBMH,CAK,SWMIN,SWMAX
    
    ! Read in one day's worth of data at a time.
    DO IHR = 1,KHRS
        READ (UMET,*,IOSTAT = IOERROR) (DATAIN(IHR,I), I = 1,NOMETCOLS)
        IF (IOERROR.NE.0) CALL SUBERROR('ERROR READING MET DATA',IFATAL,IOERROR)
    END DO

    ! Read date array, for checking (RAD).
    IF(METCOLS(MHDATE).NE.MISSING) THEN
        DO IHR = 1,KHRS
            IDATES(IHR) = DATAIN(IHR, METCOLS(MHDATE))
        END DO
        DO IHR = 2,KHRS
            IF(IDATES(IHR).NE.IDATES(IHR-1))THEN
                CALL SUBERROR('ERROR: NUMBER OF MET LINES NOT EQUAL',IFATAL,-1)    
            ENDIF    
        END DO
    ENDIF

    ! Set up pressure array
    IF (METCOLS(MHPRESS).NE.MISSING) THEN
        DO IHR = 1,KHRS
            PRESS(IHR) = DATAIN(IHR,METCOLS(MHPRESS))
        END DO
    ELSE
        DO IHR = 1,KHRS
            PRESS(IHR) = PRESSK
        END DO
    END IF

    ! Set up wind speed array
    DO IHR = 1,KHRS
        IF (METCOLS(MHWIND).EQ.MISSING) THEN
            WINDAH(IHR) = DEFWIND     ! Very default value!!
        ELSE IF (DATAIN(IHR,METCOLS(MHWIND)).LE.0) THEN
            WINDAH(IHR) = DEFWIND
        ELSE
            WINDAH(IHR) = DATAIN(IHR,METCOLS(MHWIND))
        END IF
    END DO

    ! Set hourly air temperatures
    TMIN = DATAIN(1,METCOLS(MHTAIR))
    TMAX = DATAIN(1,METCOLS(MHTAIR))
    DO IHR = 1,KHRS
        TAIR(IHR) = DATAIN(IHR,METCOLS(MHTAIR))
        IF (TAIR(IHR).LT.TMIN) TMIN = TAIR(IHR)
        IF (TAIR(IHR).GT.TMAX) TMAX = TAIR(IHR)
    END DO

    ! Set hourly soil temperatures
    IF (METCOLS(MHTSOIL).EQ.MISSING) THEN
        TSOILDATA = 0
        CALL CALCTSOIL(TAIR,TSOIL)
    ELSE
        DO IHR = 1,KHRS
            TSOIL(IHR) = DATAIN(IHR,METCOLS(MHTSOIL))
            TSOILDATA = 1
        END DO
    END IF

    ! Read in RH, VPD, VMFD, Tdew
    IF (METCOLS(MHRH).NE.MISSING) THEN
        DO IHR = 1,KHRS
            RH(IHR) = DATAIN(IHR,METCOLS(MHRH))
        END DO
    ELSE IF (METCOLS(MHRHP).NE.MISSING) THEN
        DO IHR = 1,KHRS
            RH(IHR) = DATAIN(IHR,METCOLS(MHRHP))/100.0
        END DO
    END IF
    IF (METCOLS(MHVPD).NE.MISSING) THEN
        DO IHR = 1,KHRS
            VPD(IHR) = DATAIN(IHR,METCOLS(MHVPD))
        END DO
    END IF
    IF (METCOLS(MHTDEW).NE.MISSING) THEN
        DO IHR = 1,KHRS
            TDEW(IHR) = DATAIN(IHR,METCOLS(MHTDEW))
        END DO
    END IF
    IF (METCOLS(MHMFD).NE.MISSING) THEN
        DO IHR = 1,KHRS
            VMFD(IHR) = DATAIN(IHR,METCOLS(MHMFD))
        END DO
    END IF

    ! Calculate RH if not in file
    IF ((METCOLS(MHRH).EQ.MISSING).AND.(METCOLS(MHRHP).EQ.MISSING)) THEN
        IF (METCOLS(MHVPD).NE.MISSING) THEN
            CALL VPDTORH(VPD,TAIR,RH)
        ELSEIF (METCOLS(MHTDEW).NE.MISSING) THEN
            CALL TDEWTORH(TDEW,TAIR,RH)
        ELSEIF (METCOLS(MHMFD).NE.MISSING) THEN
            CALL MFDTORH(VMFD,PRESS,TAIR,RH)
        ELSE
            CALL CALCRH(TMIN,TAIR,RH)
        END IF
    END IF
    ! Calculate VPD if not in file
    IF (METCOLS(MHVPD).EQ.MISSING) THEN
        CALL RHTOVPD(RH,TAIR,VPD)
    END IF
    ! Calculate VMFD if not in file
    IF (METCOLS(MHMFD).EQ.MISSING) THEN
        CALL VPDTOMFD(VPD,PRESS,VMFD)
    END IF

    ! Read in rainfall if present.
    ! Don't convert to mmol m-2 as before (RAD)
	DAYPPT = 0.0
    IF (METCOLS(MHPPT).NE.MISSING) THEN 
        DO IHR = 1,KHRS
            PPT(IHR) = DATAIN(IHR,METCOLS(MHPPT))   !*1E6/18.
            DAYPPT = DAYPPT+PPT(IHR)
        END DO
    ELSE
        DO IHR = 1,KHRS
            PPT(IHR) = 0.0
        END DO
    END IF

    ! Read in measured ET, if present.
    IF(METCOLS(MHET).NE.MISSING) THEN
        DO IHR = 1,KHRS
            ETMEAS(IHR) = DATAIN(IHR,METCOLS(MHET))
        END DO
    ENDIF

    ! Must have either PAR (umol m-2 s-1) or global radiation (W m-2)
    IF (METCOLS(MHPAR).NE.MISSING) THEN
        DO IHR = 1,KHRS
            RADABV(IHR,1) = DATAIN(IHR,METCOLS(MHPAR)) / UMOLPERJ
        END DO
    ELSE IF(METCOLS(MHRAD).NE.MISSING) THEN
        DO IHR = 1,KHRS
            RADABV(IHR,1) = DATAIN(IHR,METCOLS(MHRAD)) * FPAR
        END DO
    ELSE
        CALL BRISTO(IDATE,TMAX,TMIN,DAYPPT,DELTAT,ALAT,DEC,DAYL,PAR)
        CALL CALCFBMD(IDATE,ZEN,PAR,FBM)
        RADBM = PAR*FBM
        RADDF = PAR*(1.-FBM)
        CALL CALCPARHRLY(RADBM,RADDF,ZEN,RADABV,FBEAM)
    END IF

    ! Calculate beam fractions
    IF (METCOLS(MHFBEAM).NE.MISSING) THEN
        DO IHR = 1,KHRS
            FBEAM(IHR,1) = DATAIN(IHR,METCOLS(MHFBEAM))
        END DO
    ELSE
        DO IHR = 1,KHRS
            FBEAM(IHR,1) = CALCFBMH(IDATE,ZEN(IHR),RADABV(IHR,1))
            ! 29/3     FBEAM(IHR,1) = CALCFBMH(IDATE,IHR,ALAT,DEC,RADABV(IHR,1))
            !          FBEAM(IHR,1) = CALCFBMWN(IDATE,IHR,ZEN(IHR),RADABV(IHR,1))
        END DO
    END IF

    ! Calculate NIR
    CALL CALCNIR(RADABV,FBEAM)

    ! Calculate FSUN
    CALL CALCFSUN(FBEAM,FSUN)

    ! Calculate thermal radiation
    CALL THERMAL(TAIR,VPD,FSUN,RADABV)

    ! Read in values of CA if present
    IF (METCOLS(MHCA).NE.MISSING) THEN
        DO IHR = 1,KHRS
            CA(IHR) = DATAIN(IHR,METCOLS(MHCA))
        END DO
    ELSE
        DO IHR = 1,KHRS
            CA(IHR) = CAK
        END DO
    END IF

    ! Get soil moisture - either soil water potential or soil moisture deficit

    ! Initialize
    SOILMOIST = 0.0  
    SOILDATA = NONE  
     
    IF (METCOLS(MHSWP).NE.MISSING) THEN
        SOILDATA = POTENTIAL
        DO IHR = 1,KHRS
            SOILMOIST(IHR) = DATAIN(IHR,METCOLS(MHSWP))
        END DO
	ELSE IF (METCOLS(MHSMD).NE.MISSING) THEN
        SOILDATA = DEFICIT
        DO IHR = 1,KHRS
            SOILMOIST(IHR) = DATAIN(IHR,METCOLS(MHSMD))
            ! No conversion as was done before: assume that data is already deficit.
            !	    SOILMOIST(IHR) = (SWMAX - DATAIN(IHR,METCOLS(MHSW)))
            !     &    /(SWMAX - SWMIN)
        END DO
    ELSE IF (METCOLS(MHSWC).NE.MISSING) THEN
        SOILDATA = CONTENT
        DO IHR = 1,KHRS
            SOILMOIST(IHR) = DATAIN(IHR,METCOLS(MHSWC))
        END DO
    END IF

    RETURN
END SUBROUTINE GETMETHR


!**********************************************************************
SUBROUTINE CALCTHRLY(TMAX,TMIN,DAYL,TAIR)
! Calculate a daily variation in temperature from max & min temperatures.
! Temp varies linearly between sunset & sunrise, and sinusoidally during the day.
! INPUTS:
! TMIN, TMAX - minimum and maximum daily temperatures, °C
! DAYL - daylength, hours
! OUTPUTS:
! TAIR - array of hourly air temperatures, °C
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL TAIR(MAXHRS),TMAX,TMIN,DAYL
    REAL TAV,TAMPL,HRTIME,TIME
    
    TAV = (TMAX+TMIN)/2.0
    TAMPL = (TMAX-TMIN)/2.0

    DO I=1,KHRS
        HRTIME = REAL(I) - 0.5
        TIME = HRTIME + DAYL*0.5 - REAL(KHRS)/2.0
        IF (TIME.LT.0.0.OR.TIME.GT.DAYL) THEN
            IF (TIME.LT.0.0) HRTIME = HRTIME + KHRS
            TAIR(I) = TAV - (TAV-TMIN)*(HRTIME-DAYL*0.5-(REAL(KHRS)/2.0))/(REAL(KHRS)-DAYL)
        ELSE
            TAIR(I) = TAV - TAMPL*COS(1.5*PI*TIME/DAYL)
        ENDIF
    END DO    
    RETURN
END SUBROUTINE CALCTHRLY


!**********************************************************************
SUBROUTINE CALCTSOIL(TAIR,TSOIL)
! Calculate soil temperatures.
! Set equal to average daily air temperature.
! INPUTS:
! TAIR - array of hourly air temperatures, °C
! OUTPUTS:
! TSOIL - array of hourly soil temperatures, °C
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL TSOIL(MAXHRS),TAIR(MAXHRS)
    REAL TAV

    TAV = 0.0
    DO I = 1,KHRS
        TAV = TAV + TAIR(I)
    END DO
    TAV = TAV / REAL(KHRS)

    DO I = 1,KHRS
        TSOIL(I) = TAV
    END DO
    RETURN  
END SUBROUTINE CALCTSOIL


!**********************************************************************
SUBROUTINE ASSIGNRAIN(TOTAL,PPT)
! Given daily total of PPT, assign to hours through the day.
! Use algorithm from GRAECO (model of D. Loustau).
!**********************************************************************

!	USE MSFLIB        ! Library required by COMPAQ VISUAL FORTRAN - superseded
!  USE IFPORT        ! Library required by Intel Visual Fortran 
    USE maestcom
    
    IMPLICIT NONE
    INTEGER IHR,IRAIN,IHRSWITHRAIN,I
    INTEGER(4)  IFLAG
    REAL(4) RANVAL
    REAL PPT(MAXHRS)    ! Hourly rainfall: assume already set to zero
    REAL TOTAL,RAIN,RATE

    ! Reset rain array (RAD).
    PPT = 0.

    !      CALL SEED(1995) ! For testing - same random numbers - change to 0 later. 
    IFLAG = 0
      
    ! 1. All rain falls in one hour for light storms (<2 mm)
    IF (TOTAL.LE.2.) THEN
        RANVAL = RAND(IFLAG)   ! randomly select hour
        IRAIN = INT(RANVAL*KHRS)+1
        PPT(IRAIN) = TOTAL
    ! 2. All rain falls in 24 hours for storms >46 mm 
    ELSE IF (TOTAL.GT.46.) THEN
        RAIN = TOTAL/REAL(KHRS)
        DO IHR = 1,KHRS
            PPT(IHR) = RAIN
        END DO
    ! 3. All rain falls at 2mm/hour at a random time of the day */
    ELSE
        IHRSWITHRAIN = MIN(INT((TOTAL/2)*KHRS/24),KHRS)
        RATE = TOTAL/REAL(IHRSWITHRAIN)
        DO I = 1, IHRSWITHRAIN
            RANVAL = RAND(IFLAG)   ! randomly select hours
            IRAIN = INT(RANVAL*KHRS)+1
            PPT(IRAIN) = PPT(IRAIN)+RATE
        END DO
    END IF

    ! Don't do that (RAD).
    !! Convert from mm to mmol m-2
    !      DO 30 IHR = 1, KHRS
    !	  PPT(IHR) = PPT(IHR)*1E6/18.
    !30    CONTINUE

    RETURN
END SUBROUTINE ASSIGNRAIN


!**********************************************************************
 SUBROUTINE CALCRH(TMIN,TAIR,RH)
! Calculate hourly relative humidity from air temperature and minimum
! daily temperature (assumed to be the dewpoint).
! INPUTS:
! TMIN - daily minimum temperature, °C
! TAIR - array of hourly air temperatures, °C
! OUTPUTS:
! RH - array of hourly relative humidity, fraction
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL TAIR(MAXHRS),RH(MAXHRS),TMIN,HUMABS
    REAL, EXTERNAL :: SATUR

    HUMABS = SATUR(TMIN)
    DO I = 1,KHRS
        RH(I) = HUMABS/SATUR(TAIR(I))
    END DO

    RETURN
END SUBROUTINE CALCRH


!**********************************************************************
SUBROUTINE RHTOVPD(RH,TAIR,VPD)
! Convert from RH to VPD.
! Inputs: RH(KHRS) - hourly relative humidity, fraction
!   TAIR(KHRS) - hourly air temperature, degrees !
! Outputs: VPD(KHRS) - hourly VPD, Pa
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL TAIR(MAXHRS),RH(MAXHRS),VPD(MAXHRS)
    REAL, EXTERNAL :: SATUR

    DO I = 1,KHRS
        VPD(I)=(1.0-RH(I))*SATUR(TAIR(I))
    END DO

    RETURN
END SUBROUTINE RHTOVPD

!**********************************************************************
SUBROUTINE VPDTORH(VPD,TAIR,RH)
! Convert from VPD to RH.
! Inputs:
!   VPD(KHRS) - hourly VPD, Pa
!   TAIR(KHRS) - hourly air temperature, degrees !
! Outputs:
!   RH(KHRS) - hourly relative humidity, fraction
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL TAIR(MAXHRS),RH(MAXHRS),VPD(MAXHRS)
    REAL, EXTERNAL :: SATUR

    DO I = 1,KHRS
        RH(I) = 1.0 - VPD(I)/SATUR(TAIR(I))
    END DO

    RETURN
END SUBROUTINE VPDTORH

!**********************************************************************
SUBROUTINE TDEWTORH(TDEW,TAIR,RH)
! Convert from TDEW to RH.
! Inputs:
!   TDEW(KHRS) - hourly dewpoint temperature, degrees !
!   TAIR(KHRS) - hourly air temperature, degrees !
! Outputs:
!   RH(KHRS) - hourly relative humidity, fraction
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL, EXTERNAL :: SATUR

    REAL TAIR(MAXHRS),RH(MAXHRS),TDEW(MAXHRS)

    DO I = 1,KHRS
        RH(I) = SATUR(TDEW(I))/SATUR(TAIR(I))
    END DO

    RETURN
END SUBROUTINE TDEWTORH


!**********************************************************************
SUBROUTINE MFDTORH(VMFD,PRESS,TAIR,RH)
! Convert from mole fraction deficit to relative humidity.
! Inputs:
!   VMFD(KHRS) - hourly mole fraction deficit, mmol mol-1
!   TAIR(KHRS) - hourly air temperature, degrees !
!   PRESS(KHRS) - hourly pressure, Pa
! Outputs:
!   RH(KHRS) - hourly relative humidity, fraction
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL VMFD(MAXHRS),PRESS(MAXHRS),RH(MAXHRS),TAIR(MAXHRS)
    REAL, EXTERNAL :: SATUR

    DO I = 1,KHRS
        RH(I) = 1 - VMFD(I)*PRESS(I)*1E-3/SATUR(TAIR(I))
    END DO

    RETURN
END SUBROUTINE MFDTORH


!**********************************************************************
SUBROUTINE VPDTOMFD(VPD,PRESS,VMFD)
! Convert from VPD to VMFD.
! Inputs:
!   VPD(KHRS) - hourly VPD, kPa
!   PRESS(KHRS) - hourly pressure, Pa
! Outputs:
!   VMFD(KHRS) - hourly mole fraction deficit, mmol mol-1
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I

    REAL VPD(MAXHRS),PRESS(MAXHRS),VMFD(MAXHRS)

    DO I = 1,KHRS
        VMFD(I)= VPD(I)/PRESS(I)*1E3
    END DO

    RETURN
END SUBROUTINE VPDTOMFD


!*********************************************************************
SUBROUTINE BRISTO(IDAY,TMAX,TMIN,PRECIP,DELTAT,ALAT,DEC,DAYL,PAR)
! Subroutine implements Bristow & Campbell (1984) Ag For Met 31:159-166
! Calculates incident daily radiation from temperature data
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDAY,IMON,MONTH
    REAL DELTAT(12)
    REAL TMAX,TMIN,PRECIP,ALAT,DEC,DAYL,PAR,ANGDL
    REAL DELT,BRISTOK,TRANSM,RADCLR,RADTOT

    DELT = TMAX - TMIN					!Daily temperature amplitude - should include prev. day but ..
    IF (PRECIP.GT.10.) DELT = DELT*0.75	!Reduce if rainy
    IMON = MONTH(IDAY)					!Month
    BRISTOK = 0.036*EXP(-0.154*DELTAT(IMON))	
    TRANSM = 1 - EXP(-BRISTOK*(DELT)**2.4)	!Transmittance

    ANGDL=DAYL*PI/KHRS					! Angular daylength
    !Incident radiation on clear day
    RADCLR=0.0864*SOLARC/PI*TAU* (ANGDL*SIN(ALAT)*SIN(DEC)+ COS(ALAT)*COS(DEC)*SIN(ANGDL))	
    RADTOT = RADCLR*TRANSM				!Estimated solar radiation in MJ m-2 d-1??
    PAR = RADTOT*FPAR					!Estimated PAR in MJ m-2 d-1

    RETURN
END SUBROUTINE BRISTO


!**********************************************************************
SUBROUTINE CALCPARHRLY(RADBM,RADDF,ZEN,RADABV,FBEAM)
! Calculate daily course of incident PAR from daily totals of
! beam and diffuse PAR (in MJ m-2 d-1).
! INPUTS:
! RADBM - daily total beam PAR (MJ m-2 d-1)
! RADDF - daily total diffuse PAR (MJ m-2 d-1)
! ZEN - zenith angle (radians)
! OUTPUTS:
! RADABV - array of hourly incident PAR (J m-2 s-1)
! FBEAM - array of hourly beam fractions
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL ZEN(MAXHRS),RADABV(MAXHRS,3),FBEAM(MAXHRS,3)
    REAL COSBM(MAXHRS),COSDF(MAXHRS)
    REAL SUMBM,SUMDF,COSZEN,HRTIME,RDBM,RADBM,RDDF,RADDF

    SUMBM = 0.0
    SUMDF = 0.0

    DO I=1,KHRS
        COSBM(I) = 0.0
        COSDF(I) = 0.0
        HRTIME = REAL(I) - 0.5
        COSZEN = COS(ZEN(I))
        IF (COSZEN.GT.0.0) THEN
            IF (ZEN(I).LT.80*PID180) THEN  !Set FBM = 0.0 for ZEN > 80 degrees
                COSBM(I)=COSZEN*TAU**(1.0/COSZEN)
            ELSE
                COSBM(I)=0.0
            END IF
            COSDF(I)=COSZEN
            SUMBM=SUMBM+COSBM(I)
            SUMDF=SUMDF+COSDF(I)
        ENDIF
    END DO

    DO I=1,KHRS
        IF (SUMBM.GT.0.0) THEN
            RDBM = RADBM*COSBM(I)/SUMBM
        ELSE
            RDBM = 0.0
        END IF
        IF (SUMDF.GT.0.0) THEN
            RDDF = RADDF*COSDF(I)/SUMDF
        ELSE
            RDDF = 0.0
        END IF
    
        ! Convert from MJ hr-1 to J s-1
        RADABV(I,1)=(RDDF+RDBM)*1e6/SPERHR
        IF ((RDBM+RDDF).GT.0.0) THEN
            FBEAM(I,1) = RDBM/(RDBM+RDDF)
        ELSE
            FBEAM(I,1)= 0.00
        END IF
    END DO

      RETURN
END SUBROUTINE CALCPARHRLY


!**********************************************************************
SUBROUTINE CALCFBMD(IDATE,ZEN,PAR,FBM)
! Calculate the beam fraction from the total daily incident radiation.
! Use the formula of Spitters et al. (1986) Agric For Met 38:217-229.
! INPUTS:
! IDATE - date in days-since-1950 format
! ZEN - array of hourly sun zenith angle (radians)
! PAR - daily total incident PAR (MJ m-2 d-1)
! OUTPUTS:
! FBM - daily beam fraction
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDATE,IHR
	REAL ZEN(MAXHRS),PAR,FBM
	REAL S0,SINB,TRANS,ETRAD,FDIF

    ! Calculate extra-terrestrial radiation
    S0 = 0.0
    DO IHR = 1,KHRS
        SINB = SIN(PID2-ZEN(IHR))
        S0 = S0 + ETRAD(IDATE,SINB)*SPERHR/1E6
    END DO

    ! Spitter's formula
    TRANS = (PAR/FPAR) / S0
    IF (TRANS.LT.0.07) THEN
        FDIF = 1.
    ELSE IF (TRANS.LT.0.35) THEN
        FDIF = 1. - 2.3*(TRANS-0.07)**2
    ELSE IF (TRANS.LT.0.75) THEN
        FDIF = 1.33 - 1.46*TRANS
    ELSE
        FDIF = 0.23   
    END IF   
    FBM = 1. - FDIF

    RETURN
END SUBROUTINE CALCFBMD

!**********************************************************************
FUNCTION CALCFBMH(IDATE,ZEN,RADABV)
! Calculate the beam fraction from the hourly incident radiation.
! Use the formula of Spitters et al. (1986) Agric For Met 38:217-229.
! INPUTS:
! IDATE - date in days-since-1950 format
! ZEN - hourly sun zenith angle (radians)
! RADABV - total incident PAR (J m-2 s-1)
! RETURNS:
! CALCFBMH - hourly beam fraction
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDATE
    REAL SINB,ZEN,RADABV
    REAL SPITR,SPITK,TRANS,CALCFBMH,S0,ETRAD,FDIF

    ! SINB is the sine of the elevation of the sun above the horizon
    ! (= SIN(90-ZEN)). For zenith angles > 80 degrees will put FBEAM = 0.
    SINB = SIN(PID2-ZEN)

    ! Calculate extra-terrestrial radiation
    S0 = ETRAD(IDATE,SINB)

    IF (SINB.GT.0.17) THEN
        ! Spitter's formula
        TRANS = (RADABV/FPAR) / S0
        SPITR = 0.847 - 1.61*SINB + 1.04*SINB**2
        SPITK = (1.47 - SPITR)/1.66
        IF (TRANS.LE.0.22) THEN
            FDIF = 1.
        ELSE IF (TRANS.LE.0.35) THEN
            FDIF = 1. - 6.4*(TRANS-0.22)**2
        ELSE IF (TRANS.LE.SPITK) THEN
            FDIF = 1.47 - 1.66*TRANS
        ELSE
            FDIF = SPITR
        END IF
        CALCFBMH = 1. - FDIF
    ELSE
        CALCFBMH = 0.0
    END IF
    RETURN
END FUNCTION CALCFBMH


!**********************************************************************
SUBROUTINE CALCNIR(RADABV,FBEAM)
! From incident PAR, calculate incident NIR.
! Parameters:
! RADABV - hourly incident radiation in 3 wavelengths (J m-2 s-1)
! FBEAM - beam fraction of incident radiation for 3 wavelengths
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL RADABV(MAXHRS,3),FBEAM(MAXHRS,3)

    ! Fbeam is assumed the same for NIR as PAR. Amount of NIR is calculated from
    ! amount of PAR & the fraction of total solar that is PAR (FPAR).
    DO I = 1,KHRS
        RADABV(I,2) = RADABV(I,1)*(1./FPAR - 1.0)
        FBEAM(I,2) = FBEAM(I,1)
    END DO

    RETURN
END SUBROUTINE CALCNIR


!**********************************************************************
SUBROUTINE THERMAL(TAIR,VPD,FSUN,RADABV)
! Calculate incident thermal radiation, if it has not been measured. 
! Several different formulae exist:
! Ying-Ping had the following formula, from MH Unsworth & JL Monteith
! (1975) Quart. J. R. Met. Soc. 101. pp13-24, pp25-34
!        RADABV(I,3) = 1.06*SIGMA* (TK(TAIR(I))**4) - 119.
! I replaced this with a formula from Brutsaert (1975) Water Resources
! Research 11: 742-744. Leuning et al (1995) PC&E 18: 1183-1200 say that
! Hatfield et al (1983) Water Resources Research 19: 285-288 tested
! a range of alternative formulae & found Brutsaert's to be the best.
! BM 22/5/98.
! This formula is appropriate for clear skies only. The thermal radiation
! from a cloudy sky will be larger since clouds radiate more effectively. 
! Changed to the formula taking this into account, given by Monteith & Unsworth
! 1990 Principles of Environmental Physics 2nd ed p53. 
!**********************************************************************

    USE maestcom
    USE switches
    IMPLICIT NONE
    INTEGER I
    REAL RADABV(MAXHRS,3),TAIR(MAXHRS),VPD(MAXHRS),FSUN(MAXHRS)
    REAL TMP,EA,EMCLEAR,EMSKY
    REAL, EXTERNAL :: TK
    REAL, EXTERNAL :: SATUR
    
    DO I = 1,KHRS
        ! Brutsaert et al formula
        !        EA = SATUR(TAIR(I)) - VPD(I) !Old formula - see comments
        !        EMSKY = 0.642*(EA/TK(TAIR(I)))**(1./7.)
        !        RADABV(I,3) = EMSKY*SIGMA*(TK(TAIR(I))**4)

        ! Monteith and Unsworth formula
        !	  SIGT4 = SIGMA*(TK(TAIR(I))**4)
        !	  EMCLEAR = 1.06*SIGT4 - 119.
        !	  RADABV(I,3) = FSUN(I)*EMCLEAR + (1.-FSUN(I))*(0.84+0.16*EMCLEAR)

        ! BM 12/05: Am hybridising. Use Brutsaert et al for clear sky emission but
        ! Monteith and Unsworth correction for cloudy sky. 
        ! The 1.06 SIGT4 - 119 formula seems problematic in some conditions. 
        ! Also suspect I had correction incorrect - check this!!
        tmp = satur(tair(I))
        EA = SATUR(TAIR(I)) - VPD(I) !Old formula - see comments
!        IF (IOFORMAT .EQ. 0) THEN
!            WRITE(UWATTEST,*)EA
!        ELSE IF(IOFORMAT .EQ. 1) THEN
!            WRITE(UWATTEST) EA 
!        ENDIF
        EMCLEAR = 0.642*(EA/TK(TAIR(I)))**(1./7.)
        EMSKY = FSUN(I)*EMCLEAR + (1.-FSUN(I))*(0.84+0.16*EMCLEAR)
        RADABV(I,3) = EMSKY*SIGMA*(TK(TAIR(I))**4)
    END DO

    RETURN
END SUBROUTINE THERMAL


!**********************************************************************
SUBROUTINE CALCFSUN(FBEAM,FSUN)
! Calculates the sunlit time from the beam fraction.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL FBEAM(MAXHRS,3),FSUN(MAXHRS)

    DO I = 1,KHRS
        FSUN(I) = 1.16*FBEAM(I,1)-0.01
        IF (FSUN(I).LT.0.00) FSUN(I)=0.00
        IF (FSUN(I).GT.1.00) FSUN(I)=1.00
    END DO

    RETURN
END SUBROUTINE CALCFSUN


!**********************************************************************
REAL FUNCTION ETRAD(IDATE,SINB)
! Calculate the radiation incident on the atmosphere.
! Using formulae from Spitters et al (1986) Agric For Met 38:217
! Returns value in J m-2 s-1.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDATE,JDAY
    INTEGER, EXTERNAL :: JDATE
    REAL SINB
    
    JDAY = JDATE(IDATE)

    ! Spitters' formula
    IF (SINB.GT.0.0) THEN
        ETRAD = SOLARC * (1 + 0.033*COS(REAL(JDAY)/365.0*TWOPI)) * SINB
    ELSE
        ETRAD = 0.0
    END IF

    RETURN
END FUNCTION ETRAD


!**********************************************************************
SUBROUTINE SKIPMET(MFLAG,NSTEP)
! Read data from the met file until in the right position
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    INTEGER MFLAG,NSTEP,LINESTOSKIP
    CHARACTER MTITLE

    IF (MFLAG.EQ.0) THEN
        LINESTOSKIP = NSTEP - 1
    ELSE
        LINESTOSKIP = KHRS*(NSTEP - 1)
    END IF

    DO I = 1,LINESTOSKIP
        READ (UMET, *) MTITLE
    END DO

    RETURN
END SUBROUTINE SKIPMET


!**********************************************************************
SUBROUTINE READLAT(UFILE, ALAT, TTIMD)
! Read in Lat & Long details - from met file:
! ALAT - latitude, in radians (-ve for S hemisphere)
! TTIMD - time difference between longitude of plot & longitude of time
!   zone, in hours
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
      
    INTEGER UFILE,IOERROR
    REAL LAT(3),LONG(3)
    CHARACTER*1 LATHEM,LONHEM
    REAL TZLONG,ALAT,ALONG,TIMDIF,TTIMD
    NAMELIST /LATLONG/ LATHEM,LAT,LONHEM,LONG,TZLONG

    REWIND (UFILE)
    READ (UFILE, LATLONG, IOSTAT = IOERROR)
    IF (IOERROR.NE.0) CALL SUBERROR('ERROR READING LATITUDE: IT MUST BE IN MET.DAT',IFATAL,IOERROR)
    ALAT = (LAT(1)+LAT(2)/60.0+LAT(3)/3600.0)*PID180
    IF (LATHEM.EQ.'S') ALAT = -ALAT
    ALONG = LONG(1) + LONG(2)/60.0 + LONG(3)/3600.0
    IF (LONHEM.EQ.'E') THEN
        ALONG = 360.0 - ALONG
        TZLONG = 360.0 - TZLONG
    END IF
    ALONG = ALONG*PID180
    TZLONG = TZLONG*PID180
    TIMDIF = KHRS/2.0*ALONG/PI
    TTIMD = KHRS/2.0/PI*(ALONG - TZLONG)

    RETURN
END SUBROUTINE READLAT


!**********************************************************************
SUBROUTINE GETWIND(FOLLAY,FOLNTR,TOTLAI,EXTWIND,WINDLAY)
! Calculate decrease in wind speed with increasing canopy depth
! Uses exponential decline. Does not take into account different sizes of trees.
! INPUTS:
! FOLLAY - amount of foliage in each layer of target tree
! FOLNTR - total foliage of target tree
! TOTLAI - total LAI of forest
! EXTWIND - extinction coefficient of wind
! OUTPUTS:
! WINDLAY - windspeed in each layer of the canopy (m s-1)
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER ILAY
    REAL FOLLAY(MAXLAY),FOLABV(MAXLAY),WINDLAY(MAXLAY)
    REAL TOTLAI,FOLNTR,EXTWIND
    
    ! Calculate approximately LAI above each point (for wind speed)
    FOLABV(1) = FOLLAY(1)/2.
    DO ILAY = 2,MAXLAY
        FOLABV(ILAY) = FOLABV(ILAY-1) + (FOLLAY(ILAY-1) + FOLLAY(ILAY))/2.
    END DO

    ! Calculate decrease in wind speed with canopy depth
    DO ILAY = 1,MAXLAY
        IF (FOLNTR.GT.0.0) THEN
            FOLABV(ILAY) = FOLABV(ILAY)/FOLNTR * TOTLAI
        ELSE
            FOLABV(ILAY) = 0.0
        END IF
        WINDLAY(ILAY) = EXP(-EXTWIND*FOLABV(ILAY))
    END DO

    RETURN
END SUBROUTINE GETWIND


!**********************************************************************
SUBROUTINE ALTERMETCC(CA,TAIR,TSOIL,RH,VPD,VMFD,PRESS,CO2INC,TINC)
! Subroutine to change met data according to climate change scenario.
! Currently can increase CO2 by fixed amount (CO2INC, ppm)
! and/or T by a fixed amount (TINC, deg !).
! The absolute humidity (in g m-3) is maintained constant.
! VPD, RH and VMFD are re-calculated if temperature is changed.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IHR
    
    REAL CA(MAXHRS),TAIR(MAXHRS),TSOIL(MAXHRS)
    REAL VPD(MAXHRS),RH(MAXHRS),VMFD(MAXHRS),PRESS(MAXHRS)
    REAL CO2INC,ABSHUM,TINC
    REAL, EXTERNAL :: AHTORH,RHTOVPD,VPDTOMFD,CALCAH

    DO IHR = 1,KHRS
        CA(IHR) = CA(IHR) + CO2INC
        ABSHUM = CALCAH(RH(IHR),TAIR(IHR))
        TAIR(IHR) = TAIR(IHR) + TINC
        TSOIL(IHR) = TSOIL(IHR) + TINC
        RH(IHR) = AHTORH(ABSHUM,TAIR(IHR))
    END DO
    IF (TINC.NE.0.0) THEN
!        CALL RHTOVPD(RH,TAIR,VPD)
!        CALL VPDTOMFD(VPD,PRESS,VMFD)
      CALL SUBERROR('TINC NOT SUPPORTED : FAULTY CODE DISABLED, CHECK ALTERMETCC IN GETMET.F90',IFATAL,-1)

    END IF

    RETURN
END SUBROUTINE ALTERMETCC


!**********************************************************************
REAL FUNCTION CALCAH(RH,TAIR)
! Calculate absolute humidity (g m-3) from relative humidity & T.
! Conversion from Jones (1992) pp. 109-110.
!**********************************************************************
    
    IMPLICIT NONE
    REAL, EXTERNAL :: SATUR,TK
    REAL TAIR,RH
    
    CALCAH = RH*SATUR(TAIR)*2.17/TK(TAIR)

    RETURN
END FUNCTION CALCAH


!**********************************************************************
REAL FUNCTION AHTORH(ABSHUM,TAIR)
! Calculate relative humidity from absolute humidity & T.
! Conversion from Jones (1992) pp. 109-110.
!**********************************************************************
    IMPLICIT NONE
    REAL TAIR,ABSHUM
    REAL, EXTERNAL :: TK,SATUR
    
    AHTORH = ABSHUM*TK(TAIR)/2.17/SATUR(TAIR)

    RETURN
END FUNCTION AHTORH


!**********************************************************************
SUBROUTINE ALTERMETOTC(TOTC,WINDOTC,PAROTC,FBEAMOTC,TAIR,TSOIL,WINDAH,RADABV,&
                        FBEAM,RH,VPD,VMFD,PRESS)
! Subroutine to change met data according to effects of OTC.
! Currently can increase T by fixed amount (deg !),
! set wind speed to constant amount (m s-1),
! reduce radiation by a constant fraction,
! and/or alter the beam fraction of radiation.
! The absolute humidity (in g m-3) is maintained constant.
! VPD, RH and VMFD are re-calculated if temperature is changed.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IHR
    
    REAL windah(MAXHRS),tsoil(MAXHRS),tair(MAXHRS)
    REAL radabv(MAXHRS,3),fbeam(MAXHRS,3)
    REAL VPD(MAXHRS),RH(MAXHRS),VMFD(MAXHRS),PRESS(MAXHRS)
    REAL ABSHUM,TOTC,WINDOTC,PAROTC,FBEAMOTC
    REAL ZEN
    REAL, EXTERNAL :: CALCAH,AHTORH,CALCFBMWN

    DO IHR = 1,KHRS
        ABSHUM = CALCAH(RH(IHR),TAIR(IHR))
        TAIR(IHR) = TAIR(IHR) + TOTC
        TSOIL(IHR) = TSOIL(IHR) + TOTC
        RH(IHR) = AHTORH(ABSHUM,TAIR(IHR))
        IF (WINDOTC.GT.0.0) WINDAH(IHR) = WINDOTC
        RADABV(IHR,1) = RADABV(IHR,1)*PAROTC
        FBEAM(IHR,1) = FBEAM(IHR,1)*FBEAMOTC
    END DO
    IF (TOTC.NE.0.0) THEN
        CALL RHTOVPD(RH,TAIR,VPD)
        CALL VPDTOMFD(VPD,PRESS,VMFD)
    END IF

    RETURN
END SUBROUTINE ALTERMETOTC


!**********************************************************************
REAL FUNCTION CALCFBMWN(IDATE,IHR,ZEN,RADABV)
! Calculate the beam fraction of PAR using the formula of Weiss &
! Norman (1985) Agric For Met 34:205-214.
! An alternative expression to that of Spitters et al (1986)
! (see CALCFBMD and CALCFBMH).
! Assume that pressure = sea level pressure (101.325kPa)
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER NIGHT,IDATE,IHR
    REAL FBEAM(MAXHRS,3),ZEN,RADABV,RBV,RDV,AXLOG
    REAL PARMA,PARMB,PARMC,PARMD,COSZEN,AIRMASS,ADUM,WATABS
    REAL RBN,RDN,RTOTV,RTOTN,RMEAS,RATIO,FSUN

    PARMA = 0.90
    PARMB = 0.70
    PARMC = 0.88
    PARMD = 0.68

    !      DO 10 IWAVE = 1,3
    !       FBEAM(IHR,IWAVE) = 0.00
    !10    CONTINUE
    CALCFBMWN = 0.0

    IF (NIGHT(ZEN,RADABV).EQ.1) THEN
        RETURN
    END IF

    COSZEN = COS(ZEN)
    AIRMASS = 1./COSZEN

    ! Potential direct-beam PAR on a horizontal surface can be approx. by
    RBV = 600.*EXP(-0.185*AIRMASS)*COSZEN
    ! Potential visible diffuse radiation is approx. by
    ! BK 2/95 Changed after advice by F. Villalobos, Cordoba
    RDV = 0.4*(600.*COSZEN - RBV)

    AXLOG = ALOG10(AIRMASS)
    ADUM = (-1.195+.4495*AXLOG-.0345*AXLOG*AXLOG)
    WATABS = 1320.*(10.**(ADUM))
    ! Potential direct-beam NIR
    RBN = (720.*EXP(-0.06*AIRMASS)-WATABS)*COSZEN
    ! Potential diffuse NIR
    ! BK 2/95 Changed after advice by F. Villalobos, Cordoba
    RDN = 0.6* (720.*COSZEN-RBN-WATABS*COSZEN)

    ! Total potential visible radiation (estimated)
    RTOTV = RBV + RDV
    ! Total potential NIR radiation (estimated)
    RTOTN = RBN + RDN
    ! Total incoming radiation (measured)
    RMEAS = RADABV/FPAR
    ! Ratio of measured to potential solar radiation
    RATIO = RMEAS/(RTOTV + RTOTN)

    IF (RATIO.GT.PARMA) THEN
        ! FBEAM(IHR,1) = RBV/RTOTV
        CALCFBMWN = RBV/RTOTV
    ELSE
        ! FBEAM(IHR,1) = RBV/RTOTV* (1.- ((PARMA-RATIO)/PARMB)** (2./3.))
        CALCFBMWN = RBV/RTOTV* (1.- ((PARMA-RATIO)/PARMB)** (2./3.))
    END IF

    !      IF (RATIO.GT.PARMC) THEN
    !         FBEAM(IHR,2) = RBN/RTOTN
    !      ELSE
    !         FBEAM(IHR,2) = RBN/RTOTNN* (1.- ((PARMC-RATIO)/PARMD)** (2./3.))
    !      END IF

    IF (CALCFBMWN.LT.0.000) CALCFBMWN = 0.00
    !      IF (FBEAM(IHR,1).LT.0.000) FBEAM(IHR,1) = 0.00
    !      IF (FBEAM(IHR,2).LT.0.000) FBEAM(IHR,2) = 0.00

    !    *****************************************
    ! The following equation is from the Ph.D thesis of M.H. Steven,
    ! Univ. of Nottingham, 1977, entitled  'Angular distribution and
    ! interception of diffuse solar radiation'.
    !   diffuse/global=-0.86*sunshine duration+0.99
    FSUN = 1.16*CALCFBMWN - 0.01
    !         FSUN = 1.16*FBEAM(IHR,1) - 0.01
    ! Add the circumsolar diffuse components to the direct radiation
    ! for the sake of geometric simplicity.
    ! According to M.H. Steven this component is 51.0 per cent of the total
    ! diffuse.

    IF (FSUN.GT.0.8) THEN
        CALCFBMWN = CALCFBMWN + 0.51* (1.0-CALCFBMWN)*FSUN
    !         FBEAM(IHR,1) = FBEAM(IHR,1) + 0.51* (1.0-FBEAM(IHR,1))*FSUN
    !         FBEAM(IHR,2) = FBEAM(IHR,2) + 0.51* (1.0-FBEAM(IHR,2))*FSUN
    END IF
    
END FUNCTION CALCFBMWN

!***************************************************************************

!**********************************************************************
    Subroutine Ephrat(TMAX, TMIN, taun, p1, p2, Ik, DAYL, TAIR)
! Calculate a daily variation in temperature from max & min temperatures.
! using the equations for daytime and nighttime air temperature described in Ephrat et al. (1996).
! INPUTS:
! TMIN, TMAX - minimum and maximum daily temperatures
! DAYL - daylength, hours
! taun, p1, p2 and Ik - parameters of the model by Ephrat et al. (1996) with delay of minimum temperature. Ik is inverse of Tk
! OUTPUTS:
! TAIR - array of hourly air temperatures

! Added by A. Morales on Feb 2014
!**********************************************************************
Use maestcom
Implicit None
Integer :: I
Real, intent(in) :: TMAX,TMIN,DAYL,taun, p1, p2, Ik
Real, intent(out) :: TAIR(MAXHRS)
Real :: HRTIME,TIME, S, sunset, sunrise, Tss
! Time of sunrise and sunset (h)
sunrise = 12.0 - DAYL/2.0
sunset = 12.0 + DAYL/2.0

! Loop through daytime and calculate air temperature (Tair, °C) at each time of the day
Day: do I = 1,KHRS
! Calculates the time of the day (h).
   time = real(I)/KHRS*24
! Calculate the function S (relative variation of temperature during daytime)
   S = sin(pi*(time  - p1 - 12.0 + DAYL/2.0)/(DAYL + 2.0*p2))
! And now calculate temperature during the day
   if(time > sunrise +  p1 .AND. time < sunset) then
! When we do not consider the effect of buoyancy (I¬k, °C–1). In theory Ik = 0, but to avoid floating-point errors is better to give a negative value. 
     if(Ik < 0) then
! Use the simple version from Goudriaan & van Laar (1994) 
       Tair(I) = Tmin + (Tmax - Tmin)*S
! When we consider the effect of buoyancy (i.e. Ik > 0)
     else 
       Tair(I) = Tmin - 0.5/Ik + 0.5*sqrt(1/Ik**2 + 4*(Tmax - Tmin)*(1/Ik + Tmax - Tmin)*S)
     endif
   endif  
! Get the temperature the closest to sunset (Tss, °C). This variable will be replaced with the temperature until we go over the sunset   
     if(time - 24/KHRS < sunset) then
       Tss = Tair(I - 1)
     endif     
enddo Day

! Loop through nightime and calculate air temperature (Tair, °C) at each time of the night
Night: do I = 1,KHRS
! Calculates the time of the day (h).
   time = real(I)/KHRS*24.0
! Calculate the function S (relative variation of temperature during daytime)
   S = sin(pi*(time  - p1 - 12 + DAYL/2.0)/(DAYL + 2.0*p2))
! And now calculate temperature during the night
   if(time > sunset .OR. time < sunrise +  p1) then
! To make sure that it always counts the hours since sunset, especially after midnight
     if(time < sunrise +  p1) then 
       time = time + 24.0
     endif   
! Exponential decay of air temperature. Parameter taun (h) is the time coefficient. 
     Tair(I) = (Tmin - Tss*exp(-(24.0 - DAYL)/taun) + (Tss - Tmin)*exp(-(time - sunset)/taun))/(1.0 - exp(-(24.0 - DAYL)/taun))  
   endif 
enddo Night
RETURN
END

!**********************************************************************
      SUBROUTINE Spitters(IDATE,SGD,ZEN,RADABV,FBEAM,OpMie,OpRayleigh,DAYL,c,Lat)
! Calculate daily course of incident PAR from daily totals of
! beam and diffuse PAR (in MJ m-2 d-1).
! INPUTS:
! RADBM - daily total beam PAR (MJ m-2 d-1)
! RADDF - daily total diffuse PAR (MJ m-2 d-1)
! ZEN - zenith angle (radians)
! IDATE - date in days-since-1950 format
! OpMie & OpRayleigh - Options to use corrections in Spitters et al. (1986)
! DAYL - Daytime length in hours
! c - Parameter controlling the disaggregation of daily solar radiation
! Lat - Latitude of the simulated plot in radians
! OUTPUTS:
! RADABV - array of hourly incident PAR (J m-2 s-1)
! FBEAM - array of hourly beam fractions for PAR, NIR and Thermal

! Added by A. Morales on Feb 2014
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I
    REAL, intent(in) :: ZEN(MAXHRS)
    REAL, intent(in) :: SGD,DAYL,c,Lat
    Integer, intent(in) :: OpMie,OpRayleigh, IDATE
    Real, intent(out) :: RADABV(MAXHRS,3), FBEAM(MAXHRS,3)
    Real :: Beta(maxhrs), SGT
    Real :: SOD,SINB,Fdif,DailyBeta,FdifPrime,FPARdif,TauDifPAR, dec, TauDifNIR, Integral, SG, PARDF, NIRDF, PARBM, NIRBM, SO, a, b
    Real, external :: ETRAD
    INTEGER, EXTERNAL :: JDATE
! Daily solar radiation as input (MJ m-2 day-1)

! Calculate extra-terrestrial radiation (copied from CALCFBMD) in MJ m-2 day-1
    SOD = 0.0
    DO I = 1,KHRS
      SINB = SIN(PID2-ZEN(I))
! Calculate the extraterrestrial solar radiation corrected by the solar angle. SPERHR is the no. of seconds in each ?t and 1E6 is to convert from (W m-2) as generated by ETRAD into MJ m-2 day-1 through the integration
      SOD = SOD + ETRAD(IDATE,SINB)*SPERHR/1E6
    ENDDO

! Uncorrected daily diffuse fraction of solar radiation
    if(SGD/SOD < 0.07) then
     Fdif = 1
    else if (SGD/SOD > 0.07 .AND. SGD/SOD < 0.35) then 
      Fdif =  1.0 - 2.3*(SGD/SOD - 0.07)**2
    else if (SGD/SOD > 0.35 .AND. SGD/SOD < 0.75) then
      Fdif =  1.33 - 1.46*SGD/SOD
    else
      Fdif = 0.23
    endif

! Mie correction (a.k.a. circumsolar correction)
    if(OpMie == 1) then
! Calculates the solar elevation angle, beta (rad). Note that beta is a vector
      Beta = PID2 - ZEN 
! Averages the daytime beta values. Conditional sum function. Conditional counting function and vectorized maximum function (Fortran 95). This calculates the average solar elevation during daytime. 
      DailyBeta = sum(Beta, Beta>0)/(max(1,count(Beta>0)))
! Diffuse fraction of solar radiation corrected by Mie scattering (Spitters et al., 1986)      
      FdifPrime = Fdif/(1 + (1 - Fdif**2)*cos(PID2 - DailyBeta)**2*cos(DailyBeta)**3) 
    else 
! If Mie correction is not used, Fdif does not change, though the change of variable is needed for consistency in the code.    
      FdifPrime = Fdif
    endif   
! Rayleigh correction due to effect of scattering on solar spectrum
    if(OpRayleigh == 1) then
! Correction depends on both fdif and fdifprime calculated in the above   
       FPARdif = (1.0 + 0.3*(1.0 - Fdif**2.0))*FdifPrime
    else
! Even if the corection is not applied, the change in the variable is needed
      FPARdif = FdifPrime
    endif    
! Calculate tau (atmospheric transmission) for diffuse PAR with both corrections
    TauDifPAR = FPARdif*SGD/SOD 
! Calculate tau for diffuse NIR. The effect of Rayleigh is calculated indirectly by assuming equal share of energy in PAR and NIR bands for global and diffuse radiations. 
    TauDifNIR = (FdifPrime*2.0 - FPARdif)*SGD/SOD
! Calculate declination angle (dec, radian) and the parameter a and b (Goudriaan & van Laar, 1994)
   Dec = asin(-sin(pi*23.45/180.0)*cos(2.0*pi*(JDATE(IDATE) + 10.0)/365.0))             
   a = sin(Lat)*sin(Dec)                                           
   b = cos(Lat)*cos(Dec)    
! Analytical solution to integration of Sg for disaggregation (Spitters et al., 1986) in units of s
    Integral =  3600.0*(DayL*(a + c*(a**2.0 + 0.5*b**2.0)) + 24.0/pi*b*(1.0 + 1.5*c*a)*sqrt(1 - a**2.0/b**2.0))
! Calculate the FBEAM during the day (this could be vectorized...)
    SGT = 0
    DO I=1,KHRS
! Calculate instantaneous extraterrestrial solar radiation during the day (W m-2)
      SINB = SIN(PID2-ZEN(I))
      SO = ETRAD(IDATE,SINB)
! Calculate disaggregated solar radiation (W m-2)
      SG = max(cos(ZEN(I))*(1.0 + c*cos(ZEN(I)))*SGD*1e6/Integral, 0.0)
      SGT = SGT + SG
! Calculate diffuse PAR assuming constant TauDif throughout the day (W m-2)
      PARDF = min(SO*TauDifPAR*0.5, SG*0.5)
! Calculate diffuse NIR assuming constant TauDif throughout the day (W m-2)
      NIRDF = min(SO*TauDifNIR*0.5, SG*0.5)
! Calculate direct PAR from difference (W m-2)
      PARBM = SG*0.5 - PARDF
! Calculate direct NIR from difference (W m-2)
      NIRBM = SG*0.5 - NIRDF
! Calculate FBEAM for PAR
      FBEAM(I,1) =  PARDF/(PARDF + PARBM)
! Calculate FBEAM for NIR
      FBEAM(I,2) =  NIRDF/(NIRDF + NIRBM)
! Calculate total incident PAR (W m-2)
      RADABV(I,1) = PARDF + PARBM   
! Calculate total incident NIR (W m-2)
      RADABV(I,2) = NIRDF + NIRBM         
    ENDDO
    write(*,*) SGT*SPERHR/1E6/SGD
    RETURN
    END 


