module clm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CL
  use shr_sys_mod , only: shr_sys_abort ! cannot use endrun here due to circular dependency
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: clm_varctl_set    ! Set variables
  public :: cnallocate_carbon_only_set
  public :: cnallocate_carbon_only
  public :: cnallocate_carbonnitrogen_only_set
  public :: cnallocate_carbonnitrogen_only
  public :: cnallocate_carbonphosphorus_only_set
  public :: cnallocate_carbonphosphorus_only
  !
  private
  save
  !
  ! !PUBLIC TYPES:
  !
  integer , parameter, public ::  iundef = -9999999
  real(r8), parameter, public ::  rundef = -9999999._r8
  integer , parameter, public ::  fname_len = SHR_KIND_CL   ! max length of file names in this module
  !----------------------------------------------------------
  !
  ! Run control variables
  !
  ! case id
  character(len=256), public :: caseid  = ' '                            

  ! case title
  character(len=256), public :: ctitle  = ' '                            

  ! Type of run
  integer, public :: nsrest             = iundef                         

  ! Startup from initial conditions
  integer, public, parameter :: nsrStartup  = 0                          

  ! Continue from restart files
  integer, public, parameter :: nsrContinue = 1                          

  ! Branch from restart files
  integer, public, parameter :: nsrBranch   = 2                          

  ! true => allow case name to remain the same for branch run
  ! by default this is not allowed
  logical, public :: brnch_retain_casename = .false.                     

  !true => no valid land points -- do NOT run
  logical, public :: noland = .false.                                    

  ! Hostname of machine running on
  character(len=256), public :: hostname = ' '                           

  ! username of user running program
  character(len=256), public :: username = ' '                           

  ! description of this source
  character(len=256), public :: source   = "Community Land Model CLM4.0" 

  ! version of program
  character(len=256), public :: version  = " "                           

  ! dataset conventions
  character(len=256), public :: conventions = "CF-1.0"                   

  !----------------------------------------------------------
  ! Unit Numbers
  !----------------------------------------------------------
  !
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6

  !----------------------------------------------------------
  ! Output NetCDF files
  !----------------------------------------------------------

  logical, public :: outnc_large_files = .true.         ! large file support for output NetCDF files

  !----------------------------------------------------------
  ! Run input files
  !----------------------------------------------------------

  character(len=fname_len), public :: finidat    = ' '        ! initial conditions file name
  character(len=fname_len), public :: fsurdat    = ' '        ! surface data file name
  character(len=fname_len), public :: fatmgrid   = ' '        ! atm grid file name
  character(len=fname_len), public :: fatmlndfrc = ' '        ! lnd frac file on atm grid
  character(len=fname_len), public :: fatmtopo   = ' '        ! topography on atm grid
  character(len=fname_len), public :: flndtopo   = ' '        ! topography on lnd grid
  character(len=fname_len), public :: flanduse_timeseries    = ' '        ! dynamic landuse dataset
  character(len=fname_len), public :: paramfile  = ' '        ! ASCII data file with PFT physiological constants
  character(len=fname_len), public :: nrevsn     = ' '        ! restart data file name for branch run
  character(len=fname_len), public :: fsnowoptics  = ' '      ! snow optical properties file name
  character(len=fname_len), public :: fsnowaging   = ' '      ! snow aging parameters file name

 !! X. YANG  : add soil order dependent parameter file
  character(len=fname_len), public :: fsoilordercon    = ' '  ! ASCII data file with soil order dependent  constants

  !----------------------------------------------------------
  ! Flag to turn on MEGAN VOC's
  !----------------------------------------------------------

  logical, public :: use_voc = .true. 

  !----------------------------------------------------------
  ! Interpolation of finidat if requested
  !----------------------------------------------------------

  logical, public :: bound_h2osoi = .true. ! for debugging 

  ! If finidat_interp_source is non-blank and finidat is blank then interpolation will be done from
  ! finidat_interp_source to finidat_interp_dest

  character(len=fname_len), public :: finidat_interp_source = ' '
  character(len=fname_len), public :: finidat_interp_dest   = 'finidat_interp_dest.nc'     

  !----------------------------------------------------------
  ! Irrigate logic
  !----------------------------------------------------------

  ! do not irrigate by default
  logical, public :: irrigate = .false.            

  !----------------------------------------------------------
  ! Landunit logic
  !----------------------------------------------------------

  ! true => separate crop landunit is not created by default
  logical, public :: create_crop_landunit = .false.     
  
  !----------------------------------------------------------
  ! Other subgrid logic
  !----------------------------------------------------------

  ! true => make ALL patches, cols & landunits active (even if weight is 0)
  logical, public :: all_active = .false.          

  !----------------------------------------------------------
  ! BGC logic and datasets
  !----------------------------------------------------------

  ! values of 'prognostic','diagnostic','constant'
  character(len=16), public :: co2_type = 'constant'    

  ! State of the model for the accelerated decomposition (AD) spinup. 
  ! 0 (default) = normal model; 1 = AD SPINUP
  integer, public :: spinup_state = 0 

  ! true => anoxia is applied to heterotrophic respiration also considered in CH4 model
  ! default value reset in controlMod
  logical, public :: anoxia  = .true. 

  ! used to override an error check on reading in restart files
  logical, public :: override_bgc_restart_mismatch_dump = .false. 

  ! Set in CNAllocationInit (TODO - had to move it here to avoid circular dependency)
  logical, private:: carbon_only      
  logical, private:: carbonnitrogen_only      
  logical, private:: carbonphosphorus_only      

  !----------------------------------------------------------
  ! Physics
  !----------------------------------------------------------

  ! use subgrid fluxes
  integer,  public :: subgridflag = 1                   

  ! true => write global average diagnostics to std out
  logical,  public :: wrtdia       = .false.            

  ! atmospheric CO2 molar ratio (by volume) (umol/mol)
  real(r8), public :: co2_ppmv     = 355._r8            !

  !----------------------------------------------------------
  ! C isotopes
  !----------------------------------------------------------

  logical, public :: use_c13 = .false.                  ! true => use C-13 model
  logical, public :: use_c14 = .false.                  ! true => use C-14 model

  !----------------------------------------------------------
  !  ED switches
  !----------------------------------------------------------

  logical, public :: use_ed = .false.            ! true => use  ED
  logical, public :: use_ed_spit_fire = .false.  ! true => use spitfire model

  !----------------------------------------------------------
  ! lai streams switch for Sat. Phenology
  !----------------------------------------------------------

  logical, public :: use_lai_streams = .false. ! true => use lai streams in SatellitePhenologyMod.F90

  !----------------------------------------------------------
  ! glacier_mec control variables: default values (may be overwritten by namelist)
  ! NOTE: glc_smb must have the same values for CLM and GLC
  !----------------------------------------------------------

  ! glacier_mec landunit is not created (set in controlMod)
  logical , public :: create_glacier_mec_landunit = .false. 

  ! if true, pass surface mass balance info to GLC
  logical , public :: glc_smb = .true.                      

  ! if false, pass positive-degree-day info to GLC

  ! true => CLM glacier area & topography changes dynamically 
  logical , public :: glc_do_dynglacier = .false.           

  ! true => downscale precip division into rain & snow
  logical , public :: glcmec_downscale_rain_snow_convert = .false.     

  ! true => downscale longwave radiation
  logical , public :: glcmec_downscale_longwave = .true.    

  ! number of days before one considers the perennially snow-covered point 'land ice'
  integer , public :: glc_snow_persistence_max_days = 7300  

  ! glc_grid used to determine fglcmask  
  character(len=256), public :: glc_grid = ' '              

  ! glacier mask file name (based on glc_grid)
  character(len=fname_len), public :: fglcmask = ' '        
  !
  !----------------------------------------------------------
  ! single column control variables
  !----------------------------------------------------------

  logical,  public :: single_column = .false. ! true => single column mode
  real(r8), public :: scmlat        = rundef  ! single column lat
  real(r8), public :: scmlon        = rundef  ! single column lon

  !----------------------------------------------------------
  ! instance control
  !----------------------------------------------------------

  integer, public :: inst_index
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix

  !----------------------------------------------------------
  ! Decomp control variables
  !----------------------------------------------------------

  ! number of segments per clump for decomp
  integer, public :: nsegspc = 20                       

  !----------------------------------------------------------
  ! Derived variables (run, history and restart file)
  !----------------------------------------------------------

  ! directory name for local restart pointer file
  character(len=256), public :: rpntdir = '.'            

  ! file name for local restart pointer file
  character(len=256), public :: rpntfil = 'rpointer.lnd' 

  ! moved hist_wrtch4diag from histFileMod.F90 to here - caused compiler error with intel
  ! namelist: write CH4 extra diagnostic output
  logical, public :: hist_wrtch4diag = .false.         

  !----------------------------------------------------------
  ! Migration of CPP variables
  !----------------------------------------------------------

  logical, public :: use_nofire          = .false.
  logical, public :: use_lch4            = .false.
  logical, public :: use_nitrif_denitrif = .false.
  logical, public :: use_vertsoilc       = .false.
  logical, public :: use_extralakelayers = .false.
  logical, public :: use_vichydro        = .false.
  logical, public :: use_century_decomp  = .false.
  logical, public :: use_cn              = .false.
  logical, public :: use_cndv            = .false.
  logical, public :: use_crop            = .false.
  logical, public :: use_snicar_frc      = .false.
  logical, public :: use_vancouver       = .false.
  logical, public :: use_mexicocity      = .false.
  logical, public :: use_noio            = .false.

  !----------------------------------------------------------
  ! To retrieve namelist
  !----------------------------------------------------------
  character(len=SHR_KIND_CL), public :: NLFilename_in ! Namelist filename
  !
  logical, private :: clmvarctl_isset = .false.
 !-----------------------------------------------------------------------

contains

  !---------------------------------------------------------------------------
  subroutine clm_varctl_set( caseid_in, ctitle_in, brnch_retain_casename_in,    &
       single_column_in, scmlat_in, scmlon_in, nsrest_in, &
       version_in, hostname_in, username_in)
    !
    ! !DESCRIPTION:
    ! Set input control variables.
    !
    ! !ARGUMENTS:
    character(len=256), optional, intent(IN) :: caseid_in                ! case id
    character(len=256), optional, intent(IN) :: ctitle_in                ! case title
    logical,            optional, intent(IN) :: brnch_retain_casename_in ! true => allow case name to remain the 
                                                                         ! same for branch run
    logical,            optional, intent(IN) :: single_column_in         ! true => single column mode
    real(r8),           optional, intent(IN) :: scmlat_in                ! single column lat
    real(r8),           optional, intent(IN) :: scmlon_in                ! single column lon
    integer,            optional, intent(IN) :: nsrest_in                ! 0: initial run. 1: restart: 3: branch
    character(len=256), optional, intent(IN) :: version_in               ! model version
    character(len=256), optional, intent(IN) :: hostname_in              ! hostname running on
    character(len=256), optional, intent(IN) :: username_in              ! username running job
    !-----------------------------------------------------------------------

    if ( clmvarctl_isset )then
       call shr_sys_abort(' ERROR:: control variables already set, cannot call this routine')
    end if

    if ( present(caseid_in       ) ) caseid        = caseid_in
    if ( present(ctitle_in       ) ) ctitle        = ctitle_in
    if ( present(single_column_in) ) single_column = single_column_in
    if ( present(scmlat_in       ) ) scmlat        = scmlat_in
    if ( present(scmlon_in       ) ) scmlon        = scmlon_in
    if ( present(nsrest_in       ) ) nsrest        = nsrest_in
    if ( present(brnch_retain_casename_in) ) brnch_retain_casename = brnch_retain_casename_in
    if ( present(version_in      ) ) version       = version_in
    if ( present(username_in     ) ) username      = username_in
    if ( present(hostname_in     ) ) hostname      = hostname_in

  end subroutine clm_varctl_set

  ! Set module carbon_only flag
  subroutine cnallocate_carbon_only_set(carbon_only_in)
    logical, intent(in) :: carbon_only_in
    carbon_only = carbon_only_in
  end subroutine cnallocate_carbon_only_set

  ! Get module carbon_only flag
  logical function CNAllocate_Carbon_only()
    cnallocate_carbon_only = carbon_only
  end function CNAllocate_Carbon_only


  ! Set module carbonnitrogen_only flag
  subroutine cnallocate_carbonnitrogen_only_set(carbonnitrogen_only_in)
    logical, intent(in) :: carbonnitrogen_only_in
    carbonnitrogen_only = carbonnitrogen_only_in
  end subroutine cnallocate_carbonnitrogen_only_set

  ! Get module carbonnitrogen_only flag
  logical function CNAllocate_CarbonNitrogen_only()
    cnallocate_carbonnitrogen_only = carbonnitrogen_only
  end function CNAllocate_CarbonNitrogen_only


  ! Set module carbonphosphorus_only flag
  subroutine cnallocate_carbonphosphorus_only_set(carbonphosphorus_only_in)
    logical, intent(in) :: carbonphosphorus_only_in
    carbonphosphorus_only = carbonphosphorus_only_in
  end subroutine cnallocate_carbonphosphorus_only_set

  ! Get module carbonphosphorus_only flag
  logical function CNAllocate_CarbonPhosphorus_only()
    cnallocate_carbonphosphorus_only = carbonphosphorus_only
  end function CNAllocate_CarbonPhosphorus_only



end module clm_varctl
