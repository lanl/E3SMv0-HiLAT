!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module tracegas_mod

!BOP
! !MODULE: tracegas_mod
!
! !DESCRIPTION:
!
 !------------------------------------------------------------------------------
 !   Multispecies trace gas processing routine based on Chu, Elliott, Maltrud,
 !   Erickson and company papers > 2000, on references therein,
 !   and also on assorted supplementary materials as indicated in comments.
 !   Designed for insertion into POP and CCSM along with DML ecosys_mod
 !   Trace gases defined roughly as those of concentration order nanomolar or less
 !   This excludes CO2 and O2 which are handled inside NCAR ecodynamics
 !------------------------------------------------------------------------------
 !------------------------------------------------------------------------------
 !   variables/subroutines/functions used from other modules
 !   The following are called upon extensively in tracegas, and so appear at
 !   the module level. The use statements for variables that are only needed
 !   locally are dealt with at the module subprogram level.
 !------------------------------------------------------------------------------
 !------------------------------------------------------------------------------
 !   More recent update in Elliott, 2009, JGR; Wang et al., 2015, JGR-Biogeosci
 !------------------------------------------------------------------------------

! !REVISION HISTORY:

!  SVN:$Id:  $

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! !USES:

!maltrud debug
   use POP_CommMod

   use POP_KindsMod
   use POP_ErrorMod

   use kinds_mod
   use constants
   use communicate
   use broadcast
   use global_reductions
   use blocks
   use domain_size
   use domain
   use exit_mod
   use prognostic
   use grid
   use io
   use io_types
   use io_tools
   use tavg
   use timers
   use passive_tracer_tools
   use named_field_mod
   use forcing_tools
   use time_management
   use ecosys_parms
!   use tracegas_parms (SW:not used)
   use registry
   use named_field_mod
   use ecosys_mod, only: &
!       po4_ind,  & ! dissolved inorganic phosphate
!       no3_ind,  & ! dissolved inorganic nitrate
!       sio3_ind,  & ! dissolved inorganic silicate
!       nh4_ind,  & ! dissolved ammonia
!       fe_ind,  & ! dissolved inorganic iron
!       o2_ind,  & ! dissolved oxygen
!       dic_ind,  & ! dissolved inorganic carbon
!       alk_ind,  & ! alkalinity
       doc_ind,  & ! dissolved organic carbon
       zooC_ind,  & ! zooplankton carbon
!       don_ind,  & ! dissolved organic nitrogen
       ecosys_tracer_cnt, &! size of total eco array
       auto_ind, &! autotrophs index (SW)
       PAR_out, PAR_avg, comp_surf_avg
!       spCk,spChlk,diatCk,diazCk,spCaCO3k

! !INPUT PARAMETERS:
  !-----------------------------------------------------------------------------
  !   include tracegas and ecosystem parameters
  !   all variables from these modules have a parm_ prefix
  !-----------------------------------------------------------------------------
!(swang: not used)
!   use tracegas_parms
!   use ecosys_parms, only : parm_kappa_nitrif, parm_nitrif_par_lim

   implicit none
   save
   private

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
       tracegas_tracer_cnt, &
       tracegas_init, &
       tracegas_init_sflux,  &
!       tracegas_init_interior_restore, &
       tracegas_tracer_ref_val,        &
       tracegas_set_interior, &
       tracegas_set_sflux, &
       tracegas_write_restart, &
       tracegas_tavg_forcing

!-----------------------------------------------------------------------
!  module variables required by forcing_passive_tracer
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      tracegas_tracer_cnt = 2

!-----------------------------------------------------------------------
!  flags controlling which portion of code are executed
!  usefull for debugging
!-----------------------------------------------------------------------
   
  logical (log_kind) :: &
     tracegas_diurnal_cycle, &
     lsource_sink, &
     lflux_gas_dms
   
  logical (log_kind), dimension(:,:,:), allocatable :: &
     LAND_MASK
   
!----------------------------------------------------------------------------
!  relative tracegas tracer indices
!------------------------------------------------------------------------------

   integer(int_kind), parameter :: &
      dms_ind          = 1,  & ! dms
      dmsp_ind         = 2     ! dmsp

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(tracegas_tracer_cnt) :: &
      ind_name_table

!-----------------------------------------------------------------------
!  options for forcing of gas fluxes
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      gas_flux_forcing_iopt_model = 1,   &
      gas_flux_forcing_iopt_file  = 2

   integer (int_kind) :: &
      gas_flux_forcing_iopt

   character(char_len) :: &
      gas_flux_forcing_file    ! file containing gas flux forcing fields

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   type(tracer_read) :: &
      gas_flux_fice,       & ! ice fraction for gas fluxes
      gas_flux_ws,         & ! wind speed for gas fluxes
      gas_flux_ap            ! atmospheric pressure for gas fluxes

   type(forcing_monthly_every_ts) :: &
      fice_file,                 & ! ice fraction, if read from file
      xkw_file,                  & ! a * wind-speed ** 2, if read from file
      ap_file                      ! atmoshperic pressure, if read from file

!-----------------------------------------------------------------------
!  restoring climatologies for dms
!-----------------------------------------------------------------------

!   character(char_len) :: &
!      dms_rest_file               ! file containing dms field

!maltrud variable restoring
!   logical (log_kind) :: &
!      ldms_variable_restore       ! geographically varying dms restoring

!   character(char_len) :: &
!      dms_variable_rest_file,   & ! file containing variable restoring info
!      dms_variable_rest_file_fmt  ! format of file containing variable restoring info

!   real (r8), dimension(:,:,:), allocatable, target :: &
!      DMS_RESTORE_RTAU            ! inverse restoring timescale for variable
                                   ! interior restoring

!   integer (int_kind), dimension(:,:,:), allocatable :: &
!      DMS_RESTORE_MAX_LEVEL       ! maximum level for applying variable
                                   ! interior restoring

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK                  ! temp array for interpolate_forcing output

!-----------------------------------------------------------------------------
!   define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_DMS_SURF,      &! tavg id for DMS at surface
      tavg_DMS_WS,           &! tavg id for wind speed
      tavg_DMS_IFRAC,         &! tavg id for ice fraction
      tavg_DMS_XKW,          &! tavg id for xkw
      tavg_DMS_ATM_PRESS,    &! tavg id for atmospheric pressure
      tavg_DMS_SCHMIDT,   &! tavg id for DMS schmidt number
      tavg_DMS_SAT,        &! tavg id for DMS saturation
      tavg_DMS_PV          ! tavg id for DMS piston velocity

!-----------------------------------------------------------------------------
!   define tavg id for source/sink terms
!-----------------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_DMS_S_TOTAL,  &
      tavg_DMS_S_DMSP,  &
      tavg_DMS_R_TOTAL,  &
      tavg_DMS_R_B,  &
      tavg_DMS_R_PHOT,  &
      tavg_DMS_R_BKGND,  &
      tavg_DMSP_S_PHAEO,  &
      tavg_DMSP_S_NONPHAEO,  &
      tavg_DMSP_S_ZOO,  &
      tavg_DMSP_S_TOTAL,  &
      tavg_DMSP_R_B,  &
      tavg_DMSP_R_BKGND,  &
      tavg_DMSP_R_TOTAL,  &
      tavg_Cyano_frac,  &
      tavg_Cocco_frac,  &
      tavg_Eukar_frac,  &
!      tavg_Phaeo_frac,  &		(swang: not used with explicit phaeo)
      tavg_diatS,  &
      tavg_diatN,  &
      tavg_phytoN,  &
      tavg_coccoS,  &
      tavg_cyanoS,  &
      tavg_eukarS,  &
      tavg_diazS,  &
      tavg_phaeoS,  &
! (swang) to remove output of phaeoC generated by tracegas module & add phaeo
!      tavg_phaeoC,  &	
	  tavg_phaeoN,	 &		
      tavg_phaeonS,  &
	  tavg_phaeonN,	 &		
      tavg_zooCC,  &
      tavg_zooS,  &
      tavg_RSNzoo

!-----------------------------------------------------------------------
!  define array for holding flux-related quantities that need to be time-averaged
!  this is necessary since the forcing routines are called before tavg flags
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable :: &
      TRACEGAS_SFLUX_TAVG

!-----------------------------------------------------------------------
!  average surface tracer value related variables
!  used as reference value for virtual flux computations
!-----------------------------------------------------------------------

   logical (log_kind), dimension(tracegas_tracer_cnt) :: &
      vflux_flag                ! which tracers get virtual fluxes applied

   integer (int_kind) :: &
      comp_surf_avg_flag        ! time flag id for computing average
                                ! surface tracer values

   real (r8), dimension(tracegas_tracer_cnt) :: &
      surf_avg                  ! average surface tracer values

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tracegas_interior_timer,                   &
      tracegas_sflux_timer

!-----------------------------------------------------------------------
!  named field indices
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      surf_dms_nf_ind   = 0,    & ! surface DMS
      surf_dmsp_nf_ind  = 0,    & ! surface DMSP
      sflux_dms_nf_ind  = 0       ! dms flux

!-----------------------------------------------------------------------

! (swang: not used!)
!   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
!      UV_out           ! generic UV radiation (W/m^2)

!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: tracegas_init
! !INTERFACE:

 subroutine tracegas_init(init_ts_file_fmt, read_restart_filename, &
                        tracer_d_module, TRACER_MODULE, tadvect_ctype, &
                        errorCode)

! !DESCRIPTION:
!  Initialize tracegas tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(tracegas_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,tracegas_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   character (char_len), dimension(:), intent(out) :: &
      tadvect_ctype     ! advection method for tracegas tracers

   integer (POP_i4), intent(out) :: &
      errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'tracegas_mod:tracegas_init'

   character(char_len) :: &
      init_tracegas_option,        & ! option for initialization of bgc
      init_tracegas_init_file,     & ! filename for option 'file'
      init_tracegas_init_file_fmt, & ! file format for option 'file'
      comp_surf_avg_freq_opt,      & ! choice for freq of comp_surf_avg
      gas_flux_forcing_opt,        & ! option for forcing gas fluxes
      tracegas_tadvect_ctype         ! advection method for tracegas tracers

   type(tracer_read), dimension(tracegas_tracer_cnt) :: &
      tracer_init_ext              ! namelist variable for initializing tracers

   logical (log_kind) :: &
      default,                   & ! arg to init_time_flag
      lnml_found,                & ! Was tracegas_nml found ?
      lmarginal_seas               ! Is tracegastem active in marginal seas ?

   integer (int_kind) :: &
      n,                         & ! index for looping over tracers
      k,                         & ! index for looping over depth levels
      l,                         & ! index for looping over time levels
      ind,                       & ! tracer index for tracer name from namelist
      iblock,                    & ! index for looping over blocks
      nml_error                    ! namelist i/o error flag

   integer (int_kind) :: &
      freq_opt, freq,            & ! args for init_time_flag
      comp_surf_avg_freq_iopt,   & ! choice for freq of comp_surf_avg
      comp_surf_avg_freq           ! choice for freq of comp_surf_avg

   logical (log_kind) :: &
      use_nml_surf_vals            ! do namelist surf values override values from restart file

!-----------------------------------------------------------------------
!  values to be used when comp_surf_avg_freq_opt==never
!-----------------------------------------------------------------------

   namelist /tracegas_nml/ &
      init_tracegas_option, init_tracegas_init_file, tracer_init_ext, &
      init_tracegas_init_file_fmt, &
      gas_flux_forcing_opt, &
      gas_flux_fice, gas_flux_ws, gas_flux_ap, &
!      rest_time_inv_surf, rest_time_inv_deep, rest_z0, rest_z1, &
      comp_surf_avg_freq_opt, comp_surf_avg_freq,  &
      use_nml_surf_vals,   &
      lmarginal_seas, tracegas_diurnal_cycle, &
      lsource_sink, lflux_gas_dms, &
!      ldms_variable_restore, dms_variable_rest_file,  &
!      dms_variable_rest_file_fmt, &
      tracegas_tadvect_ctype

   character (char_len) :: &
      tracegas_restart_filename  ! modified file name for restart file

   real (r8), dimension (nx_block,ny_block) :: WORK

!-----------------------------------------------------------------------
!  initialize name table 
!-----------------------------------------------------------------------

   ind_name_table( 1) = ind_name_pair(dms_ind,     'DMS')
   ind_name_table( 2) = ind_name_pair(dmsp_ind,    'DMSP')

!-----------------------------------------------------------------------
!  initialize forcing_monthly_every_ts variables
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_forcing_monthly_every_ts(fice_file)
   call init_forcing_monthly_every_ts(xkw_file)
   call init_forcing_monthly_every_ts(ap_file)

!-----------------------------------------------------------------------
!  initialize tracegas parameters
!-----------------------------------------------------------------------

!   call tracegas_parms_init	(SW:not used)

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   do n = 1, tracegas_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'mmol/m^3'
      tracer_d_module(n)%tend_units = 'mmol/m^3/s'
      tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
   end do
   tracer_d_module(dms_ind    )%long_name='DiMethyl Sulfide'
   tracer_d_module(dmsp_ind   )%long_name='DiMethyl SulfiPropitionate'

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_tracegas_option = 'unknown'
   init_tracegas_init_file = 'unknown'
   init_tracegas_init_file_fmt = 'bin'

   gas_flux_forcing_opt  = 'drv'
   gas_flux_forcing_file = 'unknown'

   gas_flux_fice%filename     = 'unknown'
   gas_flux_fice%file_varname = 'FICE'
   gas_flux_fice%scale_factor = c1
   gas_flux_fice%default_val  = c0
   gas_flux_fice%file_fmt     = 'bin'

   gas_flux_ws%filename     = 'unknown'
   gas_flux_ws%file_varname = 'XKW'
   gas_flux_ws%scale_factor = c1
   gas_flux_ws%default_val  = c0
   gas_flux_ws%file_fmt     = 'bin'

   gas_flux_ap%filename     = 'unknown'
   gas_flux_ap%file_varname = 'P'
   gas_flux_ap%scale_factor = c1
   gas_flux_ap%default_val  = c0
   gas_flux_ap%file_fmt     = 'bin'

!   rest_time_inv_surf = c0
!   rest_time_inv_deep = c0
!   rest_z0            = c1000
!   rest_z1            = c2 * c1000

!maltrud variable restoring
!   ldms_variable_restore      = .false.
!   dms_variable_rest_file     = 'unknown'
!   dms_variable_rest_file_fmt = 'bin'

   do n = 1,tracegas_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   lmarginal_seas        = .true.
   tracegas_diurnal_cycle  = .false.
   lsource_sink          = .true.
   lflux_gas_dms         = .true.

   comp_surf_avg_freq_opt        = 'never'
   comp_surf_avg_freq            = 1
   use_nml_surf_vals             = .false.

   tracegas_tadvect_ctype = 'base_model'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=tracegas_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'tracegas_nml not found')
      call exit_POP(sigAbort, 'ERROR : stopping in '/&
                           &/ subname)
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' tracegas:'
      write(stdout,blank_fmt)
      write(stdout,*) ' tracegas_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,tracegas_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_tracegas_option, master_task)
   call broadcast_scalar(init_tracegas_init_file, master_task)
   call broadcast_scalar(init_tracegas_init_file_fmt, master_task)

   call broadcast_scalar(gas_flux_forcing_opt, master_task)
   if (trim(gas_flux_forcing_opt) == 'drv') then
      gas_flux_forcing_iopt = gas_flux_forcing_iopt_model
   else if (trim(gas_flux_forcing_opt) == 'file') then
      gas_flux_forcing_iopt = gas_flux_forcing_iopt_file
   else
      call document(subname, 'gas_flux_forcing_opt', gas_flux_forcing_opt)
      call exit_POP(sigAbort, 'unknown gas_flux_forcing_opt')
   endif

   call broadcast_scalar(gas_flux_forcing_file, master_task)

   call broadcast_scalar(gas_flux_fice%filename, master_task)
   call broadcast_scalar(gas_flux_fice%file_varname, master_task)
   call broadcast_scalar(gas_flux_fice%scale_factor, master_task)
   call broadcast_scalar(gas_flux_fice%default_val, master_task)
   call broadcast_scalar(gas_flux_fice%file_fmt, master_task)

   fice_file%input = gas_flux_fice

   call broadcast_scalar(gas_flux_ws%filename, master_task)
   call broadcast_scalar(gas_flux_ws%file_varname, master_task)
   call broadcast_scalar(gas_flux_ws%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ws%default_val, master_task)
   call broadcast_scalar(gas_flux_ws%file_fmt, master_task)

   xkw_file%input = gas_flux_ws

   call broadcast_scalar(gas_flux_ap%filename, master_task)
   call broadcast_scalar(gas_flux_ap%file_varname, master_task)
   call broadcast_scalar(gas_flux_ap%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ap%default_val, master_task)
   call broadcast_scalar(gas_flux_ap%file_fmt, master_task)

   ap_file%input = gas_flux_ap

!   call broadcast_scalar(rest_time_inv_surf, master_task)
!   call broadcast_scalar(rest_time_inv_deep, master_task)
!   call broadcast_scalar(rest_z0, master_task)
!   call broadcast_scalar(rest_z1, master_task)

!maltrud variable restoring
!   call broadcast_scalar(ldms_variable_restore, master_task)
!   call broadcast_scalar(dms_variable_rest_file, master_task)
!   call broadcast_scalar(dms_variable_rest_file_fmt, master_task)

   do n = 1,tracegas_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(comp_surf_avg_freq_opt, master_task)
   call broadcast_scalar(comp_surf_avg_freq, master_task)
   call broadcast_scalar(use_nml_surf_vals, master_task)

   call broadcast_scalar(lmarginal_seas, master_task)
   call broadcast_scalar(tracegas_diurnal_cycle, master_task)
   call broadcast_scalar(lsource_sink, master_task)
   call broadcast_scalar(lflux_gas_dms, master_task)

   call broadcast_scalar(tracegas_tadvect_ctype, master_task)
   tadvect_ctype = tracegas_tadvect_ctype

!-----------------------------------------------------------------------
!  set variables immediately dependent on namelist variables
!-----------------------------------------------------------------------

   select case (comp_surf_avg_freq_opt)
   case ('never')
      comp_surf_avg_freq_iopt = freq_opt_never
   case ('nyear')
      comp_surf_avg_freq_iopt = freq_opt_nyear
   case ('nmonth')
      comp_surf_avg_freq_iopt = freq_opt_nmonth
   case default
      call document(subname, 'comp_surf_avg_freq_opt', comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'unknown comp_surf_avg_freq_opt')
   end select

   call init_time_flag('tracegas_comp_surf_avg', comp_surf_avg_flag, &
      default=.false., freq_opt=comp_surf_avg_freq_iopt,  &
      freq=comp_surf_avg_freq, owner='tracegas_init')

!-----------------------------------------------------------------------
!  namelist consistency checking
!-----------------------------------------------------------------------

   if (use_nml_surf_vals .and. comp_surf_avg_freq_iopt /= freq_opt_never) then
      call document(subname, 'use_nml_surf_vals', use_nml_surf_vals)
      call document(subname, 'comp_surf_avg_freq_opt', comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'use_nml_surf_vals can only be .true. if ' /&
                           &/ ' comp_surf_avg_freq_opt is never')
   endif

!-----------------------------------------------------------------------
!  initialize virtual flux flag array
!-----------------------------------------------------------------------

   vflux_flag = .false.

!-----------------------------------------------------------------------
!  allocate and initialize LAND_MASK
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,nblocks_clinic) )

   if (lmarginal_seas) then
!maltrud debug - swang: not for CESM runs
      LAND_MASK = REGION_MASK /= 0
!      LAND_MASK = REGION_MASK > 0
   else
      LAND_MASK = REGION_MASK > 0
   endif

!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (init_tracegas_option)

!   case ('restart', 'ccsm_continue', 'branch', 'ccsm_hybrid' )
   case ('restart', 'ccsm_continue', 'branch', 'ccsm_branch', 'ccsm_hybrid' ) !PJC: added ccsm_branch. Didn't replace 'branch' since it is a bit different in stand-alone POP (according to Mat Maltrud).

      tracegas_restart_filename = char_blank

      if (init_tracegas_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read tracegas from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         tracegas_restart_filename = read_restart_filename
         init_tracegas_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         tracegas_restart_filename = trim(init_tracegas_init_file)

      endif

      call rest_read_tracer_block(init_tracegas_init_file_fmt, &
                                  tracegas_restart_filename,   &
                                  tracer_d_module,           &
                                  TRACER_MODULE)

      if (use_nml_surf_vals) then
         surf_avg = c0
      else
         call extract_surf_avg(init_tracegas_init_file_fmt, &
                               tracegas_restart_filename)
      endif

      call eval_time_flag(comp_surf_avg_flag) ! evaluates time_flag(comp_surf_avg_flag)%value via time_to_do

      if (check_time_flag(comp_surf_avg_flag)) &
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:),TRACER_MODULE(:,:,1,:,curtime,:))

   case ('file', 'startup')
      call document(subname, 'tracegas vars being read from separate files')

      call file_read_tracer_block(init_tracegas_init_file_fmt, &
                                  init_tracegas_init_file,     &
                                  tracer_d_module,           &
                                  ind_name_table,            &
                                  tracer_init_ext,           &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, tracegas_tracer_cnt
            do k=1,km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'ecosys_init: error in fill points for tracers')
                  return
               endif

            enddo
         enddo
      endif

      if (use_nml_surf_vals) then
         surf_avg = c0
      else
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:),TRACER_MODULE(:,:,1,:,curtime,:))
      endif

   case ('default_value', 'ccsm_startup')
      call document(subname, 'tracegas vars being set to default value')

      do n = 1,tracegas_tracer_cnt
         do k=1,km
            TRACER_MODULE(:,:,k,n,curtime,:) = &
               tracer_init_ext(n)%scale_factor * tracer_init_ext(n)%default_val
         enddo
      end do

      if (n_topo_smooth > 0) then
         do n = 1, tracegas_tracer_cnt
            do k=1,km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'ecosys_init: error in fill points for tracers')
                  return
               endif

            enddo
         enddo
      endif

      if (use_nml_surf_vals) then
         surf_avg = c0
      else
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:),TRACER_MODULE(:,:,1,:,curtime,:))
      endif

   case default
      call document(subname, 'init_tracegas_option', init_tracegas_option)
      call exit_POP(sigAbort, 'unknown init_tracegas_option')

   end select

!-----------------------------------------------------------------------
!  register surf_dms field for passing to ice module; set surf_dms field 
!-----------------------------------------------------------------------

   call named_field_register('oceanSurfaceDMS', surf_dms_nf_ind)

   !$OMP PARALLEL DO PRIVATE(iblock,n,k,WORK)
   do iblock=1,nblocks_clinic
      do n = 1,tracegas_tracer_cnt
         do k = 1,km
            where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do

      WORK = c0
      WORK = max(c0,p5*(TRACER_MODULE(:,:,1,dms_ind,oldtime,iblock) + &
                        TRACER_MODULE(:,:,1,dms_ind,curtime,iblock)))
                         
      call named_field_set(surf_dms_nf_ind, iblock, WORK)
   enddo
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  register surf_dmsp field for passing to ice module; set surf_dmsp field 
!-----------------------------------------------------------------------

   call named_field_register('oceanSurfaceDMSP', surf_dmsp_nf_ind)

   !$OMP PARALLEL DO PRIVATE(iblock,n,k,WORK)
   do iblock=1,nblocks_clinic
      do n = 1,tracegas_tracer_cnt
         do k = 1,km
            where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do

      WORK = c0
      WORK = max(c0,p5*(TRACER_MODULE(:,:,1,dmsp_ind,oldtime,iblock) + &
                        TRACER_MODULE(:,:,1,dmsp_ind,curtime,iblock)))
                         
      call named_field_set(surf_dmsp_nf_ind, iblock, WORK)
   enddo
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  timer init
!-----------------------------------------------------------------------

   call get_timer(tracegas_interior_timer, 'TRACEGAS_INTERIOR', &
                  nblocks_clinic, distrb_clinic%nprocs)
   call get_timer(tracegas_sflux_timer, 'TRACEGAS_SFLUX',1, &
                  distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call tracegas_init_tavg
   call tracegas_init_sflux
!   call tracegas_init_interior_restore

!-----------------------------------------------------------------------
!EOC

 end subroutine tracegas_init

!***********************************************************************
!BOP
! !IROUTINE: extract_surf_avg
! !INTERFACE:

 subroutine extract_surf_avg(init_tracegas_init_file_fmt, &
                             tracegas_restart_filename)

! !DESCRIPTION:
!  Extract average surface values from restart file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_tracegas_init_file_fmt, & ! file format (bin or nc)
      tracegas_restart_filename      ! file name for restart file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type (datafile) ::&
      restart_file    ! io file descriptor

   integer (int_kind) :: &
      n               ! tracer index

   character (char_len) :: &
      short_name      ! tracer name temporaries

!-----------------------------------------------------------------------

   surf_avg = c0

   restart_file = construct_file(init_tracegas_init_file_fmt, &
                                 full_name=trim(tracegas_restart_filename), &
                                 record_length=rec_type_dbl, &
                                 recl_words=nx_global*ny_global)

   do n = 1, tracegas_tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call add_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set(restart_file, 'open_read')

   do n = 1, tracegas_tracer_cnt
      if (vflux_flag(n)) then
         short_name = 'surf_avg_' /&
                   &/ ind_name_table(n)%name
         call extract_attrib_file(restart_file, trim(short_name), surf_avg(n))
      endif
   end do

   call data_set (restart_file, 'close')

   call destroy_file (restart_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine extract_surf_avg

!***********************************************************************
!BOP
! !IROUTINE: tracegas_init_tavg
! !INTERFACE:

 subroutine tracegas_init_tavg

! !DESCRIPTION:
!  call define_tavg_field for nonstandard tavg fields
!
! !REVISION HISTORY:
!  same as module
!

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      var_cnt             ! how many tavg variables are defined

!-----------------------------------------------------------------------
!  2D fields related to surface fluxes
!-----------------------------------------------------------------------

   var_cnt = 0

   call define_tavg_field(tavg_DMS_IFRAC,'DMS_IFRAC',2,          &
                          long_name='Ice Fraction for DMS fluxes',  &
                          units='fraction', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_XKW,'DMS_XKW',2,              &
                          long_name='XKW for DMS fluxes',           &
                          units='cm/s', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_ATM_PRESS,'DMS_ATM_PRESS',2,  &
                          long_name='Atmospheric Pressure for DMS fluxes', &
                          units='atmospheres', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_PV,'DMS_PV',2,                        &
                          long_name='Piston Velocity for DMS fluxes',    &
                          units='cm/s', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_SCHMIDT,'DMS_SCHMIDT',2,              &
                          long_name='DMS Schmidt Number',               &
                          units='none', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_SAT,'DMS_SAT',2,                        &
                          long_name='DMS Saturation',                   &
                          units='mmol/m^3', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_SURF,'DMS_SURF',2,                        &
                          long_name='DMS Surface Values',                   &
                          units='mmol/m^3', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_WS,'DMS_WS',2,                        &
                          long_name='Windspeed used for DMS Flux',      &
                          units='m/s', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_DMS_S_TOTAL,'DMS_S_TOTAL',3,                        &
                          long_name='Total DMS source',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMS_S_DMSP,'DMS_S_DMSP',3,                        &
                          long_name='DMS source from DMSP',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMS_R_TOTAL,'DMS_R_TOTAL',3,                        &
                          long_name='Total DMS removal',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMS_R_B,'DMS_R_B',3,                        &
                          long_name='DMS removal from bacteria',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMS_R_PHOT,'DMS_R_PHOT',3,                        &
                          long_name='DMS removal from photolysis',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMS_R_BKGND,'DMS_R_BKGND',3,                        &
                          long_name='DMS removal from background',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMSP_S_TOTAL,'DMSP_S_TOTAL',3,                        &
                          long_name='Total DMSP source',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMSP_S_PHAEO,'DMSP_S_PHAEO',3,                        &
                          long_name='DMSP source from Phaeocystis',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMSP_S_NONPHAEO,'DMSP_S_NONPHAEO',3,                        &
                          long_name='DMSP source from non-Phaeo Phyto',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMSP_S_ZOO,'DMSP_S_ZOO',3,                        &
                          long_name='DMSP source from Zooplankton',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMSP_R_TOTAL,'DMSP_R_TOTAL',3,                        &
                          long_name='Total DMSP removal',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMSP_R_B,'DMSP_R_B',3,                        &
                          long_name='DMSP removal from bacteria',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_DMSP_R_BKGND,'DMSP_R_BKGND',3,                        &
                          long_name='DMSP removal from background',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_Cyano_frac,'Cyano_frac',3,                        &
                          long_name='Cyanobacteria Fraction',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_Cocco_frac,'Cocco_frac',3,                        &
                          long_name='Coccolithophore Fraction',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_Eukar_frac,'Eukar_frac',3,                        &
                          long_name='Eukaryote Fraction',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_diatS,'diatS',3,                        &
                          long_name='Diatom Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_diatN,'diatN',3,                        &
                          long_name='Diatom Nitrogen Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_phytoN,'phytoN',3,                        &
                          long_name='Total Phyto Nitrogen Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_coccoS,'coccoS',3,                        &
                          long_name='Coccolithophore Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_cyanoS,'cyanoS',3,                        &
                          long_name='Cyanobacteria Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_eukarS,'eukarS',3,                        &
                          long_name='Eukaryote Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_diazS,'diazS',3,                        &
                          long_name='Diazotroph Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_phaeoS,'phaeoS',3,                        &
                          long_name='Phaeocystis Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_phaeoN,'phaeoN',3,                        &
                          long_name='Phaeocystis Nitrogen Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_phaeonS,'phaeonS',3,                        &
                          long_name='Phaeocystis NH Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_phaeonN,'phaeonN',3,                        &
                          long_name='Phaeocystis NH Nitrogen Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

! (swang) to remove output of phaeoC generated by tracegas module
!   call define_tavg_field(tavg_phaeoC,'phaeoC',3,                        &
!                         long_name='Phaeocyctis Carbon Concentration',      &
!                          units='mmol/m^3/s', grid_loc='3111',           &
!                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_zooS,'zooS',3,                        &
                          long_name='Zooplankton Sulfur Concentration',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_zooCC,'zooCC',3,                        &
                          long_name='Zooplankton Carbon Concentration in tracegas',      &
                          units='mmol/m^3/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_RSNzoo,'RSNzoo',3,                        &
                          long_name='Zooplankton Sulfur-Nitrogen Redfield Ratio',      &
                          units='none', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!
   allocate(TRACEGAS_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   TRACEGAS_SFLUX_TAVG = c0
!
!-----------------------------------------------------------------------

!EOC

 end subroutine tracegas_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: tracegas_set_interior
! !INTERFACE:

 subroutine tracegas_set_interior(k, TEMP_OLD, TEMP_CUR, &
    TRACER_ECO_OLD, TRACER_ECO_CUR, &
    TRACER_MODULE_OLD, TRACER_MODULE_CUR, DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  Compute time derivatives for tracegas state variables
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! vertical level index

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      TEMP_OLD,          &! old potential temperature (C)
      TEMP_CUR            ! current potential temperature (C)

  !real (r8), dimension(nx_block,ny_block,km,tracegas_tracer_cnt), intent(in) :: &
   real (r8), dimension(:,:,:,:), intent(in) :: &
      TRACER_ECO_OLD, &! old tracer values
      TRACER_ECO_CUR,  &! current tracer values
      TRACER_MODULE_OLD, &! old tracer values
      TRACER_MODULE_CUR   ! current tracer values

   type (block), intent(in) :: &
      this_block          ! block info for the current block

! !OUTPUT PARAMETERS:

!  real (r8), dimension(nx_block,ny_block,tracegas_tracer_cnt), intent(out) :: &
   real (r8), dimension(:,:,:), intent(out) :: &
      DTRACER_MODULE      ! computed source/sink terms

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    character(*), parameter :: &
         subname = 'tracegas_mod:tracegas_set_interior'

    !---------------------------------------------------------------------------
    !   begin with those recreated from ecosystem module
    !---------------------------------------------------------------------------

    real(r8), parameter :: &
         epsC      = 1.00e-8 ! small C concentration (mmol C/m^3)

!(swang:not used!)
!    real(r8), dimension(nx_block,ny_block) :: &
!         UV_in,          & ! generic UV radiation (W/m^2)
!         KUVdz,          & ! UV adsorption coefficient (non-dim)
!         UV_avg            ! average UV (W/m^2)

    real(r8), dimension(nx_block,ny_block) :: & !mmol base element/m^3, Moore et al. 2002
!        NO3_loc,        & ! local copy of model NO3, can determine cyano fraction (swang: not used)
!        DOC_loc,        & ! local copy of model DOC, source for OCS, NMHC, others
         spC_loc,        & ! local copy of model spC, emission of several gases
         spChl_loc,      & ! local copy of model spChl, Beers Law and stress calculations.
         spCaCO3_loc,    & ! local copy of model spCaCO3, characterize the prym.
         diatC_loc,      & ! local copy of model diatC, bloom specialists
!        diatChl_loc,    & ! local copy of model diatChl, need for Beers Law (swang: not used)
         phaeoC_loc,     & ! local copy of model phaeoC, bloom specialists
         phaeonC_loc,    & ! local copy of model phaeonC, bloom specialists (for N.H.)
         zooC_loc,       & ! local copy of model zooC, relate to DMSP and others`
         diazC_loc         ! local copy of model diazC, gyre specialists
!        diazChl_loc       ! local copy of model diazChl, need for Beers Law(swang: not used)

    real(r8), dimension(nx_block,ny_block) :: &
         work              ! temporary for summed quantities to be averaged

    real(r8) :: &
         R                 ! N/C mole ratio for all phyto and zoopl

    real(r8), dimension(nx_block,ny_block) :: &
         Fcocco            ! Fraction of small carbon associated with calcite

    !--------------------------------------------------------------------------------
    !   Conversion to nitrogen currency as the base for initial sulfur simulations.
    !   Links will most often be made to models such as Vallina et al. 2008,
    !   which tend to track organisms as N.
    !   Thus the CCSM standard carbon quantities are transposed.
    !   The trace gas system really begins here.
    !   Note that the DML phyto-classes are often subdivided,
    !   e.g. into specialists such as phaeocystis or nonS producers cyanobacteria.
    !   Ordering of the organisms in the code is intended to reflect
    !   prioritization of this decomposition process.
    !   E.g. coccolithophorids are partitioned prior to the cyano
    !   but after phaeocystis due to its traditional local dominance.
    !   Diatom, diazotroph and zooplanktonic bins are already sulfur appropriate
    !   and so are not directly altered here.
    !   Heterotrophic, recycling bacteria will be decoupled from phaeocystis
    !   which is colonial and generates antibiotics.
    !--------------------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block) :: & !all mmol N/m^3
         diatN_loc,       & ! analog to diatC_loc
         phaeoN_loc,      & ! nitrogen associated with phaeocystis
         phaeonN_loc,     & ! nitrogen associated with phaeocystis NH
         coccoN_loc,      & ! nitrogen associated with coccolithophores
         cyanoN_loc,      & ! nitrogen associated with cyanobacteria
         eukarN_loc,      & ! analog to spC_loc but remove specialists
         diazN_loc,       & ! analog to diazC_loc
         phytoN_loc,      & ! total phytoplanktonic nitrogen excluding phaeo
         zooN_loc           ! analog to zooC_loc

    !--------------------------------------------------------------------------------
    !   Further conversion to sulfur quantities.
    !   These may be thought of as particulate DMSP by class
    !   and so may be zero as for example in the case of prokaryotes.
    !   A simple renaming is included at the end of the list
    !   to align with DMS community convention.
    !   Phaeocystis is pulled out of the total again,
    !   but this time because it is not subject to grazing.
    !   An even, steady state leakage from colonies is assumed.
    !--------------------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block) :: & !all mmol S/m^3
         diatS_loc,       & ! DMSP in diatoms
         phaeoS_loc,      & ! DMSP in phaeocystis
         phaeonS_loc,     & ! DMSP in phaeocystis NH
         coccoS_loc,      & ! DMSP in coccolithophores
         cyanoS_loc,      & ! DMSP in cyanobacteria
         eukarS_loc,      & ! DMSP in smalls remaining
         diazS_loc,       & ! DMSP in diazotrophs
         phytoS_loc,      & ! total phytoplanktonic sulfur excluding phaeo
         zooS_loc           ! total zooplanktonic sulfur

    real(r8), dimension(nx_block,ny_block) :: & !all mmol S/m^3
         DMSPp_loc          ! renaming total per sulfur community convention

    !--------------------------------------------------------------------------------
    !   The mechanism traces to reviews such as Kiene et al. 2000 and Simo 2004
    !   supplemented by detail taken from low dimensionality models,
    !   examples being Lefevre et al. 2002, Vallina et al. 2008, Toole et al. 2008.
    !   Best descriptions to date appear in our reports on piston velocity testing.
    !   Channels are similar but not identical to those incorporated into
    !   the first (SciDAC) coupled ocean-atmosphere chemistry runs in CCSM.
    !   Key points include:
    !   Phytoplanktonic intracellular contents consistent with Stefels 2000.
    !   Zooplanktonic dependence built into the release constant k
    !   which relaxes coupling to chlorophyll and works toward a summer phase lag.
    !   Sunda et al. (2002) cell internal oxidant stress under high irradiance
    !   treated ad hoc by enhancing picoeukaryotic cell content
    !   in proportion to an arbitrary chlorophyll decrement.
    !   This situation arose because ultraviolet calculations were running late
    !   relative to the recent CODiM intercomparison.
    !   Real UV penetration will soon be computed per Chu et al. CO simulations.
    !   Direct summer exudation remains a speculative explanation for summer peaks.
    !   It may be thought of as closely related to the decrement approach
    !   but is not actually simulated here in a direct sense.
    !   DMS yield from bacterial processing of DMSP is given a T dependence.
    !   The concept is that microbial sulfur demand will be higher
    !   in cold, nutrient rich waters. See Kiene et al. 2000 for concepts involved.
    !   True sulfur utilization must clearly be incorporated very soon.
    !   This can begin with application of a cycling time to diagnosed densities,
    !   but may eventually entail addition of a dynamic bacterial module.
    !   Microbial consumption is rendered 2nd order per the above reviews,
    !   but in particular based on the equatorial data of Kiene and Bates 1990.
    !   Densities depend on free (non-colonial) phytoplankton distributions.
    !   Injection scale is a dial for dealing simultaneously with uncertainties in
    !   disruption rate and intracellular DMSP content.
    !   Optimization issues are thus focused upon average release.
    !   This collapses their dimensionality considerably.
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   DMS and DMSP are the dissolved forms carried as tracers by POP.
    !   The community would refer to them as DMS/Pd.
    !   The analytically particulate form DMSPp is for the moment
    !   defined as computationally local in order to minimize expense.
    !   Averaging and output will be added as necessary.
    !--------------------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block) :: & !mmol molecule/m^3
         DMS_loc,         & ! local copy of model DMS
         DMSP_loc           ! local copy of model DMSP

    !--------------------------------------------------------------------------------
    !   Overall uncertainties are enfolded in the injection scaling factor.
    !   The analyst is then free to morph the shape of the DMS distribution
    !   through mechanism aspects with better understood variation.
    !   Examples include the yield and uptake dynamics
    !   plus variability in cell internal sulfur metabolism whether these
    !   happen to be inter- or intraspecies.
    !--------------------------------------------------------------------------------

    real(r8):: &
         k_S_p_base,     & ! baseline phyto release to be augmented by grazing (1/sec)
         zooC_avg,       & ! median zooplankton mass (mmol C/m^3)
         mort,           & ! background of mortality and nongrazing (non-dim)
         k_conv,         & ! 1st order constant for DMSP conversion (1/sec)
         k_S_z,          & ! 1st order constant for release from zooplankton (1/sec)
         B_preexp,       & ! preexponential for diagnosed bacterial density (mmol N/m^3)
         B_exp,          & ! exponent for diagnosed bacterial density (non-dim)
         k_S_B,          & ! 2nd order constant for bacterial consumption (1/sec/B)
         k_bkgnd,        & ! Sulfur removal at low level in thermocline (1/sec)
         j_dms_perI        ! coefficient for photolytic removal of DMS (1/sec/W/m^2)

    real(r8):: &
         inject_scale      ! S release and content semi-arbitrary so can be scaled

    !--------------------------------------------------------------------------------
    !   These kinetically important quantities have geographic dependence.
    !--------------------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block) :: &
         k_S_p,        & ! 1st order constant for DMSP release from phyto (1/sec)
         yield,        & ! fraction of DMSP converion to DMS
         B_diagnosed,  & ! local bacterial density (mmol N/m^3)
         j_dms           ! overall dms photolysis rate (1/sec)

    !--------------------------------------------------------------------------------
    !   Phaeocystis are very high in cell internal reduced sulfur content.
    !   Cyanobacteria by contrast are low.
    !   If DMS patterns are to be realistically simulated,
    !   the two classes are must be segregated relative to the incoming ecosystem.
    !   Further, bacterial yield effects will be accounted at latitudes 30-60
    !   along with upregulation stress in the gyres.
    !   Fuzzy logic is now applied to achieve this fractionation/regulation.
    !   It should of course be relaxed in favor of real ecology in the medium term.
    !   Temperature or concentration references have standard model units.
    !   All indices and membership ramps are dimensionless.
    !   Phaeocystis is handled not via a ramp but rather as a temperature switch.
    !   Witholding from general phytoplankton totals due to colonization
    !   is simulated in some module versions.
    !   Presently a single T band handles both speciation and yield functions.
    !   Needs for membership overlap and blending are obviated somewhat.
    !   Stratification stress is signaled by low DML picoplanktonic densities.
    !   Hence the symbol "sp" for small phytoplankton reappears.
    !   This method was criticized at CODiM as redundant-
    !   chlorophyll predicts production etc.
    !   Indeed it must soon be modernized to reflect the true photobiology.
    !--------------------------------------------------------------------------------

    real(r8) :: &          ! (C)
! applying explicit phaeocystis. No need for a temperautre switch (swang)
!         T_phaeo           ! phaeocystis triggered at a single temperature
! Temperature limits for cryoprotection effect (S content increase at cold temperature)
		 T_cryo_hi,		 & ! upper temperature limit
		 T_cryo_lo		   ! lower temperature limit
		 
    real(r8) :: &          ! (C)
         T_lo,           & ! lower reference temperature
         T_hi              ! upper reference temperature

    real(r8), dimension(nx_block,ny_block) :: &
         T_ind             ! upward temperature index

    real(r8) :: &
         Min_cyano_frac,   & ! lower limit to fraction of smalls as cyanos
         Max_cyano_frac      ! upper limit to fraction of smalls as cyanos

    real(r8), dimension(nx_block,ny_block) :: &
!         Phaeo_frac,       & ! local fraction of smalls as phaeocystis		(swang: not used)
         Cocco_frac,       & ! local fraction of smalls as coccos
         Cyano_frac,       & ! local fraction of smalls as cyanos
         Eukar_frac          ! local fraction of remaining smalls

    real(r8) :: &
         Min_yld,        & ! lower limit on bacterial sulfur demand
         Max_yld           ! upper limit on bacterial sulfur demand

    real(r8) :: &          ! (mg/m^3)
         Sp_ref            ! reference chlorophyll for decrement calculation

    real(r8), dimension(nx_block,ny_block) :: &
         Sp_dec            ! small phytoplanktonic decrement

    real(r8) :: &
         Stress_mult       ! slope of the chlorophyll ramps

    real(r8), dimension(nx_block,ny_block) :: &
         Stress_fac        ! local up regulation

    !--------------------------------------------------------------------------------
    !   Surface temperature fixed moving down the column
    !   for several purposes within the fuzzy and binary logic 
    !--------------------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block) :: &
         TEMP_loc           ! local surface temperature (C)

    !--------------------------------------------------------------------------------
    !   Sulfur to nitrogen ratios will be set per Stefels 2000 and related,
    !   with local adjustments under the present fuzzy logic.
    !--------------------------------------------------------------------------------

    real(r8) :: &
         Rs2n_diat,    & ! S/N mole ratio for diatoms
         Rs2n_phaeo,   & ! S/N mole ratio for phaeocystis
         Rs2n_phaeon,  & ! S/N mole ratio for phaeocystis NH
         Rs2n_cocco,   & ! S/N mole ratio for coccolithophorids
         Rs2n_cyano,   & ! S/N mole ratio for cyanobacteria
         Rs2n_eukar,   & ! S/N mole ratio for remaining smalls
         Rs2n_diaz       ! S/N mole ratio for diazotrophs

    !--------------------------------------------------------------------------------
    !   Average cell contents have geographic dependence for zooplantkon
    !   because they are weighted averages over food items.
    !--------------------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block) :: &
         Rs2n_zoo        ! S/N assuming weighted average of consumable content

    !--------------------------------------------------------------------------------
    !   Begin declaration of source sink terms.
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   By and large the structure of the mechanism is reflected.
    !   For example in the present configuration,
    !   most sulfur release from phytoplankton takes the form DMSP.
    !   In anticipation of increased mechanistic complexity
    !   certain likely future permutations are included, e.g. exudation.
    !   The concepts here are clarity and motivation.
    !   The ports are simply readied for later connections.
    !   A few omissions will be obvious.
    !   There are to date no indications that DMSP photolyzes in the column.
    !   Note that sea-air transfer is time split into other routines.
    !--------------------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block) :: &
         dms_s_exu,    & ! DMS source from exudation (mmol S/m^3/sec)
         dms_s_dmsp,   & ! DMS source by conversion of DMSP (mmol S/m^3/sec)
         dms_s           ! DMS source total (mmol S/m^3/sec)

    real(r8), dimension(nx_block,ny_block) :: &
         dms_r_B,      & ! DMS removal by bacteria (mmol S/m^3/sec)
         dms_r_phot,   & ! DMS removal by photolysis (mmol S/m^3/sec)
         dms_r_bkgnd,  & ! DMS removal low level in thermocline
         dms_r           ! DMS removal total (mmol S/m^3/sec)

    real(r8), dimension(nx_block,ny_block) :: &
         dmsp_s_phaeo, & ! DMSP source from phaeo (mmol S/m^3/sec)
         dmsp_s_nonphaeo,   & ! DMSP source from other phytoplankton (mmol S/m^3/sec)
         dmsp_s_zoo,   & ! DMSP source from zooplankton (mmol S/m^3/sec)
         dmsp_s          ! DMSP source total (mmol S/m^3/sec)

   real(r8), dimension(nx_block,ny_block) :: &
         dmsp_r_B,       & ! DMSP removal by bacteria (mmol S/m^3/sec)
         dmsp_r_bkgnd,   & ! DMSP removal low level in thermocline
         dmsp_r            ! DMSP removal total (mmol S/m^3/sec)

   integer (int_kind) :: &
      n, auto_ind,       &
      bid                 ! local_block id

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

   bid = this_block%local_id

   call timer_start(tracegas_interior_timer, block_id=bid)

   DTRACER_MODULE = c0

!-----------------------------------------------------------------------
!  exit immediately if computations are not to be performed
!-----------------------------------------------------------------------

!maltrud debug
!if (my_task == master_task) write(stdout,*)k,' starting tracegas_set_interior'

   if (.not. lsource_sink) then
      call timer_stop(tracegas_interior_timer, block_id=bid)
      return
   endif

    !---------------------------------------------------------------------------
    !   create local copies of requisite ecotracers
    !   treat negative values as zero and apply mask to locals
    !---------------------------------------------------------------------------
!(swang: NO3_loc not used)
!   NO3_loc      = max(c0, p5*(TRACER_ECO_OLD(:,:,k,no3_ind) + &
!                              TRACER_ECO_CUR(:,:,k,no3_ind)))
!   DOC_loc      = max(c0, p5*(TRACER_ECO_OLD(:,:,k,doc_ind) + &
!                              TRACER_ECO_CUR(:,:,k,doc_ind)))
   zooC_loc     = max(c0, p5*(TRACER_ECO_OLD(:,:,k,zooC_ind) + &
                              TRACER_ECO_CUR(:,:,k,zooC_ind)))

!maltrud debug
!BUG -- should have used TRACER_ECO_* for below (fixed - swang)
  auto_ind = 1 ! small phyto
     n = autotrophs(auto_ind)%Chl_ind
     spChl_loc = max(c0, &
        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

     n = autotrophs(auto_ind)%C_ind
     spC_loc = max(c0, &
        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

     n = autotrophs(auto_ind)%CaCO3_ind
     if (n > 0) then
        spCaCO3_loc = max(c0, &
           p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))
     endif

  auto_ind = 2 ! diatoms
!     n = autotrophs(auto_ind)%Chl_ind
!     diatChl_loc = max(c0, &
!        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

     n = autotrophs(auto_ind)%C_ind
     diatC_loc = max(c0, &
        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

  auto_ind = 3 ! diazotrophs
!     n = autotrophs(auto_ind)%Chl_ind
!     diazChl_loc = max(c0, &
!        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

     n = autotrophs(auto_ind)%C_ind
     diazC_loc = max(c0, &
        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

  auto_ind = 4 ! phaeocystis
     n = autotrophs(auto_ind)%C_ind
     phaeoC_loc = max(c0, &
        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

  auto_ind = 5 ! phaeocystis NH
     n = autotrophs(auto_ind)%C_ind
     phaeonC_loc = max(c0, &
        p5*(TRACER_ECO_OLD(:,:,k,n) + TRACER_ECO_CUR(:,:,k,n)))

    where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid))
!       DOC_loc      = c0
       spC_loc      = c0
       spChl_loc    = c0
       spCaCO3_loc  = c0 
       diatC_loc    = c0 
!       diatChl_loc  = c0 
       phaeoC_loc   = c0 
       phaeonC_loc  = c0 
       zooC_loc     = c0 
       diazC_loc    = c0
!       diazChl_loc  = c0
    end where

    !---------------------------------------------------------------------------
    !   create local copies of trace gas fields in direct analogy with above
    !   treat negative values as zero and apply mask to locals
    !---------------------------------------------------------------------------

    DMS_loc      = MAX(c0, p5*(TRACER_MODULE_OLD(:,:,k,dms_ind) + &
                               TRACER_MODULE_CUR(:,:,k,dms_ind)))
    DMSP_loc     = MAX(c0, p5*(TRACER_MODULE_OLD(:,:,k,dmsp_ind) + &
                               TRACER_MODULE_CUR(:,:,k,dmsp_ind)))

    where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid))
       DMS_loc     = c0
       DMSP_loc    = c0
    end where

    !---------------------------------------------------------------------------
    !   Set parameter values not included in namelists.
    !   Conversion of time from days to seconds takes place when necessary,
    !   to maintain consistency with DML ecosystem settings.
    !   Zooplankton reference level adapted from Moore et al. 2002
    !   and general experience with the DML model.
    !   Since it is used as a normalization the currency is irrelevant.
    !   Mortality is zeroed for the present in order to maximize decoupling
    !   which can be derived from grazing.
    !   The absolute rate of DMSP loss is unimportant in setting DMS levels
    !   since they will be determined by the overall sulfur flow and yield.
    !   In other words, DMSP is a relatively short lived intermediate
    !   lying along a sequential kinetic equilibrium sequence.
    !   Kiene argues that it is removed quickly and further,
    !   that most measurements in fact constitute overestimates.
    !   As a compromise the conversion rate is set at a round one day.
    !   It could be demonstrated by sensitivity testing faster values
    !   that the major effect is suppression of DMSPd,
    !   for which the distributions are poorly understood in any case.
    !   Microbial N is log linear to free phyto density with a slope one half.
    !   In some versions phaeocystis is not included in this calculation.
    !   It is considered decoupled from bulk organic material.
    !   Lower case j usually indicates photolysis in the trace gas modules.
    !---------------------------------------------------------------------------

    k_S_p_base  = 0.1_r8 * dps
    zooC_avg    = 0.3_r8
    mort        = 0.0_r8
    k_conv      = 1.0_r8 * dps
    k_S_z       = 0.1_r8 * dps
    B_preexp    = 0.1_r8
    B_exp       = 0.5_r8
    k_S_B       = 30.0_r8 * dps
    k_bkgnd     = 0.01_r8 * dps
    j_dms_perI  = 0.005_r8 * dps

    !---------------------------------------------------------------------------
    !  Exploratory scaling appended for disruption rate/DMSP content.
    !---------------------------------------------------------------------------

!maltrud global avg of surface DMS from Lana is 2.35
!   inject_scale = 1.00_r8  ! globl surface avg about 4.5
    inject_scale = 2.00_r8  

    !---------------------------------------------------------------------------
    !  Baseline phytoplanktonic sulfur release rate constant is here adjusted
    !  by local zooplanktonic densities normalized to an average.
    !  This has the effect of partially decoupling sulfur from chlorophyll.
    !  Carbon is used in this case since conversions would simply cancel.
    !  Mortality modulates the decoupling but may be zeroed.
    !  Note that phaeocystis release is treated as an independent leak
    !  adjusted in fact to control grazing.
    !  This is normally given the base rate but could be as slow as senescence.
    !---------------------------------------------------------------------------

    k_S_p = k_S_p_base * (mort + (zooC_loc/0.3_r8))

    !---------------------------------------------------------------------------
    !   Fuzzy logical parameters
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !   The fuzzy logic style settings may be readily altered.
    !   Proportion of low S prokaryote constrained by Gregg et al. 2003 review.
    !   Stress functions are designed to align gyres with the Kettle data base.
    !   They may be compensating in part for a lack of exudation.
    !   Phaeocystis is currently a hard switch.
    !----------------------------------------------------------------------------

    T_cryo_hi      =  1.0_r8
    T_cryo_lo      =  -1.0_r8
    T_lo           =  15.0_r8
    T_hi           =  20.0_r8
    Min_cyano_frac =  c0
    Max_cyano_frac =  0.5_r8
    Min_yld        =  0.2_r8
    Max_yld        =  0.7_r8
    Sp_ref         =  0.1_r8
    Stress_mult    = 10.0_r8

    !---------------------------------------------------------------------------
    !   Organism compositions
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !   N/C ratio taken directly from Moore et al. 2002.
    !   S/N usually rounded to powers of ten from classic Keller et al. 1989
    !   and are consistent with the intermediate review Stefels 2000.
    !   Phaeocystis is exceptional based on Liss et al. 1994 and Schoemann et al. 2005.
    !   Note that with the exception of gyre stresses,
    !   ratios  between the intracellular relationships remain fixed.
    !   The scale quantity thus operates in part to raise/lower them together.
    !---------------------------------------------------------------------------

    R             = 0.137_r8
    Rs2n_diat     = 0.01_r8
    Rs2n_phaeo    = 0.3_r8
    Rs2n_phaeon   = 0.3_r8
    Rs2n_cocco    = 0.1_r8
    Rs2n_cyano    = c0
    Rs2n_eukar    = 0.1_r8
    Rs2n_diaz     = c0

   TEMP_loc = p5*(TEMP_OLD + TEMP_CUR)

!(swang: not used!)
   !UV_in  = UV_out(:,:,bid)

   !where (.not. LAND_MASK(:,:,bid) .or. k > KMT(:,:,bid)) UV_in = c0
    
   !if (partial_bottom_cells) then
   !   KUVdz = merge((0.01e-2_r8 * DOC_loc+0.04e-4_r8)*DZT(:,:,k,bid), c0, k<=KMT(:,:,bid))
   !else
   !   KUVdz = merge((0.01e-2_r8 * DOC_loc + 0.04e-4_r8) * dz(k), c0, k<=KMT(:,:,bid))
   !endif

   !UV_out(:,:,bid) = UV_in * EXP(-KUVdz)
   !UV_avg = UV_in * (c1 - EXP(-KUVdz)) / KUVdz
    
    !---------------------------------------------------------------------------
    !   Momentarily we follow Gabric et al. 1993 who speculate implicitly that
    !   DMS photolytic wavelengths attenuate as PAR.
    !   This is confirmed in part by the review Mopper and Kieber 2002,
    !   But modernization will be required per Toole et al. 2004 and related.
    !   Photolysis per unit intensity is based by Gabric on
    !   Brimblecombe and Shooter 1986 results for wavelengths unknown.
    !---------------------------------------------------------------------------

    j_dms = j_dms_perI * PAR_avg(:,:,bid)

    !---------------------------------------------------------------------------
    !   Fcocco is the fraction of sp organic matter in coccolithophores
    !   as imported from the driver DML ecology.
    !   In the S model this takes precedence over all but phaeocystis.
    !---------------------------------------------------------------------------

    Fcocco = spCaCO3_loc / (spC_loc + epsC)
    where (Fcocco > 0.4_r8) Fcocco = 0.4_r8

    Cocco_frac = Fcocco

    !--------------------------------------------------------------------------
    !   Fuzzy logical segregation into relevant subclasses
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !   Further distribution of the DML small phytoplanktonic biomass
    !   is now undertaken per requirements of the sulfur cycle.
    !   Phaeocystis is considered dominant and ultimately trumps other autotrophs.
    !--------------------------------------------------------------------------

    T_ind  = (TEMP_loc - T_lo)/(T_hi - T_lo)

    where (T_ind <= c0) T_ind = c0
    where (T_ind >= c1) T_ind = c1

    Cyano_frac  = (T_ind * (Max_cyano_frac - Min_cyano_frac)) + Min_cyano_frac
    Cyano_frac  = (c1 - Cocco_frac) * Cyano_frac

    Eukar_frac = c1 - Cocco_frac - Cyano_frac


    !--------------------------------------------------------------------------
    !  Convert to nitrogen currency plus distribute into new classes
    !--------------------------------------------------------------------------

!maltrud debug
! diatC_loc = 1.e-5
! spC_loc   = 1.e-5
! diazC_loc = 1.e-5
    diatN_loc    = R * diatC_loc
    phaeoN_loc   = R * phaeoC_loc
    phaeonN_loc  = R * phaeonC_loc
    coccoN_loc   = Cocco_frac * R * spC_loc
    cyanoN_loc   = Cyano_frac * R * spC_loc
    eukarN_loc   = Eukar_frac * R * spC_loc
    diazN_loc    = R*diazC_loc
    zooN_loc     = R*zooC_loc

    !--------------------------------------------------------------------------
    !  Collect noncolonial nitrogen
    !  (swang) since phaeocystis is grazed in BEC, include it in grazing-dmsp production too
    !--------------------------------------------------------------------------

    phytoN_loc = diatN_loc + coccoN_loc + cyanoN_loc + eukarN_loc + diazN_loc  &
    			+ phaeoN_loc + phaeonN_loc

    !--------------------------------------------------------------------------
    !  Gyre oxidant stress now exerts upregulation per chlorophyll decrement.
    !  Functional form maintains round figure parameters and powers,
    !  but was essentially determined by offline, postprocessing trial and error.
    !  Real physiology must be incorporated in short order.
    !--------------------------------------------------------------------------

    Sp_dec = (Sp_ref - spChl_loc)/Sp_ref

    where (Sp_dec <= c0) Sp_dec = c0
    where (Sp_dec >= c1) Sp_dec = c1

    Stress_fac = c1 + Stress_mult* Sp_dec * Sp_dec

    where (Stress_fac >= 10.0_r8) Stress_fac = 10.0_r8

    !--------------------------------------------------------------------------
    !  Yields for bacterial conversion of DMSP also determined by fuzzy logic.
    !  Since Phaeocystis is colonial and produces antimicrobials,
    !  its habitat constitutes an exception at very low S demand.
    !--------------------------------------------------------------------------

    yield = (T_ind * (Max_yld - Min_yld)) + Min_yld

!	(swang) assume an optimal temperautre for cryoprotention 
    where (TEMP_loc < T_cryo_hi .and. TEMP_loc > T_cryo_lo) yield = 0.5_r8 
    where (TEMP_loc < -1.0_r8) yield = 0.25_r8 

    !--------------------------------------------------------------------------
    !  Phytoplanktonic sulfur content determined.
    !  This is where oxidant stress is in fact applied.
    !--------------------------------------------------------------------------

    diatS_loc  = Rs2n_diat  * diatN_loc
    phaeoS_loc = Rs2n_phaeo * phaeoN_loc
    phaeonS_loc= Rs2n_phaeon * phaeonN_loc
    coccoS_loc = Rs2n_cocco * coccoN_loc
    cyanoS_loc = Rs2n_cyano * cyanoN_loc
    eukarS_loc = Rs2n_eukar * eukarN_loc * Stress_fac
    diazS_loc  = Rs2n_diaz   * diazN_loc

    !--------------------------------------------------------------------------
    !  Collect noncolonials
    ! (swang) assume only a fraction (35%) of phaeo S contribute to dsmp release when grazed
    !--------------------------------------------------------------------------

    phytoS_loc = diatS_loc + coccoS_loc + cyanoS_loc + eukarS_loc + diazS_loc &
                + 0.4_r8 * (phaeoS_loc + phaeonS_loc) 

    !--------------------------------------------------------------------------
    !  Weight the zooplanktonic reduced sulfur content per the various
    !  phytoplankton types available unprotected in the column as food.
    !  Observe that plant nitrogen may total zero under some circumstances
    !  and in this case it is assumed that all concentrations are equally small.
    !  Note that phaeocystis has been segregated entirely.
    !  In colonial form it exhibits multiple grazing inhibition strategies.
    !--------------------------------------------------------------------------

    where (phytoN_loc > c0)
    Rs2n_zoo = (Rs2n_diat   * diatN_loc  + &
    			 0.4_r8 * Rs2n_phaeo   * phaeoN_loc  + &
    			 0.4_r8 * Rs2n_phaeon   * phaeonN_loc  + &
                 Rs2n_cocco * coccoN_loc + &
                 Rs2n_cyano * cyanoN_loc + &
                 Rs2n_eukar * eukarN_loc * Stress_fac + &
                 Rs2n_diaz   * diazN_loc)/phytoN_loc
    elsewhere
! (swang) since phaeo and phaeon in 2 hemi., considered as 1 group here for averaging
    Rs2n_zoo = (Rs2n_diat + Rs2n_cocco + Rs2n_cyano + Rs2n_eukar + Rs2n_diaz &
    			+ Rs2n_phaeo)/6.0
    end where

    zooS_loc = Rs2n_zoo * zooN_loc

    !---------------------------------------------------------------------------
    !   Given a nitrogen distribution diagnose B which is bacterial N.
    !   Phaeocystis excluded because it is colonial and exudes antibiotics.
    !   Given a turnover time B could serve as the basis for S demand calcs.
    !   However our tendency is to move directly to bacterial populations.
    !   Second order DMS removal is based on reports of biotic uptake presented
    !   by Kiene and Bates 1991 for the eastern equatorial Pacific.
    !---------------------------------------------------------------------------

    B_diagnosed = B_preexp*(phytoN_loc**B_exp)

    !-------------------------------------------------------------------------
    !   Construction of kinetic terms for the sulfur cycle
    !-------------------------------------------------------------------------

    dms_s_dmsp   = yield * k_conv * DMSP_loc
    dms_s        = dms_s_dmsp

    dms_r_B      = k_S_B * B_diagnosed * DMS_loc
    dms_r_phot   = j_dms * DMS_loc
    dms_r_bkgnd  = k_bkgnd * DMS_loc
    dms_r        = dms_r_B + dms_r_phot + dms_r_bkgnd

    dmsp_s_phaeo = inject_scale * k_S_p_base * (phaeoS_loc + phaeonS_loc)
    dmsp_s_nonphaeo   = inject_scale * k_S_p * phytoS_loc
    dmsp_s_zoo   = inject_scale * k_S_z * zooS_loc
    dmsp_s       = dmsp_s_phaeo + dmsp_s_nonphaeo + dmsp_s_zoo

    dmsp_r_B     = k_conv * DMSP_loc
    dmsp_r_bkgnd = k_bkgnd * DMSP_loc
    dmsp_r       = dmsp_r_B + dmsp_r_bkgnd

    DTRACER_MODULE(:,:,dms_ind)  = dms_s - dms_r
    DTRACER_MODULE(:,:,dmsp_ind) = dmsp_s - dmsp_r

! DMS source terms
       call accumulate_tavg_field(dms_s_dmsp,tavg_DMS_S_DMSP,bid,k)
       call accumulate_tavg_field(dms_s,tavg_DMS_S_TOTAL,bid,k)

! DMS removal terms
       call accumulate_tavg_field(dms_r_B,tavg_DMS_R_B,bid,k)
       call accumulate_tavg_field(dms_r_phot,tavg_DMS_R_PHOT,bid,k)
       call accumulate_tavg_field(dms_r_bkgnd,tavg_DMS_R_BKGND,bid,k)
       call accumulate_tavg_field(dms_r,tavg_DMS_R_TOTAL,bid,k)

! DMSP source terms
       call accumulate_tavg_field(dmsp_s_phaeo,tavg_DMSP_S_PHAEO,bid,k)
       call accumulate_tavg_field(dmsp_s_nonphaeo,tavg_DMSP_S_NONPHAEO,bid,k)
       call accumulate_tavg_field(dmsp_s_zoo,tavg_DMSP_S_ZOO,bid,k)
       call accumulate_tavg_field(dmsp_s,tavg_DMSP_S_TOTAL,bid,k)

! DMSP removal terms
       call accumulate_tavg_field(dmsp_r_B,tavg_DMSP_R_B,bid,k)
       call accumulate_tavg_field(dmsp_r_bkgnd,tavg_DMSP_R_BKGND,bid,k)
       call accumulate_tavg_field(dmsp_r,tavg_DMSP_R_TOTAL,bid,k)

! fractional compositions
       call accumulate_tavg_field(Cyano_frac,tavg_Cyano_frac,bid,k)
       call accumulate_tavg_field(Cocco_frac,tavg_Cocco_frac,bid,k)
       call accumulate_tavg_field(Eukar_frac,tavg_Eukar_frac,bid,k)

! sulfur content
       call accumulate_tavg_field(diatS_loc,tavg_diatS,bid,k)
       call accumulate_tavg_field(diatN_loc,tavg_diatN,bid,k)
       call accumulate_tavg_field(phytoN_loc,tavg_phytoN,bid,k)
       call accumulate_tavg_field(coccoS_loc,tavg_coccoS,bid,k)
       call accumulate_tavg_field(cyanoS_loc,tavg_cyanoS,bid,k)
       call accumulate_tavg_field(eukarS_loc,tavg_eukarS,bid,k)
       call accumulate_tavg_field(diazS_loc,tavg_diazS,bid,k)
       call accumulate_tavg_field(phaeoS_loc,tavg_phaeoS,bid,k)
       call accumulate_tavg_field(phaeonS_loc,tavg_phaeonS,bid,k)
       call accumulate_tavg_field(zooS_loc,tavg_zooS,bid,k)


! other
       call accumulate_tavg_field(zooC_loc,tavg_zooCC,bid,k)
       call accumulate_tavg_field(Rs2n_zoo,tavg_RSNzoo,bid,k)

   call timer_stop(tracegas_interior_timer, block_id=bid)

!-----------------------------------------------------------------------
!EOC

 end subroutine tracegas_set_interior

!***********************************************************************
!BOP
! !IROUTINE: tracegas_init_sflux
! !INTERFACE:

 subroutine tracegas_init_sflux

! !DESCRIPTION:
!  Initialize surface flux computations for tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      luse_INTERP_WORK     ! does INTERP_WORK need to be allocated

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block) :: WORK

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields

!-----------------------------------------------------------------------

   luse_INTERP_WORK = .false.

!-----------------------------------------------------------------------
!  read gas flux forcing (if required)
!-----------------------------------------------------------------------

   if (lflux_gas_dms .and. &
       gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then

      luse_INTERP_WORK = .true.

!-----------------------------------------------------------------------
!  first, read ice file
!-----------------------------------------------------------------------

      allocate(fice_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(fice_file%input%filename) == 'unknown') &
         fice_file%input%filename = gas_flux_forcing_file

      call read_field(fice_file%input%file_fmt, &
                      fice_file%input%filename, &
                      fice_file%input%file_varname, &
                      WORK_READ)

      do iblock=1,nblocks_clinic
      do n=1,12
         fice_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            fice_file%DATA(:,:,iblock,1,n) = c0
         fice_file%DATA(:,:,iblock,1,n) = &
            fice_file%DATA(:,:,iblock,1,n) * fice_file%input%scale_factor
      end do
      end do

      call find_forcing_times(fice_file%data_time, &
                              fice_file%data_inc, fice_file%interp_type, &
                              fice_file%data_next, fice_file%data_time_min_loc, &
                              fice_file%data_update, fice_file%data_type)

!-----------------------------------------------------------------------
!  next, read piston velocity file
!-----------------------------------------------------------------------

      allocate(xkw_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(xkw_file%input%filename) == 'unknown') &
         xkw_file%input%filename = gas_flux_forcing_file

      call read_field(xkw_file%input%file_fmt, &
                      xkw_file%input%filename, &
                      xkw_file%input%file_varname, &
                      WORK_READ)

      do iblock=1,nblocks_clinic
      do n=1,12
         xkw_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            xkw_file%DATA(:,:,iblock,1,n) = c0
         xkw_file%DATA(:,:,iblock,1,n) = &
            xkw_file%DATA(:,:,iblock,1,n) * xkw_file%input%scale_factor
      end do
      end do

      call find_forcing_times(xkw_file%data_time, &
                              xkw_file%data_inc, xkw_file%interp_type, &
                              xkw_file%data_next, xkw_file%data_time_min_loc, &
                              xkw_file%data_update, xkw_file%data_type)

!-----------------------------------------------------------------------
!  last, read atmospheric pressure file
!-----------------------------------------------------------------------

      allocate(ap_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
      if (trim(ap_file%input%filename) == 'unknown') &
         ap_file%input%filename = gas_flux_forcing_file

      call read_field(ap_file%input%file_fmt, &
                      ap_file%input%filename, &
                      ap_file%input%file_varname, &
                      WORK_READ)

      do iblock=1,nblocks_clinic
      do n=1,12
         ap_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            ap_file%DATA(:,:,iblock,1,n) = c0
         ap_file%DATA(:,:,iblock,1,n) = &
            ap_file%DATA(:,:,iblock,1,n) * ap_file%input%scale_factor
      end do
      end do

      call find_forcing_times(ap_file%data_time, &
                              ap_file%data_inc, ap_file%interp_type, &
                              ap_file%data_next, ap_file%data_time_min_loc, &
                              ap_file%data_update, ap_file%data_type)

    endif

!-----------------------------------------------------------------------
!  allocate space for interpolate_forcing
!-----------------------------------------------------------------------

   if (luse_INTERP_WORK) &
      allocate(INTERP_WORK(nx_block,ny_block,max_blocks_clinic,1))

!-----------------------------------------------------------------------
!  register and set DMS flux
!-----------------------------------------------------------------------

   call named_field_register('SFLUX_DMS', sflux_dms_nf_ind)
   WORK = c0
   do iblock=1,nblocks_clinic
      call named_field_set(sflux_dms_nf_ind, iblock, WORK)
   end do

!-----------------------------------------------------------------------
!  verify running coupled if gas fluxes use coupler forcing
!-----------------------------------------------------------------------

   if (lflux_gas_dms .and. &
       (gas_flux_forcing_iopt == gas_flux_forcing_iopt_model) .and. &
       .not. registry_match('lcoupled')) then
      call exit_POP(sigAbort, 'tracegas_init: tracegas module requires the ' /&
                           &/ 'flux coupler when gas_flux_forcing_opt=drv')
   endif

!-----------------------------------------------------------------------
!  verify running coupled if diurnal cycle being used
!-----------------------------------------------------------------------

   if (tracegas_diurnal_cycle .and. .not. registry_match('lcoupled')) then
      call exit_POP(sigAbort, 'tracegas_init: tracegas module requires the ' /&
                           &/ 'flux coupler when tracegas_diurnal_cycle=.true.')
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine tracegas_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: tracegas_set_sflux
! !INTERFACE:

 subroutine tracegas_set_sflux(U10_SQR,IFRAC,PRESS,SST,SSS,  &
                               SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute surface fluxes for tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,&! 10m wind speed squared (cm/s)**2
      IFRAC, & ! sea ice fraction (non-dimensional)
      PRESS, & ! sea level atmospheric pressure (dyne/cm**2)
      SST,   & ! sea surface temperature (C)
      SSS      ! sea surface salinity (psu)

!  real (r8), dimension(nx_block,ny_block,tracegas_tracer_cnt,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:), &
      intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,tracegas_tracer_cnt,max_blocks_clinic), &
      intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'tracegas_mod:tracegas_set_sflux'

   integer (int_kind) :: &
      i,j, iblock     ! loop indices

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,   & ! used ice fraction (non-dimensional)
      XKW_USED,     & ! portion of piston velocity (cm/s)
      XKW,          & ! compromise between the two per H04
      AP_USED         ! used atm pressure (converted from dyne/cm**2 to atm)

   real (r8), dimension(nx_block,ny_block) :: &
      XKW_ICE,      & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      SCHMIDT_USED, & ! used Schmidt number
      PV,           & ! piston velocity (cm/s)
      FW92,         & ! apportionment to Wanninkhof 1992
      FLM86,        & ! apportionment to Liss and Merlivat 1986
      XKW_W92,      & ! the Wanninkhof limit
      XKW_LM86,     & ! the Liss and Merlivat limit
      WIND_SPEED,   & ! optional Huebert restriction
      DMSSAT_1atm,  & ! DMS saturation @ 1 atm (mmol/m^3)
      DMSSAT_USED,  & ! used DMS saturation (mmol/m^3)
      FLUX            ! tracer flux (nmol/cm^2/s)

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,          &! location and field type for ghost
      tracer_bndy_type           !    cell updates

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, WORK2 ! temporaries for averages

   real (r8) :: scalar_temp

!-----------------------------------------------------------------------
!  local parameters
!-----------------------------------------------------------------------

    real(r8), parameter :: &
         a  = 0.31_r8,     &    ! W92
         e1 = 0.17_r8,     &    ! LM86 from here 
         e2 = 2.85_r8,     &
         e3 = 0.612_r8,    &    
         e4 = 5.9_r8,      &    
         e5 = 26.79_r8,    &
         e6 = 0.612_r8           

   call timer_start(tracegas_sflux_timer)

!-----------------------------------------------------------------------

   if (check_time_flag(comp_surf_avg_flag))  &
      call comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR)

!-----------------------------------------------------------------------
!  set DMS field for passing to ice(SW) 
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,WORK1)
   do iblock = 1, nblocks_clinic

      WORK1 = max(c0,p5*(SURF_VALS_OLD(:,:,dms_ind,iblock) + &
                         SURF_VALS_CUR(:,:,dms_ind,iblock))) 

      call named_field_set(surf_dms_nf_ind, iblock, WORK1)

   enddo
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  set DMSP field for passing to ice(SW) 
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,WORK1)
   do iblock = 1, nblocks_clinic

      WORK1 = max(c0,p5*(SURF_VALS_OLD(:,:,dmsp_ind,iblock) + &
                         SURF_VALS_CUR(:,:,dmsp_ind,iblock))) 

      call named_field_set(surf_dmsp_nf_ind, iblock, WORK1)

   enddo
   !$OMP END PARALLEL DO

!---------------------------------------------------------------------------
!   Some k==1 initializations occur in this area.
!   0.45 is the fraction of incoming SW converted to PAR (non-dimensional)
!   Early models and Kieber references indicate that DMS is degraded
!   through photosentitization driven by wavelengths near 450 nm.
!   UV = 1% PAR currently represents an arbitrary portion of the B bands.
!   Indications from recent quantum yield studies for other volatiles
!   are that UVA must also be resolved in the medium term.
!---------------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
      STF_MODULE(:,:,:,iblock) = c0
! (swang: not used)
!      where (LAND_MASK(:,:,iblock))
!         UV_out(:,:,iblock) = 0.01_r8 * PAR_out(:,:,iblock)
!      elsewhere
!         UV_out(:,:,iblock)  = c0
!      end where
   enddo

!-----------------------------------------------------------------------
!  Interpolate gas flux forcing data if necessary
!-----------------------------------------------------------------------

   IFRAC_USED = c0
   XKW_USED   = c0
   AP_USED    = c0

!maltrud logic is messed for using file option
   if (lflux_gas_dms .and. &
        gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then
      if (thour00 >= fice_file%data_update) then
         tracer_data_names = fice_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Ice Fraction'
         call update_forcing_data(fice_file%data_time,      &
            fice_file%data_time_min_loc,  fice_file%interp_type,    &
            fice_file%data_next,          fice_file%data_update,    &
            fice_file%data_type,          fice_file%data_inc,       &
            fice_file%DATA(:,:,:,:,1:12), fice_file%data_renorm,    &
            tracer_data_label,            tracer_data_names,        &
            tracer_bndy_loc,              tracer_bndy_type,         &
            fice_file%filename,           fice_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
         fice_file%DATA(:,:,:,:,1:12), &
         fice_file%data_time,         fice_file%interp_type, &
         fice_file%data_time_min_loc, fice_file%interp_freq, &
         fice_file%interp_inc,        fice_file%interp_next, &
         fice_file%interp_last,       0)
      IFRAC_USED = INTERP_WORK(:,:,:,1)

      if (thour00 >= xkw_file%data_update) then
         tracer_data_names = xkw_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Piston Velocity'
         call update_forcing_data(xkw_file%data_time,      &
            xkw_file%data_time_min_loc,  xkw_file%interp_type,    &
            xkw_file%data_next,          xkw_file%data_update,    &
            xkw_file%data_type,          xkw_file%data_inc,       &
            xkw_file%DATA(:,:,:,:,1:12), xkw_file%data_renorm,    &
            tracer_data_label,           tracer_data_names,       &
            tracer_bndy_loc,             tracer_bndy_type,        &
            xkw_file%filename,           xkw_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK,     &
         xkw_file%DATA(:,:,:,:,1:12), &
         xkw_file%data_time,         xkw_file%interp_type, &
         xkw_file%data_time_min_loc, xkw_file%interp_freq, &
         xkw_file%interp_inc,        xkw_file%interp_next, &
         xkw_file%interp_last,       0)
      XKW_USED = INTERP_WORK(:,:,:,1)

      if (thour00 >= ap_file%data_update) then
         tracer_data_names = ap_file%input%file_varname
         tracer_bndy_loc   = field_loc_center
         tracer_bndy_type  = field_type_scalar
         tracer_data_label = 'Atmospheric Pressure'
         call update_forcing_data(ap_file%data_time,    &
            ap_file%data_time_min_loc,  ap_file%interp_type,    &
            ap_file%data_next,          ap_file%data_update,    &
            ap_file%data_type,          ap_file%data_inc,       &
            ap_file%DATA(:,:,:,:,1:12), ap_file%data_renorm,    &
            tracer_data_label,          tracer_data_names,      &
            tracer_bndy_loc,            tracer_bndy_type,       &
            ap_file%filename,           ap_file%input%file_fmt)
      endif
      call interpolate_forcing(INTERP_WORK, &
         ap_file%DATA(:,:,:,:,1:12), &
         ap_file%data_time,         ap_file%interp_type, &
         ap_file%data_time_min_loc, ap_file%interp_freq, &
         ap_file%interp_inc,        ap_file%interp_next, &
         ap_file%interp_last,       0)
      AP_USED = INTERP_WORK(:,:,:,1)

   endif

!-----------------------------------------------------------------------
!  calculate gas flux quantities if necessary
!-----------------------------------------------------------------------

   if (lflux_gas_dms) then

      do iblock = 1, nblocks_clinic

         if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then
            where (IFRAC_USED(:,:,iblock) < 0.2000_r8) &
               IFRAC_USED(:,:,iblock) = 0.2000_r8
            where (IFRAC_USED(:,:,iblock) > 0.9999_r8) &
               IFRAC_USED(:,:,iblock) = 0.9999_r8
            XKW(:,:,iblock) = XKW_USED(:,:,iblock)
         endif

         if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_model) then
            IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)
            where (IFRAC_USED(:,:,iblock) < c0) IFRAC_USED(:,:,iblock) = c0
            where (IFRAC_USED(:,:,iblock) > c1) IFRAC_USED(:,:,iblock) = c1
            AP_USED(:,:,iblock) = PRESS(:,:,iblock)
            XKW(:,:,iblock) = c0
         endif

!maltrud u10_sqr is cm/s, not m/s
         SCHMIDT_USED = SCHMIDT_DMS(SST(:,:,iblock), LAND_MASK(:,:,iblock))
!maltrud debug
SCHMIDT_USED = min(SCHMIDT_USED,1.e20)
SCHMIDT_USED = max(SCHMIDT_USED,1.)

!maltrud debug
!        WIND_SPEED = SQRT(U10_SQR(:,:,iblock)) * 0.01_r8
         WIND_SPEED = SQRT(abs(U10_SQR(:,:,iblock))) * 0.01_r8

         XKW_W92 = &
            a*((660.0_r8/SCHMIDT_USED)**0.500_r8)*U10_SQR(:,:,iblock)*0.0001_r8
         XKW_LM86 = &
            e2*((600.0_r8/SCHMIDT_USED)**0.500_r8)*(WIND_SPEED - 3.6_r8) &
          + e3*((600.0_r8/SCHMIDT_USED)**0.667_r8)

!maltrud original Scott coding...
!      HERE eWIND_SPEED(:,:,iblock) < 3.6_r8) XKW(:,:,iblock) = XKW_W92(:,:,iblock)
!      WHERE (eWIND_SPEED(:,:,iblock) >= 3.6_r8).and.(WIND_SPEED(:,:,iblock) < 5.6_r8))
!            FLM86(:,:,iblock) = (WIND_SPEED(:,:,iblock) - 3.6_r8)/2.0_r8
!            FW92(:,:,iblock)  = c1 - FLM86(:,:,iblock)
!            XKW(:,:,iblock) = FW92(:,:,iblock)*XKW_W92(:,:,iblock) &
!                            + FLM86(:,:,iblock)*XKW_LM86(:,:,iblock)
!      END WHERE
!      WHERE (WIND_SPEED(:,:,iblock) >= 5.6_r8) XKW(:,:,iblock) = XKW_LM86(:,:,iblock)
!maltrud now my coding...
         where (WIND_SPEED < 3.6_r8) XKW(:,:,iblock) = XKW_W92
         where ((WIND_SPEED >= 3.6_r8) .and. (WIND_SPEED < 5.6_r8))
            FLM86 = p5*(WIND_SPEED - 3.6_r8)
            FW92  = c1 - FLM86
            XKW(:,:,iblock) = FW92*XKW_W92 + FLM86*XKW_LM86
         end where
         where (WIND_SPEED >= 5.6_r8) XKW(:,:,iblock) = XKW_LM86

!maltrud what is this?
         XKW(:,:,iblock) = XKW(:,:,iblock)/3600.0_r8 ! conversion to cm/s

         TRACEGAS_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
!maltrud debug
!write(*,*)'min,max IFRAC, SFLUX = ',minval(IFRAC_USED(:,:,iblock)), maxval(IFRAC_USED(:,:,iblock)),minval(TRACEGAS_SFLUX_TAVG(:,:,1,iblock)), maxval(TRACEGAS_SFLUX_TAVG(:,:,1,iblock))

         TRACEGAS_SFLUX_TAVG(:,:,5,iblock) = SCHMIDT_USED

         TRACEGAS_SFLUX_TAVG(:,:,8,iblock) = WIND_SPEED

       !------------------------------------------------------------------------
       !   assume PRESS is in cgs units (dyne/cm**2) since that is what is
       !     required for pressure forcing in barotropic
       !   want units to be atmospheres
       !   conversion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
       !   conversion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
       !
       !   Set bad AP values to 1. This is necessary for runs restarting off
       !   a run in which the flux coupler didnot restart on AP correctly.
       !------------------------------------------------------------------------

       AP_USED(:,:,iblock) = AP_USED(:,:,iblock) / 101.325e+4_r8

       AP_USED(:,:,iblock) = merge( c1, AP_USED(:,:,iblock),   &
                                    (AP_USED(:,:,iblock) > 1.5_r8 .or.   &
                                     AP_USED(:,:,iblock) < 0.5_r8) )

       !------------------------------------------------------------------------
       !   Compute XKW_Ice. XKW is zero over land, so XKW_ICE is too.
       !------------------------------------------------------------------------

          XKW_ICE = XKW(:,:,iblock)
          where (IFRAC_USED(:,:,iblock) > 0.2_r8   &
           .and. IFRAC_USED(:,:,iblock) < 0.9999_r8)
             XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_ICE
          end where
          where (IFRAC_USED(:,:,iblock) >= 0.9999_r8)
             XKW_ICE = c0
          end where

    !---------------------------------------------------------------------------
    !   compute DMS flux
    !---------------------------------------------------------------------------

          DMSSAT_1atm = DMSSAT(SST(:,:,iblock),SSS(:,:,iblock), LAND_MASK(:,:,iblock))
          FLUX = c0
          where (LAND_MASK(:,:,iblock))
             PV = XKW_ICE 
             DMSSAT_USED = AP_USED(:,:,iblock) * DMSSAT_1atm
             FLUX = PV * (DMSSAT_USED - p5*(SURF_VALS_OLD(:,:,dms_ind,iblock) +  &
                                            SURF_VALS_CUR(:,:,dms_ind,iblock)) )
             STF_MODULE(:,:,dms_ind,iblock) = FLUX
          elsewhere
             STF_MODULE(:,:,dms_ind,iblock) = c0
          end where

          call named_field_set(sflux_dms_nf_ind, iblock, FLUX)
   
!maltrud debug
!           TRACEGAS_SFLUX_TAVG(:,:,2,iblock) = XKW_ICE
            where (LAND_MASK(:,:,iblock))
              TRACEGAS_SFLUX_TAVG(:,:,2,iblock) = c0
            elsewhere
              TRACEGAS_SFLUX_TAVG(:,:,2,iblock) = c1
            endwhere

            TRACEGAS_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)

            TRACEGAS_SFLUX_TAVG(:,:,4,iblock) = PV

            TRACEGAS_SFLUX_TAVG(:,:,6,iblock) = DMSSAT_USED

            TRACEGAS_SFLUX_TAVG(:,:,7,iblock) = p5*(SURF_VALS_OLD(:,:,dms_ind,iblock) +  &
                                                    SURF_VALS_CUR(:,:,dms_ind,iblock))

       enddo

    endif  ! lflux_gas_dms

   call timer_stop(tracegas_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine tracegas_set_sflux

!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_DMS
! !INTERFACE:

 function SCHMIDT_DMS(SST, LAND_MASK)

! !DESCRIPTION:
!---------------------------------------------------------------------------
!   Compute Schmidt number in seawater as function of SST
!   where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!   ref : Kettle and Andreae 2000
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: SST

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: SCHMIDT_DMS

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   coefficients in expansion
    !---------------------------------------------------------------------------

    real(r8), parameter :: &
         a = 2674.0_r8, &
         b = 147.12_r8, &
         c = 3.726_r8, &
         d = 0.038_r8

!-----------------------------------------------------------------------

    where (LAND_MASK)
       SCHMIDT_DMS = a + SST * (-b + SST * (c + SST * (-d)))
    elsewhere
       SCHMIDT_DMS = c0
    end where

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_DMS

!*****************************************************************************
!BOP
! !IROUTINE: DMSSAT
! !INTERFACE:

 function DMSSAT(SST, SSS, LAND_MASK)

! !DESCRIPTION:
!
!   Sat functions normally compute for a given molecule
!   a sea surface saturation concentration estimate at 1 atm total pressure
!   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
!   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
!   For DMS the assumption is made per Kettle and Andreae 2000
!   that the atmospheric concentration is negligible.
!   Henrys Law is thus not accounted here but passes for temperature and
!   salinity influences are preserved pending.
!
! !REVISION HISTORY:
!  same as module

! !USES:

    use constants, only : T0_Kelvin,radian

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST, & ! sea surface temperature (C)
      SSS    ! sea surface salinity (psu)

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK

! !OUTPUT PARAMETERS:

    real (r8), dimension(nx_block,ny_block) :: DMSSAT

!EOP
!BOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!---------------------------------------------------------------------------
!   Units should lead to mmol/m^3 for saturation level
!---------------------------------------------------------------------------

    where (LAND_MASK)
       DMSSAT = c0
    elsewhere
       DMSSAT = c0
    end where

!-----------------------------------------------------------------------
!EOC

 end function DMSSAT

!*****************************************************************************
!BOP
! !IROUTINE: tracegas_tavg_forcing
! !INTERFACE:

 subroutine tracegas_tavg_forcing(STF_MODULE)

! !DESCRIPTION:
!  Accumulate non-standard forcing related tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real (r8), dimension(nx_block,ny_block,tracegas_tracer_cnt,max_blocks_clinic), &
     intent(in) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------
!  do components of flux calculations saved in hack array
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic

!maltrud debug
!write(*,*)'min,max accum SFLUX = ',minval(TRACEGAS_SFLUX_TAVG(:,:,1,iblock)), maxval(TRACEGAS_SFLUX_TAVG(:,:,1,iblock))
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,1,iblock),  &
                                    tavg_DMS_IFRAC,iblock,1)
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,2,iblock),  &
                                    tavg_DMS_XKW,iblock,1)
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,3,iblock),  &
                                    tavg_DMS_ATM_PRESS,iblock,1)
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,4,iblock),  &
                                    tavg_DMS_PV,iblock,1)
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,5,iblock),  &
                                    tavg_DMS_SCHMIDT,iblock,1)
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,6,iblock),  &
                                    tavg_DMS_SAT,iblock,1)
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,7,iblock),  &
                                    tavg_DMS_SURF,iblock,1)
         call accumulate_tavg_field(TRACEGAS_SFLUX_TAVG(:,:,8,iblock),  &
                                    tavg_DMS_WS,iblock,1)

   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine tracegas_tavg_forcing

!*****************************************************************************
!BOP
! !IROUTINE: tracegas_write_restart
! !INTERFACE:

 subroutine tracegas_write_restart(restart_file, action)

! !DESCRIPTION:
!  write auxiliary fields & scalars to restart files
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!  NO DMS ADDITIONAL RESTART FIELDS
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!EOC

 end subroutine tracegas_write_restart

!*****************************************************************************
!BOP
! !IROUTINE: tracegas_tracer_ref_val
! !INTERFACE:

 function tracegas_tracer_ref_val(ind)

! !DESCRIPTION:
!  return reference value for tracers using virtual fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: ind

! !OUTPUT PARAMETERS:

   real (r8) :: tracegas_tracer_ref_val

!EOP
!BOC
!-----------------------------------------------------------------------

   if (vflux_flag(ind)) then
      tracegas_tracer_ref_val = surf_avg(ind)
   else
      tracegas_tracer_ref_val = c0
   endif

!-----------------------------------------------------------------------
!EOC

 end function tracegas_tracer_ref_val

!***********************************************************************

 end module tracegas_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
