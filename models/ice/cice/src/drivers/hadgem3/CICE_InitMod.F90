!  SVN:$Id: CICE_InitMod.F90 1115 2016-04-08 17:10:40Z eclare $
!=======================================================================
!
!  This module contains the CICE initialization routine that sets model
!  parameters and initializes the grid and CICE state variables.
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
! 2008: E. Hunke moved ESMF code to its own driver

      module CICE_InitMod

      use ice_kinds_mod

      implicit none
      private
      public :: CICE_Initialize, cice_init
      save

!=======================================================================

      contains

!=======================================================================

!  Initialize the basic state, grid and all necessary parameters for
!  running the CICE model.  Return the initial state in routine
!  export state.
!  Note: This initialization driver is designed for standalone and
!        CCSM-coupled applications.  For other
!        applications (e.g., standalone CAM), this driver would be
!        replaced by a different driver that calls subroutine cice_init,
!        where most of the work is done.

      subroutine CICE_Initialize

   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      call cice_init

      end subroutine CICE_Initialize

!=======================================================================
!
!  Initialize CICE model.

      subroutine cice_init

      use ice_arrays_column, only: hin_max, c_hi_range, zfswin, trcrn_sw, &
          ocean_bio_all, ice_bio_net, snow_bio_net
      use ice_calendar, only: dt, dt_dyn, time, istep, istep1, write_ic, &
          init_calendar, calendar
      use ice_colpkg, only: colpkg_init_itd, colpkg_init_itd_hist
      use ice_communicate, only: init_communicate, my_task, master_task
      use ice_diagnostics, only: init_diags
      use ice_domain, only: init_domain_blocks
      use ice_domain_size, only: ncat
      use ice_dyn_eap, only: init_eap
      use ice_dyn_shared, only: kdyn, init_evp
      use ice_fileunits, only: init_fileunits, nu_diag
      use ice_flux, only: init_coupler_flux, init_history_therm, &
          init_history_dyn, init_flux_atm, init_flux_ocn
      use ice_forcing, only: init_forcing_ocn, init_forcing_atmo, &
          get_forcing_atmo, get_forcing_ocn
      use ice_forcing_bgc, only: get_forcing_bgc, get_atm_bgc, &
          faero_data, faero_default, faero_optics
      use ice_grid, only: init_grid1, init_grid2
      use ice_history, only: init_hist, accum_hist
      use ice_restart_shared, only: restart, runid, runtype
      use ice_init, only: input_data, init_state
      use ice_init_column, only: init_thermo_vertical, init_shortwave, init_zbgc
      use ice_kinds_mod
      use ice_restoring, only: ice_HaloRestore_init
      use ice_colpkg_tracers, only: tr_aero, tr_zaero
      use ice_timers, only: timer_total, init_ice_timers, ice_timer_start
      use ice_transport_driver, only: init_transport
      use ice_colpkg_shared, only: skl_bgc, z_tracers
#ifdef popcice
      use drv_forcing, only: sst_sss
#endif

      call init_communicate     ! initial setup for message passing
      call init_fileunits       ! unit numbers
      call input_data           ! namelist variables
      if (trim(runid) == 'bering') call check_finished_file
      call init_zbgc            ! vertical biogeochemistry namelist

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file

      if (kdyn == 2) then
         call init_eap (dt_dyn) ! define eap dynamics parameters, variables
      else                      ! for both kdyn = 0 or 1
         call init_evp (dt_dyn) ! define evp dynamics parameters, variables
      endif

      call init_coupler_flux    ! initialize fluxes exchanged with coupler
#ifdef popcice
      call sst_sss              ! POP data for CICE initialization
#endif 
      call init_thermo_vertical ! initialize vertical thermodynamics
      call colpkg_init_itd(ncat, hin_max, nu_diag) ! ice thickness distribution
      if (my_task == master_task) &
         call colpkg_init_itd_hist(ncat, hin_max, nu_diag, c_hi_range) ! output
      call calendar(time)       ! determine the initial date

#ifndef CICE_IN_NEMO
      call init_forcing_ocn(dt) ! initialize sss and sst from data
#endif
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport
      call ice_HaloRestore_init ! restored boundary conditions

      call init_restart         ! initialize restart variables

      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      if (tr_aero .or. tr_zaero) call faero_optics !initialize aerosol optical 
                                                   !property tables

      ! Initialize shortwave components using swdn from previous timestep 
      ! if restarting. These components will be scaled to current forcing 
      ! in prep_radiation.
      if (trim(runtype) == 'continue' .or. restart) &
         call init_shortwave    ! initialize radiative transfer

      istep  = istep  + 1    ! update time step counters
      istep1 = istep1 + 1
      time = time + dt       ! determine the time and date
      call calendar(time)    ! at the end of the first timestep

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

#ifndef CICE_IN_NEMO
      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)
#endif

#ifndef coupled
#ifndef CCSMCOUPLED
      call get_forcing_atmo     ! atmospheric forcing from data
      call get_forcing_ocn(dt)  ! ocean forcing from data

      ! aerosols
      ! if (tr_aero)  call faero_data                   ! data file
      ! if (tr_zaero) call fzaero_data                  ! data file (gx1)
      if (tr_aero .or. tr_zaero)  call faero_default    ! default values

      if (skl_bgc .or. z_tracers) call get_forcing_bgc  ! biogeochemistry
#endif
#endif
      if (z_tracers) call get_atm_bgc                   ! biogeochemistry

      if (runtype == 'initial' .and. .not. restart) &
         call init_shortwave    ! initialize radiative transfer using current swdn

      call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler

      if (write_ic) call accum_hist(dt) ! write initial conditions 

      end subroutine cice_init

!=======================================================================

      subroutine init_restart

      use ice_arrays_column, only: dhsn
      use ice_blocks, only: nx_block, ny_block
      use ice_calendar, only: time, calendar
      use ice_colpkg, only: colpkg_aggregate
      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat, max_ntrcr, n_aero
      use ice_dyn_eap, only: read_restart_eap
      use ice_dyn_shared, only: kdyn
      use ice_flux, only: sss
      use ice_grid, only: tmask
      use ice_init, only: ice_ic
      use ice_init_column, only: init_age, init_FY, init_lvl, &
          init_meltponds_cesm,  init_meltponds_lvl, init_meltponds_topo, &
          init_aerosol, init_hbrine, init_bgc
      use ice_restart_column, only: restart_age, read_restart_age, &
          restart_FY, read_restart_FY, restart_lvl, read_restart_lvl, &
          restart_pond_cesm, read_restart_pond_cesm, &
          restart_pond_lvl, read_restart_pond_lvl, &
          restart_pond_topo, read_restart_pond_topo, &
          restart_aero, read_restart_aero, &
          restart_hbrine, read_restart_hbrine, &
          restart_zsal, restart_bgc
      use ice_restart_driver, only: restartfile, restartfile_v4
      use ice_restart_shared, only: runtype, restart
      use ice_state ! almost everything
      use ice_colpkg_tracers, only: tr_iage, tr_FY, tr_lvl, nt_alvl, nt_vlvl, &
          tr_pond_cesm, nt_apnd, nt_hpnd, tr_pond_lvl, nt_ipnd, &
          tr_pond_topo, tr_aero, tr_brine, nt_iage, nt_FY, nt_aero
      use ice_colpkg_shared, only: skl_bgc, z_tracers, solve_zsal

      integer(kind=int_kind) :: &
         i, j        , & ! horizontal indices
         iblk            ! block index

      if (trim(runtype) == 'continue') then 
         ! start from core restart file
         call restartfile()           ! given by pointer in ice_in
         call calendar(time)          ! update time parameters
         if (kdyn == 2) call read_restart_eap ! EAP
      else if (restart) then          ! ice_ic = core restart file
         call restartfile (ice_ic)    !  or 'default' or 'none'
         !!! uncomment to create netcdf
         ! call restartfile_v4 (ice_ic)  ! CICE v4.1 binary restart file
         !!! uncomment if EAP restart data exists
         ! if (kdyn == 2) call read_restart_eap
      endif         

      ! tracers
      ! ice age tracer   
      if (tr_iage) then 
         if (trim(runtype) == 'continue') &
              restart_age = .true.
         if (restart_age) then
            call read_restart_age
         else
            do iblk = 1, nblocks 
               call init_age(trcrn(:,:,nt_iage,:,iblk))
            enddo ! iblk
         endif
      endif
      ! first-year area tracer
      if (tr_FY) then
         if (trim(runtype) == 'continue') restart_FY = .true.
         if (restart_FY) then
            call read_restart_FY
         else
            do iblk = 1, nblocks 
               call init_FY(trcrn(:,:,nt_FY,:,iblk))
            enddo ! iblk
         endif
      endif
      ! level ice tracer
      if (tr_lvl) then
         if (trim(runtype) == 'continue') restart_lvl = .true.
         if (restart_lvl) then
            call read_restart_lvl
         else
            do iblk = 1, nblocks 
               call init_lvl(trcrn(:,:,nt_alvl,:,iblk), &
                             trcrn(:,:,nt_vlvl,:,iblk))
            enddo ! iblk
         endif
      endif
      ! CESM melt ponds
      if (tr_pond_cesm) then
         if (trim(runtype) == 'continue') &
              restart_pond_cesm = .true.
         if (restart_pond_cesm) then
            call read_restart_pond_cesm
         else
            do iblk = 1, nblocks 
               call init_meltponds_cesm(trcrn(:,:,nt_apnd,:,iblk), &
                                        trcrn(:,:,nt_hpnd,:,iblk))
            enddo ! iblk
         endif
      endif
      ! level-ice melt ponds
      if (tr_pond_lvl) then
         if (trim(runtype) == 'continue') &
              restart_pond_lvl = .true.
         if (restart_pond_lvl) then
            call read_restart_pond_lvl
         else
            do iblk = 1, nblocks 
               call init_meltponds_lvl(trcrn(:,:,nt_apnd,:,iblk), &
                                       trcrn(:,:,nt_hpnd,:,iblk), &
                                       trcrn(:,:,nt_ipnd,:,iblk), &
                                       dhsn(:,:,:,iblk))
            enddo ! iblk
         endif
      endif
      ! topographic melt ponds
      if (tr_pond_topo) then
         if (trim(runtype) == 'continue') &
              restart_pond_topo = .true.
         if (restart_pond_topo) then
            call read_restart_pond_topo
         else
            do iblk = 1, nblocks 
               call init_meltponds_topo(trcrn(:,:,nt_apnd,:,iblk), &
                                        trcrn(:,:,nt_hpnd,:,iblk), &
                                        trcrn(:,:,nt_ipnd,:,iblk))
            enddo ! iblk
         endif ! .not. restart_pond
      endif
      if (tr_aero) then ! ice aerosol
         if (trim(runtype) == 'continue') restart_aero = .true.
         if (restart_aero) then
            call read_restart_aero
         else
            do iblk = 1, nblocks 
               call init_aerosol(trcrn(:,:,nt_aero:nt_aero+4*n_aero-1,:,iblk))
            enddo ! iblk
         endif ! .not. restart_aero
      endif

      if (trim(runtype) == 'continue') then
         if (tr_brine) &
             restart_hbrine = .true.
         if (solve_zsal) &
             restart_zsal = .true.
         if (skl_bgc .or. z_tracers) &
             restart_bgc = .true.
      endif

      if (tr_brine .or. skl_bgc) then ! brine height tracer
         call init_hbrine
         if (tr_brine .and. restart_hbrine) call read_restart_hbrine
      endif

      if (solve_zsal .or. skl_bgc .or. z_tracers) call init_bgc ! biogeochemistry

      !-----------------------------------------------------------------
      ! aggregate tracers
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j,iblk)) &
         call colpkg_aggregate (ncat,               &
                                aicen(i,j,:,iblk),  &
                                trcrn(i,j,:,:,iblk),&
                                vicen(i,j,:,iblk),  &
                                vsnon(i,j,:,iblk),  &
                                aice (i,j,  iblk),  &
                                trcr (i,j,:,iblk),  &
                                vice (i,j,  iblk),  &
                                vsno (i,j,  iblk),  &
                                aice0(i,j,  iblk),  &
                                max_ntrcr,          &
                                trcr_depend,        &
                                trcr_base,          &
                                n_trcr_strata,      &
                                nt_strata)
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine init_restart

!=======================================================================
!
! Check whether a file indicating that the previous run finished cleanly
! If so, then do not continue the current restart.  This is needed only 
! for runs on machine 'bering' (set using runid = 'bering').
!
!  author: Adrian Turner, LANL

      subroutine check_finished_file()

      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_restart_shared, only: restart_dir

      character(len=char_len_long) :: filename
      logical :: lexist = .false.

      if (my_task == master_task) then
           
         filename = trim(restart_dir)//"finished"
         inquire(file=filename, exist=lexist)
         if (lexist) then
            call abort_ice("Found already finished file - quitting")
         end if

      endif

      end subroutine check_finished_file

!=======================================================================

      end module CICE_InitMod

!=======================================================================
