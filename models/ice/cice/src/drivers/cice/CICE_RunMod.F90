!  SVN:$Id: CICE_RunMod.F90 1150 2016-09-08 16:41:28Z eclare $
!=======================================================================
!
!  Main driver for time stepping of CICE.
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2006 ECH: moved exit timeLoop to prevent execution of unnecessary timestep
! 2006 ECH: Streamlined for efficiency 
! 2006 ECH: Converted to free source form (F90)
! 2007 BPB: Modified Delta-Eddington shortwave interface
! 2008 ECH: moved ESMF code to its own driver

      module CICE_RunMod

      use ice_kinds_mod

      implicit none
      private
      public :: CICE_Run, ice_step
      save

!=======================================================================

      contains

!=======================================================================
!
!  This is the main driver routine for advancing CICE forward in time.
!
!  author Elizabeth C. Hunke, LANL
!         Philip W. Jones, LANL
!         William H. Lipscomb, LANL

      subroutine CICE_Run

      use ice_calendar, only: istep, istep1, time, dt, stop_now, calendar
      use ice_forcing, only: get_forcing_atmo, get_forcing_ocn
      use ice_forcing_bgc, only: get_forcing_bgc, get_atm_bgc, fzaero_data, & 
          faero_default
      use ice_flux, only: init_flux_atm, init_flux_ocn
      use ice_colpkg_tracers, only: tr_aero, tr_zaero
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_couple, timer_step
      use ice_colpkg_shared, only: skl_bgc, z_tracers

   !--------------------------------------------------------------------
   !  initialize error code and step timer
   !--------------------------------------------------------------------

      call ice_timer_start(timer_step)   ! start timing entire run

#ifndef CICE_IN_NEMO
   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      timeLoop: do
#endif

         call ice_step

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date

         call calendar(time)    ! at the end of the timestep

#ifndef CICE_IN_NEMO
         if (stop_now >= 1) exit timeLoop
#endif

         call ice_timer_start(timer_couple)  ! atm/ocn coupling

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

         call init_flux_atm     ! initialize atmosphere fluxes sent to coupler
         call init_flux_ocn     ! initialize ocean fluxes sent to coupler

         call ice_timer_stop(timer_couple)    ! atm/ocn coupling

#ifndef CICE_IN_NEMO
      enddo timeLoop
#endif

   !--------------------------------------------------------------------
   ! end of timestep loop
   !--------------------------------------------------------------------

      call ice_timer_stop(timer_step)   ! end timestepping loop timer     

      end subroutine CICE_Run

!=======================================================================
!
!  Calls drivers for physics components, some initialization, and output
!
!  author Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL

      subroutine ice_step

      use ice_boundary, only: ice_HaloUpdate
      use ice_calendar, only: dt, dt_dyn, ndtd, diagfreq, write_restart, istep
      use ice_constants, only: field_loc_center, field_type_scalar, c0
      use ice_diagnostics, only: init_mass_diags, runtime_diags
      use ice_diagnostics_bgc, only: hbrine_diags, zsal_diags, bgc_diags
      use ice_domain, only: halo_info, nblocks
      use ice_domain_size, only: nslyr
      use ice_dyn_eap, only: write_restart_eap
      use ice_dyn_shared, only: kdyn
      use ice_flux, only: scale_factor, init_history_therm, &
          daidtt, daidtd, dvidtt, dvidtd, dagedtt, dagedtd
      use ice_history, only: accum_hist
      use ice_history_bgc, only: init_history_bgc
      use ice_restart, only: final_restart
      use ice_restart_column, only: write_restart_age, write_restart_FY, &
          write_restart_lvl, write_restart_pond_cesm, write_restart_pond_lvl, &
          write_restart_pond_topo, write_restart_aero, &
          write_restart_bgc, write_restart_hbrine
      use ice_restart_driver, only: dumpfile
      use ice_restoring, only: restore_ice, ice_HaloRestore
      use ice_state, only: trcrn
      use ice_colpkg_tracers, only: tr_iage, tr_FY, tr_lvl, &
          tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_brine, tr_aero
      use ice_step_mod, only: prep_radiation, step_therm1, step_therm2, &
          update_state, step_dyn_horiz, step_dyn_ridge, step_radiation, &
          biogeochemistry
      use ice_colpkg_shared, only: calc_Tsfc, skl_bgc, solve_zsal, z_tracers
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_diags, timer_column, timer_thermo, timer_bound, &
          timer_hist, timer_readwrite

      integer (kind=int_kind) :: &
         iblk        , & ! block index 
         k               ! dynamics supercycling index

      real (kind=dbl_kind) :: &
         offset          ! d(age)/dt time offset

      !-----------------------------------------------------------------
      ! restoring on grid boundaries
      !-----------------------------------------------------------------

         if (restore_ice) call ice_HaloRestore

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics/history
         call init_mass_diags   ! diagnostics per timestep
         call init_history_therm
         call init_history_bgc
         call ice_timer_stop(timer_diags)   ! diagnostics/history

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Scale radiation fields
      !-----------------------------------------------------------------

            if (calc_Tsfc) call prep_radiation (dt, iblk)

      !-----------------------------------------------------------------
      ! thermodynamics and biogeochemistry
      !-----------------------------------------------------------------
            
            call step_therm1     (dt, iblk) ! vertical thermodynamics
            call biogeochemistry (dt, iblk) ! biogeochemistry
            call step_therm2     (dt, iblk) ! ice thickness distribution thermo

         enddo ! iblk
         !$OMP END PARALLEL DO

         ! clean up, update tendency diagnostics
         offset = dt
         call update_state (dt, daidtt, dvidtt, dagedtt, offset)

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         do k = 1, ndtd

            ! momentum, stress, transport
            call step_dyn_horiz (dt_dyn)

            ! ridging
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call step_dyn_ridge (dt_dyn, ndtd, iblk)
            enddo
            !$OMP END PARALLEL DO

            ! clean up, update tendency diagnostics
            offset = c0
            call update_state (dt_dyn, daidtd, dvidtd, dagedtd, offset)

         enddo

      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

            call step_radiation (dt, iblk)

      !-----------------------------------------------------------------
      ! get ready for coupling and the next time step
      !-----------------------------------------------------------------

            call coupling_prep (iblk)

         enddo ! iblk
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (scale_factor,     halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics
         if (mod(istep,diagfreq) == 0) then
            call runtime_diags(dt)          ! log file
            if (solve_zsal) call zsal_diags(dt)  
            if (skl_bgc .or. z_tracers)  call bgc_diags (dt)
            if (tr_brine) call hbrine_diags(dt)
         endif
         call ice_timer_stop(timer_diags)   ! diagnostics

         call ice_timer_start(timer_hist)   ! history
         call accum_hist (dt)               ! history file
         call ice_timer_stop(timer_hist)    ! history

         call ice_timer_start(timer_readwrite)  ! reading/writing
         if (write_restart == 1) then
            call dumpfile     ! core variables for restarting
            if (tr_iage)      call write_restart_age
            if (tr_FY)        call write_restart_FY
            if (tr_lvl)       call write_restart_lvl
            if (tr_pond_cesm) call write_restart_pond_cesm
            if (tr_pond_lvl)  call write_restart_pond_lvl
            if (tr_pond_topo) call write_restart_pond_topo
            if (tr_aero)      call write_restart_aero
            if (solve_zsal .or. skl_bgc .or. z_tracers) &
                              call write_restart_bgc 
            if (tr_brine)     call write_restart_hbrine
            if (kdyn == 2)    call write_restart_eap
            call final_restart
         endif

         call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_step
    
!=======================================================================
!
! Prepare for coupling
!
! authors: Elizabeth C. Hunke, LANL

      subroutine coupling_prep (iblk)

      use ice_arrays_column, only: alvdfn, alidfn, alvdrn, alidrn, &
          albicen, albsnon, albpndn, apeffn, fzsal_g, fzsal, snowfracn
      use ice_blocks, only: block, nx_block, ny_block
      use ice_calendar, only: dt, nstreams
      use ice_colpkg_shared, only: calc_Tsfc, oceanmixed_ice, max_aero
      use ice_colpkg_tracers, only: nbtrcr
      use ice_constants, only: c0, c1, puny, rhofresh
      use ice_domain_size, only: ncat
      use ice_flux, only: alvdf, alidf, alvdr, alidr, albice, albsno, &
          albpnd, albcnt, apeff_ai, coszen, fpond, fresh, l_mpond_fresh, &
          alvdf_ai, alidf_ai, alvdr_ai, alidr_ai, fhocn_ai, &
          fresh_ai, fsalt_ai, fsalt, &
          fswthru_ai, fhocn, fswthru, scale_factor, snowfrac, &
          swvdr, swidr, swvdf, swidf, Tf, Tair, Qa, strairxT, strairyt, &
          fsens, flat, fswabs, flwout, evap, Tref, Qref, &
          fsurfn_f, flatn_f, scale_fluxes, frzmlt_init, frzmlt
      use ice_flux_bgc, only: faero_ocn, fzsal_ai, fzsal_g_ai, flux_bio, flux_bio_ai
      use ice_grid, only: tmask
      use ice_state, only: aicen, aice, aice_init
      use ice_step_mod, only: ocean_mixed_layer
      use ice_timers, only: timer_couple, ice_timer_start, ice_timer_stop

      integer (kind=int_kind), intent(in) :: & 
         iblk            ! block index 

      ! local variables

      integer (kind=int_kind) :: & 
         n           , & ! thickness category index
         i,j         , & ! horizontal indices
         k               ! tracer index

      real (kind=dbl_kind) :: &
         cszn        , & ! counter for history averaging
         netsw           ! flag for shortwave radiation presence

      !-----------------------------------------------------------------
      ! Save current value of frzmlt for diagnostics.
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            frzmlt_init  (i,j,iblk) = frzmlt(i,j,iblk)
         enddo
         enddo

         call ice_timer_start(timer_couple)   ! atm/ocn coupling

         if (oceanmixed_ice) &
         call ocean_mixed_layer (dt,iblk) ! ocean surface fluxes and sst

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0

            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
            apeff_ai(i,j,iblk) = c0
            snowfrac(i,j,iblk) = c0

            ! for history averaging
            cszn = c0
            netsw = swvdr(i,j,iblk)+swidr(i,j,iblk)+swvdf(i,j,iblk)+swidf(i,j,iblk)
            if (netsw > puny) cszn = c1
            do n = 1, nstreams
               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
            enddo
         enddo
         enddo
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (aicen(i,j,n,iblk) > puny) then
                  
            alvdf(i,j,iblk) = alvdf(i,j,iblk) &
               + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidf(i,j,iblk) = alidf(i,j,iblk) &
               + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alvdr(i,j,iblk) = alvdr(i,j,iblk) &
               + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidr(i,j,iblk) = alidr(i,j,iblk) &
               + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)

            netsw = swvdr(i,j,iblk) + swidr(i,j,iblk) &
                  + swvdf(i,j,iblk) + swidf(i,j,iblk)
            if (netsw > puny) then ! sun above horizon
            albice(i,j,iblk) = albice(i,j,iblk) &
               + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
            albsno(i,j,iblk) = albsno(i,j,iblk) &
               + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
            albpnd(i,j,iblk) = albpnd(i,j,iblk) &
               + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif

            apeff_ai(i,j,iblk) = apeff_ai(i,j,iblk) &       ! for history
               + apeffn(i,j,n,iblk)*aicen(i,j,n,iblk)
            snowfrac(i,j,iblk) = snowfrac(i,j,iblk) &       ! for history
               + snowfracn(i,j,n,iblk)*aicen(i,j,n,iblk)
               
            endif ! aicen > puny
         enddo
         enddo
         enddo

         do j = 1, ny_block
         do i = 1, nx_block

      !-----------------------------------------------------------------
      ! reduce fresh by fpond for coupling
      !-----------------------------------------------------------------

            if (l_mpond_fresh) then
               fpond(i,j,iblk) = fpond(i,j,iblk) * rhofresh/dt
               fresh(i,j,iblk) = fresh(i,j,iblk) - fpond(i,j,iblk)
            endif

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_ai  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_ai  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_ai  (i,j,iblk) = alidr  (i,j,iblk)
            fresh_ai  (i,j,iblk) = fresh  (i,j,iblk)
            fsalt_ai  (i,j,iblk) = fsalt  (i,j,iblk)
            fhocn_ai  (i,j,iblk) = fhocn  (i,j,iblk)
            fswthru_ai(i,j,iblk) = fswthru(i,j,iblk)
            fzsal_ai  (i,j,iblk) = fzsal  (i,j,iblk) 
            fzsal_g_ai(i,j,iblk) = fzsal_g(i,j,iblk)  

            if (nbtrcr > 0) then
            do k = 1, nbtrcr
              flux_bio_ai  (i,j,k,iblk) = flux_bio  (i,j,k,iblk)
            enddo
            endif

      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            scale_factor(i,j,iblk) = &
                       swvdr(i,j,iblk)*(c1 - alvdr_ai(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_ai(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_ai(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_ai(i,j,iblk))

         enddo
         enddo

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area 
      !  - the CCSM coupler assumes fluxes are per unit ice area
      !  - also needed for global budget in diagnostics
      !-----------------------------------------------------------------

         call scale_fluxes (nx_block,            ny_block,           &
                            tmask    (:,:,iblk), nbtrcr, max_aero,   &
                            aice     (:,:,iblk), Tf      (:,:,iblk), &
                            Tair     (:,:,iblk), Qa      (:,:,iblk), &
                            strairxT (:,:,iblk), strairyT(:,:,iblk), &
                            fsens    (:,:,iblk), flat    (:,:,iblk), &
                            fswabs   (:,:,iblk), flwout  (:,:,iblk), &
                            evap     (:,:,iblk),                     &
                            Tref     (:,:,iblk), Qref    (:,:,iblk), &
                            fresh    (:,:,iblk), fsalt   (:,:,iblk), &
                            fhocn    (:,:,iblk), fswthru (:,:,iblk), &
                            faero_ocn(:,:,:,iblk),                   &
                            alvdr    (:,:,iblk), alidr   (:,:,iblk), &
                            alvdf    (:,:,iblk), alidf   (:,:,iblk), &
                            fzsal    (:,:,iblk), fzsal_g (:,:,iblk), &
                            flux_bio(:,:,1:nbtrcr,iblk))
 
!echmod - comment this out for efficiency, if .not. calc_Tsfc
         if (.not. calc_Tsfc) then

       !---------------------------------------------------------------
       ! If surface fluxes were provided, conserve these fluxes at ice 
       ! free points by passing to ocean. 
       !---------------------------------------------------------------

            call sfcflux_to_ocn & 
                         (nx_block,              ny_block,             &
                          tmask   (:,:,iblk),    aice_init(:,:,iblk),  &
                          fsurfn_f (:,:,:,iblk), flatn_f(:,:,:,iblk),  &
                          fresh    (:,:,iblk),   fhocn    (:,:,iblk))
         endif                 
!echmod

         call ice_timer_stop(timer_couple)   ! atm/ocn coupling

      end subroutine coupling_prep

!=======================================================================
!
! If surface heat fluxes are provided to CICE instead of CICE calculating
! them internally (i.e. .not. calc_Tsfc), then these heat fluxes can 
! be provided at points which do not have ice.  (This is could be due to
! the heat fluxes being calculated on a lower resolution grid or the
! heat fluxes not recalculated at every CICE timestep.)  At ice free points, 
! conserve energy and water by passing these fluxes to the ocean.
!
! author: A. McLaren, Met Office

      subroutine sfcflux_to_ocn(nx_block,   ny_block,     &
                                tmask,      aice,         &
                                fsurfn_f,   flatn_f,      &
                                fresh,      fhocn)

      use ice_domain_size, only: ncat

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block  ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: &
          aice        ! initial ice concentration

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
          intent(in) :: &
          fsurfn_f, & ! net surface heat flux (provided as forcing)
          flatn_f     ! latent heat flux (provided as forcing)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          fresh        , & ! fresh water flux to ocean         (kg/m2/s)
          fhocn            ! actual ocn/ice heat flx           (W/m**2)

#ifdef CICE_IN_NEMO

      ! local variables
      integer (kind=int_kind) :: &
          i, j, n    ! horizontal indices
      
      real (kind=dbl_kind)    :: &
          rLsub            ! 1/Lsub

      rLsub = c1 / Lsub

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j) .and. aice(i,j) <= puny) then
               fhocn(i,j)      = fhocn(i,j)              &
                            + fsurfn_f(i,j,n) + flatn_f(i,j,n)
               fresh(i,j)      = fresh(i,j)              &
                                 + flatn_f(i,j,n) * rLsub
            endif
         enddo   ! i
         enddo   ! j
      enddo      ! n

#endif 

      end subroutine sfcflux_to_ocn

!=======================================================================

      end module CICE_RunMod

!=======================================================================
