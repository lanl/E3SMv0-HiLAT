!  SVN:$Id: ice_init_column.F90 1110 2016-03-08 21:22:20Z njeffery $
!=========================================================================
!
! Initialization routines for the column package.
!
! author: Elizabeth C. Hunke, LANL
!
! 2014: Moved subroutines from column package modules

      module ice_init_column

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: ncat, nilyr, nslyr, max_blocks 

      implicit none
      save

      private
      public :: init_thermo_vertical, init_shortwave, &
                init_age, init_FY, init_lvl, &
                init_meltponds_cesm, init_meltponds_lvl, init_meltponds_topo, &
                init_aerosol, init_bgc, init_hbrine, init_zbgc

!=======================================================================

      contains

!=======================================================================
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine init_thermo_vertical

      use ice_blocks, only: nx_block, ny_block
      use ice_colpkg, only: colpkg_init_thermo
      use ice_flux, only: salinz, Tmltz

      integer (kind=int_kind) :: &
         i, j, iblk, &  ! horizontal indices
         k              ! ice layer index

      real (kind=dbl_kind), dimension(nilyr+1) :: &
         sprofile                         ! vertical salinity profile

      !-----------------------------------------------------------------
      ! initialize heat_capacity, l_brine, and salinity profile
      !-----------------------------------------------------------------

      call colpkg_init_thermo(nilyr, sprofile)

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,k)
      do iblk = 1,max_blocks
      do j = 1, ny_block
      do i = 1, nx_block
         do k = 1, nilyr+1
            salinz(i,j,k,iblk) = sprofile(k)
            Tmltz (i,j,k,iblk) = -salinz(i,j,k,iblk)*depressT
         enddo ! k
      enddo    ! i
      enddo    ! j
      enddo    ! iblk
      !$OMP END PARALLEL DO

      end subroutine init_thermo_vertical

!=======================================================================
!
!  Initialize shortwave

      subroutine init_shortwave

      use ice_arrays_column, only: fswpenln, Iswabsn, Sswabsn, albicen, &
          albsnon, alvdrn, alidrn, alvdfn, alidfn, fswsfcn, fswthrun, &
          fswintn, albpndn, apeffn, trcrn_sw, dhsn, ffracn, snowfracn, &
          kaer_tab, waer_tab, gaer_tab, kaer_bc_tab, waer_bc_tab, gaer_bc_tab, bcenh, &
          swgrid, igrid
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_calendar, only: nstreams, istep1, dt, calendar_type, &
          days_per_year, nextsw_cday, yday, sec
      use ice_communicate, only: my_task, master_task
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, &
          diagnostic_abort
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: n_aero, n_zaero, ncat, nilyr, nslyr, n_algae, nblyr
      use ice_fileunits, only: nu_diag
      use ice_flux, only: alvdf, alidf, alvdr, alidr, &
                          alvdr_ai, alidr_ai, alvdf_ai, alidf_ai, &
                          swvdr, swvdf, swidr, swidf, scale_factor, snowfrac, &
                          albice, albsno, albpnd, apeff_ai, albcnt, coszen, fsnow
      use ice_grid, only: tlat, tlon, tmask
      use ice_restart_shared, only: restart, runtype
      use ice_state, only: aicen, vicen, vsnon, trcrn

      ! column package includes
      use ice_colpkg, only: colpkg_step_radiation, colpkg_init_orbit
      use ice_colpkg_shared, only: shortwave, dEdd_algae, modal_aero
      use ice_colpkg_tracers, only: nt_Tsfc, &
          nt_alvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero, tr_bgc_N, &
          tr_zaero, nlt_chl_sw, nlt_zaero_sw, ntrcr, nbtrcr, nbtrcr_sw, nt_fbri, tr_brine, &
          nt_zaero

      integer (kind=int_kind) :: &
         i, j , k    , & ! horizontal indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n                  ! thickness category index

      real (kind=dbl_kind) :: &
         cszn        , & ! counter for history averaging
         netsw           ! flag for shortwave radiation presence

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop       , & ! if true, abort the model
         l_print_point, & ! flag to print designated grid point diagnostics
         debug            ! if true, print diagnostics

      character (char_len) :: stop_label

      integer (kind=int_kind) :: &
         ipoint

      real (kind=dbl_kind), dimension(ncat) :: &
         fbri                 ! brine height to ice thickness

      real(kind= dbl_kind), dimension(ntrcr, ncat) :: &
         ztrcr

      real(kind= dbl_kind), dimension(nbtrcr_sw, ncat) :: &
         ztrcr_sw

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     cszn,l_print_point,debug,ipoint)
      do iblk=1,nblocks

         ! Initialize
         fswpenln(:,:,:,:,iblk) = c0
         Iswabsn(:,:,:,:,iblk) = c0
         Sswabsn(:,:,:,:,iblk) = c0

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j = 1, ny_block ! can be jlo, jhi
         do i = 1, nx_block ! can be ilo, ihi

            l_print_point = .false.
            debug = .false.
            if (debug .and. print_points) then
               do ipoint = 1, npnt
                  if (my_task == pmloc(ipoint) .and. &
                       i == piloc(ipoint) .and. &
                       j == pjloc(ipoint)) &
                       l_print_point = .true.
                       write (nu_diag, *) 'my_task = ',my_task
               enddo ! n
            endif

            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0
            alvdr_ai(i,j,iblk) = c0
            alidr_ai(i,j,iblk) = c0
            alvdf_ai(i,j,iblk) = c0
            alidf_ai(i,j,iblk) = c0
            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
            snowfrac(i,j,iblk) = c0
            apeff_ai(i,j,iblk) = c0

            do n = 1, ncat
               alvdrn(i,j,n,iblk) = c0
               alidrn(i,j,n,iblk) = c0
               alvdfn(i,j,n,iblk) = c0
               alidfn(i,j,n,iblk) = c0
               fswsfcn(i,j,n,iblk) = c0
               fswintn(i,j,n,iblk) = c0
               fswthrun(i,j,n,iblk) = c0
            enddo   ! ncat

         enddo
         enddo
         do j = jlo, jhi
         do i = ilo, ihi

            if (trim(shortwave) == 'dEdd') then ! delta Eddington

#ifndef CCSMCOUPLED
               ! initialize orbital parameters
               ! These come from the driver in the coupled model.
               call colpkg_init_orbit(nu_diag, l_stop, stop_label) 

               if (l_stop) then
                  call diagnostic_abort(i, j, iblk, istep1, stop_label)
               endif
#endif
            endif

         fbri(:) = c0
         ztrcr_sw(:,:) = c0
         do n = 1, ncat
           do k = 1, ntrcr
             ztrcr(k,n) = trcrn(i,j,k,n,iblk)
           enddo
           if (tr_brine)  fbri(n) = trcrn(i,j,nt_fbri,n,iblk)
         enddo

         if (tmask(i,j,iblk))&
         call colpkg_step_radiation (dt,         ncat,                    &
                          n_algae,   tr_zaero, nblyr,                     &
                          ntrcr,     nbtrcr,   nbtrcr_sw,                 &
                          nilyr,    nslyr,       n_aero,                  &
                          n_zaero,  dEdd_algae,  nlt_chl_sw,              &
                          nlt_zaero_sw(:),                                &
                          swgrid(:),           igrid(:),                  &
                          fbri(:),                                        &
                          aicen(i,j,:,iblk),     vicen(i,j,:,iblk),       &
                          vsnon(i,j,:,iblk),                              &
                          trcrn(i,j,nt_Tsfc,:,iblk),                      &
                          trcrn(i,j,nt_alvl,:,iblk),                      &
                          trcrn(i,j,nt_apnd,:,iblk),                      &
                          trcrn(i,j,nt_hpnd,:,iblk),                      &
                          trcrn(i,j,nt_ipnd,:,iblk),                      &
                          trcrn(i,j,nt_aero:nt_aero+4*n_aero-1,:,iblk),   &
                          ztrcr_sw,                                       &
                          ztrcr,                                          &
                          TLAT(i,j,iblk),        TLON(i,j,iblk),          &
                          calendar_type,         days_per_year,           &
                          nextsw_cday,           yday,                    &
                          sec,                                            &
                          kaer_tab, waer_tab,                             &
                          gaer_tab,                                       &
                          kaer_bc_tab(:,:),      waer_bc_tab(:,:),        &
                          gaer_bc_tab(:,:),      bcenh(:,:,:),            &
                          modal_aero,                                     &
                          swvdr(i,j,iblk),       swvdf(i,j,iblk),         &
                          swidr(i,j,iblk),       swidf(i,j,iblk),         &
                          coszen(i,j,iblk),      fsnow(i,j,iblk),         &
                          alvdrn(i,j,:,iblk),    alvdfn(i,j,:,iblk),      &
                          alidrn(i,j,:,iblk),    alidfn(i,j,:,iblk),      &
                          fswsfcn(i,j,:,iblk),   fswintn(i,j,:,iblk),     &
                          fswthrun(i,j,:,iblk),  fswpenln(i,j,:,:,iblk),  &
                          Sswabsn(i,j,:,:,iblk), Iswabsn(i,j,:,:,iblk),   &
                          albicen(i,j,:,iblk),   albsnon(i,j,:,iblk),     &
                          albpndn(i,j,:,iblk),   apeffn(i,j,:,iblk),      &
                          snowfracn(i,j,:,iblk),                          &
                          dhsn(i,j,:,iblk),      ffracn(i,j,:,iblk),      &
                          nu_diag,               l_print_point,           &
                          initonly = .true.)

      !-----------------------------------------------------------------
      ! Define aerosol tracer on shortwave grid
      !-----------------------------------------------------------------

      if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
        do n = 1, ncat
           do k = 1, nbtrcr_sw
              trcrn_sw(i,j,k,n,iblk) = ztrcr_sw(k,n)
           enddo
        enddo
      endif

      !-----------------------------------------------------------------
      ! Aggregate albedos 
      !-----------------------------------------------------------------

            do n = 1, ncat
               
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
                  
                  apeff_ai(i,j,iblk) = apeff_ai(i,j,iblk) &
                       + apeffn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  snowfrac(i,j,iblk) = snowfrac(i,j,iblk) &
                       + snowfracn(i,j,n,iblk)*aicen(i,j,n,iblk)
               
               endif ! aicen > puny

            enddo  ! ncat

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_ai  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_ai  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_ai  (i,j,iblk) = alidr  (i,j,iblk)
            
            ! for history averaging
!echmod?            cszn = c0
!echmod            if (coszen(i,j,iblk) > puny) cszn = c1
!echmod            do n = 1, nstreams
!echmod               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
!echmod            enddo
            
      !----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !----------------------------------------------------------------
            if (runtype == 'initial' .and. .not. restart) then
               scale_factor(i,j,iblk) = &
	 	      swvdr(i,j,iblk)*(c1 - alvdr_ai(i,j,iblk)) &
	 	    + swvdf(i,j,iblk)*(c1 - alvdf_ai(i,j,iblk)) &
 	            + swidr(i,j,iblk)*(c1 - alidr_ai(i,j,iblk)) &
	 	    + swidf(i,j,iblk)*(c1 - alidf_ai(i,j,iblk))
	    endif

         enddo ! i
         enddo ! j
      enddo ! iblk

      end subroutine init_shortwave

!=======================================================================

!  Initialize ice age tracer (call prior to reading restart data)

      subroutine init_age(iage)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: iage

      iage(:,:,:) = c0

      end subroutine init_age

!=======================================================================

!  Initialize ice FY tracer (call prior to reading restart data)

      subroutine init_FY(firstyear)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: firstyear

      firstyear(:,:,:) = c0

      end subroutine init_FY

!=======================================================================

!  Initialize ice lvl tracers (call prior to reading restart data)

      subroutine init_lvl(alvl, vlvl) 

      use ice_constants, only: c1

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         alvl , & ! level ice area fraction
         vlvl     ! level ice volume

      alvl(:,:,:) = c1 ! level ice area fraction
      vlvl(:,:,:) = c1 ! level ice volume

      end subroutine init_lvl

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_cesm(apnd, hpnd)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd     ! melt pond depth

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0

      end subroutine init_meltponds_cesm

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_lvl(apnd, hpnd, ipnd, dhsn)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd , & ! melt pond depth
         ipnd , & ! melt pond refrozen lid thickness
         dhsn     ! depth difference for snow on sea ice and pond ice

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0
      ipnd(:,:,:) = c0
      dhsn(:,:,:) = c0

      end subroutine init_meltponds_lvl

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_topo(apnd, hpnd, ipnd)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd , & ! melt pond depth
         ipnd     ! melt pond refrozen lid thickness

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0
      ipnd(:,:,:) = c0
        
      end subroutine init_meltponds_topo

!=======================================================================

!  Initialize ice aerosol tracer (call prior to reading restart data)

      subroutine init_aerosol(aero)

      real(kind=dbl_kind), dimension(:,:,:,:), intent(out) :: &
         aero ! aerosol tracers

      aero(:,:,:,:) = c0

      end subroutine init_aerosol

!=======================================================================

!  Initialize vertical profile for biogeochemistry

      subroutine init_bgc() 

      use ice_arrays_column, only: zfswin, trcrn_sw, &
          ocean_bio_all, ice_bio_net, snow_bio_net, &
          cgrid, igrid, bphi, iDi, bTiz, iki, &
          Rayleigh_criteria, Rayleigh_real
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_calendar, only: dt, istep1
      use ice_communicate, only: my_task
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, &
          diagnostic_abort
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: nblyr, nilyr
      use ice_fileunits, only: nu_diag
      use ice_flux, only: sss
      use ice_flux_bgc, only: nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use ice_forcing_bgc, only: init_bgc_data, get_forcing_bgc
      use ice_restart_column, only: restart_zsal, &
          read_restart_bgc, restart_bgc
      use ice_state, only: trcrn, aicen, vicen, vsnon
      use ice_colpkg_shared, only: solve_zsal, &
         max_algae, max_don, max_doc, max_dic, max_aero, max_fe, &
         max_nbtrcr

      ! column package includes
      use ice_colpkg, only: colpkg_init_bgc, colpkg_init_zsalinity, colpkg_init_ocean_conc 
      use ice_colpkg_tracers, only: nbtrcr, ntrcr, nt_bgc_S, ntrcr_o, &
           nt_sice, nt_fbri

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk       , & ! horizontal indices
         ilo,ihi,jlo,jhi  , & ! beginning and end of physical domain
         k,m              , & ! vertical index
         n                , & ! category index
         ipoint

      logical (kind=log_kind) :: &
         l_print_point, & ! flag to print designated grid point diagnostics
         debug        , & ! prints debugging output if true
         l_stop           ! if true, print diagnostics and abort on return
        
      character (char_len) :: stop_label

      type (block) :: &
         this_block      ! block information for current block

      real(kind=dbl_kind), dimension(ntrcr,ncat) :: &
         trcrn_bgc 
      
      real(kind=dbl_kind), dimension(nilyr,ncat) :: &
         sicen    

      ! Initialize

      l_stop = .false.

      bphi(:,:,:,:,:) = c0   ! initial porosity for no ice 
      iDi (:,:,:,:,:) = c0   ! interface diffusivity
      bTiz(:,:,:,:,:) = c0   ! initial bio grid ice temperature
      iki (:,:,:,:,:) = c0   ! permeability

      ocean_bio_all(:,:,:,:)   = c0
      ice_bio_net  (:,:,:,:)   = c0 ! integrated ice tracer conc (mmol/m^2 or mg/m^2) 
      snow_bio_net (:,:,:,:)   = c0 ! integrated snow tracer conc (mmol/m^2 or mg/m^2)
      zfswin       (:,:,:,:,:) = c0 ! shortwave flux on bio grid
      trcrn_sw     (:,:,:,:,:) = c0 ! tracers active in the shortwave calculation
      trcrn_bgc    (:,:) = c0

      if (.not. solve_zsal) restart_zsal = .false.
      if (restart_zsal .or. restart_bgc) call read_restart_bgc  

      !-----------------------------------------------------------------
      ! zsalinity initialization
      !-----------------------------------------------------------------
      
      if (solve_zsal) then

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi  
               call colpkg_init_zsalinity(nblyr, ntrcr_o, restart_zsal, Rayleigh_criteria(i,j,iblk), &
                      Rayleigh_real(i,j,iblk), trcrn_bgc, nt_bgc_S, ncat, sss(i,j,iblk))
               if (.not. restart_zsal) then
               do n = 1,ncat
                 do k  = 1, nblyr
                   trcrn(i,j,nt_bgc_S+k-1,n,iblk) = trcrn_bgc(nt_bgc_S-1+k-ntrcr_o,n)
                 enddo
               enddo
               endif
            enddo      ! i
            enddo      ! j  
         enddo         ! iblk
      endif ! solve_zsal

      !-----------------------------------------------------------------
      ! biogeochemistry initialization
      !-----------------------------------------------------------------

      if (restart_bgc) then       
     
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi  
                  l_print_point = .false.
                  debug = .false.
                  if (debug .and. print_points) then
                     do ipoint = 1, npnt
                        if (my_task == pmloc(ipoint) .and. &
                                  i == piloc(ipoint) .and. &
                                  j == pjloc(ipoint)) &
                           l_print_point = .true.
                           write (nu_diag, *) 'my_task = ',my_task
                     enddo ! ipoint
                  endif
            enddo  ! i
            enddo  ! j
         enddo     ! iblk

      else  ! not restarting

      !-----------------------------------------------------------------
      ! Initial Ocean Values if not coupled to the ocean bgc
      !-----------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi  
               call colpkg_init_ocean_conc (amm(i,j,  iblk), dmsp(i,j,  iblk), dms(i,j,  iblk), &
                    algalN(i,j, :, iblk), doc(i,j,:,  iblk), dic(i,j,:,  iblk), don(i,j,:,  iblk), &
                    fed(i,j,:,  iblk), fep(i,j,:,  iblk), hum(i,j,  iblk),  nit(i,j,  iblk), &
                    sil(i,j,  iblk), zaeros(i,j, :, iblk), max_dic, max_don, max_fe, max_aero)

            enddo  ! i
            enddo  ! j

         enddo     ! iblk

         call init_bgc_data(fed(:,:,1,:),fep(:,:,1,:)) ! input dFe from file
         call get_forcing_bgc                          ! defines nit and sil

      endif     ! restart

      !-----------------------------------------------------------------
      ! Complete bgc initialization
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi  

            do n = 1, ncat
            do k = 1, nilyr
               sicen(k,n) = trcrn(i,j,nt_sice+k-1,n,iblk)
            enddo
            do k = ntrcr_o+1, ntrcr
               trcrn_bgc(k-ntrcr_o,n) = trcrn(i,j,k,n,iblk)
            enddo
            enddo

            call colpkg_init_bgc(dt, ncat, nblyr, nilyr, &
               ntrcr_o, &
               cgrid, igrid, &
               restart_bgc, ntrcr, nbtrcr, &
               sicen(:,:), &
               trcrn_bgc(:,:), &
               sss(i,j,  iblk), &
               nit(i,j,  iblk), amm   (i,j,  iblk), &
               sil(i,j,  iblk), dmsp  (i,j,  iblk), &
               dms(i,j,  iblk), algalN(i,j,:,iblk), &
               doc(i,j,:,iblk), don   (i,j,:,iblk), &
               dic(i,j,:,iblk), fed   (i,j,:,iblk), &
               fep(i,j,:,iblk), zaeros(i,j,:,iblk), &
               hum(i,j,  iblk),           &
               ocean_bio_all(i,j,:,iblk), &
               max_algae, max_doc, max_dic, max_don,  max_fe, max_nbtrcr, max_aero, &
               l_stop, stop_label)

            if (l_stop) then
               call diagnostic_abort(i, j, iblk, istep1, stop_label)
            endif
            do n = 1, ncat
            do k = ntrcr_o+1, ntrcr
               trcrn(i,j,k,n,iblk) = trcrn_bgc(k-ntrcr_o,n)
            enddo
            enddo

         enddo  ! i
         enddo  ! j

      enddo     ! iblk

      end subroutine init_bgc

!=======================================================================

!  Initialize brine height tracer

      subroutine init_hbrine()

      use ice_arrays_column, only: first_ice, bgrid, igrid, cgrid, &
          icgrid, swgrid
      use ice_domain_size, only: nblyr
      use ice_state, only: trcrn
      use ice_colpkg, only: colpkg_init_hbrine
      use ice_colpkg_tracers, only: nt_fbri, tr_brine
      use ice_colpkg_shared, only: phi_snow

      call colpkg_init_hbrine(bgrid, igrid, cgrid, icgrid, &
            swgrid, nblyr, nilyr, phi_snow)

      first_ice(:,:,:,:) = .true.            
      if (tr_brine) trcrn(:,:,nt_fbri,:,:) = c1

      end subroutine init_hbrine

!=======================================================================

! Namelist variables, set to default values; may be altered at run time
! 
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL

      subroutine init_zbgc

      use ice_broadcast, only: broadcast_scalar
      use ice_communicate, only: my_task, master_task
      use ice_constants_colpkg, only: c1, p5, c0, c5, rhos, rhoi
      use ice_domain_size, only: max_ntrcr, nblyr, nilyr, nslyr, &
                           n_algae, n_zaero, n_doc, n_dic, n_don, &
                           n_fed, n_fep, max_nsw, n_bgc
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nml_filename, get_fileunit, &
                               release_fileunit, nu_diag
      use ice_restart_column, only: restart_bgc, restart_zsal, &
          restart_hbrine
      use ice_state, only: trcr_base, trcr_depend, n_trcr_strata, &
          nt_strata      

      use ice_colpkg, only:colpkg_init_bgc_trcr,  colpkg_init_zbgc

      use ice_colpkg_tracers, only: tr_brine, &
          ntrcr,         nbtrcr,       nbtrcr_sw,    &
          ntrcr_o, &
          tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil,   &
          tr_bgc_DMS,    tr_bgc_PON,   tr_bgc_S,     &
          tr_bgc_N,      tr_bgc_C,     tr_bgc_chl,   &
          tr_bgc_DON,    tr_bgc_Fe,    tr_zaero,     &
          tr_bgc_hum,    tr_aero,      nt_fbri,      &  
          nt_bgc_Nit,    nt_bgc_Am,    nt_bgc_Sil,   &
          nt_bgc_DMS,    nt_bgc_PON,   nt_bgc_S,     &
          nt_bgc_N,      nt_bgc_C,     nt_bgc_chl,   &
          nt_bgc_DOC,    nt_bgc_DON,   nt_bgc_DIC,   &
          nt_zaero,      nt_bgc_DMSPp, nt_bgc_DMSPd, &
          nt_bgc_Fed,    nt_bgc_Fep,  nt_zbgc_frac, &
          nlt_zaero_sw,  nlt_chl_sw, &
          nlt_bgc_N, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, &
          nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, nlt_bgc_C, nlt_bgc_chl, &
          nlt_bgc_DIC, nlt_bgc_DOC, nlt_bgc_PON, &
          nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, nlt_zaero, &
          nt_bgc_hum,  nlt_bgc_hum, bio_index_o, bio_index
 
      use ice_colpkg_shared, only: ktherm, shortwave, solve_zsal, &
          skl_bgc, z_tracers, scale_bgc, dEdd_algae, solve_zbgc, &
          bgc_data_dir, sil_data_type, nit_data_type, fe_data_type, &
          bgc_flux_type, grid_o, l_sk, grid_o_t, initbio_frac, &
          frazil_scav, grid_oS, l_skS, max_nbtrcr, max_algae, max_aero, &
          max_doc, max_dic, max_don, max_fe, restore_bgc, phi_snow, &
          modal_aero

      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        k            ! loop index  

      !------------------------------------------------------------
      !        Tracers have mobile and  stationary phases. 
      ! ice growth allows for retention, ice melt facilitates mobility
      ! bgc_tracer_type defines the exchange timescales between these phases
      ! -1 : entirely in the mobile phase, no exchange  (this is the default)
      !  0 : retention time scale is tau_min, release time scale is tau_max
      !  1 : retention time scale is tau_max, release time scale is tau_min
      ! 0.5: retention time scale is tau_min, release time scale is tau_min
      !  2 : retention time scale is tau_max, release time scale is tau_max
      ! tau_min and tau_max are defined in ice_zbgc_shared.f90
      !------------------------------------------------------------

      !-----------------------------------------------------------------
      ! namelist variables
      !-----------------------------------------------------------------

      namelist /zbgc_nml/  &
        tr_brine, restart_hbrine, tr_zaero, modal_aero, skl_bgc, z_tracers, &
        dEdd_algae, solve_zbgc, bgc_flux_type, &
        restore_bgc, restart_bgc, scale_bgc, solve_zsal, restart_zsal, &
        bgc_data_dir, sil_data_type, nit_data_type,  fe_data_type, &
        tr_bgc_Nit, tr_bgc_C, tr_bgc_chl, tr_bgc_Am, tr_bgc_Sil, &
        tr_bgc_DMS, tr_bgc_PON, tr_bgc_hum, tr_bgc_DON, tr_bgc_Fe, &
        grid_o, grid_o_t, l_sk, grid_oS, &   
        l_skS, phi_snow,  initbio_frac, frazil_scav

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------
      tr_brine        = .false.  ! brine height differs from ice height
      tr_zaero        = .false.  ! z aerosol tracers
      modal_aero      = .false.  ! use modal aerosol treatment of aerosols
      restore_bgc     = .false.  ! restore bgc if true
      solve_zsal      = .false.  ! update salinity tracer profile from solve_S_dt
      bgc_data_dir    = 'unknown_bgc_data_dir'
      sil_data_type   = 'default'
      nit_data_type   = 'default'
      fe_data_type    = 'default'
      restart_bgc     = .false.  ! biogeochemistry restart
      restart_zsal    = .false.  ! salinity restart
      restart_hbrine  = .false.  ! hbrine restart
      scale_bgc       = .false.  ! initial bgc tracers proportional to S  
      skl_bgc         = .false.  ! solve skeletal biochemistry 
      z_tracers       = .false.  ! solve vertically resolved tracers
      dEdd_algae      = .false.  ! dynamic algae contributes to shortwave absorption
                                 ! in delta-Eddington calculation
      solve_zbgc      = .false.  ! turn on z layer biochemistry 
      tr_bgc_PON      = .false.  !---------------------------------------------   
      tr_bgc_Nit      = .false.  ! biogeochemistry (skl or zbgc)
      tr_bgc_C        = .false.  ! if skl_bgc = .true. then skl
      tr_bgc_chl      = .false.  ! if z_tracers = .true. then vertically resolved
      tr_bgc_Sil      = .false.  ! if z_tracers + solve_zbgc = .true. then
      tr_bgc_Am       = .false.  ! vertically resolved with reactions  
      tr_bgc_DMS      = .false.  !------------------------------------------------
      tr_bgc_DON      = .false.  ! 
      tr_bgc_hum      = .false.  !
      tr_bgc_Fe       = .false.  ! 
      tr_bgc_N        = .true.   !

      ! brine height parameter
      phi_snow        = p5       ! snow porosity

      ! skl biology parameters
      bgc_flux_type   = 'Jin2006'! type of ocean-ice poston velocity ('constant')

      ! z biology parameters  
      grid_o          = c5           ! for bottom flux        
      grid_o_t        = c5           ! for top flux        
      l_sk            = 7.0_dbl_kind ! characteristic diffusive scale (m)   
      initbio_frac    = c1           ! fraction of ocean trcr concentration in bio trcrs
      frazil_scav     = c1           ! increase in initial bio tracer from ocean scavenging 

      ! z salinity  parameters
      grid_oS         = c5            ! for bottom flux         
      l_skS           = 7.0_dbl_kind  ! characteristic diffusive scale (m)  

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)

      if (my_task == master_task) then
         open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif 

         print*,'Reading zbgc_nml'
         do while (nml_error > 0)
            read(nu_nml, nml=zbgc_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice('error reading zbgc namelist')
      endif
      call release_fileunit(nu_nml)

      !-----------------------------------------------------------------
      ! zsalinity and brine
      !-----------------------------------------------------------------
      if (solve_zsal .and. TRZS == 0) then
         write(nu_diag,*) 'WARNING: solve_zsal=T but 0 zsalinity tracers'
         write(nu_diag,*) 'WARNING: setting solve_zsal = F'
         solve_zsal = .false.      
      elseif (solve_zsal .and. nblyr < 1)  then
         write(nu_diag,*) 'WARNING: solve_zsal=T but 0 zsalinity tracers'
         write(nu_diag,*) 'WARNING: setting solve_zsal = F'
         solve_zsal = .false.     
      endif 

      if (solve_zsal .and. ((.not. tr_brine) .or. (ktherm /= 1))) then
         write(nu_diag,*) 'WARNING: solve_zsal needs tr_brine=T and ktherm=1'
         write(nu_diag,*) 'WARNING: setting tr_brine=T and ktherm=1'
         tr_brine = .true.
         ktherm = 1
      endif

      if (tr_brine .and. TRBRI == 0 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no brine height compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      elseif (tr_brine .and. nblyr < 1 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no biology layers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      endif 

      call broadcast_scalar(solve_zsal,         master_task)  
      call broadcast_scalar(restart_zsal,       master_task)  
      call broadcast_scalar(tr_brine,           master_task)
      call broadcast_scalar(restart_hbrine,     master_task) 

      call broadcast_scalar(phi_snow,           master_task)
      call broadcast_scalar(grid_oS,            master_task)
      call broadcast_scalar(l_skS,              master_task)

      if (my_task == master_task) then
         write(nu_diag,1010) ' tr_brine                  = ', tr_brine
         if (tr_brine) then
         write(nu_diag,1010) ' restart_hbrine            = ', restart_hbrine
         write(nu_diag,1005) ' phi_snow                  = ', phi_snow
         endif
         if (solve_zsal) then
         write(nu_diag,1010) ' solve_zsal                = ', solve_zsal
         write(nu_diag,1010) ' restart_zsal              = ', restart_zsal
         write(nu_diag,1000) ' grid_oS                   = ', grid_oS
         write(nu_diag,1005) ' l_skS                     = ', l_skS
         endif
      endif

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      if (.not. tr_brine) then
         if (solve_zbgc) then
            write(nu_diag,*) 'WARNING: tr_brine = F and solve_zbgc = T'
            write(nu_diag,*) 'WARNING: setting solve_zbgc = F'
            solve_zbgc = .false.
         endif
         if (tr_zaero) then
            write(nu_diag,*) 'WARNING: tr_brine = F and tr_zaero = T'
            write(nu_diag,*) 'WARNING: setting tr_zaero = F'
            tr_zaero = .false.
         endif
      endif

      if ((skl_bgc .AND. solve_zbgc) .or. (skl_bgc .AND. z_tracers)) &
              call abort_ice('error:skl_bgc &
              and solve_zbgc or z_tracers are both true')

      if (skl_bgc .AND. tr_zaero) then
         write(nu_diag,*) 'WARNING: skl bgc does not use vertical tracers'
         write(nu_diag,*) 'WARNING: setting tr_zaero = F'
         tr_zaero = .false.
      endif

      if (dEdd_algae .AND. trim(shortwave) /= 'dEdd') then 
         write(nu_diag,*) 'WARNING: dEdd_algae = T but shortwave /= dEdd'
         write(nu_diag,*) 'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif

      if (dEdd_algae .AND. (.NOT. tr_bgc_N) .AND. (.NOT. tr_zaero)) then 
         write(nu_diag,*) 'WARNING: need tr_bgc_N or tr_zaero for dEdd_algae'
         write(nu_diag,*) 'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif

      if (modal_aero .AND. (.NOT. tr_zaero) .AND. (.NOT. tr_aero)) then
         modal_aero = .false.
      endif
         
      if (modal_aero .AND. trim(shortwave) /= 'dEdd') then 
         write(nu_diag,*) 'WARNING: modal_aero = T but shortwave /= dEdd'
         write(nu_diag,*) 'WARNING: setting modal_aero = F'
         modal_aero = .false.
      endif
      if (n_algae > max_algae) call abort_ice('error:number of algal &
            types exceeds max_algae')
      if (n_doc > max_doc) call abort_ice('error:number of doc &
            types exceeds max_doc')
      if (n_dic > max_doc) call abort_ice('error:number of dic &
            types exceeds max_dic')
      if (n_don > max_don) call abort_ice('error:number of don &
            types exceeds max_don')
      if (n_fed  > max_fe ) call abort_ice('error:number of dissolved fe &
            types exceeds max_fe ')
      if (n_fep  > max_fe ) call abort_ice('error:number of particulate fe &
            types exceeds max_fe ')
      if ((TRBGCS == 0 .and. skl_bgc) .or. (TRALG == 0 .and. skl_bgc)) then
         write(nu_diag,*) &
            'WARNING: skl_bgc=T but 0 bgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting skl_bgc = F'
         skl_bgc = .false.
      endif

      if ((TRBGCZ == 0 .and. solve_zbgc) .or. (TRALG == 0 .and. solve_zbgc)) then
         write(nu_diag,*) &
            'WARNING: solve_zbgc=T but 0 zbgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zbgc = F'
         solve_zbgc = .false.
      endif

      if (solve_zbgc .and. .not. z_tracers) z_tracers = .true.
      if (skl_bgc .or. solve_zbgc) then
         tr_bgc_N         = .true.   ! minimum NP biogeochemistry
         tr_bgc_Nit       = .true.
      else
         tr_bgc_N         = .false.
         tr_bgc_C         = .false.
         tr_bgc_chl       = .false.
         tr_bgc_Nit       = .false.
         tr_bgc_Am        = .false.
         tr_bgc_Sil       = .false.
         tr_bgc_hum       = .false.
         tr_bgc_DMS       = .false.
         tr_bgc_PON       = .false.
         tr_bgc_DON       = .false.
         tr_bgc_Fe        = .false.
      endif

      call broadcast_scalar(solve_zbgc,         master_task)
      call broadcast_scalar(skl_bgc,            master_task)
      call broadcast_scalar(restart_bgc,        master_task)
      call broadcast_scalar(bgc_flux_type,      master_task)
      call broadcast_scalar(restore_bgc,        master_task)
      call broadcast_scalar(bgc_data_dir,       master_task)
      call broadcast_scalar(sil_data_type,      master_task)
      call broadcast_scalar(nit_data_type,      master_task)
      call broadcast_scalar(fe_data_type,       master_task)
      call broadcast_scalar(tr_bgc_N,           master_task)
      call broadcast_scalar(tr_bgc_C,           master_task)
      call broadcast_scalar(tr_bgc_chl,         master_task)
      call broadcast_scalar(tr_bgc_Nit,         master_task)
      call broadcast_scalar(tr_bgc_Am,          master_task)
      call broadcast_scalar(tr_bgc_Sil,         master_task)
      call broadcast_scalar(tr_bgc_hum,         master_task)
      call broadcast_scalar(tr_bgc_DMS,         master_task) 
      call broadcast_scalar(tr_bgc_PON,         master_task) 
      call broadcast_scalar(tr_bgc_DON,         master_task) 
      call broadcast_scalar(tr_bgc_Fe,          master_task) 

      !-----------------------------------------------------------------
      ! z layer aerosols
      !-----------------------------------------------------------------
      if (tr_zaero .and. .not. z_tracers) z_tracers = .true.

      if (n_zaero > max_aero) call abort_ice('error:number of z aerosols &
            exceeds max_aero')
         
      call broadcast_scalar(z_tracers,          master_task)
      call broadcast_scalar(tr_zaero,           master_task)
      call broadcast_scalar(dEdd_algae,         master_task) 
      call broadcast_scalar(modal_aero,         master_task)
      call broadcast_scalar(grid_o,             master_task)
      call broadcast_scalar(grid_o_t,           master_task)
      call broadcast_scalar(l_sk,               master_task)
      call broadcast_scalar(scale_bgc,          master_task)
      call broadcast_scalar(initbio_frac,       master_task)
      call broadcast_scalar(frazil_scav,        master_task)

      if (skl_bgc .and. n_bgc < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of bgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',n_bgc
         call abort_ice ('ice_zbgc error: skl_bgc and n_bgc < 2')
      endif

      if (solve_zbgc .and. n_bgc < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of zbgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',n_bgc
         call abort_ice ('ice_zbgc error: solve_zbgc and n_bgc < 2')
      endif

      if (tr_zaero .and. TRZAERO <  1) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of TRZAERO > 0'
         write (nu_diag,*) 'in order to solve z aerosols:',TRZAERO
         call abort_ice ('ice_zbgc error: tr_zaero and tr zaero < 1')
      endif

      !-----------------------------------------------------------------
      ! initialize tracers etc in the column package
      !-----------------------------------------------------------------
      call colpkg_init_zbgc (nblyr, nilyr, nslyr, &
                 n_algae, n_zaero, n_doc, n_dic, n_don, n_fed, n_fep, &
                 trcr_base, trcr_depend, n_trcr_strata, nt_strata, nbtrcr_sw, &
                 tr_brine, nt_fbri, ntrcr, nbtrcr, nt_bgc_Nit, nt_bgc_Am, &
                 nt_bgc_Sil, nt_bgc_DMS, nt_bgc_PON, nt_bgc_S, nt_bgc_N, &
                 nt_bgc_C, nt_bgc_chl, nt_bgc_DOC, nt_bgc_DON, nt_bgc_DIC, & 
                 nt_zaero, nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_Fed, nt_bgc_Fep, &
                 nt_zbgc_frac, tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS, &
                 tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C, tr_bgc_chl, &
                 tr_bgc_DON, tr_bgc_Fe, tr_zaero, nlt_zaero_sw, nlt_chl_sw, &
                 nlt_bgc_N, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, &
                 nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
                 nlt_bgc_C, nlt_bgc_chl, nlt_bgc_DIC, nlt_bgc_DOC, &
                 nlt_bgc_PON, nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, &
                 nlt_zaero, &
                 nt_bgc_hum, nlt_bgc_hum, tr_bgc_hum, solve_zsal, &
                 skl_bgc, z_tracers, dEdd_algae, solve_zbgc, &
                 frazil_scav, initbio_frac, bio_index_o, bio_index, ntrcr_o, &
                 max_algae, max_doc, max_dic, max_don, max_fe)

      !-----------------------------------------------------------------
      ! final consistency checks
      !----------------------------------------------------------------- 
      if (nbtrcr > max_nbtrcr) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr > max_nbtrcr'
         write (nu_diag,*) 'nbtrcr, max_nbtrcr:',nbtrcr, max_nbtrcr
         call abort_ice ('ice_zbgc error: nbtrcr > max_nbtrcr')
      endif
      if (.NOT. dEdd_algae) nbtrcr_sw = 1

      if (nbtrcr_sw > max_nsw) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr_sw > max_nsw'
         write (nu_diag,*) 'nbtrcr_sw, max_nsw:',nbtrcr_sw, max_nsw
         call abort_ice ('ice_zbgc error: nbtrcr_sw > max_nsw')
      endif

      if (ntrcr > max_ntrcr) then
         write(nu_diag,*) 'max_ntrcr < number of namelist tracers'
         write(nu_diag,*) 'max_ntrcr = ',max_ntrcr,' ntrcr = ',ntrcr
         call abort_ice('max_ntrcr < number of namelist tracers')
      endif                               

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------
      if (my_task == master_task) then
      if (skl_bgc) then

         write(nu_diag,1010) ' skl_bgc                   = ', skl_bgc
         write(nu_diag,1030) ' bgc_flux_type             = ', bgc_flux_type
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' restore_bgc               = ', restore_bgc
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,*)    ' fe_data_type              = ', &
                               trim(fe_data_type)
         write(nu_diag,1020) ' number of bio tracers     = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc          = ', n_doc
         write(nu_diag,1020) ' number of dic          = ', n_dic
         write(nu_diag,1020) ' number of don          = ', n_don
         write(nu_diag,1020) ' number of fed          = ', n_fed
         write(nu_diag,1020) ' number of fep          = ', n_fep
         write(nu_diag,1010) ' tr_bgc_N               = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_C               = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_chl             = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_Nit             = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_Am              = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_Sil             = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_hum             = ', tr_bgc_hum
         write(nu_diag,1010) ' tr_bgc_DMS             = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON             = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON             = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe              = ', tr_bgc_Fe 
        
      elseif (z_tracers) then

         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,*)    ' fe_data_type              = ', &
                               trim(fe_data_type)
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' dEdd_algae                = ', dEdd_algae  
         write(nu_diag,1010) ' modal_aero                = ', modal_aero  
         write(nu_diag,1010) ' scale_bgc                 = ', scale_bgc
         write(nu_diag,1010) ' solve_zbgc                = ', solve_zbgc
         write(nu_diag,1020) ' number of ztracers        = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc             = ', n_doc
         write(nu_diag,1020) ' number of dic             = ', n_dic
         write(nu_diag,1020) ' number of fed             = ', n_fed
         write(nu_diag,1020) ' number of fep             = ', n_fep
         write(nu_diag,1020) ' number of aerosols        = ', n_zaero
         write(nu_diag,1010) ' tr_zaero                  = ', tr_zaero
         write(nu_diag,1010) ' tr_bgc_Nit                = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_N                  = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_Am                 = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_C                  = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_Sil                = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_hum                = ', tr_bgc_hum
         write(nu_diag,1010) ' tr_bgc_chl                = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_DMS                = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON                = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON                = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe                 = ', tr_bgc_Fe 
         !bio parameters
         write(nu_diag,1000) ' grid_o                    = ', grid_o
         write(nu_diag,1000) ' grid_o_t                  = ', grid_o_t
         write(nu_diag,1005) ' l_sk                      = ', l_sk
         write(nu_diag,1000) ' initbio_frac              = ', initbio_frac
         write(nu_diag,1000) ' frazil_scav               = ', frazil_scav  

      endif  ! skl_bgc or solve_bgc
      endif  ! master_task

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character

      end subroutine init_zbgc

!=======================================================================

      end module ice_init_column

!=======================================================================

