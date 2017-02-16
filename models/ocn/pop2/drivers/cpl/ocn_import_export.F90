module ocn_import_export

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod
   use POP_IOUnitsMod
   use POP_MCT_vars_mod
   use POP_CplIndices

   use seq_flds_mod
   use seq_timemgr_mod
   use shr_file_mod 
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod

   use perf_mod
   use kinds_mod,         only: int_kind, r8
   use communicate,       only: my_task, master_task
   use constants
   use blocks
   use domain,            only: distrb_clinic, POP_haloClinic
   use exit_mod
   use forcing_shf,       only: SHF_QSW
   use forcing_sfwf,      only: lsend_precip_fact, precip_fact
   use forcing_fields
   use forcing_coupled,   only: ncouple_per_day,  &
                                update_ghost_cells_coupler_fluxes, &
                                rotate_wind_stress, pop_set_coupled_forcing, &
                                pop_init_coupled,  &
                                orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp
   use ice,               only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, &
                                tlast_ice
   use global_reductions, only: global_sum_prod
   use io_tools,          only: document
   use named_field_mod,   only: named_field_register, named_field_get_index, &
                                named_field_set, named_field_get
   use prognostic
   use time_management
   use registry

   implicit none
   public
   save

   ! accumulated sum of send buffer quantities for averaging before being sent
   real (r8), dimension(:,:,:,:), allocatable ::  SBUFF_SUM 

   real (r8) :: tlast_coupled 

contains

!***********************************************************************
!BOP
! !IROUTINE: ocn_import
! !INTERFACE:

  subroutine ocn_import(x2o, ldiag_cpl, errorCode)

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
! 
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction
!    o  rofl   -- river runoff flux                        (kg/m2/s)
!    o  rofi   -- ice runoff flux                          (kg/m2/s)
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  pslv   -- sea-level pressure                       (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
! 
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    real(r8)           , intent(inout) :: x2o(:,:)
    logical (log_kind) , intent(in)    :: ldiag_cpl
    integer (POP_i4)   , intent(out)   :: errorCode  ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) ::   &
      label,                 &
      message
 
   integer (int_kind) ::  &
      i,j,k,n,iblock

   real (r8), dimension(nx_block,ny_block) ::  &
      WORKB

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
      WORK1, WORK2        ! local work space

   real (r8) ::  &
      m2percm2,  &
      gsum

   type (block) :: this_block ! local block info

   integer (int_kind) :: nrecv

!-----------------------------------------------------------------------
!
!  zero out padded cells 
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   WORK1 = c0
   WORK2 = c0

!-----------------------------------------------------------------------
!
!  unpack and distribute wind stress, then convert to correct units
!  and rotate components to local coordinates
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         WORK1(i,j,iblock) = x2o(index_x2o_Foxx_taux,n)
         WORK2(i,j,iblock) = x2o(index_x2o_Foxx_tauy,n)
      enddo
      enddo
   enddo ! iblock

   !***
   !*** do NOT perform halo updates now, because vector updates must
   !***   be done after the rotation is completed.
   !***

!-----------------------------------------------------------------------
!
!  rotate true zonal/meridional wind stress into local coordinates,
!  convert to dyne/cm**2, and shift SMFT to U grid
!
!  halo updates are performed in subroutine rotate_wind_stress, 
!  following the rotation
!
!-----------------------------------------------------------------------

      call rotate_wind_stress(WORK1, WORK2)

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

!-----------------------------------------------------------------------
!
!  unpack and distribute fresh water flux and salt flux
!
!  NOTE: if there are code changes associated with changing the names or
!        the number of fluxes received from the coupler, then subroutine
!        update_ghost_cells_coupler_fluxes will need to be modified also
!
!-----------------------------------------------------------------------


      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         SNOW_F(i,j,iblock) = x2o(index_x2o_Faxa_snow,n)
         WORKB (i,j       ) = x2o(index_x2o_Faxa_rain,n)
         EVAP_F(i,j,iblock) = x2o(index_x2o_Foxx_evap,n)
         MELT_F(i,j,iblock) = x2o(index_x2o_Fioi_meltw,n)
         ROFF_F(i,j,iblock) = x2o(index_x2o_Foxx_rofl,n)
         IOFF_F(i,j,iblock) = x2o(index_x2o_Foxx_rofi,n)
         SALT_F(i,j,iblock) = x2o(index_x2o_Fioi_salt,n)

         PREC_F(i,j,iblock) = WORKB(i,j) + SNOW_F(i,j,iblock)    ! rain + snow

         WORKB(i,j        ) = x2o(index_x2o_Foxx_swnet,n)
         SHF_QSW(i,j,iblock) = WORKB(i,j)*  &
            RCALCT(i,j,iblock)*hflux_factor  !  convert from W/m**2
         SENH_F(i,j,iblock)  = x2o(index_x2o_Foxx_sen,n)
         LWUP_F(i,j,iblock)  = x2o(index_x2o_Foxx_lwup,n)
         LWDN_F(i,j,iblock)  = x2o(index_x2o_Faxa_lwdn,n)
         MELTH_F(i,j,iblock) = x2o(index_x2o_Fioi_melth,n)

         WORKB(i,j       ) = x2o(index_x2o_Si_ifrac,n)
         IFRAC(i,j,iblock) = WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from Pa to dynes/cm**2
         WORKB(i,j       ) = x2o(index_x2o_Sa_pslv,n)
         ATM_PRESS(i,j,iblock) = c10 * WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from m**2/s**2 to cm**2/s**2
         WORKB(i,j       ) = x2o(index_x2o_So_duu10n,n)
         U10_SQR(i,j,iblock) = cmperm * cmperm * WORKB(i,j) * RCALCT(i,j,iblock)

      enddo
      enddo

   enddo

!-----------------------------------------------------------------------
!
!  incoming data quality control
!
!-----------------------------------------------------------------------
#ifdef CCSMCOUPLED
      if ( any(IOFF_F < c0) ) then
        write(message, "(A,1x,e10.3,A)") 'Error: incoming IOFF_F has min value', &
                                      minval(IOFF_F), '; value can not be negative.'
        ! call shr_sys_abort ('Error: incoming IOFF_F is negative')
        call shr_sys_abort (trim(message))
      endif
#endif


!-----------------------------------------------------------------------
!
!  update ghost cells for fluxes received from the coupler
!
!-----------------------------------------------------------------------

   call update_ghost_cells_coupler_fluxes(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ocn_import: error in update_ghost_cells_coupler_fluxes')
      return
   endif

!-----------------------------------------------------------------------
!
!  unpack atmospheric CO2
!
!-----------------------------------------------------------------------

   if (index_x2o_Sa_co2prog > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Sa_co2prog,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating PROG CO2 halo')
         return
      endif

      call named_field_set(ATM_CO2_PROG_nf_ind, WORK1)
   endif

   if (index_x2o_Sa_co2diag > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Sa_co2diag,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating DIAG CO2 halo')
         return
      endif

      call named_field_set(ATM_CO2_DIAG_nf_ind, WORK1)
   endif

!-----------------------------------------------------------------------
!
!  unpack sea ice fluxes (swang)
!
!-----------------------------------------------------------------------

   if (index_x2o_Fioi_diat > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_diat,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_diat halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s, as POP needs the latter
      call named_field_set(sflux_diat_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_sp > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_sp,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_sp halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_sp_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_phaeo > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_phaeo,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_phaeo halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_phaeo_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_fed > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_fed,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_dFe halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_dFe_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_no3 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_no3,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_NO3 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_NO3_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_nh4 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_nh4,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_NH4 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_NH4_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_sio3 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_sio3,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_SiO3 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_SiO3_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_doc > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_doc,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_DOC halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_DOC_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_don > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_don,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_DON halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_DON_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_donr > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_donr,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_DONr halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_DONr_nf_ind, WORK1 * cmperm)
   endif


   if (index_x2o_Fioi_dms > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_dms,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_idms halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_idms_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_dmsp > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_dmsp,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_idmsp halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_idmsp_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_dic1 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_dic1,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_dic1 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_dic1_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_doc2 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_doc2,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_doc2 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_doc2_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_doc3 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_doc3,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_doc3 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_doc3_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_dmspp > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_dmspp,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_dmspp halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_dmspp_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_fed2 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_fed2,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_fed2 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_fed2_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_fep1 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_fep1,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_fep1 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_fep1_nf_ind, WORK1 * cmperm)
   endif
   
      if (index_x2o_Fioi_fep2 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_fep2,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_fep2 halo')
         return
      endif
! convert flux unit from mmol/m^2/s to mmol/m^3 cm/s
      call named_field_set(sflux_fep2_nf_ind, WORK1 * cmperm)
   endif

   if (index_x2o_Fioi_dust > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Fioi_dust,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating sflux_dust halo')
         return
      endif
! NOTE: Current dust is in kg/m^2/s, 'cmperm' should be changed to the correct
! factor (swang)
      call named_field_set(sflux_dust_nf_ind, WORK1 * cmperm)
   endif
 
!-----------------------------------------------------------------------
!
!  diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then

     write(message,'(6a,1x,5a)')  &
         ' Global averages of fluxes received from cpl at ',  &
           cyear,'/',cmonth ,'/',cday,  chour,':',cminute,':',csecond
     call document ('pop_recv_from_coupler', trim(message))
 
     m2percm2  = mpercm*mpercm
     nrecv = size(x2o, dim=1)
     do k = 1,nrecv

         n = 0
         do iblock = 1, nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)

            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               n = n + 1
               WORK1(i,j,iblock) = x2o(k,n)  ! mult. by TAREA in global_sum_prod
            enddo
            enddo
         enddo

         gsum = global_sum_prod(WORK1 , TAREA, distrb_clinic, &
                                 field_loc_center, RCALCT)*m2percm2
         if (my_task == master_task) then
            call seq_flds_getField(label,k,seq_flds_x2o_fields)
            write(stdout,1100)'ocn','recv', label ,gsum
            call shr_sys_flush(stdout)
         endif
      enddo
   endif


1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import

!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export(o2x, ldiag_cpl, errorCode)   

! !DESCRIPTION:
!  This routine calls the routines necessary to send pop fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

   real(r8)           , intent(inout) :: o2x(:,:)
   logical (log_kind) , intent(in)    :: ldiag_cpl
   integer (POP_i4)   , intent(out)   :: errorCode  ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n, iblock
           
   character (char_len)    :: label
 
   integer (int_kind) ::  &
        i,j,k

   real (r8), dimension(nx_block,ny_block) ::   &
        WORK1, WORK2,      &! local work space
        WORK3, WORK4

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
        WORKA               ! local work space with full block dimension

   real (r8) ::   &
      m2percm2,   &
      gsum

   type (block) :: this_block ! local block info

   integer (int_kind) :: nsend

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  initialize control buffer
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points and rotate on T grid
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)

     call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2x_So_u),iblock)
     call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2x_So_v),iblock)

     WORK1 = (WORK3*cos(ANGLET(:,:,iblock))+WORK4*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled
     WORK2 = (WORK4*cos(ANGLET(:,:,iblock))-WORK3*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled

     do j=this_block%jb,this_block%je
     do i=this_block%ib,this_block%ie
        n = n + 1
        o2x(index_o2x_So_u,n) = WORK1(i,j)
        o2x(index_o2x_So_v,n) = WORK2(i,j)
     enddo
     enddo
  enddo

!-----------------------------------------------------------------------
!
!     convert and pack surface temperature
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_So_t,n) =   &
             SBUFF_SUM(i,j,iblock,index_o2x_So_t)/tlast_coupled + T0_Kelvin
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     convert and pack salinity
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_So_s,n) =   &
             SBUFF_SUM(i,j,iblock,index_o2x_So_s)*salt_to_ppt/tlast_coupled
!----- if ocean salinity is negative then send zero:
         o2x(index_o2x_So_s,n) = max(o2x(index_o2x_So_s,n), c0)
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points, then rotate on T grid
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx),iblock)
      call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy),iblock)
 
      WORK1 = (WORK3*cos(ANGLET(:,:,iblock)) + WORK4*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled
      WORK2 = (WORK4*cos(ANGLET(:,:,iblock)) - WORK3*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_So_dhdx,n) = WORK1(i,j)
         o2x(index_o2x_So_dhdy,n) = WORK2(i,j)
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     pack heat flux due to freezing/melting (W/m^2)
!     QFLUX computation and units conversion occurs in ice.F
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_Fioo_q,n) = QFLUX(i,j,iblock)
      enddo
      enddo
   enddo

   tlast_ice = c0
   AQICE     = c0
   QICE      = c0

!-----------------------------------------------------------------------
!
!     pack co2 flux, if requested (kg CO2/m^2/s)
!     units conversion occurs where co2 flux is computed
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fco2_ocn > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_Faoo_fco2_ocn,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_Faoo_fco2_ocn)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!
!     pack dms flux, if requested (kg dms/m^2/s)
!     units conversion occurs where dms flux is computed
!
!-----------------------------------------------------------------------

!maltrud debug
!   write(stdout,*) 'pjc: ocn_export_mct: ', my_task, index_o2x_Faoo_fdms_ocn

   if (index_o2x_Faoo_fdms_ocn > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_Faoo_fdms_ocn,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_Faoo_fdms_ocn)/tlast_coupled * &
               62.132400_r8 *1D-8     !pjc, convert from (mmol/m3)(cm/s) to kg(DMS)/m2/s, using M(DMS)=62.1324 g/mol
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_dms, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

!   write(stdout,*) 'swang: ocn_export_mct: ', my_task, index_o2x_So_dms

   if (index_o2x_So_dms > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_dms,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_dms)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_dmsp, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_dmsp > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_dmsp,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_dmsp)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_diat, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_diat > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_diat,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_diat)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_sp, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_sp > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_sp,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_sp)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_phaeo, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_phaeo > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_phaeo,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_phaeo)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_dFe, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_fed > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_fed,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_fed)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_NH4, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_nh4 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_nh4,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_nh4)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_NO3, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_no3 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_no3,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_no3)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_SiO3, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_sio3 > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_sio3,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_sio3)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_DIC, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_dic > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_dic,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_dic)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_DOC, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_doc > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_doc,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_doc)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_DON, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_don > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_don,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_don)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!     pack surf_DONr, NOTE: unit is mmol/m^3 (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_donr > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_So_donr,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_So_donr)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!
!     diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then
      call ccsm_char_date_and_time
      !DEBUG      write(message,'(6a,1x,5a)')' Global averages of fluxes sent to cpl at ', &
      !DEBUG           cyear,'/',cmonth, '/',cday,  chour,':',cminute,':',csecond
      !DEBUG      call document ('pop_send_to_coupler', message)
      write(stdout,*)'pop_send_to_coupler'

      m2percm2  = mpercm*mpercm
      nsend = size(o2x,dim=1)
      do k = 1,nsend
        n = 0
        do iblock = 1, nblocks_clinic
           this_block = get_block(blocks_clinic(iblock),iblock)
           do j=this_block%jb,this_block%je
           do i=this_block%ib,this_block%ie
              n = n + 1
              WORKA(i,j,iblock) = o2x(k,n)
           enddo
           enddo
        enddo

        call POP_HaloUpdate(WORKA,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)
       
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'ocn_export_mct: error updating halo for state')
            return
         endif

        gsum = global_sum_prod(WORKA , TAREA, distrb_clinic, &
                                   field_loc_center, RCALCT)*m2percm2
        if (my_task == master_task) then
           call seq_flds_getField(label,k,seq_flds_o2x_fields)
           write(stdout,1100)'ocn','send', label ,gsum
        endif
      enddo ! k
      if (my_task == master_task) call shr_sys_flush(stdout)
   endif

1100 format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)
    
    tlast_coupled = c0

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_export

!***********************************************************************

!BOP
! !IROUTINE: POP_sum_buffer
! !INTERFACE:

 subroutine POP_sum_buffer

! !DESCRIPTION:
!  This routine accumulates sums for averaging fields to
!  be sent to the coupler
!
! !REVISION HISTORY:
!  same as module
! 
!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK                ! local work arrays

   real (r8) ::   &
      delt,             & ! time interval since last step
      delt_last           ! time interval for previous step

   integer (int_kind) :: &
      iblock,           & ! block index
      sflux_co2_nf_ind = 0, & ! named field index of fco2
      sflux_dms_nf_ind = 0, & ! named field index of fDMS
      surf_dms_nf_ind  = 0, & ! named field index of DMS
      surf_dmsp_nf_ind = 0, & ! named field index of DMSP
      surf_dFe_nf_ind  = 0, & ! named field index of dissolved iron
      surf_NH4_nf_ind  = 0, & ! named field index of NH4
      surf_NO3_nf_ind  = 0, & ! named field index of NO3
      surf_SiO3_nf_ind = 0, & ! named field index of SiO3
      surf_DIC_nf_ind  = 0, & ! named field index of DIC
      surf_DOC_nf_ind  = 0, & ! named field index of DOC
      surf_DON_nf_ind  = 0, & ! named field index of DON      
      surf_DONr_nf_ind = 0, & ! named field index of DONr 
      surf_diat_nf_ind = 0, & ! named field index of diatom C
      surf_sp_nf_ind   = 0, & ! named field index of small phyto C
      surf_phaeo_nf_ind= 0    ! named field index of phaeocystis C      

   logical (log_kind) :: &
      first = .true.      ! only true for first call

   save first

!-----------------------------------------------------------------------
!
!  zero buffer if this is the first time after a coupling interval
!
!-----------------------------------------------------------------------

   if (tlast_coupled == c0) SBUFF_SUM = c0
   WORK = c0

!-----------------------------------------------------------------------
!
!  update time since last coupling
!
!-----------------------------------------------------------------------

   if (avg_ts .or. back_to_back) then
      delt = p5*dtt
   else
      delt =    dtt
   endif
   tlast_coupled = tlast_coupled + delt

!-----------------------------------------------------------------------
!
!  allow for fco2 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fco2_ocn > 0) then
      if (sflux_co2_nf_ind == 0) then
         call named_field_get_index('SFLUX_CO2', sflux_co2_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!-----------------------------------------------------------------------
!
!  same as above for DMS
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fdms_ocn > 0) then
      if (sflux_dms_nf_ind == 0) then
         call named_field_get_index('SFLUX_DMS', sflux_dms_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!maltrud debug
!  write(stdout,*) 'pjc: pop_sum_buffer1: ', my_task, &
!    index_o2x_Faoo_fdms_ocn, sflux_dms_nf_ind, delt_last

!-----------------------------------------------------------------------
!  allow for surf_dms field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_dms > 0) then
      if (surf_dms_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceDMS', surf_dms_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!maltrud debug
!  write(stdout,*) 'swang: pop_sum_buffer1: ', my_task, &
!    index_o2x_So_dms, surf_dms_nf_ind, delt_last

!-----------------------------------------------------------------------
!  same as above for DMSP
!-----------------------------------------------------------------------

   if (index_o2x_So_dmsp > 0) then
      if (surf_dmsp_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceDMSP', surf_dmsp_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!-----------------------------------------------------------------------
!  allow for surf_diat field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_diat > 0) then
      if (surf_diat_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceDiatomCarbon',surf_diat_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_sp field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_sp > 0) then
      if (surf_sp_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceSmallPhytoCarbon', surf_sp_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_phaeo field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_phaeo > 0) then
      if (surf_phaeo_nf_ind == 0) then
         call named_field_get_index('oceanSurfacePhaeoCarbon', surf_phaeo_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!
!  allow for surf_Fe field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep
!       (SW, copied index_o2x_Faoo_fco2 and sflux_co2_nf_ind
!-----------------------------------------------------------------------

   if (index_o2x_So_fed > 0) then
      if (surf_dFe_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceFeBioavailable', surf_dFe_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_NH4 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_nh4 > 0) then
      if (surf_NH4_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceNH4', surf_NH4_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_NO3 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_no3 > 0) then
      if (surf_NO3_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceNO3', surf_NO3_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_SiO3 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_sio3 > 0) then
      if (surf_SiO3_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceSiO3', surf_SiO3_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_DIC field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_dic > 0) then
      if (surf_DIC_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceDIC', surf_DIC_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_DOC field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_doc > 0) then
      if (surf_DOC_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceDOC', surf_DOC_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif
!-----------------------------------------------------------------------
!  allow for surf_DON field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep (SW)
!-----------------------------------------------------------------------

   if (index_o2x_So_don > 0) then
      if (surf_DON_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceDON', surf_DON_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!-----------------------------------------------------------------------
!  same as above 
!-----------------------------------------------------------------------

   if (index_o2x_So_donr > 0) then
      if (surf_DONr_nf_ind == 0) then
         call named_field_get_index('oceanSurfaceDONr', surf_DONr_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!-----------------------------------------------------------------------
!
!  accumulate sums of U,V,T,S and GRADP
!  accumulate sum of co2 flux, if requested
!     implicitly use zero flux if fco2 field not registered yet
!  ice formation flux is handled separately in ice routine
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
   SBUFF_SUM(:,:,iblock,index_o2x_So_u) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_u) + delt*  &
                                   UVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_v) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_v) + delt*  &
                                   VVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_t ) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_t ) + delt*  &
                                   TRACER(:,:,1,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_s ) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_s ) + delt*  &
                                   TRACER(:,:,1,2,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx) + delt*  &
                                   GRADPX(:,:,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy) + delt*  &
                                   GRADPY(:,:,curtime,iblock)

   if (index_o2x_Faoo_fco2_ocn > 0 .and. sflux_co2_nf_ind > 0) then
      call named_field_get(sflux_co2_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fco2_ocn) = &
         SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fco2_ocn) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_Faoo_fdms_ocn > 0 .and. sflux_dms_nf_ind > 0) then
      call named_field_get(sflux_dms_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fdms_ocn) = &
         SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fdms_ocn) + delt_last*WORK(:,:,iblock)

!maltrud debug
!  write(stdout,*) 'pjc: pop_sum_buffer2: ', my_task,  &
!    minval(WORK(:,:,iblock)), maxval(WORK(:,:,iblock))

   endif

   if (index_o2x_So_fed > 0 .and. surf_dFe_nf_ind > 0) then
      print *,'surf_dFe_nf_ind',surf_dFe_nf_ind
      call named_field_get(surf_dFe_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_fed) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_fed) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_nh4 > 0 .and. surf_NH4_nf_ind > 0) then
      print *,'surf_NH4_nf_ind',surf_NH4_nf_ind
      call named_field_get(surf_NH4_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_nh4) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_nh4) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_no3 > 0 .and. surf_NO3_nf_ind > 0) then
      print *,'surf_NO3_nf_ind',surf_NO3_nf_ind
      call named_field_get(surf_NO3_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_no3) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_no3) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_sio3 > 0 .and. surf_SiO3_nf_ind > 0) then
      print *,'surf_SiO3_nf_ind',surf_SiO3_nf_ind
      call named_field_get(surf_SiO3_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_sio3) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_sio3) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_dic > 0 .and. surf_DIC_nf_ind > 0) then
      print *,'surf_DIC_nf_ind',surf_DIC_nf_ind
      call named_field_get(surf_DIC_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_dic) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_dic) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_doc > 0 .and. surf_DOC_nf_ind > 0) then
      print *,'surf_DOC_nf_ind',surf_DOC_nf_ind
      call named_field_get(surf_DOC_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_doc) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_doc) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_don > 0 .and. surf_DON_nf_ind > 0) then
      print *,'surf_DON_nf_ind',surf_DON_nf_ind
      call named_field_get(surf_DON_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_don) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_don) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_donr > 0 .and. surf_DONr_nf_ind > 0) then
      print *,'surf_DONr_nf_ind',surf_DONr_nf_ind
      call named_field_get(surf_DONr_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_donr) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_donr) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_dms > 0 .and. surf_dms_nf_ind > 0) then
      print *,'surf_dms_nf_ind',surf_dms_nf_ind
      call named_field_get(surf_dms_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_dms) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_dms) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_dmsp > 0 .and. surf_dmsp_nf_ind > 0) then
      print *,'surf_dmsp_nf_ind',surf_dmsp_nf_ind
      call named_field_get(surf_dmsp_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_dmsp) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_dmsp) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_diat > 0 .and. surf_diat_nf_ind > 0) then
      print *,'surf_diat_nf_ind',surf_diat_nf_ind
      call named_field_get(surf_diat_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_diat) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_diat) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_sp > 0 .and. surf_sp_nf_ind > 0) then
      print *,'surf_sp_nf_ind',surf_sp_nf_ind
      call named_field_get(surf_sp_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_sp) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_sp) + delt_last*WORK(:,:,iblock)
   endif

   if (index_o2x_So_phaeo > 0 .and. surf_phaeo_nf_ind > 0) then
      print *,'surf_phaeo_nf_ind',surf_phaeo_nf_ind
      call named_field_get(surf_phaeo_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_So_phaeo) = &
         SBUFF_SUM(:,:,iblock,index_o2x_So_phaeo) + delt_last*WORK(:,:,iblock)
   endif


   enddo
   !$OMP END PARALLEL DO

   first = .false.

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_sum_buffer
 
!***********************************************************************

end module ocn_import_export
