!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module forcing_fields

!BOP
! !MODULE: forcing_fields

! !DESCRIPTION:
!  Contains the forcing fields necessary for supporting high-level coupling.
!  These fields originally resided in modules forcing and forcing_coupled.

! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use kinds_mod
   use blocks,      only: nx_block, ny_block
   use constants,   only: c0
   use domain_size, only: max_blocks_clinic,nt
      
   implicit none
   save

!EOP
!BOC
! !PUBLIC DATA MEMBERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),public ::  &
      EVAP_F = c0,       &! evaporation   flux    from cpl (kg/m2/s)
      PREC_F = c0,       &! precipitation flux    from cpl (kg/m2/s)
                          ! (rain + snow)
      SNOW_F = c0,       &! snow          flux    from cpl (kg/m2/s)
      MELT_F = c0,       &! melt          flux    from cpl (kg/m2/s)
      ROFF_F = c0,       &! river runoff  flux    from cpl (kg/m2/s)
      IOFF_F = c0,       &! ice   runoff  flux    from cpl (kg/m2/s)
      SALT_F = c0,       &! salt          flux    from cpl (kg(salt)/m2/s)
      SENH_F = c0,       &! sensible heat flux    from cpl (W/m2   )
      LWUP_F = c0,       &! longwave heat flux up from cpl (W/m2   )
      LWDN_F = c0,       &! longwave heat flux dn from cpl (W/m2   )
      MELTH_F= c0         ! melt     heat flux    from cpl (W/m2   )


   integer(kind=int_kind), public :: &
      ATM_CO2_PROG_nf_ind = 0, & ! bottom atm level prognostic co2
      ATM_CO2_DIAG_nf_ind = 0, & ! bottom atm level diagnostic co2
      sflux_diat_nf_ind   = 0, & ! ice -> ocn diat flux
      sflux_sp_nf_ind     = 0, & ! ice -> ocn sp flux
      sflux_phaeo_nf_ind  = 0, & ! ice -> ocn phaeo flux
      sflux_dFe_nf_ind    = 0, & ! ice -> ocn fed flux
      sflux_NO3_nf_ind    = 0, & ! ice -> ocn NO3 flux
      sflux_NH4_nf_ind    = 0, & ! ice -> ocn NH4 flux
      sflux_SiO3_nf_ind   = 0, & ! ice -> ocn SiO3 flux
      sflux_DOC_nf_ind    = 0, & ! ice -> ocn DOC flux
      sflux_DON_nf_ind    = 0, & ! ice -> ocn DON flux
      sflux_DONr_nf_ind   = 0, & ! ice -> ocn DONr flux
      sflux_idms_nf_ind   = 0, & ! ice -> ocn idms flux
      sflux_idmsp_nf_ind  = 0, & ! ice -> ocn idmsp flux
      sflux_dic1_nf_ind   = 0, & ! ice -> ocn dic1 flux
      sflux_doc2_nf_ind   = 0, & ! ice -> ocn doc2 flux
      sflux_doc3_nf_ind   = 0, & ! ice -> ocn doc3 flux
      sflux_fed2_nf_ind   = 0, & ! ice -> ocn fed2 flux
      sflux_fep1_nf_ind   = 0, & ! ice -> ocn fep1 flux
      sflux_fep2_nf_ind   = 0, & ! ice -> ocn fep2 flux
      sflux_dust_nf_ind   = 0, & ! ice -> ocn dust flux
      sflux_dmspp_nf_ind  = 0    ! ice -> ocn dmspp flux


   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
      public, target :: &
      SMF,  &!  surface momentum fluxes (wind stress)
      SMFT   !  surface momentum fluxes on T points if avail

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      public, target :: &
      STF,  &!  surface tracer fluxes
      TFW    ! tracer content in freshwater flux


   logical (log_kind), public :: &
      lsmft_avail   ! true if SMFT is an available field

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      public, target ::  &
      IFRAC,             &! ice fraction; not initialized in this routine
      U10_SQR,           &! 10m wind speed squared; not initialized in this routine
      ATM_PRESS,         &! atmospheric pressure forcing
      FW,FW_OLD           ! freshwater flux at T points (cm/s)
                          ! FW_OLD is at time n-1


!***********************************************************************

 end module forcing_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
