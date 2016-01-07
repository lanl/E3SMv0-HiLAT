!  SVN:$Id: $
!=======================================================================

! Flux variable declarations for biogeochemistry
!
! author Elizabeth C. Hunke, LANL
!
      module ice_flux_bgc

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_colpkg_shared, only: max_aero, max_nbtrcr, &
                max_algae, max_doc, max_don, max_dic, max_fe
      use ice_domain_size, only: max_blocks, ncat

      implicit none
      private

      save

      ! in from atmosphere

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_aero,max_blocks), public :: &
         faero_atm   ! aerosol deposition rate (kg/m^2 s)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_nbtrcr,max_blocks), public :: &
         flux_bio_atm  ! all bio fluxes to ice from atmosphere

      ! in from ocean

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_aero,max_blocks), public :: &
         faero_ocn   ! aerosol flux to ocean  (kg/m^2/s)

      ! out to ocean 

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_nbtrcr,max_blocks), public :: &
         flux_bio   , & ! all bio fluxes to ocean
         flux_bio_ai    ! all bio fluxes to ocean, averaged over grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         fzsal_ai, & ! salt flux to ocean from zsalinity (kg/m^2/s) 
         fzsal_g_ai  ! gravity drainage salt flux to ocean (kg/m^2/s) 

      ! internal

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         hin_old     , & ! old ice thickness
         dsnown          ! change in snow thickness in category n (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         nit        , & ! ocean nitrate (mmol/m^3)          
         amm        , & ! ammonia/um (mmol/m^3)
         sil        , & ! silicate (mmol/m^3)
         dmsp       , & ! dmsp (mmol/m^3)
         dms        , & ! dms (mmol/m^3)
         hum            ! humic material (mmol/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_algae, max_blocks), public :: &
         algalN         ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_doc, max_blocks), public :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_don, max_blocks), public :: &
         don             ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_dic, max_blocks), public :: &
         dic             ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_fe, max_blocks), public :: &
         fed, fep        ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_aero, max_blocks), public :: &
         zaeros          ! ocean aerosols (mmol/m^3) 

!=======================================================================

      end module ice_flux_bgc

!=======================================================================
