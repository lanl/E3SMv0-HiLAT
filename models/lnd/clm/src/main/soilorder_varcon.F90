module soilorder_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use abortutils  , only : endrun
  use clm_varpar  , only : nsoilorder
  use clm_varctl  , only : iulog, use_cndv, use_vertsoilc
  !
  ! !PUBLIC TYPES:
  implicit none
  save

!
! Soil order  constants
!
  character(len=40) soilordername(0:nsoilorder)    ! soil order description
  integer :: Water
  integer :: Andisols
  integer :: Gelisols
  integer :: Histosols
  integer :: Entisols
  integer :: Inceptisols
  integer :: Aridsols
  integer :: Vertisols
  integer :: Mollisols
  integer :: Alfisols
  integer :: Spodosols
  integer :: Ultisols
  integer :: Oxisols
  integer :: Shifting_sand
  integer :: rock_land
  integer :: Ice_Glacier

  real(r8), allocatable :: smax(:)
  real(r8), allocatable :: ks_sorption(:)
  real(r8), allocatable :: r_weather(:)
  real(r8), allocatable :: r_adsorp(:)
  real(r8), allocatable :: r_desorp(:)
  real(r8), allocatable :: r_occlude(:)
  real(r8), allocatable :: k_s1_biochem(:)
  real(r8), allocatable :: k_s2_biochem(:)
  real(r8), allocatable :: k_s3_biochem(:)
  real(r8), allocatable :: k_s4_biochem(:)



  ! !PUBLIC MEMBER FUNCTIONS:
  public :: soilorder_conrd ! Read and initialize soil order dependent  constants
  !
  ! !REVISION HISTORY:
  ! Created by X.YANG 01/12/2015
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine soilorder_conrd
    !
    ! !DESCRIPTION:
    ! Read and initialize soil order dependent constants
    !
    ! !USES:
    use fileutils ,  only : getfil
    use ncdio_pio ,  only : ncd_io, ncd_pio_closefile, ncd_pio_openfile,file_desc_t, &
                            ncd_inqdid, ncd_inqdlen
    use clm_varctl,  only : fsoilordercon
    use spmdMod   ,  only : masterproc

    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    ! Created by Gordon Bonan
    ! F. Li and S. Levis (11/06/12)
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn ! local file name
    integer :: i,n              ! loop indices
    integer :: ier              ! error code
    type(file_desc_t) :: ncid   ! pio netCDF file id
    integer :: dimid            ! netCDF dimension id
    integer :: nsoil             ! number of pfts on pft-physiology file
    logical :: readv            ! read variable in or not
    character(len=32) :: subname = 'soilorder_conrd'              ! subroutine name

    ! Expected soil order names: The names expected on the  file and the order
    ! they are expected to be in.
    !
    character(len=40), parameter :: expected_soilnames(0:nsoilorder) = (/ &
                 'Water                      '  &
               , 'Andisols                   '  &
               , 'Gelisols                   '  &
               , 'Histosols                  '  &
               , 'Entisols                   '  &
               , 'Inceptisols                '  &
               , 'Aridsols                   '  &
               , 'Vertisols                  '  &
               , 'Mollisols                  '  &
               , 'Alfisols                   '  &
               , 'Spodosols                  '  &
               , 'Ultisols                   '  &
               , 'Oxisols                    '  &
               , 'shiftingsand               '  &
               , 'rockland                   '  &
               , 'iceglacier                 '  &
    /)
!-----------------------------------------------------------------------

    allocate( smax        (0:nsoilorder) )
    allocate( ks_sorption (0:nsoilorder) )
    allocate( r_weather   (0:nsoilorder) )
    allocate( r_adsorp    (0:nsoilorder) )
    allocate( r_desorp    (0:nsoilorder) )
    allocate( r_occlude   (0:nsoilorder) )
    allocate(k_s1_biochem (0:nsoilorder) )
    allocate(k_s2_biochem (0:nsoilorder) )
    allocate(k_s3_biochem (0:nsoilorder) )
    allocate(k_s4_biochem (0:nsoilorder) )

   ! Set specific soil order values

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil order dependent parameters .....'
    end if
    call getfil (fsoilordercon, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'soilorder',dimid)
    call ncd_inqdlen(ncid,dimid,nsoil)

    call ncd_io('soilordername',soilordername, 'read', ncid, readvar=readv,posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('r_weather',r_weather, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('r_adsorp',r_adsorp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('r_desorp',r_desorp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('r_occlude',r_occlude, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('k_s1_biochem',k_s1_biochem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('k_s2_biochem',k_s2_biochem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('k_s3_biochem',k_s3_biochem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('k_s4_biochem',k_s4_biochem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('smax',smax, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))
    call ncd_io('ks_sorption',ks_sorption, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in soil order parameter'//errMsg(__FILE__, __LINE__))

    call ncd_pio_closefile(ncid)



    do i = 1,nsoilorder

       if ( trim(adjustl(soilordername(i))) /= trim(expected_soilnames(i)) )then
          write(iulog,*)'soilorder_conrd: soil order name is NOT what is expected, name = ', &
                        trim(soilordername(i)), ', expected name = ',trim(expected_soilnames(i))
          call endrun( 'soilorder_conrd: bad name for soil order onfsoilordercon dataset' )
       end if

    enddo

    if (masterproc) then
       write(iulog,*) 'Successfully read PFT physiological data'
       write(iulog,*)
    end if

   end subroutine soilorder_conrd


end module soilorder_varcon

