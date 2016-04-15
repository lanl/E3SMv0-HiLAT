!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module surfdyes_mod

!BOP
! !MODULE: surfdyes_mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: surfdyes_mod.F90 26603 2011-01-28 23:09:02Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use blocks
   use domain_size, only: max_blocks_clinic, km, nx_global
   use domain, only: nblocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use prognostic, only: tracer_field
   use kinds_mod
   use constants, only: c0, c1, char_blank, &
       delim_fmt,ocn_ref_salinity,ppt_to_salt
   use io, only: data_set
   use io_types, only: stdout, nml_in, nml_filename
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   use passive_tracer_tools, only: ind_name_pair, tracer_read, &
       rest_read_tracer_block, file_read_tracer_block
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: surfdyes_tracer_cnt,        &
             surfdyes_init,              &
             surfdyes_set_sflux

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      surfdyes_tracer_cnt     = 4    ! All dye tracers released at the surface

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      dyePREC_ind = 1,   &! Precipitation
      dyeMELT_ind = 2,   &! Sea ice melt
      dyeROFF_ind = 3,   &! River run-off
      dyeIOFF_ind = 4     ! Ice run-off

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(surfdyes_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(dyePREC_ind, 'DYE_PREC'), &
      ind_name_pair(dyeMELT_ind, 'DYE_MELT'), &
      ind_name_pair(dyeROFF_ind, 'DYE_ROFF'), &
      ind_name_pair(dyeIOFF_ind, 'DYE_IOFF') /)

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: surfdyes_init
! !INTERFACE:

 subroutine surfdyes_init(init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, tadvect_ctype, errorCode)

! !DESCRIPTION:
!  Initialize surfdyes tracer module. This involves setting metadata, reading
!  the module namelist and setting initial conditions.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use broadcast, only: broadcast_scalar
   use prognostic, only: curtime, oldtime
   use grid, only: KMT, n_topo_smooth, fill_points

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(surfdyes_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,surfdyes_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   character (char_len), dimension(:), intent(out) :: &
      tadvect_ctype     ! advection method for surfdyes tracers

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'surfdyes_mod:surfdyes_init'

   character(char_len) :: &
      init_surfdyes_option,           & ! option for initialization of surfdyes
      init_surfdyes_init_file,        & ! filename for option 'file'
      init_surfdyes_init_file_fmt,    & ! file format for option 'file'
      surfdyes_tadvect_ctype            ! advection method for surfdyes tracers

   logical(log_kind) :: &
      lnml_found             ! Was surfdyes_nml found ?

   integer(int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over depth levels
      iblock,              & ! index for looping over blocks
      nml_error              ! namelist i/o error flag

!     l,                   & ! index for looping over time levels

   type(tracer_read), dimension(surfdyes_tracer_cnt) :: &
      tracer_init_ext        ! namelist variable for initializing tracers

   namelist /surfdyes_nml/ &
      init_surfdyes_option, init_surfdyes_init_file, tracer_init_ext, &
      init_surfdyes_init_file_fmt, surfdyes_tadvect_ctype

   character (char_len) ::  &
      surfdyes_restart_filename  ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   tracer_d_module(dyePREC_ind)%short_name = 'DYE_PREC'
   tracer_d_module(dyePREC_ind)%long_name  = 'Dye Tracer Precipitation'
   tracer_d_module(dyePREC_ind)%units      = 'g/kg'
   tracer_d_module(dyePREC_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dyePREC_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dyeMELT_ind)%short_name = 'DYE_MELT'
   tracer_d_module(dyeMELT_ind)%long_name  = 'Dye Tracer Sea Ice Melt'
   tracer_d_module(dyeMELT_ind)%units      = 'g/kg'
   tracer_d_module(dyeMELT_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dyeMELT_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dyeROFF_ind)%short_name = 'DYE_ROFF'
   tracer_d_module(dyeROFF_ind)%long_name  = 'Dye Tracer River Run-Off'
   tracer_d_module(dyeROFF_ind)%units      = 'g/kg'
   tracer_d_module(dyeROFF_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dyeROFF_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dyeIOFF_ind)%short_name = 'DYE_IOFF'
   tracer_d_module(dyeIOFF_ind)%long_name  = 'Dye Tracer Ice Run-Off'
   tracer_d_module(dyeIOFF_ind)%units      = 'g/kg'
   tracer_d_module(dyeIOFF_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dyeIOFF_ind)%flux_units = 'g/kg cm/s'

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_surfdyes_option = 'unknown'
   init_surfdyes_init_file = 'unknown'
   init_surfdyes_init_file_fmt = 'bin'

   surfdyes_tadvect_ctype = 'base_model'

   do n = 1,surfdyes_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then  
         nml_error = -1
      else
         nml_error =  1      
      endif
      do while (nml_error > 0)
         read(nml_in, nml=surfdyes_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'surfdyes_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_surfdyes_option , master_task)
   call broadcast_scalar(init_surfdyes_init_file, master_task)
   call broadcast_scalar(init_surfdyes_init_file_fmt, master_task)

   call broadcast_scalar(surfdyes_tadvect_ctype, master_task)
   tadvect_ctype = surfdyes_tadvect_ctype

   do n = 1,surfdyes_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (init_surfdyes_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d Dye Tracers set to all zeros' 
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
       
   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      surfdyes_restart_filename = char_blank

      if (init_surfdyes_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read surfdyes from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         surfdyes_restart_filename = read_restart_filename
         init_surfdyes_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         surfdyes_restart_filename = trim(init_surfdyes_init_file)

      endif

      call rest_read_tracer_block(init_surfdyes_init_file_fmt, &
                                  surfdyes_restart_filename,   &
                                  tracer_d_module,         &
                                  TRACER_MODULE)

   case ('file')
      call document(subname, 'surfdyes being read from separate file')

      call file_read_tracer_block(init_surfdyes_init_file_fmt, &
                                  init_surfdyes_init_file,     &
                                  tracer_d_module,         &
                                  ind_name_table,          &
                                  tracer_init_ext,         &
                                  TRACER_MODULE)
 
      if (n_topo_smooth > 0) then
         do n = 1, surfdyes_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'surfdyes_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

   case default
      call document(subname, 'unknown init_surfdyes_option = ', init_surfdyes_option)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,surfdyes_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine surfdyes_init

!***********************************************************************
!BOP
! !IROUTINE: surfdyes_set_sflux
! !INTERFACE:

 subroutine surfdyes_set_sflux(STF_MODULE,PREC_F,MELT_F,ROFF_F,IOFF_F)

! !DESCRIPTION:
!  Compute surface fluxes for surfdyes tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      PREC_F,   & ! precipitation flux kg/m^2/s  fw
      MELT_F,   & ! snow&ice melt flux kg/m^2/s  fw
      ROFF_F,   & ! river runoff  flux kg/m^2/s  fw
      IOFF_F      ! ice   runoff  flux kg/m^2/s  fw

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,surfdyes_tracer_cnt,max_blocks_clinic), &
      intent(inout) :: STF_MODULE    ! Virtual surface fluxes

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n, iblock     ! loop indices

   real (r8), dimension(nx_block,ny_block) :: &
      FPER

   do iblock = 1,nblocks_clinic

      STF_MODULE(:,:,:,iblock) = c0

      do n=1,surfdyes_tracer_cnt

! Source distributions
! Assuming unit of dye concentration is g/kg
! Flux wanted by POP is g/kg cm/s
! FW [cm/s] = FW [kg/m^2/s] / rho [kg/m^3] * 100 cm/m

         FPER = c0
         if (n == dyePREC_ind) then
            FPER = c1 * PREC_F(:,:,iblock) * 1e-1_r8
         elseif (n == dyeMELT_ind) then
            where (MELT_F(:,:,iblock) > 0)
               FPER = c1 * MELT_F(:,:,iblock) * 1e-1_r8
            end where
         elseif (n == dyeROFF_ind) then
            FPER = c1 * ROFF_F(:,:,iblock) * 1e-1_r8
         elseif (n == dyeIOFF_ind) then
            FPER = c1 * IOFF_F(:,:,iblock) * 1e-1_r8
         endif
         STF_MODULE(:,:,n,iblock) = STF_MODULE(:,:,n,iblock) + FPER

      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine surfdyes_set_sflux

!*****************************************************************************

end module surfdyes_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
