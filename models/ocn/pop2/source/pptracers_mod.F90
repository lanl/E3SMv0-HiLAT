!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module pptracers_mod

!BOP
! !MODULE: pptracers_mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: pptracers_mod.F90 26603 2011-01-28 23:09:02Z njn01 $

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
   use passive_tracer_tools, only: ind_name_pair, tracer_read, &
       rest_read_tracer_block, file_read_tracer_block
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: pptracers_tracer_cnt,        &
             pptracers_init,              &
             pptracers_set_sflux

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by pseudo-pseudo-tracers
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      pptracers_tracer_cnt     = 2    ! All pseudo-pseudo-tracers

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      ppTEMP_ind = 1,   &! pseudo-pseudo-temperature
      ppSALT_ind = 2     ! pseudo-pseudo-salinity

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(pptracers_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(ppTEMP_ind, 'ppTEMP'), &
      ind_name_pair(ppSALT_ind, 'ppSALT') /)

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: pptracers_init
! !INTERFACE:

 subroutine pptracers_init(init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, TEMP, SALT, tadvect_ctype, &
                      errorCode)

! !DESCRIPTION:
!  Initialize pptracers tracer module. This involves setting metadata, reading
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

   real(r8), dimension(nx_block,ny_block,km,3,max_blocks_clinic), &
      intent(in) :: TEMP, SALT

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(pptracers_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,pptracers_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   character (char_len), dimension(:), intent(out) :: &
      tadvect_ctype     ! advection method for pseudo-pseudo-tracers tracers

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'pptracers_mod:pptracers_init'

   character(char_len) :: &
      init_pptracers_option,           & ! option for initialization of pseudo-pseudo-tracers
      init_pptracers_init_file,        & ! filename for option 'file'
      init_pptracers_init_file_fmt,    & ! file format for option 'file'
      pptracers_tadvect_ctype            ! advection method for pseudo-pseudo-tracers

   logical(log_kind) :: &
      lnml_found             ! Was pptracers_nml found ?

   integer(int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over depth levels
      iblock,              & ! index for looping over blocks
      nml_error              ! namelist i/o error flag

!     l,                   & ! index for looping over time levels

   type(tracer_read), dimension(pptracers_tracer_cnt) :: &
      tracer_init_ext        ! namelist variable for initializing tracers

   namelist /pptracers_nml/ &
      init_pptracers_option, init_pptracers_init_file, tracer_init_ext, &
      init_pptracers_init_file_fmt, pptracers_tadvect_ctype

   character (char_len) ::  &
      pptracers_restart_filename  ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   tracer_d_module(ppTEMP_ind)%short_name = 'ppTEMP'
   tracer_d_module(ppTEMP_ind)%long_name  = 'Pseudo-pseudo-temperature'
   tracer_d_module(ppTEMP_ind)%units      = 'degC'
   tracer_d_module(ppTEMP_ind)%tend_units = 'degC/s'
   tracer_d_module(ppTEMP_ind)%flux_units = 'degC cm/s'
   tracer_d_module(ppTEMP_ind)%scale_factor = 1.0_rtavg

   tracer_d_module(ppSALT_ind)%short_name = 'ppSALT'
   tracer_d_module(ppSALT_ind)%long_name  = 'Pseudo-pseudo-salinity'
   tracer_d_module(ppSALT_ind)%units      = 'g/kg'
   tracer_d_module(ppSALT_ind)%tend_units = 'g/kg/s'
   tracer_d_module(ppSALT_ind)%flux_units = 'g/kg cm/s'
   tracer_d_module(ppSALT_ind)%scale_factor = 1000.0_rtavg

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_pptracers_option = 'unknown'
   init_pptracers_init_file = 'unknown'
   init_pptracers_init_file_fmt = 'bin'

   pptracers_tadvect_ctype = 'base_model'

   do n = 1,pptracers_tracer_cnt
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
         read(nml_in, nml=pptracers_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'pptracers_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_pptracers_option , master_task)
   call broadcast_scalar(init_pptracers_init_file, master_task)
   call broadcast_scalar(init_pptracers_init_file_fmt, master_task)

   call broadcast_scalar(pptracers_tadvect_ctype, master_task)
   tadvect_ctype = pptracers_tadvect_ctype

   do n = 1,pptracers_tracer_cnt
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

   select case (init_pptracers_option)

   case ('ccsm_startup', 'ccsm_startup_spunup')
      TRACER_MODULE(:,:,:,1,:,:) = TEMP
      TRACER_MODULE(:,:,:,2,:,:) = SALT

      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d pseudo-pseudo-tracers same as T and S' 
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif

   case ('zero')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d pseudo-pseudo-tracers set to zero'
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
       
   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      pptracers_restart_filename = char_blank

      if (init_pptracers_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read pseudo-pseudo-tracers from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         pptracers_restart_filename = read_restart_filename
         init_pptracers_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         pptracers_restart_filename = trim(init_pptracers_init_file)

      endif

      call rest_read_tracer_block(init_pptracers_init_file_fmt, &
                                  pptracers_restart_filename,   &
                                  tracer_d_module,         &
                                  TRACER_MODULE)

! Temporary hack to test step-by-step changes

!        TRACER_MODULE(:,:,:,1,:,:) = TEMP
!        TRACER_MODULE(:,:,:,2,:,:) = SALT

!        if (my_task == master_task) then
!           write(stdout,delim_fmt)
!           write(stdout,*) 'Init Max of T and pT ', &
!                           maxval(abs(TEMP)), &
!                           maxval(abs(TRACER_MODULE(:,:,:,1,:,:)))
!           write(stdout,*) 'Init Max of S and pS ', &
!                           maxval(abs(SALT)), &
!                           maxval(abs(TRACER_MODULE(:,:,:,2,:,:)))
!        endif

   case ('file')
      call document(subname, 'pseudo-pseudo-tracers being read from separate file')

      call file_read_tracer_block(init_pptracers_init_file_fmt, &
                                  init_pptracers_init_file,     &
                                  tracer_d_module,         &
                                  ind_name_table,          &
                                  tracer_init_ext,         &
                                  TRACER_MODULE)
 
      if (n_topo_smooth > 0) then
         do n = 1, pptracers_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'pptracers_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

   case default
      call document(subname, 'unknown init_pptracers_option = ', init_pptracers_option)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,pptracers_tracer_cnt
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

 end subroutine pptracers_init

!BOP
! !IROUTINE: pptracers_set_sflux
! !INTERFACE:

 subroutine pptracers_set_sflux(STF_MODULE,SHF,SFWF)

! !DESCRIPTION:
!  Compute surface fluxes for pptracers tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      SHF,   & ! Surface heat flux
      SFWF     ! Surface freshwater flux

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,pptracers_tracer_cnt,max_blocks_clinic), &
      intent(inout) :: STF_MODULE    ! Surface fluxes

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

      STF_MODULE(:,:,1,:) = SHF
      STF_MODULE(:,:,2,:) = SFWF

!-----------------------------------------------------------------------
!EOC

 end subroutine pptracers_set_sflux

!*****************************************************************************

end module pptracers_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
