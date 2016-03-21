!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module sectdyes_mod

!BOP
! !MODULE: sectdyes_mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: sectdyes_mod.F90 26603 2011-01-28 23:09:02Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use blocks
   use domain_size, only: max_blocks_clinic, km, nx_global
   use domain, only: nblocks_clinic, blocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use prognostic, only: tracer_field
   use kinds_mod
   use constants, only: c0, c1, char_blank, &
       delim_fmt,ocn_ref_salinity,ppt_to_salt
   use io, only: data_set
   use io_types, only: stdout, nml_in, nml_filename, get_unit, release_unit
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   use passive_tracer_tools, only: ind_name_pair, tracer_read, &
       rest_read_tracer_block, file_read_tracer_block
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: sectdyes_tracer_cnt,        &
             sectdyes_init,              &
             sectdyes_reset

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      sectdyes_tracer_cnt  = 6    ! All dye tracers released on interior sections

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      dye01_ind =  1,  & 
      dye02_ind =  2,  &
      dye03_ind =  3,  &
      dye04_ind =  4,  &
      dye05_ind =  5,  &
      dye06_ind =  6

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(sectdyes_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(dye01_ind, 'DYE01'), &
      ind_name_pair(dye02_ind, 'DYE02'), &
      ind_name_pair(dye03_ind, 'DYE03'), &
      ind_name_pair(dye04_ind, 'DYE04'), &
      ind_name_pair(dye05_ind, 'DYE05'), &
      ind_name_pair(dye06_ind, 'DYE06') /)

!-----------------------------------------------------------------------
!  Section definitions
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      sectdyes_section_cnt          ! number of source sections read in from file

   type :: section
      character(char_len) :: name   ! name for section
      integer (int_kind)  :: dyeID  ! tracer number this section will be assigned to
      integer (int_kind), dimension(6,max_blocks_clinic) :: &
                             add    ! address range for computing transp
   end type

   type (section), dimension(:), allocatable :: &
      sections                      ! section info for all potential sections

   integer (r8), dimension(:,:), allocatable :: &
                     dyeVal         ! Dye values assigned to all dyes on each section


!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: sectdyes_init
! !INTERFACE:

 subroutine sectdyes_init(init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, tadvect_ctype, errorCode)

! !DESCRIPTION:
!  Initialize sectdyes tracer module. This involves setting metadata, reading
!  the module namelist and setting initial conditions.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use broadcast, only: broadcast_scalar, broadcast_array
   use prognostic, only: curtime, oldtime
   use grid, only: KMT, n_topo_smooth, fill_points

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(sectdyes_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,sectdyes_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   character (char_len), dimension(:), intent(out) :: &
      tadvect_ctype     ! advection method for sectdyes tracers

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'sectdyes_mod:sectdyes_init'

   character(char_len) :: &
      init_sectdyes_option,           & ! option for initialization of sectdyes
      init_sectdyes_init_file,        & ! filename for option 'file'
      init_sectdyes_init_file_fmt,    & ! file format for option 'file'
      sectdyes_tadvect_ctype,         & ! advection method for sectdyes tracers
      sectdyes_sections_file            ! filename for section definitions

   logical(log_kind) :: &
      lnml_found             ! Was sectdyes_nml found ?

   integer(int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over depth levels
      iblock,              & ! index for looping over blocks
      nml_error,           & ! namelist i/o error flag
      ib,ie,jb,je,         & ! beg,end indices for block physical domain
      nu

   type(tracer_read), dimension(sectdyes_tracer_cnt) :: &
      tracer_init_ext        ! namelist variable for initializing tracers

   namelist /sectdyes_nml/ &
      init_sectdyes_option, init_sectdyes_init_file, tracer_init_ext, &
      init_sectdyes_init_file_fmt, sectdyes_tadvect_ctype, &
      sectdyes_sections_file

   character (char_len) ::  &
      sectdyes_restart_filename  ! modified file name for restart file

   integer (int_kind), dimension(6) :: &
      tmp_add             ! section global addresses
                          ! (ibeg,iend,jbeg,jend,kbeg,kend)

   type (block) ::         &
      this_block           ! block information for current block

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   tracer_d_module(dye01_ind)%short_name = 'DYE01'
   tracer_d_module(dye01_ind)%long_name  = 'Dye Tracer 1'
   tracer_d_module(dye01_ind)%units      = 'g/kg'
   tracer_d_module(dye01_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dye01_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dye02_ind)%short_name = 'DYE02'
   tracer_d_module(dye02_ind)%long_name  = 'Dye Tracer 2'
   tracer_d_module(dye02_ind)%units      = 'g/kg'
   tracer_d_module(dye02_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dye02_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dye03_ind)%short_name = 'DYE03'
   tracer_d_module(dye03_ind)%long_name  = 'Dye Tracer 3'
   tracer_d_module(dye03_ind)%units      = 'g/kg'
   tracer_d_module(dye03_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dye03_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dye04_ind)%short_name = 'DYE04'
   tracer_d_module(dye04_ind)%long_name  = 'Dye Tracer 4'
   tracer_d_module(dye04_ind)%units      = 'g/kg'
   tracer_d_module(dye04_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dye04_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dye05_ind)%short_name = 'DYE05'
   tracer_d_module(dye05_ind)%long_name  = 'Dye Tracer 5'
   tracer_d_module(dye05_ind)%units      = 'g/kg'
   tracer_d_module(dye05_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dye05_ind)%flux_units = 'g/kg cm/s'

   tracer_d_module(dye06_ind)%short_name = 'DYE06'
   tracer_d_module(dye06_ind)%long_name  = 'Dye Tracer 6'
   tracer_d_module(dye06_ind)%units      = 'g/kg'
   tracer_d_module(dye06_ind)%tend_units = 'g/kg/s'
   tracer_d_module(dye06_ind)%flux_units = 'g/kg cm/s'

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_sectdyes_option = 'unknown'
   init_sectdyes_init_file = 'unknown'
   init_sectdyes_init_file_fmt = 'bin'

   sectdyes_tadvect_ctype = 'base_model'

   sectdyes_sections_file = 'unknown'

   do n = 1,sectdyes_tracer_cnt
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
         read(nml_in, nml=sectdyes_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'sectdyes_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_sectdyes_option , master_task)
   call broadcast_scalar(init_sectdyes_init_file, master_task)
   call broadcast_scalar(init_sectdyes_init_file_fmt, master_task)

   call broadcast_scalar(sectdyes_tadvect_ctype, master_task)
   tadvect_ctype = sectdyes_tadvect_ctype

   call broadcast_scalar(sectdyes_sections_file, master_task)

   do n = 1,sectdyes_tracer_cnt
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

   select case (init_sectdyes_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d Dye Tracers set to all zeros' 
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
       
   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      sectdyes_restart_filename = char_blank

      if (init_sectdyes_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read sectdyes from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         sectdyes_restart_filename = read_restart_filename
         init_sectdyes_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         sectdyes_restart_filename = trim(init_sectdyes_init_file)

      endif

      call rest_read_tracer_block(init_sectdyes_init_file_fmt, &
                                  sectdyes_restart_filename,   &
                                  tracer_d_module,         &
                                  TRACER_MODULE)

   case ('file')
      call document(subname, 'sectdyes being read from separate file')

      call file_read_tracer_block(init_sectdyes_init_file_fmt, &
                                  init_sectdyes_init_file,     &
                                  tracer_d_module,         &
                                  ind_name_table,          &
                                  tracer_init_ext,         &
                                  TRACER_MODULE)
 
      if (n_topo_smooth > 0) then
         do n = 1, sectdyes_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'sectdyes_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

   case default
      call document(subname, 'unknown init_sectdyes_option = ', init_sectdyes_option)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,sectdyes_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do
   enddo

!-----------------------------------------------------------------------
!  Read in definition of dye release sections
!-----------------------------------------------------------------------

   call get_unit(nu)
   if (my_task == master_task) then

      open(nu, file=sectdyes_sections_file, status='old')
      read(nu,*) sectdyes_section_cnt

      write(stdout,'(a14,i4,1x,a27)') 'The following ',sectdyes_section_cnt, &
                                    'sections will be read in, and assigned  &
                                     to dye tracer number'
   endif

   call broadcast_scalar(sectdyes_section_cnt, master_task)

   allocate( sections(sectdyes_section_cnt) )
   allocate( dyeVal(sectdyes_section_cnt,sectdyes_tracer_cnt) )

   dyeVal = c0

   do n=1,sectdyes_section_cnt
      if (my_task == master_task) then

         read(nu,'(6(i4,1x),i2,1x,a)') tmp_add, &
                     sections(n)%dyeID, sections(n)%name
         write(stdout,'(a2,a,x,i2)') '  ',trim(sections(n)%name),sections(n)%dyeID

         if (sections(n)%dyeID > sectdyes_tracer_cnt) then
            write(stdout,'(3(x,i2))') n,sections(n)%dyeID,sectdyes_tracer_cnt
            call exit_POP(sigAbort,'ERROR: Requested dye exceeds available tracers')
         endif

      endif

      call broadcast_array(tmp_add,master_task)
      call broadcast_scalar(sections(n)%dyeID,master_task)

      dyeVal(n,sections(n)%dyeID) = c1

!-----------------------------------------------------------------------
!  Translate global section specs to local addresses
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,ib,ie,jb,je)

      do iblock=1,nblocks_clinic

         this_block = get_block(blocks_clinic(iblock),iblock)

         ib = this_block%ib
         ie = this_block%ie
         jb = this_block%jb
         je = this_block%je

         sections(n)%add(1,iblock) = nx_block
         sections(n)%add(2,iblock) = 1
         sections(n)%add(3,iblock) = ny_block
         sections(n)%add(4,iblock) = 1
         sections(n)%add(5,iblock) = tmp_add(5)  ! beg k index
         sections(n)%add(6,iblock) = tmp_add(6)  ! end k index

         if (tmp_add(1) <= this_block%i_glob(ie) .and. &
             tmp_add(2) >= this_block%i_glob(ib) ) then

            if (tmp_add(1) >= this_block%i_glob(ib)) then
               sections(n)%add(1,iblock) = &
                           tmp_add(1) - this_block%i_glob(ib) + ib
            else
               sections(n)%add(1,iblock) = ib
            endif
            if (tmp_add(2) <= this_block%i_glob(ie)) then
               sections(n)%add(2,iblock) = &
                           tmp_add(2) - this_block%i_glob(ib) + ib
            else
               sections(n)%add(2,iblock) = ie
            endif
         endif

         if (tmp_add(3) <= this_block%j_glob(je) .and. &
             tmp_add(4) >= this_block%j_glob(jb) ) then
            if (tmp_add(3) >= this_block%j_glob(jb)) then
               sections(n)%add(3,iblock) = &
                           tmp_add(3) - this_block%j_glob(jb) + jb
            else
               sections(n)%add(3,iblock) = jb
            endif
           if (tmp_add(4) <= this_block%j_glob(je)) then
              sections(n)%add(4,iblock) = &
                           tmp_add(4) - this_block%j_glob(jb) + jb
           else
              sections(n)%add(4,iblock) = je
           endif
         endif

      enddo ! block loop

      !$OMP END PARALLEL DO

   enddo

   if (my_task == master_task) close(nu)
   call release_unit(nu)

!-----------------------------------------------------------------------
!EOC

 end subroutine sectdyes_init

!***********************************************************************
!BOP
! !IROUTINE: sectdyes_reset
! !INTERFACE:

 subroutine sectdyes_reset(TRACER_MODULE, this_block, iblock)

! !DESCRIPTION:
!  reset internal value for dye tracers at 'seeding' sections
!  Sections are hard coded, and need to be changed for different grids or 
!  release configurations
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,km,sectdyes_tracer_cnt), intent(inout) :: &
      TRACER_MODULE      ! dye tracers

   type (block), intent(in) :: &
      this_block             ! block info for the current block

   integer(int_kind), intent(in) :: iblock


!EOP
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, m, ib, ie, jb, je, i, j, k

!-----------------------------------------------------------------------
!  Reset dye tracer values at the specified sections.
!  The dye will be reset to 1 on its assigned sections, and to 0
!  on every other section. 
!-----------------------------------------------------------------------

      ib = this_block%ib
      ie = this_block%ie
      jb = this_block%jb
      je = this_block%je

      do n=1,sectdyes_section_cnt

         do k=sections(n)%add(5,iblock), &
              sections(n)%add(6,iblock)
         do j=sections(n)%add(3,iblock), &
              sections(n)%add(4,iblock)
         do i=sections(n)%add(1,iblock), &
              sections(n)%add(2,iblock)

            do m=1,sectdyes_tracer_cnt
               TRACER_MODULE(i,j,k,m) = dyeVal(n,m)
            enddo

         end do
         end do
         end do

      end do  ! section loop


!-----------------------------------------------------------------------
!EOC

 end subroutine sectdyes_reset

!*****************************************************************************

end module sectdyes_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
