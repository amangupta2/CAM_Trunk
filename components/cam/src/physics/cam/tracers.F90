!======================================================================
! Module implements passive test tracers.
!
! Two options:
!
! 1) Specify only the number of tracers desired by setting the -nadv_tt option to configure.
!    This results in setting up the desired number of tracers making use of the tracers_suite
!    module to generate tracer names and initialize mixing ratios.
!
! 2) Specify both the number of tracers using configure's -nadv_tt option, and specify
!    the same number of tracer names using the test_tracer_names namelist variable.  This
!    specifies a set of passive tracers that are either initialized from analytic expressions
!    or by reading values from the IC file.  The tracers for which analytic expressions are
!    available are the following:
!
!    test_tracer_names     description
!    -----------------     -----------
!    TT_SLOT               Non-smooth scalar field (slotted cylinder)
!    TT_GBALL              Smooth Gaussian "ball"
!    TT_TANH               Zonally constant, tanh function of latitude.
!    TT_EM8                Constant field of size 1.e-8.
!    TT_Y2_2               Approximately Y^2_2 spherical harmonic
!    TT_Y32_16             Approximately Y^32_16 spherical harmonic
!    TT_LATP2              Zonally constant, latitude + 2.
!    TT_LONP2              Meridionally constant, cos(longitude) + 2.
!
!======================================================================

module tracers

use shr_kind_mod,    only: r8 => shr_kind_r8
use shr_sys_mod,     only: shr_sys_flush
use spmd_utils,      only: masterproc
use ppgrid,          only: pver
use physconst,       only: mwdry, cpair
use constituents,    only: cnst_add, cnst_name, cnst_longname
use tracers_suite,   only: get_tracer_name, init_cnst_tr
use cam_history,     only: addfld, add_default
use cam_logfile,     only: iulog
use cam_abortutils,  only: endrun

implicit none
private
save

public :: &
   tracers_readnl,            &! read namelist
   tracers_register,          &! register constituent
   tracers_implements_cnst,   &! true if named constituent is implemented by this package
   tracers_init_cnst,         &! initialize constituent field
   tracers_init,              &! initialize history fields, datasets
   tracers_timestep_forcing   ! ag4680@nyu.edu : new function created to force tracer value daily
                              ! Sets time dependent initial tracer values near
                              ! the
                              ! boundary every timestep


integer, parameter :: num_names_max = 30
integer, parameter :: num_analytic  = 8

! Data from namelist variables
integer           :: test_tracer_num = 0
character(len=16) :: test_tracer_names(num_names_max)

logical :: tracers_flag       = .false.  ! true => turn on test tracer code
logical :: tracers_suite_flag = .false.  ! true => test tracers provided by tracers_suite module

integer :: ixtrct=-999                   ! index of 1st constituent

character(len=16), parameter :: analytic_names(num_analytic) = &
   (/'TT_SLOT         ', 'TT_GBALL        ', 'TT_TANH         ', &
     'TT_EM8          ', 'TT_Y2_2         ', 'TT_Y32_16       ', &
     'TT_LATP2        ', 'TT_LONP2        ' /)

logical :: analytic_tracer(num_names_max)

!======================================================================
contains
!======================================================================

subroutine tracers_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_character

   ! args
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   integer :: i, j
   integer :: num_names
   character(len=*), parameter :: subname = 'tracers_readnl'

   namelist /test_tracers_nl/ test_tracer_num, test_tracer_names
   !-----------------------------------------------------------------------------

   test_tracer_names = (/ (' ', i=1,num_names_max) /)

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'test_tracers_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, test_tracers_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(test_tracer_names, len(test_tracer_names)*num_names_max, mpi_character, &
                  mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: test_tracer_names")

   call mpi_bcast(test_tracer_num,   1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: iradsw")

   ! If any tracers have been specified then turn on the tracers module
   if (test_tracer_num > 0) then
      tracers_flag = .true.
   else
      return
   end if

   ! Determine the number of tracer names supplied:
   num_names = 0
   analytic_tracer = .false.
   do i = 1, num_names_max
      if (len_trim(test_tracer_names(i)) > 0) then
         num_names = num_names + 1

         ! Does the tracer have an analytic IC?
         do j = 1, num_analytic
            if (trim(test_tracer_names(i)) == trim(analytic_names(j))) then
               analytic_tracer(i) = .true.
               exit
            end if
         end do
      else
         exit
      end if
   end do

   if (num_names > 0) then
      ! If test_tracer_names have been specified, the test_tracer_num should
      ! equal the number of names supplied.
      if (num_names /= test_tracer_num) then
         write(iulog, *) subname//' number of names, number of tracers: ', num_names, test_tracer_num
         call endrun(subname // ':: number of names does not match number of tracers')
      end if
   else
      ! If no names have been supplied then
      ! the tracers will be provided by the tracers_suite module.
      tracers_suite_flag = .true.
   end if

   ! Print summary to log file
   if (masterproc) then

      write(iulog, *) 'Test Tracers Module'
      write(iulog, *) '  Number of Test Tracers:', test_tracer_num
      if (tracers_suite_flag) then
         write(iulog, *) '  Tracers will be provided by tracers_suite module.'
      else
         do i = 1, num_names
            if (analytic_tracer(i)) then
               write(iulog, *) '  '//trim(test_tracer_names(i))//&
                  ' will be initialized from an analytic expression'
            else
               write(iulog, *) '  '//trim(test_tracer_names(i))//&
                  ' will be initialized from the IC file'
            end if
         end do
      end if
   end if

end subroutine  tracers_readnl

!======================================================================

subroutine tracers_register()

   ! Register advected tracers.

   ! Local variables
   integer  :: m, mm
   logical  :: read_from_file
   real(r8) :: minc
   character(len=16) :: name
   !-----------------------------------------------------------------------

   if (.not. tracers_flag) return

   minc = -1.e36_r8   ! min mixing ratio (disable qneg3)

   do m = 1, test_tracer_num

      read_from_file = .true.
      if (tracers_suite_flag) then
         name = get_tracer_name(m)  ! get name from suite file
         read_from_file = .false.
      else
         name = test_tracer_names(m)
         if (analytic_tracer(m)) read_from_file = .false.
      end if

      ! add constituent name to list of advected, save index number ixtrct
      call cnst_add(name, mwdry, cpair, minc, mm, &
                    readiv=read_from_file, mixtype='dry')
      if (m == 1) ixtrct = mm  ! save index number of first tracer
   end do

end subroutine tracers_register

!======================================================================

function tracers_implements_cnst(name)

   ! return true if specified constituent is implemented by this package

   ! Arguments
   character(len=*), intent(in) :: name   ! constituent name
   logical :: tracers_implements_cnst        ! return value

   ! Local variables
   integer :: m
   character(len=16) :: trc_name
   !-----------------------------------------------------------------------

   tracers_implements_cnst = .false.
   if (.not. tracers_flag) return

   do m = 1, test_tracer_num

      if (tracers_suite_flag) then
         trc_name = get_tracer_name(m)
      else
         trc_name = test_tracer_names(m)
      end if

      if (name == trc_name) then
         tracers_implements_cnst = .true.
         return
      end if
   end do

end function tracers_implements_cnst

!===============================================================================

subroutine tracers_init_cnst(name, latvals, lonvals, mask, q)

   ! Initialize test tracer mixing ratio

   character(len=*), intent(in)  :: name
   real(r8),         intent(in)  :: latvals(:) ! lat in radians (ncol)
   real(r8),         intent(in)  :: lonvals(:) ! lon in radians (ncol)
   logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev

   ! Local
   integer :: m
   logical :: found
   character(len=*), parameter :: subname = 'tracers_init_cnst'
   !----------------------------------------------------------------------------

   if (.not. tracers_flag) return

   found = .false.
   if (tracers_suite_flag) then

      do m = 1, test_tracer_num
         if (name ==  get_tracer_name(m))  then
            call init_cnst_tr(m, latvals, lonvals, mask, q)
            found = .true.
            exit
         endif
      end do

   else

      do m = 1, test_tracer_num
         if (name == test_tracer_names(m)) then
            if (analytic_tracer(m)) then
               call test_func_set(name, latvals, lonvals, mask, q)
               found = .true.
               exit
            else
               ! The initial values were supposed to be read from the IC file.  This
               ! call should not have been made in that case, so it appears that a requested
               ! tracer is not on the IC file.
               write(iulog, *) subname//': ERROR: tracer ', trim(name), ' should be on IC file'
               call endrun(subname//': ERROR: tracer missing from IC file')
            end if
         end if
      end do

   end if

   if (.not. found) then
      ! unrecognized tracer name
      write(iulog, *) subname//': ERROR: ', trim(name), ' not recognized'
      call endrun(subname//': ERROR: tracer name not recognized')
   end if

end subroutine tracers_init_cnst

!===============================================================================

subroutine tracers_init()

   ! Add tracers to history output

   ! Local
   integer m, mm
   character(len=16) :: name   ! constituent name

   if (.not. tracers_flag ) return

   do m = 1,test_tracer_num
      mm = ixtrct + m - 1
      ! ag4680@nyu.edu : Instantaneous write, not average
      call addfld(cnst_name(mm), (/ 'lev' /), 'I', 'kg/kg', cnst_longname(mm))
      !call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))
      call add_default(cnst_name(mm), 1, ' ')
   end do

end subroutine tracers_init

!=========================================================================================



!======================================================================
! Aman Gupta ag4680@nyu.edu : Begin : This function is new. Created to force the
! boundary value for
! tracer at every time step
! Function designed to compute only the clock tracer, delta tracer and pulse
! tracer. Hence num_tracer = 3
! In the main for-loop, thus, boundary values are mentioned for the 2 tracers
subroutine tracers_timestep_forcing(state, ptend, dt)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use cam_history,   only: outfld
    use time_manager,  only: get_curr_time
    use constituents,  only: pcnst
    use ppgrid,        only: pcols

    ! Arguments
    type(physics_state), intent(inout) :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(in)    :: dt                 ! timestep

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk                          ! chunk identifier
    integer :: ncol                           ! no. of column in chunk
    logical  :: lq(pcnst)
    character(len=16)  :: name                ! get tracer name
    integer :: num_tracer                     ! get number of passive tracers
    integer :: day,sec                        ! date variables
    real    :: t                              ! tracer boundary condition(time) 
    real    :: p_source                       ! tracer boundary layer
    real(r8), allocatable :: pmid(:,:)        ! mid-point pressure

    !------------------------------------------------------------------

    ! Get or manually set tracer number
    num_tracer=3          ! Set to 3 to include all tracers. 
                          ! Ideally this should be same as nadv_tt in config_files/definition.xml
    p_source = 70000.0_r8 ! p > p_source is the tracer source region

    if (.not. tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    call get_curr_time(day,sec)
    !print *,"get curr time -",day,' ',sec
    t=day + sec/86400.0 ! t is in days

    ! Do this for tracer indices
    lq(:)      = .FALSE.
    do i=1,num_tracer
        lq(ixtrct+i-1) = .TRUE.
    end do

    call physics_ptend_init(ptend,state%psetcols, 'tracers', lq=lq)

    lchnk = state%lchnk
    ncol  = state%ncol
    pmid  = state%pmid
        ! Old definition only for clock tracer and num_tracer=1 
        if (num_tracer .eq. 1) then
        do k = 1, pver
                ! For clock tracer
                if (k >= pver-2 .and. k <= pver) then
                        do i = 1, ncol
                                    ! passive tracer count begins at ixtrct to
                                    ! ixtrct + num_tracer -1
                                    ! A 300 day spinup is chosen.
                                    if (t > 300.0) then
                                            state%q(i,k,ixtrct) = t
                                    else
                                            state%q(i,k,ixtrct) = 0.
                                    end if
                        end do
                end if
        end do
        elseif (num_tracer .eq. 3) then
        ! New definition for 3 tracers - clock, delta and pulse
        do k = 1, pver
                do i = 1, ncol
                                ! passive tracer count begins in the register
                                ! from 
                                ! (ixtrct) to (ixtrct + num_tracer -1)
                                if (pmid(i,k) > 70000.) then

                                ! Clock Tracer
                                ! A 300 days spinup is chosen.
                                        if (t > 300.) then
                                            state%q( i, k, ixtrct) = t
                                        else
                                            state%q( i, k, ixtrct) = 0.
                                        end if

                                ! Delta Tracer : Force the boundary layer to
                                ! erode all at once
                                        if (t > 300. .and. t <= 330.) then
                                            state%q( i, k, ixtrct+1) = 1.
                                        else
                                            state%q( i, k, ixtrct+1) = 0.
                                        end if

                                ! Pulse Tracer : Let the bndary layer erode
                                ! naturally
                                        if (t <= 300.) then
                                            state%q( i, k, ixtrct+2) = 0.
                                        else if (t > 300. .and. t < 330.) then
                                            state%q( i, k, ixtrct+2) = 1.
                                        else if (t > 330.) then
                                            ! do nothing
                                        end if
                               end if
                end do
        end do
        end if

        ! Outputting the values to the history files
        do i=1,num_tracer
            name =  get_tracer_name(i)
            ! for simple physics the tracer order in the register is : (a) sphum(1) + (b)
            ! tracers(num_tracer) + (c) aoa_tracers(4)
            call outfld (name, state%q(:,:,i+ixtrct-1), pcols, lchnk)
        end do

end subroutine tracers_timestep_forcing

! ag4680@nyu.edu : END











subroutine test_func_set(name, latvals, lonvals, mask, q)

   ! use test_func code to set array q

   character(len=*), intent(in)  :: name
   real(r8),         intent(in)  :: latvals(:) ! lat in radians (ncol)
   real(r8),         intent(in)  :: lonvals(:) ! lon in radians (ncol)
   logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev

   ! local variables
   integer :: i, k
   !----------------------------------------------------------------------------

   do i = 1, size(mask)
      if (mask(i)) then
         do k = 1, size(q, 2)
            q(i,k) = test_func(name, latvals(i), lonvals(i), k)
         end do
      end if
   end do

end subroutine test_func_set

!=========================================================================================

function test_func(name, lat, lon, k) result(fout)

   ! Analytic test functions.

   use physconst, only: pi
   use hycoef,    only: hyam, hybm, ps0

   character(len=*), intent(in) :: name
   real(r8),         intent(in) :: lon     ! radians
   real(r8),         intent(in) :: lat     ! radians
   integer,          intent(in) :: k       ! vertical index for layer mid-points

   real(r8) :: fout

   real(r8), parameter :: psurf_moist = 100000.0_r8 ! moist surface pressure
   real(r8), parameter :: deg2rad = pi/180._r8

   real(r8) :: lon1, lat1, R0, Rg1
   real(r8) :: eta, eta_c
   real(r8) :: cos_tmp, sin_tmp
   !----------------------------------------------------------------------------

   select case(name)
   case('TT_SLOT')
      !
      !   Non-smooth scalar field (slotted cylinder)
      !
      R0 = 0.25_r8
      lon1  = 20.0_r8*deg2rad   ! longitudinal position, 20E
      lat1  = 40.0_r8*deg2rad  ! latitudinal position, 40N
      Rg1 = acos(sin(lat1)*sin(lat)+cos(lat1)*cos(lat)*cos(lon-lon1))

      if ((Rg1 <= R0) .AND. (abs(lon-lon1) >= R0/6)) then
         fout = 2.0_r8
      elseif ((Rg1 <= R0) .AND. (abs(lon-lon1) < R0/6) &
         .AND. (lat-lat1 < -5.0_r8*R0/12.0_r8)) then
         fout = 2.0_r8
      else
         fout = 1.0_r8
      endif

   case('TT_GBALL')
      !
      ! Smooth Gaussian "ball"
      !
      R0    = 10.0_r8           ! radius of the perturbation
      lon1  = 20.0_r8*deg2rad   ! longitudinal position, 20E
      lat1  = 40.0_r8*deg2rad  ! latitudinal position, 40N
      eta_c = 0.6_r8
      sin_tmp = SIN(lat1)*SIN(lat)
      cos_tmp = COS(lat1)*COS(lat)
      Rg1 = ACOS( sin_tmp + cos_tmp*COS(lon-lon1) )    ! great circle distance
      eta =  (hyam(k)*ps0 + hybm(k)*psurf_moist)/psurf_moist
      fout = EXP(- ((Rg1*R0)**2 + ((eta-eta_c)/0.1_r8)**2))
      IF (ABS(fout) < 1.0E-8_r8) fout = 0.0_r8

   case('TT_TANH')
      !
      !
      !
      fout = 0.5_r8 * ( tanh( 3.0_r8*abs(lat)-pi ) + 1.0_r8)

   case('TT_EM8')
      fout = 1.0e-8_r8

   case('TT_Y2_2')
      !
      ! approximately Y^2_2 spherical harmonic
      !
      fout = 0.5_r8 + 0.5_r8*(cos(lat)*cos(lat)*cos(2.0_r8*lon))

   case('TT_Y32_16')
      !
      ! approximately Y32_16 spherical harmonic
      !
      fout = 0.5_r8 + 0.5_r8*(cos(16*lon)*(sin(2_r8*lat)**16))

   case('TT_LATP2')
      fout = 2.0_r8 + lat

   case('TT_LONP2')
      fout = 2.0_r8 + cos(lon)

   case default
      call endrun("test_func: ERROR: name not recognized")
   end select

end function test_func

!=========================================================================================

end module tracers
