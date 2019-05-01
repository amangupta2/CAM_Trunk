#define MODHS 1
#undef MODHS
module held_suarez_cam
! ----------------------------------------------------------------------------------
! Modified to inlude an option for the Polvani and Kushner (2002) , 
! GRL, 29, 7, 10.1029/2001GL014284 (PK02) equilibrium temperature profile.
! 
! If pkstrat = .True. then Polvani and Kushner relaxation temperature profile
! is used.  
!
! Namelist parameter: vgamma, sets the vortex strength (gamma parameter in PK02)
!
! Modifications are denoted by 
! 
! !PKSTRAT
! blah blah blah
! !END-PKSTRAT
!
! Isla Simpson 8th June 2017
! ----------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8

! PKSTAT
  use cam_abortutils  ,only:endrun 
  use spmd_utils  ,only:masterproc
  use spmd_utils, only: mpicom,mstrid=>masterprocid, mpi_integer, mpi_real8, &
                        mpi_logical, mpi_character
! END-PKSTRAT



  use ppgrid,       only: pcols, pver

  implicit none
  private
  save

  public :: held_suarez_init, held_suarez_tend

  !PKSTRAT
  public :: pkstrat_readnl
  !END-PKSTRAT

  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

  !PKSTRAT
  logical :: pkstrat !.True. to use the PK02 TEQ
  logical :: relax_tropical_stratosphere ! ag4680@nyu.edu : specified wind by default
  real(r8) :: vgamma ! vortex gamma parameter (controlling vortex strength)
  integer :: lat0   ! Polar Vortex boundary
  integer :: eps    ! Hemispheric Assymetry. Seasons
  !END PKSTRAT



!======================================================================= 
contains
!======================================================================= 

  subroutine held_suarez_init(pbuf2d)
    use physics_buffer,     only: physics_buffer_desc
    use cam_history,        only: addfld, add_default
    use physconst,          only: cappa, cpair
    use ref_pres,           only: pref_mid_norm, psurf_ref
    use held_suarez,        only: held_suarez_1994_init

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! Set model constant values
    call held_suarez_1994_init(cappa, cpair, psurf_ref, pref_mid_norm)

    ! This field is added by radiation when full physics is used
    call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
         'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')

    call add_default('QRS', 1, ' ')

 end subroutine held_suarez_init

  subroutine held_suarez_tend(state, ptend, ztodt)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  algorithm 1: Held/Suarez IDEALIZED physics
    !  algorithm 2: Held/Suarez IDEALIZED physics (Williamson modified stratosphere
    !  algorithm 3: Held/Suarez IDEALIZED physics (Lin/Williamson modified strato/meso-sphere
    !
    ! Author: J. Olson
    ! 
    !-----------------------------------------------------------------------
    use physconst,          only: cpairv
    use phys_grid,          only: get_rlat_all_p
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use cam_abortutils,     only: endrun
    use cam_history,        only: outfld
    use time_manager,       only: get_nstep, get_curr_date, get_start_date, timemgr_datediff 
    use held_suarez,        only: held_suarez_1994

    !
    ! Input arguments
    !
    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: ztodt            ! Two times model timestep (2 delta-t)
                                                           !
                                                           ! Output argument
                                                           !
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
                                                           !
    !---------------------------Local workspace-----------------------------

    integer                            :: lchnk            ! chunk identifier
    integer                            :: ncol             ! number of atmospheric columns
    integer                            :: nstep            ! timestep number

    real(r8)                           :: clat(pcols)      ! latitudes(radians) for columns
    real(r8)                           :: pmid(pcols,pver) ! mid-point pressure
    integer                            :: i, k             ! Longitude, level indices
    !
    !-----------------------------------------------------------------------
    !

    lchnk = state%lchnk
    ncol  = state%ncol

    nstep   = get_nstep()

    call get_rlat_all_p(lchnk, ncol, clat)
    do k = 1, pver
      do i = 1, ncol
        pmid(i,k) = state%pmid(i,k)
      end do
    end do

    ! initialize individual parameterization tendencies
    call physics_ptend_init(ptend, state%psetcols, 'held_suarez', ls=.true., lu=.true., lv=.true.)

!    call held_suarez_1994(pcols, ncol, clat, pmid, &
!         state%u, state%v, state%t, ptend%u, ptend%v, ptend%s)
    !PKSTRAT
    !Modified call to Held_suarez_1994 to pass pkstrat and vgamma
    !call held_suarez_1994(pcols, ncol, clat, pmid, &
    !   state%u, state%v, state%t, ptend%u, ptend%v, ptend%s, pkstrat, vgamma)
    !END-PKSTRAT

    ! MODIFIED PKSTRAT
    ! ag4680@nyu.edu : added 2 arguments passed through namelists : eps and lat0
       call held_suarez_1994(pcols, ncol, clat, pmid, &
       state%u, state%v, state%t, ptend%u, ptend%v, ptend%s, pkstrat, relax_tropical_stratosphere, vgamma, eps, lat0)
    ! END MODIFIED PKSTRAT

    ! Note, we assume that there are no subcolumns in simple physics
    pmid(:ncol,:) = ptend%s(:ncol, :)/cpairv(:ncol,:,lchnk)
    if (pcols > ncol) then
      pmid(ncol+1:,:) = 0.0_r8
    end if
    call outfld('QRS', pmid, pcols, lchnk)
    
  end subroutine held_suarez_tend

!PKSTRAT
  subroutine pkstrat_readnl(nlfile)
! ---------------------------------------------
! Read in namelist parameters to control the PK02 stratosphere
! logical: pkstrat.  If (pkstrat=.True.) then use the PK02
! stratosphere
!
! Isla Simpson (6th March 2017)
!ENDPKSTRAT

  use shr_kind_mod,only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use namelist_utils  ,only:find_group_name
  use units           ,only:getunit,freeunit
  character(len=*),intent(in)::nlfile
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'pkstrat_readnl'

  ! ag4680@nyu.edu : reading eps, lat0 and relax_tropical_stratosphere via namelist too!
  namelist /pkstrat_nl/ pkstrat, relax_tropical_stratosphere, vgamma, eps, lat0

  !Set default namelist parameters
  pkstrat=.False.
  relax_tropical_stratosphere = .true.
  vgamma=2._r8
  eps=-10.      ! default eps value
  lat0=-50.     ! default lat0 value

  !Read in namelist parameters
  if (masterproc) then
   unitn = getunit()
   open(unitn,file=trim(nlfile),status='old')
   call find_group_name(unitn,'pkstrat_nl',status=ierr)
   if (ierr.eq.0) then
     read(unitn,pkstrat_nl,iostat=ierr)
     if (ierr.ne.0) then
       call endrun('pkstrat_readnl:: ERROR reading namelist')
     endif
   endif
   close(unitn)
   call freeunit(unitn)
  endif

  call mpi_bcast(pkstrat   , 1, mpi_logical, mstrid, mpicom, ierr)
  call mpi_bcast(vgamma    , 1, mpi_real8, mstrid, mpicom, ierr)
  call mpi_bcast(relax_tropical_stratosphere   , 1, mpi_logical, mstrid, mpicom, ierr)
  ! Aman Gupta : ag4680@nyu.edu : added for 2 new namelist parameters 
  call mpi_bcast(eps       , 1, mpi_integer, mstrid, mpicom, ierr)
  call mpi_bcast(lat0      , 1, mpi_integer, mstrid, mpicom, ierr)


  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: held_suarez_1994")

  return
  end subroutine

end module held_suarez_cam
