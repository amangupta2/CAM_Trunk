module ic_held_suarez

  !-----------------------------------------------------------------------
  !
  ! Purpose: Set Held-Suarez initial conditions based on input coordinates
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  implicit none
  private

  ! Public interface
  public :: hs94_set_ic

!==============================================================================
CONTAINS
!==============================================================================

  subroutine hs94_set_ic(latvals, lonvals, U, V, T, PS, PHIS,           &
       Q, m_cnst, mask, verbose)
    use const_init,    only: cnst_init_default
    use constituents,  only: cnst_name

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set Held-Suarez initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(inout) :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical,  optional, intent(in)    :: verbose    ! For internal use

    ! Local variables
    logical, allocatable              :: mask_use(:)
    logical                           :: verbose_use
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'HS94_SET_IC'

    ! ag4680@nyu.edu : Topography related parameters
    ! phi0/1: Bottom and Top latitudinal extent
    ! H is the zsurf in meters, knum = topography wavenumber
    ! === begin ===
    real(r8) :: H=3000., knum=2, phi0, lamb, phi1, phi, scale_H, zsurf
    real(r8),parameter :: pi = 3.141592, grav = 9.81, R = 287.0 ! can also use inbuilt constants
    integer :: NX,NY
    ! === end ===



    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun('cnst_init_default: input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    if (present(verbose)) then
      verbose_use = verbose
    else
      verbose_use = .true.
    end if

    ncol = size(latvals, 1)
    nlev = -1
    if (present(U)) then
      nlev = size(U, 2)
      do k = 1, nlev
        where(mask_use)
          U(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          U initialized by "',subname,'"'
      end if
    end if

    if (present(V)) then
      nlev = size(V, 2)
      do k = 1, nlev
        where(mask_use)
          V(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          V initialized by "',subname,'"'
      end if
    end if

    if (present(T)) then
      nlev = size(T, 2)
      do k = 1, nlev
        where(mask_use)
          T(:,k) = 250.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          T initialized by "',subname,'"'
      end if
    end if

    if (present(PS)) then
      !where(mask_use)
      !  PS = 100000.0_r8
      !end where

      ! ag4680@nyu.edu : defining surface pressure for the topography
       ! This pressure is in hydrostatic balance with the geopotential
       phi0 = 25.*pi/180. ! 25
       phi1 = 65.*pi/180. ! 65
       NX = size(lonvals)
       NY = size(latvals)
       scale_H = R*250.0_r8/grav ! Change for general runs

                do m=1,NX
                           phi = latvals(m)
                           lamb = lonvals(m)
                           if(mask_use(m)) then
                                if ( (phi > phi0) .and. phi < phi1 ) then

                                             ! Height of the Gerber-Polvani
                                             ! topography
                                             zsurf = H*cos(knum*lamb)*(sin(pi*((phi - phi0)/(phi1 - phi0) ) )**2 )

                                            ! Hydrostatically balanced surface pressure (T_init = 250.0 for HS94)
                                             PS(m) = 100000.0_r8*(exp(-zsurf/scale_H))

                                   else
                                             PS(m) = 100000.0_r8
                                   end if
                           end if
                 end do

      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
    end if

    if (present(PHIS)) then
      !where(mask_use)
      !  PHIS = 0.0_r8
      !end where

      ! ag4680@nyu.edu : Geopotential for the Gerber-Polvani topography
       phi0 = 25.*pi/180. ! 25
       phi1 = 65.*pi/180. ! 65
       NX = size(lonvals)
       NY = size(latvals)

                do m=1,NX
                           phi = latvals(m)
                           lamb = lonvals(m)
                           if(mask_use(m)) then
                               if ( (phi > phi0) .and. phi < phi1 ) then
                                    PHIS(m) = grav*H*cos(knum*lamb)*( sin(pi*((phi - phi0)/(phi1 - phi0) ) )**2 )
                               else
                                    PHIS(m) = 0.0_r8
                               end if
                           end if

                 end do


      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PHIS initialized by "',subname,'"'
      end if
    end if

    if (present(Q)) then
      nlev = size(Q, 2)
      ncnst = size(m_cnst, 1)
      do m = 1, ncnst
        if (m_cnst(m) == 1) then
          ! No water vapor in Held-Suarez
          do k = 1, nlev
            where(mask_use)
              Q(:,k,m_cnst(m)) = 0.0_r8
            end where
          end do
          if(masterproc .and. verbose_use) then
            write(iulog,*) '          ', trim(cnst_name(m_cnst(m))), ' initialized by "',subname,'"'
          end if
        else
          call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
               mask=mask_use, verbose=verbose_use, notfound=.false.)
        end if
      end do
    end if

    deallocate(mask_use)

  end subroutine hs94_set_ic

end module ic_held_suarez
