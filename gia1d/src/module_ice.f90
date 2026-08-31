module module_ice

  !=========================================================================!
  !  Reads a deglaciation history from the netCDF files as they are         !
  !  distributed, and holds the properties of the ice model in use.         !
  !                                                                         !
  !  Three models are supported:                                            !
  !                                                                         !
  !     ice5g   ICE-5G v1.2,   21 ka to present, 39 times                   !
  !     ice6g   ICE-6G_C VM5a, 26 ka to present, 48 times                   !
  !     ice7g   ICE-7G_NA VM7, 26 ka to present, 48 times                   !
  !                                                                         !
  !  All three are on a one degree grid and are tabulated every 1000 years  !
  !  at the oldest times and every 500 years thereafter, but they differ in !
  !  their file naming, in the names of their variables, and in whether     !
  !  their longitudes sit on cell edges or cell centres.  The differences   !
  !  are confined to the table in set_ice_model_type below.                 !
  !                                                                         !
  !  The files are looked for in a directory named after the model, so      !
  !  ./ice5g, ./ice6g or ./ice7g; setting ICE_DIR moves that search to      !
  !  $ICE_DIR/ice5g and so on.                                             !
  !  The build links data/ into build/work, so programs run from there find !
  !  them with no further setup.  The model itself is chosen with an        !
  !  optional second command line argument, or with the ICE_MODEL           !
  !  environment variable, and defaults to ice5g.                           !
  !=========================================================================!

  use nrtype
  use module_util
  use netcdf

  implicit none

  ! the ice grid, common to all three models
  integer(i4b), protected :: nlat_ice  = 180
  integer(i4b), protected :: nlong_ice = 360

  ! properties of the model in use, set by set_ice_model_type
  character(len=8), protected  :: ice_model = 'ice5g'
  integer(i4b), protected :: nt_ice = 39     ! number of times
  integer(i4b), protected :: ice_ib = 5      ! index at which 1000 year
                                             ! spacing gives way to 500

  ! how this model's files are named and what is in them
  character(len=64)  :: ice_prefix = 'ice5g_v1.2_'
  character(len=64)  :: ice_suffix = 'k_1deg.nc'
  logical(lgt)       :: ice_pad_age = .true.  ! ages as 21.0, 09.5, 00.0
  character(len=16)  :: ice_var_thickness = 'sftgit'
  character(len=16)  :: ice_var_topography = 'orog'
  logical(lgt)       :: ice_lon_centred = .false. ! longitudes at 0.5, 1.5, ...

  private :: ice_prefix,ice_suffix,ice_pad_age,ice_var_thickness, &
             ice_var_topography,ice_lon_centred,check

contains


  subroutine init_ice_model
    !  Chooses the ice model from the second command line argument if there is
    !  one, otherwise from ICE_MODEL, otherwise leaves the default in place.
    implicit none
    character(len=256) :: name
    integer(i4b) :: length,status

    if(command_argument_count() >= 2) then
       call get_command_argument(2,name)
       call set_ice_model_type(name)
       return
    end if

    call get_environment_variable('ICE_MODEL',name,length,status)
    if(status == 0 .and. length > 0) call set_ice_model_type(name)

  end subroutine init_ice_model


  function is_ice_model(name) result(known)
    !  Is this the name of one of the ice models?
    implicit none
    character(len=*), intent(in) :: name
    logical(lgt) :: known
    select case(trim(adjustl(name)))
    case('ice5g','ice6g','ice7g')
       known = .true.
    case default
       known = .false.
    end select
  end function is_ice_model


  subroutine set_ice_model_type(name)
    !  Everything that differs between the three models lives here.
    implicit none
    character(len=*), intent(in) :: name

    select case(trim(adjustl(name)))

    case('ice5g')
       ice_model = 'ice5g'
       nt_ice = 39                     ! 21, 20, 19, 18, 17, 16.5, ... 0 ka
       ice_ib = 5
       ice_prefix = 'ice5g_v1.2_'
       ice_suffix = 'k_1deg.nc'
       ice_pad_age = .true.
       ice_var_thickness = 'sftgit'
       ice_var_topography = 'orog'
       ice_lon_centred = .false.

    case('ice6g')
       ice_model = 'ice6g'
       nt_ice = 48                     ! 26, 25, 24, 23, 22, 21, 20.5, ... 0 ka
       ice_ib = 6
       ice_prefix = 'I6_C.VM5a_1deg.'
       ice_suffix = '.nc'
       ice_pad_age = .false.
       ice_var_thickness = 'stgit'
       ! ICE6G also carries an orog, but it is clipped at zero and so has no
       ! bathymetry; Topo is the field matching ICE5G's orog.
       ice_var_topography = 'Topo'
       ice_lon_centred = .true.

    case('ice7g')
       ice_model = 'ice7g'
       nt_ice = 48                     ! 26, 25, 24, 23, 22, 21, 20.5, ... 0 ka
       ice_ib = 6
       ice_prefix = 'I7G_NA.VM7_1deg.'
       ice_suffix = '.nc'
       ice_pad_age = .false.
       ice_var_thickness = 'stgit'
       ice_var_topography = 'Topo'
       ice_lon_centred = .true.

    case default
       print *, 'unknown ice model: ',trim(adjustl(name))
       print *, 'known models are ice5g, ice6g and ice7g'
       stop

    end select

  end subroutine set_ice_model_type


  function ice_coarse_time() result(t)
    !  Length of the leading, 1000 year spaced, part of the history in years.
    implicit none
    real(dp) :: t
    t = (ice_ib-1)*1000.0_dp
  end function ice_coarse_time


  function ice_total_time() result(t)
    !  Length of the whole history in years: 21000 for ICE5G, 26000 otherwise.
    implicit none
    real(dp) :: t
    t = ice_coarse_time() + (nt_ice-ice_ib)*500.0_dp
  end function ice_total_time


  function ice_time(it) result(t)
    !  Time in years of index it, measured forward from the start of the run.
    implicit none
    integer(i4b), intent(in) :: it
    real(dp) :: t
    if(it <= ice_ib) then
       t = (it-1)*1000.0_dp
    else
       t = ice_coarse_time() + (it-ice_ib)*500.0_dp
    end if
  end function ice_time


  function ice_file(it) result(fname)
    !  The file holding time index it, where it = 1 is the oldest time.
    implicit none
    integer(i4b), intent(in) :: it
    character(len=256) :: fname

    integer(i4b) :: age_half,whole
    logical(lgt) :: half
    character(len=256) :: stem,dir

    ! age of this index, in half thousands of years
    age_half = nint(2.0_dp*(ice_total_time() - ice_time(it))/1000.0_dp)
    whole = age_half/2
    half = (mod(age_half,2) /= 0)

    if(ice_pad_age) then
       ! ICE5G writes every age with two digits and one decimal: 21.0, 09.5
       if(whole < 10) then
          call string_cat_int(trim(ice_prefix) // '0',whole,stem)
       else
          call string_cat_int(trim(ice_prefix),whole,stem)
       end if
       if(half) then
          stem = trim(stem) // '.5'
       else
          stem = trim(stem) // '.0'
       end if
    else
       ! ICE6G and ICE7G write 26, 21, 20.5, 0
       call string_cat_int(trim(ice_prefix),whole,stem)
       if(half) stem = trim(stem) // '.5'
    end if

    fname = trim(stem) // trim(ice_suffix)

    ! ICE_DIR names the directory that holds ice5g, ice6g and ice7g
    call get_environment_variable('ICE_DIR',dir)
    if(len_trim(dir) == 0) then
       dir = ice_model
    else
       dir = trim(dir) // '/' // trim(ice_model)
    end if

    fname = trim(dir) // '/' // trim(fname)

  end function ice_file


  subroutine read_ice_grid(it,field,grid)
    !  Reads one field for time index it: 'thickness' for the ice thickness,
    !  or 'topography' for the surface altitude, both returned in metres.
    !
    !  The files run from latitude -89.5 north, which is the codes' own
    !  convention, and from longitude 0 (ICE5G) or 0.5 (ICE6G, ICE7G) east.
    !  The grid returned is rotated by half a turn to put longitude -180 at
    !  index 1, as get_lat_long expects.  Where the file's longitudes sit at
    !  cell centres they are also averaged onto the whole degrees the codes
    !  interpolate on, which is a shift of half a cell.
    implicit none
    integer(i4b), intent(in) :: it
    character(len=*), intent(in) :: field
    real(dp), dimension(nlong_ice,nlat_ice), intent(out) :: grid

    character(len=256) :: fname
    character(len=16) :: varname
    integer(i4b) :: ncid,varid,ilong,ilonged,ilongm1
    real(dp), dimension(nlong_ice,nlat_ice) :: raw

    if(it < 1 .or. it > nt_ice) stop 'read_ice_grid: time index out of range'

    select case(trim(field))
    case('thickness')
       varname = ice_var_thickness
    case('topography')
       varname = ice_var_topography
    case default
       print *, 'read_ice_grid: unknown field ',trim(field)
       stop
    end select

    fname = ice_file(it)

    call check(nf90_open(trim(fname),nf90_nowrite,ncid),trim(fname))
    call check(nf90_inq_varid(ncid,trim(varname),varid),trim(fname))
    call check(nf90_get_var(ncid,varid,raw),trim(fname))
    call check(nf90_close(ncid),trim(fname))

    do ilong = 1,nlong_ice
       ilonged = ilong + (nlong_ice/2)
       if(ilonged > nlong_ice) ilonged = ilonged - nlong_ice
       if(ice_lon_centred) then
          ilongm1 = ilong-1
          if(ilongm1 < 1) ilongm1 = nlong_ice
          grid(ilonged,:) = 0.5_dp*(raw(ilongm1,:) + raw(ilong,:))
       else
          grid(ilonged,:) = raw(ilong,:)
       end if
    end do

  end subroutine read_ice_grid


  subroutine check(status,fname)
    !  Turns a netCDF return code into a readable error.
    implicit none
    integer(i4b), intent(in) :: status
    character(len=*), intent(in) :: fname
    if(status /= nf90_noerr) then
       print *, 'error reading ',trim(fname),': ',trim(nf90_strerror(status))
       print *, 'ice files are looked for in a directory named after the', &
                ' model; set ICE_DIR to the directory that holds them'
       stop
    end if
  end subroutine check


end module module_ice
