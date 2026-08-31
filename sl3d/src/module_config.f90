module module_config

  !=========================================================================!
  !  The parameters of a run, and where they come from.                     !
  !                                                                         !
  !  A program sets the parameters it uses to its own defaults and then     !
  !  calls read_config, which takes the command line and, if one is named   !
  !  there, a namelist file.  Anything the file does not mention keeps the  !
  !  program's default, so a program run without a config file behaves      !
  !  exactly as it did before there were any:                               !
  !                                                                         !
  !     sea_level_3d prem.simp.vm1              defaults                    !
  !     sea_level_3d prem.simp.vm1 ice6g        and another ice model       !
  !     sea_level_3d run.nml                    from a config file          !
  !                                                                         !
  !  A config file looks like                                               !
  !                                                                         !
  !     &run                                                                !
  !       earth_model  = 'prem.simp.vm1'                                    !
  !       ice_model_name = 'ice6g'                                          !
  !       lmax_run     = 64                                                 !
  !       time_step    = 50.0      ! years                                  !
  !       output_every = 500.0     ! years                                  !
  !       run_length   = -1.0      ! years; negative is the whole history   !
  !     /                                                                   !
  !                                                                         !
  !  and the kernel programs read a second group, &kernel.  Both groups are !
  !  optional, and so is every entry in them.                               !
  !                                                                         !
  !  write_config puts the parameters a run actually used beside its        !
  !  results, so that the results say what produced them.                   !
  !=========================================================================!

  use nrtype
  use module_util
  use module_ice

  implicit none

  ! the run
  character(len=256) :: earth_model    = ''        ! the earth model file
  character(len=32)  :: ice_model_name = ''        ! ice5g, ice6g or ice7g
  integer(i4b)       :: lmax_run       = 64        ! truncation degree
  real(dp)           :: time_step      = 50.0_dp   ! years
  real(dp)           :: output_every   = 500.0_dp  ! years
  real(dp)           :: run_length     = -1.0_dp   ! years; < 0 is the whole
                                                   ! of the ice history

  ! the measurement a kernel or inversion is for
  real(dp)     :: site_lat    = 0.0_dp
  real(dp)     :: site_long   = 0.0_dp
  integer(i4b) :: load_degree = 2
  integer(i4b) :: load_order  = 2

  namelist /run/ earth_model,ice_model_name,lmax_run,time_step, &
                 output_every,run_length
  namelist /kernel/ site_lat,site_long,load_degree,load_order

contains


  subroutine read_config
    !  Work out the parameters of this run from the command line and, if one
    !  is named there, a namelist file.  Anything given on the command line
    !  wins over the file.
    implicit none

    integer(i4b) :: narg,iarg,io,ios,length,status
    character(len=256) :: arg,config_file,cli_model,cli_ice
    logical(lgt) :: ltmp

    narg = command_argument_count()
    if(narg == 0) then
       print *, 'inputs: [earth model] {ice model: ice5g (default), ice6g, ice7g}'
       print *, '    or: [config file, named *.nml]'
       stop
    end if

    config_file = ''
    cli_model   = ''
    cli_ice     = ''
    do iarg = 1,narg
       call get_command_argument(iarg,arg)
       if(len_trim(arg) > 4 .and. &
            arg(max(len_trim(arg)-3,1):len_trim(arg)) == '.nml') then
          config_file = arg
       else if(is_ice_model(arg)) then
          cli_ice = arg
       else
          cli_model = arg
       end if
    end do

    ! the file first, so that the command line overrides it
    if(len_trim(config_file) > 0) then
       inquire(file=trim(config_file),exist=ltmp)
       if(.not.ltmp) then
          print *, 'cannot find config file ',trim(config_file)
          stop
       end if
       open(newunit=io,file=trim(config_file),action='read')
       ! a negative iostat is simply a group the file does not have
       read(io,nml=run,iostat=ios)
       if(ios > 0) call config_error(config_file,'run')
       rewind(io)
       read(io,nml=kernel,iostat=ios)
       if(ios > 0) call config_error(config_file,'kernel')
       close(io)
    end if

    if(len_trim(cli_model) > 0) earth_model    = cli_model
    if(len_trim(cli_ice)   > 0) ice_model_name = cli_ice

    if(len_trim(earth_model) == 0) stop 'no earth model given'
    inquire(file=trim(earth_model),exist=ltmp)
    if(.not.ltmp) then
       print *, 'cannot find earth model ',trim(earth_model)
       stop
    end if

    ! the ice model: from the command line, the file, or the environment
    if(len_trim(ice_model_name) == 0) then
       call get_environment_variable('ICE_MODEL',arg,length,status)
       if(status == 0 .and. length > 0) ice_model_name = arg
    end if
    if(len_trim(ice_model_name) > 0) call set_ice_model_type(ice_model_name)
    ice_model_name = ice_model            ! record what was actually used

  end subroutine read_config


  subroutine config_error(file,group)
    implicit none
    character(len=*), intent(in) :: file,group
    print *, 'could not read the &',trim(group),' group of ',trim(file)
    print *, 'check the parameter names against those listed in the README'
    stop
  end subroutine config_error


  subroutine write_config(file)
    !  Record the parameters this run used, beside its results.  The file is
    !  itself a valid config file, so a run can be repeated by handing it
    !  back to the program.  It is written out by hand rather than with
    !  write(io,nml=run) only because that pads every string to its full
    !  declared length, which makes the record hard to read.
    implicit none
    character(len=*), intent(in) :: file
    integer(i4b) :: io
    open(newunit=io,file=file)
    write(io,'(a)')    '! written by the run itself; usable as a config file'
    write(io,'(a)')    '&run'
    write(io,'(3a)')   "  earth_model    = '",trim(earth_model),"'"
    write(io,'(3a)')   "  ice_model_name = '",trim(ice_model_name),"'"
    write(io,'(a,i0)') '  lmax_run       = ',lmax_run
    write(io,'(a,g0)') '  time_step      = ',time_step
    write(io,'(a,g0)') '  output_every   = ',output_every
    write(io,'(a,g0)') '  run_length     = ',run_length
    write(io,'(a)')    '/'
    write(io,'(a)')    '&kernel'
    write(io,'(a,g0)') '  site_lat       = ',site_lat
    write(io,'(a,g0)') '  site_long      = ',site_long
    write(io,'(a,i0)') '  load_degree    = ',load_degree
    write(io,'(a,i0)') '  load_order     = ',load_order
    write(io,'(a)')    '/'
    close(io)
  end subroutine write_config


  function run_time() result(t)
    !  How long the calculation runs for, in years: the whole of the ice
    !  history unless a run length was asked for.
    implicit none
    real(dp) :: t
    if(run_length > 0.0_dp) then
       t = run_length
    else
       t = ice_total_time()
    end if
  end function run_time


end module module_config
