program xfem_tetra_test
  use omp_lib, only: omp_get_num_procs, omp_set_num_threads
  use types, only: ik, lk, cl, stdout, stderr, debug
  use general_routines, only: print_error, to_lower, str2i
  use read_input, only: read_data, input_file_name, inp_file
  use mesh_data, only: mesh_data_statistics, mesh_data_finish
  use volume_integral, only: test_functions
  implicit none
!*****************************************************************************80
  integer(ik) :: esta = 0
  character(len=cl) :: emsg = ''
!*****************************************************************************80
  logical(lk) :: input_file_given = .false.
!*****************************************************************************80
  integer(ik) :: cpus = 1
  logical(lk) :: cpus_given = .false.
!*****************************************************************************80
! Command line arguments
!*****************************************************************************80
  call get_command_line_args(esta,emsg)
  if ( esta /= 0 ) call print_error(esta,emsg)
  ! Command line checks
  call command_line_checks(esta,emsg)
  if ( esta /= 0 ) call print_error(esta,emsg)
  ! Open input file
  inquire(file=input_file_name,exist=input_file_given)
  if ( .not.input_file_given ) then
    esta = -1
    emsg = 'Input file does not exist.'
  end if
  if ( esta /= 0 ) call print_error(esta,emsg)
  open(newunit=inp_file,file=input_file_name,iostat=esta,iomsg=emsg,&
  &action='read',status='old')
  if ( esta /= 0 ) call print_error(esta,emsg)
!*****************************************************************************80
! Read data
!*****************************************************************************80
  call read_data(esta,emsg)
  if ( esta /= 0 ) call print_error(esta,emsg)
  write(stdout,'(a)') 'Read data complete.'
  if ( debug ) call mesh_data_statistics()
!*****************************************************************************80
! Calculate volume
!*****************************************************************************80
  call test_functions(esta,emsg)
  if ( esta /= 0 ) call print_error(esta,emsg)
  write(stdout,'(a)') 'Tests complete.'
!*****************************************************************************80
! Finish
!*****************************************************************************80
  call mesh_data_finish(esta,emsg)
  if ( esta /= 0 ) call print_error(esta,emsg)
!*****************************************************************************80
  stop 0
!*****************************************************************************80
contains
!*****************************************************************************80
! Check command line arguments
!*****************************************************************************80
  subroutine get_command_line_args(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i, n
    integer(ik) :: narg
    character(len=cl) :: argval
    character(len=cl), allocatable :: arg_list(:)
    ! Check the number of arguments
    narg = command_argument_count()
    if( narg == 0 ) then
      write(stderr,'(a)') 'Get command line arguments:&
      & no valid arguments.'
    end if
    allocate(arg_list(narg),stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    ! Set args to list and lower
    do i = 1, narg
      call get_command_argument(i,arg_list(i))
    end do
    ! Set arguments
    i = 1
    do while( i <= narg )
      call to_lower(arg_list(i))
      select case( arg_list(i) )
      ! Vhodna datoteka
      case( '-i', '--input' )
        i = i + 1
        argval = arg_list(i)
        n = 0; n = len_trim(argval)
        if( trim(argval(n-3:)) /= '.inp' ) then
          input_file_name = trim(argval) // '.inp'
        else
          input_file_name = argval
        end if
        input_file_given = .true.
      ! Set number of cpus
      case( '-c', '--cpus')
        i = i + 1
        argval = arg_list(i)
        call str2i(argval,cpus,esta,emsg)
        if  ( esta /= 0 ) return
        cpus_given = .true.
      ! Help
      case( '-h', '--help' )
        write(stdout,'(a)') 'No help.'
      case default
        write(stderr,'(a)') 'Get command line arguments:&
        & argument not supported.'
      end select
      i = i + 1
    end do
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine
!*****************************************************************************80
  subroutine command_line_checks(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: available_cpus
    ! Check input
    if ( .not. input_file_given ) then
      esta = -1
      emsg = 'Get command line arguments:&
      & no input file given'
      return
    end if
    ! Check cpus
    if ( cpus_given ) then
      available_cpus = omp_get_num_procs() ! Get num of available cpu
      if ( cpus > available_cpus ) then
        write(stdout,'(a,i0,a,i0,a)') 'The number of cpus requested (',cpus,&
        &') exceeds the number of available cpus (',available_cpus,').'
        cpus = available_cpus
      end if
    end if
    write(stdout,'(a,i0,a)') 'Using ', cpus, ' processors for solution.'
    call omp_set_num_threads(cpus)
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine command_line_checks
!*****************************************************************************80
end program

