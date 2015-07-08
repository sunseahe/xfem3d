module read_input
!*****************************************************************************80
  use types, only: ik, rk, lk, cl, stdout
  use general_routines, only: str2i, str2r, to_lower, strip_char
  use tokenize_string, only: tokenize, comma, equal
  use point, only: dom, zero_pnt, point_3d_t
  use fe_c3d10, only: nelnod, c3d10_t
  use mesh_data, only: fe_type, set_nodes, set_finite_elements
  use reinitalzation, only: set_reinitalization
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  integer(ik) :: inp_file = 0
  character(len=cl) :: input_file_name = ''
  !
  logical(lk) :: end_of_file = .false.
!*****************************************************************************80
  public :: inp_file, input_file_name, read_data
!*****************************************************************************80
  contains
!*****************************************************************************80
! Read data
!*****************************************************************************80
  subroutine read_data(esta,emsg)
    !
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    integer(ik) :: rstat
    character(len=cl) :: inp_str, keyword, sub_keyword
    character(len=cl), allocatable :: keyword_arg(:), sub_keyword_arg(:)
    logical(lk) :: fe_type_given
    !
    fe_type_given = .false.
    !
    ra: do
      if( end_of_file ) exit ra
      call read_input_string(inp_str,rstat)
      select case(rstat)
      case( 0 )
        exit ra
      case( 1 ); cycle ra
      case( 2 ); continue
      case( 3 )
        esta = -1
        emsg = 'Read data: read keyword error'
        return
      end select
      ! tokenize arguments
      call tokenize(inp_str,comma,keyword_arg,esta,emsg)
      if ( esta /= 0 ) return
      call strip_char(keyword_arg,'*')
      keyword = keyword_arg(1) ! Firstone
      select case( keyword )
!*****************************************************************************80
      case( 'node' )
        call read_nodes(esta,emsg); if ( esta /= 0 ) return
!*****************************************************************************80
! Move this to a subroutine
      case( 'element')
        do i = 2, size(keyword_arg)
          sub_keyword = keyword_arg(i)
          call tokenize(sub_keyword,equal,sub_keyword_arg,esta,emsg)
          if ( esta /= 0 ) return
          select case( sub_keyword_arg(1) )
          case( 'type' )
            if ( sub_keyword_arg(2) == fe_type ) then
              fe_type_given = .true.
              call read_connectivity(esta,emsg)
              if ( esta /= 0 ) return
            else
              esta = -1
              emsg = 'Read data: element type >' // trim(sub_keyword_arg(2)) &
              &// '< not supported'
              return
            end if
          case default
            esta = -1
            emsg = 'Read data: unknown sub keyword >' //  &
            &trim(sub_keyword_arg(1)) // '<'
            return
          end select
        end do
        if ( .not. fe_type_given ) then
          esta = -1
          emsg = 'Read data: finite element type not given'
          return
        end if
!*****************************************************************************80
      case ( 'reinitalization' )
        call read_reinitalization(esta,emsg); if ( esta /= 0 ) return
!*****************************************************************************80
      case default
        esta = -1
        emsg = 'Read data: unknown keyword >' // trim(keyword) // '<'
        return
      end select
    end do ra
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine read_data
!*****************************************************************************80
! Read nodes
!*****************************************************************************80
  subroutine read_nodes(esta,emsg)
    !
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    integer(ik) :: rstat
    real(rk) :: node_coo(dom)
    character(len=cl) :: inp_str
    character(len=cl), allocatable :: read_arg(:)
    type(point_3d_t), allocatable :: nodes(:)
    !
    write(stdout,'(a)') 'Reading nodes ...'
    allocate(nodes(0),stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    ra: do
      call read_input_string(inp_str,rstat)
      select case(rstat)
      case( 0 )
        exit ra
      case( 1 ); cycle ra
      case( 2 )
        backspace(inp_file)
        exit ra
      case( 3 ); continue
      end select
      ! tokenize arguments
      call tokenize(inp_str,comma,read_arg,esta,emsg)
      if ( esta /= 0 ) return
      if ( size(read_arg) /= dom + 1 ) then
        esta = -1
        emsg = 'Read nodes: domain not supported'
        return
      end if
      do i = 2, dom + 1
        call str2r(read_arg(i),node_coo(i-1),esta,emsg)
        if ( esta /= 0 ) return
      end do
      nodes = [ nodes, point_3d_t(x=node_coo) ]
    end do ra
    ! Add nodes
    call set_nodes(nodes,esta,emsg)
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine read_nodes
!*****************************************************************************80
! Read connectivity
!*****************************************************************************80
  subroutine read_connectivity(esta,emsg)
    !
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    integer(ik) :: rstat
    integer(ik) :: el_conn(nelnod)
    character(len=cl) :: inp_str
    character(len=cl), allocatable :: read_arg(:)
    type(c3d10_t), allocatable :: finite_elements(:)
    !
    write(stdout,'(a)') 'Reading connectivity ...'
    allocate(finite_elements(0),stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    ra: do
      call read_input_string(inp_str,rstat)
      select case(rstat)
      case( 0 )
        exit ra
      case( 1 ); cycle ra
      case( 2 )
        backspace(inp_file)
        exit ra
      case( 3 ); continue
      end select
      ! tokenize arguments
      call tokenize(inp_str,comma,read_arg,esta,emsg)
      if ( esta /= 0 ) return
      if ( size(read_arg) /= nelnod + 1 ) then
        esta = -1
        emsg = 'Read connectivity: Number of vertices is incorrect'
        return
      end if
      do i = 2, nelnod + 1
        call str2i(read_arg(i),el_conn(i-1),esta,emsg)
        if ( esta /= 0 ) return
      end do
      finite_elements = [ finite_elements, c3d10_t(nodes=zero_pnt,&
      &connectivity=el_conn) ]
      if ( esta /= 0 ) return
    end do ra
    ! Create element array
    call set_finite_elements(finite_elements,esta,emsg)
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine read_connectivity
!*****************************************************************************80
! Read reinitalization
!*****************************************************************************80
    subroutine read_reinitalization(esta,emsg)
      !
      integer(ik), intent(out) :: esta
      character(len=*), intent(out) :: emsg
      !
      integer(ik) :: i, n
      integer(ik) :: int_read
      integer(ik) :: rstat
      real(rk) :: real_read
      character(len=*), parameter :: w1 = 'Preprocess: reinitalization &
      &equation - '
      character(len=cl) :: inp_str
      character(len=cl), allocatable :: read_arg(:), arg_val(:)
      !
      write(stdout,'(a)') 'Reading reinitalization parameters ...'
      call set_reinitalization() ! default reinitalization parameters
      ra: do
        call read_input_string(inp_str,rstat)
        select case(rstat)
        case( 0 )
          exit ra
        case( 1 ); cycle ra
        case( 2 )
          backspace(inp_file)
          exit ra
        case( 3 ); continue
        end select
        ! tokenize arguments
        call tokenize(inp_str,comma,read_arg,esta,emsg)
        if ( esta /= 0 ) return
        ! Check arguments
        n = size(read_arg)
        carg: do i = 1, n
          call tokenize(read_arg(i),equal,arg_val,esta,emsg)
          if ( esta /= 0 ) return
          select case( arg_val(1) )
            case ('number of reinitalization iterations')
              call str2i(arg_val(2),int_read,esta,emsg)
              if ( esta /= 0 ) return
              call set_reinitalization(num_reinit_in=int_read)
            case ( 'number of convergence iterations')
              call str2i(arg_val(2),int_read,esta,emsg)
              if ( esta /= 0 ) return
              call set_reinitalization(num_conv_iter_in=int_read)
            case ('correction to timestep')
              call str2r(arg_val(2),real_read,esta,emsg)
              if ( esta /= 0 ) return
              call set_reinitalization(alpha_in=real_read)
            case ('diffusion coefficient')
              call str2r(arg_val(2),real_read,esta,emsg)
              if ( esta /= 0 ) return
              call set_reinitalization(c_in=real_read)
            case ('enforce dirichlet boundary coefficient')
              call str2r(arg_val(2),real_read,esta,emsg)
              if ( esta /= 0 ) return
              call set_reinitalization(rho_in=real_read)
            case ('signed distance function tolerance')
              call str2r(arg_val(2),real_read,esta,emsg)
               if ( esta /= 0 ) return
              call set_reinitalization(sign_dist_tol_in=real_read)
            case default
              esta = -1
              emsg =  w1 // 'unknown value >' // trim(arg_val(1)) // '<'
              exit carg
          end select
        end do carg
        !
      end do ra
      !
    end subroutine read_reinitalization
!*****************************************************************************80
! Input string read in keywords
!*****************************************************************************80
  subroutine read_input_string(inp_str,rstat)
    character(len=*), intent(out) :: inp_str
    integer(ik), intent(out) :: rstat
    !
    integer(ik) :: esta
    ! Read
    read(unit=inp_file,fmt='(a)',iostat=esta) inp_str
    if( is_iostat_end(esta) ) then
      rstat = 0 ! End of file
      end_of_file = .true.
      return
    end if
    ! Adjust
    call to_lower(inp_str)
    inp_str = adjustl(inp_str)
    ! Check for empty line
    if( len_trim(inp_str) == 0 ) then
      rstat = 1 ! Empty line - cycle
      return
    end if
    ! Check if comment
    if( inp_str(1:2) == '**' ) then
      rstat = 1 ! Comment - cycle
      return
    end if
    ! Keyword and keyword arguments
    if( inp_str(1:1) == '*' ) then
      rstat = 2 ! Keyword
      return
    else ! Keyword arguments
      rstat = 3
      return
    end if
    !
  end subroutine read_input_string
!*****************************************************************************80
end module read_input
