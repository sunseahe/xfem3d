module read_input
  use types
  use general_routines
  use tokenize_string, only: tokenize, comma, equal
  use tet_volume, only: dom, nver, zero_pnt, point_3d_dat, point_3d_dat_ll, &
  &  tet_dat, tet_dat_ll
  implicit none
  private
!*****************************************************************************80
  integer(ik) :: inp_file = 0
  character(len=cl) :: input_file_name = ''
  !
  character(len=*), parameter :: fe_type = 'c3d4'
  !
  logical(lk) :: end_of_file = .false.
  !
  type(point_3d_dat), allocatable :: nodes(:)
  type(tet_dat), allocatable :: finite_elements(:)
!*****************************************************************************80
  public :: inp_file, input_file_name, read_data, read_input_statistics
  public :: fe_type, nodes, finite_elements
!*****************************************************************************80
  contains
!*****************************************************************************80
! Read data
!*****************************************************************************80
  subroutine read_data(istat,emsg)
    !
    integer(ik), intent(out) :: istat
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
        istat = -1
        emsg = 'Read data: read keyword error'
        return
      end select
      ! tokenize arguments
      call tokenize(inp_str,comma,keyword_arg,istat,emsg)
      if ( istat /= 0 ) return
      keyword = keyword_arg(1)(2:) ! First one, without *
      select case( keyword )
!*****************************************************************************80
      case( 'node' )
        call read_nodes(istat,emsg); if ( istat /= 0 ) return
!*****************************************************************************80
      case( 'element')
        do i = 2, size(keyword_arg)
          sub_keyword = keyword_arg(i)
          call tokenize(sub_keyword,equal,sub_keyword_arg,istat,emsg)
          if ( istat /= 0 ) return
          select case( sub_keyword_arg(1) )
          case( 'type' )
            if ( sub_keyword_arg(2) == fe_type ) then
              fe_type_given = .true.
              call read_connectivity(istat,emsg)
              if ( istat /= 0 ) return
            else
              istat = -1
              emsg = 'Read data: element type >' // trim(sub_keyword_arg(2)) &
              &// '< not supported'
              return
            end if
          case default
            istat = -1
            emsg = 'Read data: unknown sub keyword >' //  &
            &trim(sub_keyword_arg(1)) // '<'
            return
          end select
        end do
        if ( .not. fe_type_given ) then
          istat = -1
          emsg = 'Read data: finite element type not given'
          return
        end if
!*****************************************************************************80
      case default
        istat = -1
        emsg = 'Read data: unknown keyword >' // trim(keyword) // '<'
        return
      end select
    end do ra
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine read_data
!*****************************************************************************80
! Read nodes
!*****************************************************************************80
  subroutine read_nodes(istat,emsg)
    !
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    integer(ik) :: rstat
    real(rk) :: node_coo(dom)
    character(len=cl) :: inp_str
    character(len=cl), allocatable :: read_arg(:)
    type(point_3d_dat) :: node
    type(point_3d_dat_ll) :: nodes_ll
    !
    write(stdout,'(a)') 'Reading nodes ...'
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
      call tokenize(inp_str,comma,read_arg,istat,emsg)
      if ( istat /= 0 ) return
      if ( size(read_arg) /= dom + 1 ) then
        istat = -1
        emsg = 'Read nodes: domain not supported'
        return
      end if
      do i = 2, dom + 1
        call str2r(read_arg(i),node_coo(i-1),istat,emsg)
        if ( istat /= 0 ) return
      end do
      node = point_3d_dat(node_coo)
      call nodes_ll%add(node,istat,emsg); if ( istat /= 0 ) return
    end do ra
    ! Add nodes
    call nodes_ll%fill_array(nodes,istat,emsg); if ( istat /= 0 ) return
    call nodes_ll%clean(istat,emsg); if ( istat /= 0 ) return
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine read_nodes
!*****************************************************************************80
! Read connectivity
!*****************************************************************************80
  subroutine read_connectivity(istat,emsg)
    !
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i, j
    integer(ik) :: rstat
    integer(ik) :: el_conn(nver)
    character(len=cl) :: inp_str
    character(len=cl), allocatable :: read_arg(:)
    type(tet_dat) :: fe_element
    type(tet_dat_ll) :: fe_elements_ll
    !
    write(stdout,'(a)') 'Reading connectivity ...'
    !
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
      call tokenize(inp_str,comma,read_arg,istat,emsg)
      if ( istat /= 0 ) return
      if ( size(read_arg) /= nver + 1 ) then
        istat = -1
        emsg = 'Read connectivity: Number of vertices is incorrect'
        return
      end if
      do i = 2, nver + 1
        call str2i(read_arg(i),el_conn(i-1),istat,emsg)
        if ( istat /= 0 ) return
      end do
      fe_element = tet_dat(vert=zero_pnt,connectivity=el_conn)
      call fe_elements_ll%add(fe_element,istat,emsg)
      if ( istat /= 0 ) return
      !
    end do ra
    ! Create element array
    call fe_elements_ll%fill_array(finite_elements,istat,emsg)
    if ( istat /= 0 ) return
    call fe_elements_ll%clean(istat,emsg); if ( istat /= 0 ) return
    ! Copy nodes to finite elements
    do i = 1, size(finite_elements)
      el_conn = finite_elements(i)%connectivity
      do j = 1, nver
        finite_elements(i)%vert(j) = nodes(el_conn(j))
      end do
    end do
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine read_connectivity
!*****************************************************************************80
! Input string read in keywords
!*****************************************************************************80
  subroutine read_input_string(inp_str,rstat)
    character(len=*), intent(out) :: inp_str
    integer(ik), intent(out) :: rstat
    !
    integer(ik) :: istat
    ! Read
    read(unit=inp_file,fmt='(a)',iostat=istat) inp_str
    if( is_iostat_end(istat) ) then
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
! Read input statistics
!*****************************************************************************80
  subroutine read_input_statistics()
    !
    write(stdout,'(a,i0)') 'Number of nodes is: ', size(nodes)
    write(stdout,'(a,i0)') 'Number of elements is: ', size(finite_elements)
    write(stdout,'(a)') 'Reading time: '
    !
  end subroutine read_input_statistics
!*****************************************************************************80
end module read_input
