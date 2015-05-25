module fe_c3d10
!*****************************************************************************80
  use types
  use ll
  use point, only: dom, point_3d_t
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: nelnod = 10
!*****************************************************************************80
  type :: c3d10_t
    type(point_3d_t) :: nodes(nelnod)
    integer(ik) :: connectivity(nelnod) = 0
  contains
    procedure :: write => write_c3d10
    procedure :: get_node_coo => get_node_coo_c3d10
  end type c3d10_t
  type, extends(list) :: c3d10_t_ll
    private
  contains
    procedure :: add => add_c3d10_t_ll
    procedure :: fill_array => fill_array_c3d10_t_ll
  end type c3d10_t_ll
!*****************************************************************************80
  public :: nelnod, c3d10_t, c3d10_t_ll
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine write_c3d10(self)
    class(c3d10_t), intent(in) :: self
    !
    integer(ik) :: i
    !
    write(stdout,'(a)') 'C3d10 coordinates:'
    do i = 1, nelnod
      write(stdout,'(i0,a)',advance='no') i, '. '
      call self%nodes(i)%write()
    end do
    !
  end subroutine write_c3d10
!*****************************************************************************80
  pure function get_node_coo_c3d10(self,n) result(coo)
    class(c3d10_t), intent(in) :: self
    integer(ik), intent(in) :: n
    real(rk) :: coo(dom)
    !
    if (n > nelnod ) then
      coo = 0.0_rk
      return
    end if
    coo = self%nodes(n)%x
    !
  end function get_node_coo_c3d10
!*****************************************************************************80
  pure subroutine add_c3d10_t_ll(self,arg,esta,emsg)
    class(c3d10_t_ll), intent(inout) :: self
    type(c3d10_t), intent(in) :: arg
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,esta,emsg)
    !
  end subroutine add_c3d10_t_ll
!*****************************************************************************80
  pure subroutine fill_array_c3d10_t_ll(self,arr,esta,emsg)
    class(c3d10_t_ll), intent(inout) :: self
    type(c3d10_t), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      esta = -1
      emsg = 'Fill array tet dat: list is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=esta,errmsg=emsg)
    if( esta /= 0 ) return
    call self%reset()
    do i = 1, self%get_nitem()
      call self%get_current(curr,esta,emsg)
      if( esta /= 0 ) return
      select type(curr)
      type is( c3d10_t )
        arr(i) = curr
      class default
        esta = -1
        emsg = 'Fill array tet dat: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
  end subroutine fill_array_c3d10_t_ll
!*****************************************************************************80
end module fe_c3d10
