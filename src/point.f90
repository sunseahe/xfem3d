module point
  use types
  use ll
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: dom = 3
!*****************************************************************************80
  type :: point_3d_t
    real(rk) :: x(dom) = 0.0_rk
  contains
    procedure :: subtract_pnt
    procedure :: dot_pnt
    procedure :: assign_pnt
    procedure :: norm_pnt
    generic :: operator(-) => subtract_pnt
    generic :: operator(.dot.) => dot_pnt
    generic :: operator(.norm.) => norm_pnt
    generic :: assignment(=) => assign_pnt
    procedure :: write => write_pnt
  end type point_3d_t
  type, extends(list) :: point_3d_t_ll
    private
  contains
    procedure :: add => add_point_3d_t_ll
    procedure :: fill_array => fill_array_point_3d_t_ll
  end type point_3d_t_ll
!*****************************************************************************80
  type(point_3d_t), parameter :: zero_pnt = point_3d_t([0.0_rk,0.0_rk, &
  &0.0_rk])
!*****************************************************************************80
  public :: dom, zero_pnt, point_3d_t, point_3d_t_ll
!*****************************************************************************80
contains
!*****************************************************************************80
  pure function subtract_pnt(self,other) result(res)
    class(point_3d_t), intent(in) :: self
    type(point_3d_t), intent(in) :: other
    type(point_3d_t) :: res
    real(rk) :: sub(dom)
    sub = self%x - other%x
    res = point_3d_t(sub)
  end function subtract_pnt
!*****************************************************************************80
  pure function dot_pnt(self,other) result(res)
    class(point_3d_t), intent(in) :: self
    type(point_3d_t), intent(in) :: other
    real(rk) :: res
    res = sum(self%x * other%x)
  end function dot_pnt
!*****************************************************************************80
  pure function norm_pnt(self,other) result(res)
    class(point_3d_t), intent(in) :: self
    type(point_3d_t), intent(in) :: other
    real(rk) :: res
    res = sqrt(sum(self%x - other%x))
  end function norm_pnt
!*****************************************************************************80
  pure subroutine assign_pnt(self,other)
    class(point_3d_t), intent(inout) :: self
    type(point_3d_t), intent(in) :: other
    self%x = other%x
  end subroutine assign_pnt
!*****************************************************************************80
  subroutine write_pnt(self)
    class(point_3d_t), intent(in) :: self
    character(len=*), parameter :: pw = '(a,3(' // &
    &es//',:,","))'
    write(stdout,pw,advance='no') 'Point coordinates: (', self%x
    write(stdout,'(a)') ' )'
  end subroutine write_pnt
!*****************************************************************************80
  pure subroutine add_point_3d_t_ll(self,arg,esta,emsg)
    class(point_3d_t_ll), intent(inout) :: self
    type(point_3d_t), intent(in) :: arg
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,esta,emsg)
    !
  end subroutine add_point_3d_t_ll
!*****************************************************************************80
  pure subroutine fill_array_point_3d_t_ll(self,arr,esta,emsg)
    class(point_3d_t_ll), intent(inout) :: self
    type(point_3d_t), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      esta = -1
      emsg = 'Fill array point 3d: list is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=esta,errmsg=emsg)
    if( esta /= 0 ) return
    call self%reset()
    do i = 1, self%get_nitem()
      if(allocated(curr)) then
        deallocate(curr,stat=esta,errmsg=emsg)
        if( esta /= 0 ) return
      end if
      allocate(curr,source=self%get_current(),stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      select type(curr)
      type is( point_3d_t )
        arr(i) = curr
      class default
        esta = -1
        emsg = 'Fill array point 3d: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
  end subroutine fill_array_point_3d_t_ll
!*****************************************************************************80
end module point
