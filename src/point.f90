module point
!*****************************************************************************80
  use types, only: ik, rk, cl, es, stdout
!*****************************************************************************80
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
    procedure :: norm_pnt
    generic :: operator(-) => subtract_pnt
    generic :: operator(.dot.) => dot_pnt
    generic :: operator(.norm.) => norm_pnt
    procedure :: write => write_pnt
  end type point_3d_t
!*****************************************************************************80
  type(point_3d_t), parameter :: zero_pnt = point_3d_t([0.0_rk,0.0_rk, &
  &0.0_rk])
!*****************************************************************************80
  public :: dom, zero_pnt, point_3d_t, add_point_3d
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
  subroutine write_pnt(self)
    class(point_3d_t), intent(in) :: self
    character(len=*), parameter :: pw = '(a,3(' // &
    &es//',:,","))'
    write(stdout,pw,advance='no') 'Point coordinates: (', self%x
    write(stdout,'(a)') ' )'
  end subroutine write_pnt
!*****************************************************************************80
! Add point
!*****************************************************************************80
  pure subroutine add_point_3d(vector,element,esta,emsg)
    type(point_3d_t), allocatable, intent(inout) :: vector(:)
    type(point_3d_t), intent(in) :: element
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: n
    type(point_3d_t), allocatable :: tmp_vector(:)
    !
    if ( allocated(vector) ) then
      n = size(vector,dim=1)
      allocate(tmp_vector(n+1),stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
      tmp_vector(1:n) = vector
      tmp_vector(n+1) = element
      call move_alloc(tmp_vector,vector)
    else
      allocate(vector(1),source=element,stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
    end if
    !
  end subroutine add_point_3d
!*****************************************************************************80
end module point
