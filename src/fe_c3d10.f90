module fe_c3d10
  use types
  use fe_c3d4, only: dom, point_3d_t
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: nnod = 10
!*****************************************************************************80
  type :: c3d10_t
    type(point_3d_t) :: nodes(nnod)
    integer(ik) :: connectivity(nnod) = 0
  contains
    procedure :: write => write_c3d10
  end type c3d10_t
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine write_c3d10(self)
    class(c3d10_t), intent(in) :: self
    !
    integer(ik) :: i
    !
    write(stdout,'(a)') 'C3d10 coordinates:'
    do i = 1, nnod
      write(stdout,'(i0,a)',advance='no') i, '. '
      call self%nodes(i)%write()
    end do
    !
  end subroutine write_c3d10
!*****************************************************************************80
end module fe_c3d10
