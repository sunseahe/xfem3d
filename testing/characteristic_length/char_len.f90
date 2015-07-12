module cp2d4
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: nelnod = 4
  integer(ik), parameter :: ngp = 4
!*****************************************************************************80
  real(rk), parameter :: gp(ngp,dom) = reshape() 
  real(rk), parameter :: w(ngp) = []
!*****************************************************************************80
  type :: cp2d4_t
    type(point_2d_t) :: nodes(nelnod) 
    integer(ik) :: connectivity(nelnod) = 0
  contains
    procedure :: dn_dxi_mtx => dn_dxi_mtx_cp2d4
    procedure :: gradient => gradient_cp2d4
  end type cp2d4_t
!*****************************************************************************80
contains
!*****************************************************************************80
  pure subroutine dn_dxi_mtx_cp2d4(gp_num,xi_coo_pnt,dn_dxi_mtx,esta,&
  &emsg)
    integer(ik), optional, intent(in) :: gp_num
    type(point_2d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), intent(out) :: dn_dxi_mtx(:,:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
  end subroutine dn_dxi_mtx_cp2d4
end module cp2d4
program main
  implicit none

end program main
