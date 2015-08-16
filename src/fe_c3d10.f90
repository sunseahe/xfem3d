module fe_c3d10
!*****************************************************************************80
  use blas95, only: gemm
  use types, only: ik, rk, cl, lk, stdout, debug
  use ll, only: list
  use point, only: dom, point_3d_t
  use general_routines, only:  size_mtx, exclude, real_interval, i2str, r2str
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: nelnod = 10
  integer(ik), parameter :: ngp = 15 !4
!  real(rk), parameter :: a = 0.58541020_rk, b = 0.13819660, &
!  &c = 1.0_rk / 2.4e1_rk
!  real(rk), parameter :: gp(ngp,dom) = reshape([ &
!  & a, b, b, &
!  & b, a, b, &
!  & b, b, a, &
!  & b, b, b  ], order=[2,1], shape=[ngp,dom])
!  real(rk), parameter :: w(ngp) = [ c, c, c, c ]
   real(rk), parameter :: xa(ngp) =[ 0.2500000000000000_rk, &
   & 0.0000000000000000_rk, 0.3333333333333333_rk, 0.3333333333333333_rk, &
   & 0.3333333333333333_rk, 0.7272727272727273_rk, 0.0909090909090909_rk, &
   & 0.0909090909090909_rk, 0.0909090909090909_rk, 0.4334498464263357_rk, &
   & 0.0665501535736643_rk, 0.0665501535736643_rk, 0.0665501535736643_rk, &
   & 0.4334498464263357_rk, 0.4334498464263357_rk ]
   real(rk), parameter :: ya(ngp) = [ 0.2500000000000000_rk, &
   & 0.3333333333333333_rk, 0.3333333333333333_rk, 0.3333333333333333_rk, &
   & 0.0000000000000000_rk, 0.0909090909090909_rk, 0.0909090909090909_rk, &
   & 0.0909090909090909_rk, 0.7272727272727273_rk, 0.0665501535736643_rk, &
   & 0.4334498464263357_rk, 0.0665501535736643_rk, 0.4334498464263357_rk, &
   & 0.0665501535736643_rk, 0.4334498464263357_rk ];
   real(rk), parameter :: za(ngp) =[ 0.2500000000000000_rk, &
   & 0.3333333333333333_rk, 0.3333333333333333_rk, 0.0000000000000000_rk, &
   & 0.3333333333333333_rk, 0.0909090909090909_rk, 0.0909090909090909_rk, &
   & 0.7272727272727273_rk, 0.0909090909090909_rk, 0.0665501535736643_rk, &
   & 0.0665501535736643_rk, 0.4334498464263357_rk, 0.4334498464263357_rk, &
   & 0.4334498464263357_rk, 0.0665501535736643_rk ];
   real(rk), parameter :: gp(ngp,dom) = reshape([xa,ya,za], order=[2,1], &
   & shape=[ngp,dom])
   real(rk), parameter :: w(ngp)  = [ 0.1817020685825351_rk, &
   & 0.0361607142857143_rk, 0.0361607142857143_rk, 0.0361607142857143_rk, &
   & 0.0361607142857143_rk, 0.0698714945161738_rk, 0.0698714945161738_rk, &
   & 0.0698714945161738_rk, 0.0698714945161738_rk, 0.0656948493683187_rk, &
   & 0.0656948493683187_rk, 0.0656948493683187_rk, 0.0656948493683187_rk, &
   & 0.0656948493683187_rk, 0.0656948493683187_rk ] / 6.0_rk
!*****************************************************************************80
  type :: c3d10_t
    type(point_3d_t) :: nodes(nelnod)
    integer(ik) :: connectivity(nelnod) = 0
  contains
    procedure :: write => write_c3d10
    procedure :: get_node_coo
    procedure :: main_values
  end type c3d10_t
!*****************************************************************************80
  type, extends(list) :: c3d10_t_ll
    private
  contains
    procedure :: add => add_c3d10_t_ll
    procedure :: fill_array => fill_array_c3d10_t_ll
  end type c3d10_t_ll
!*****************************************************************************80
  public :: nelnod, ngp, w, c3d10_t, c3d10_t_ll
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
  pure function get_node_coo(self,n) result(coo)
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
  end function get_node_coo
!*****************************************************************************80
! N matrix values
!*****************************************************************************80
  pure subroutine n_values(xi,n_mtx)
    real(rk), intent(in) :: xi(:)
    real(rk), intent(out) :: n_mtx(:)
    !
    real(rk) :: lambda
    !
    lambda = 1.0_rk - xi(1) - xi(2) - xi(3)
    n_mtx(1) = lambda * ( 2.0_rk * lambda - 1.0_rk )
    n_mtx(2) = xi(1) * ( 2.0_rk * xi(1) - 1.0_rk )
    n_mtx(3) = xi(2) * ( 2.0_rk * xi(2) - 1.0_rk )
    n_mtx(4) = xi(3) * ( 2.0_rk * xi(3) - 1.0_rk )
    n_mtx(5) = 4.0_rk * xi(1) * lambda
    n_mtx(6) = 4.0_rk * xi(1) * xi(2)
    n_mtx(7) = 4.0_rk * xi(2) * lambda
    n_mtx(8) = 4.0_rk * xi(3) * lambda
    n_mtx(9) = 4.0_rk * xi(1) * xi(3)
    n_mtx(10) = 4.0_rk * xi(2) * xi(3)
    !
  end subroutine n_values
!*****************************************************************************80
! dN dXi matrix values
!*****************************************************************************80
  pure subroutine dn_dxi_mtx_values(xi,dn_dxi_mtx)
    real(rk), intent(in) :: xi(:)
    real(rk), intent(out) :: dn_dxi_mtx(:,:)
    !
    real(rk) :: lambda
    !
    dn_dxi_mtx = 0.0_rk
    lambda = 1.0_rk - xi(1) - xi(2) - xi(3)
    !
    dn_dxi_mtx(1,1) = 1.0_rk-4.0_rk*lambda
    dn_dxi_mtx(2,1) = 1.0_rk-4.0_rk*lambda
    dn_dxi_mtx(3,1) = 1.0_rk-4.0_rk*lambda
    !
    dn_dxi_mtx(1,2) = -1.0_rk+4.0_rk*xi(1)
    !
    dn_dxi_mtx(2,3) = -1.0_rk+4.0_rk*xi(2)
    !
    dn_dxi_mtx(3,4) = -1.0_rk+4.0_rk*xi(3)
    !
    dn_dxi_mtx(1,5) = -4.0_rk*xi(1)+4.0_rk*lambda
    dn_dxi_mtx(2,5) = -4.0_rk*xi(1)
    dn_dxi_mtx(3,5) = -4.0_rk*xi(1)
    !
    dn_dxi_mtx(1,6) = 4.0_rk*xi(2)
    dn_dxi_mtx(2,6) = 4.0_rk*xi(1)
    !
    dn_dxi_mtx(1,7) = -4.0_rk*xi(2)
    dn_dxi_mtx(2,7) = -4.0_rk*xi(2)+4.0_rk*lambda
    dn_dxi_mtx(3,7) = -4.0_rk*xi(2)
    !
    dn_dxi_mtx(1,8) = -4.0_rk*xi(3)
    dn_dxi_mtx(2,8) = -4.0_rk*xi(3)
    dn_dxi_mtx(3,8) = 4.0_rk*lambda-4.0_rk*xi(3)
    !
    dn_dxi_mtx(1,9) = 4.0_rk*xi(3)
    dn_dxi_mtx(3,9) = 4.0_rk*xi(1)
    !
    dn_dxi_mtx(2,10) = 4.0_rk*xi(3)
    dn_dxi_mtx(3,10) = 4.0_rk*xi(2)
    !
  end subroutine dn_dxi_mtx_values
!*****************************************************************************80
! Jacobian matrix c3d10
!*****************************************************************************80
  pure subroutine calc_det_jac(j,det_jac)
    real(rk), intent(in) :: j(:,:)
    real(rk), intent(out) :: det_jac
    !
    det_jac = - j(1,3)*j(2,2)*j(3,1) + j(1,2)*j(2,3)*j(3,1) &
    &         + j(1,3)*j(2,1)*j(3,2) - j(1,1)*j(2,3)*j(3,2) &
    &         - j(1,2)*j(2,1)*j(3,3) + j(1,1)*j(2,2)*j(3,3)
    !
  end subroutine calc_det_jac
  pure subroutine calc_inv_jac_mtx(j,det_jac,inv_jac_mtx)
    real(rk), intent(in) :: j(:,:)
    real(rk), intent(in) :: det_jac
    real(rk), intent(out) :: inv_jac_mtx(:,:)
    !
    inv_jac_mtx = ( 1.0_rk / det_jac ) * reshape( [ &
    & - j(2,3) * j(3,2) + j(2,2) * j(3,3),     &
    & + j(1,3) * j(3,2) - j(1,2) * j(3,3),     &
    & - j(1,3) * j(2,2) + j(1,2) * j(2,3),     &
    & + j(2,3) * j(3,1) - j(2,1) * j(3,3),     &
    & - j(1,3) * j(3,1) + j(1,1) * j(3,3),     &
    & + j(1,3) * j(2,1) - j(1,1) * j(2,3),     &
    & - j(2,2) * j(3,1) + j(2,1) * j(3,2),     &
    & + j(1,2) * j(3,1) - j(1,1) * j(3,2),     &
    & - j(1,2) * j(2,1) + j(1,1) * j(2,2)      &
    & ], shape=[3,3],order=[2,1] )
    !
  end subroutine calc_inv_jac_mtx
  pure subroutine main_values(self,gp_num,xi_coo_pnt,n_mtx,jac_mtx,det_jac &
  &,inv_jac_mtx,b_mtx)
    class(c3d10_t), intent(in) :: self
    integer(ik), optional, intent(in) :: gp_num
    type(point_3d_t), optional,  intent(in) :: xi_coo_pnt
    real(rk), optional, intent(out) :: n_mtx(:)
    real(rk), optional, intent(out) :: jac_mtx(:,:)
    real(rk), optional, intent(out) :: det_jac
    real(rk), optional, intent(out) :: inv_jac_mtx(:,:)
    real(rk), optional, intent(out) :: b_mtx(:,:)
    !
    integer(ik) :: i
    real(rk) :: xi(dom)
    real(rk) :: dn_dxi_mtx(dom,nelnod), el_coo(nelnod,dom)
    real(rk) :: jac_mtx_int(dom,dom), det_jac_int, inv_jac_mtx_int(dom,dom)
    logical(lk) :: a1, a2, a3, a4, a5
    !
    a1 = present(n_mtx); a2 = present(jac_mtx)
    a3 = present(det_jac); a4 = present(inv_jac_mtx)
    a5 = present(b_mtx)
    !
    if ( present(gp_num) ) then
      xi = gp(gp_num,1:dom)
    else if (present (xi_coo_pnt) ) then
      xi = xi_coo_pnt%x
    else
      return
    end if
    ! N matrix
    if ( a1 ) call n_values(xi,n_mtx)
    if ( .not.a2 .and. .not.a3 .and. &
    &    .not.a4 .and. .not.a5  ) return
    ! Jacobian matrix
    call dn_dxi_mtx_values(xi,dn_dxi_mtx)
    do i = 1, nelnod
      el_coo(i,1:dom) = self%nodes(i)%x(1:dom)
    end do
    call gemm(dn_dxi_mtx,el_coo,jac_mtx_int)
    if ( a2 ) jac_mtx = jac_mtx_int
    ! Jacobian determinant
    call calc_det_jac(jac_mtx_int,det_jac_int)
    if ( a3 ) det_jac = det_jac_int
    if ( .not.a4 .and. .not.a5 ) return
    ! Jacobian inverse
    call calc_inv_jac_mtx(jac_mtx_int,det_jac_int,inv_jac_mtx_int)
    if ( a4 ) inv_jac_mtx = inv_jac_mtx_int
    if ( .not. a5 ) return
    ! B matrix
    if ( a5 ) call gemm(inv_jac_mtx_int,dn_dxi_mtx,b_mtx)
    !
  end subroutine main_values
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
      emsg = 'Fill array c3d10: list is empty'
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
        emsg = 'Fill array c3d10: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
    ! Clean
    call self%clean(esta,emsg)
    if( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine fill_array_c3d10_t_ll
!*****************************************************************************80
end module fe_c3d10
