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
  integer(ik), parameter :: ngp = 5
  real(rk), parameter :: oh = 1.0_rk / 2.0_rk, oq = 1.0_rk / 4.0_rk, &
  & os = 1.0_rk / 6.0_rk, tf = 3.0_rk / 4.0e1_rk
  real(rk), parameter :: gp(ngp,dom) = reshape([ &
  & oq, oq, oq, &
  & os, os, os, &
  & os, os, oh, &
  & os, oh, os, &
  & oh, os, os ], order=[2,1], shape=[ngp,dom])
  real(rk), parameter :: w(ngp) = [ - 2.0_rk / 1.5e1_rk, tf, tf, tf, tf ]
!*****************************************************************************80
  type :: c3d10_t
    type(point_3d_t) :: nodes(nelnod)
    integer(ik) :: connectivity(nelnod) = 0
  contains
    procedure :: write => write_c3d10
    procedure :: get_node_coo => get_node_coo_c3d10
    procedure, nopass :: n_mtx => n_mtx_c3d10
    procedure, nopass :: dn_dxi_mtx => dn_dxi_mtx_c3d10
    procedure :: gradient => gradient_c3d10
    procedure :: get_connectivity => get_connectivity_c3d10
  end type c3d10_t
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
! N matrix c3d10
!*****************************************************************************80
  pure subroutine n_mtx_c3d10(gp_num,xi_coo_pnt,n_mtx,esta,emsg)
    integer(ik), optional, intent(in) :: gp_num
    type(point_3d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), intent(out) :: n_mtx(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    real(rk) :: xi(dom), lambda
    logical(lk) :: gp_num_p, xi_coo_p
    !
    gp_num_p = present(gp_num); xi_coo_p = present(xi_coo_pnt)
    if ( debug ) then
      if ( .not. exclude(gp_num_p,xi_coo_p) ) then
        esta = -1
        emsg ='N matrix scalar tet: Gauss points and xi coordinates must&
        & be excluded'
        return
      end if
      if ( .not. size(n_mtx)==nelnod ) then
        esta = -1
        emsg ='N matrix scalar tet: N matrix size incorrect'
        return
      end if
    end if
    ! Set xi value
    if( gp_num_p ) then
      if ( debug ) then
        if( gp_num > ngp ) then
          esta = -1
          emsg ='N matrix scalar tet: max five Gauss points'
          return
        end if
      end if
      xi = gp(gp_num,1:dom)
    else if( xi_coo_p ) then
      if ( debug ) then
        do i = 1, dom
          if ( .not. real_interval(0.0e0_rk,1.0e0_rk,xi_coo_pnt%x(i))) then
            esta = -1
            emsg ='N matrix scalar tet: xi coo vector component number '&
            & // trim(i2str(i)) // ' out of bounds'
            return
          end if
        end do
      end if
      xi = xi_coo_pnt%x
    end if
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
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine n_mtx_c3d10
!*****************************************************************************80
! dN dXi matrix c3d10
!*****************************************************************************80
  pure subroutine dn_dxi_mtx_c3d10(gp_num,xi_coo_pnt,dn_dxi_mtx,esta,&
  &emsg)
    integer(ik), optional, intent(in) :: gp_num
    type(point_3d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), intent(out) :: dn_dxi_mtx(:,:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    real(rk) :: xi(dom), lambda
    logical(lk) :: gp_num_p, xi_coo_p
    !
    gp_num_p = present(gp_num); xi_coo_p = present(xi_coo_pnt)
    if ( debug ) then
      if ( .not. exclude(gp_num_p,xi_coo_p) ) then
        esta = -1
        emsg ='dN dXi matrix tet: Gauss points and xi coordinates must&
        & be excluded'
        return
      end if
      if ( .not.size_mtx(dn_dxi_mtx,dom,nelnod) ) then
        esta = -1
        emsg ='dN dXi matrix tet: matrix size incorrect'
        return
      end if
    end if
    ! Set xi value
    if( gp_num_p ) then
      if ( debug ) then
        if( gp_num > ngp ) then
          esta = -1
          emsg ='dN dXi matrix tet: max five Gauss points'
          return
        end if
      end if
      xi = gp(gp_num,1:dom)
    else if( xi_coo_p ) then
      if ( debug ) then
        do i = 1, dom
          if ( .not. real_interval(0.0e0_rk,1.0e0_rk,xi_coo_pnt%x(i))) then
            esta = -1
            emsg ='dN dXi matrix tet: xi coo vector component number '&
            & // trim(i2str(i)) // ' out of bounds'
            return
          end if
        end do
      end if
      xi = xi_coo_pnt%x
    end if
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
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine dn_dxi_mtx_c3d10
!*****************************************************************************80
! Jacobian matrix c3d10
!*****************************************************************************80
  pure subroutine gradient_c3d10(self,gp_num,xi_coo_pnt,jac_mtx,det_jac &
  &,inv_jac_mtx,b_mtx,esta,emsg)
    class(c3d10_t), intent(in) :: self
    integer(ik), optional, intent(in) :: gp_num
    type(point_3d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), optional, intent(out) :: jac_mtx(:,:)
    real(rk), optional, intent(out) :: det_jac
    real(rk), optional, intent(out) :: inv_jac_mtx(:,:)
    real(rk), optional, intent(out) :: b_mtx(:,:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    real(rk) :: dn_dxi_mtx(dom,nelnod), el_coo(nelnod,dom)
    real(rk) :: jac_mtx_int(dom,dom), det_jac_int, inv_jac_mtx_int(dom,dom)
    ! Checks
    if ( debug ) then
      if ( present(jac_mtx) ) then
        if ( .not.size_mtx(jac_mtx,dom,dom) ) then
          esta = -1
          emsg ='Jacobian: matrix jac_mtx size incorrect'
          return
        end if
      end if
      if ( present(inv_jac_mtx) ) then
        if ( .not.size_mtx(inv_jac_mtx,dom,dom) ) then
          esta = -1
          emsg ='Jacobian: matrix inv_jac_mtx size incorrect'
          return
        end if
      end if
      if ( present(b_mtx) ) then
        if ( .not.size_mtx(b_mtx,dom,nelnod) ) then
          esta = -1
          emsg ='B matrix size incorrect'
          return
        end if
      end if
    end if
    ! Jacobian matrix
    call self%dn_dxi_mtx(gp_num=gp_num,xi_coo_pnt=xi_coo_pnt,&
    &dn_dxi_mtx=dn_dxi_mtx,esta=esta,emsg=emsg)
    if( esta /= 0 ) return
    do i = 1, nelnod
      el_coo(i,1:dom) = self%nodes(i)%x(1:dom)
    end do
    call gemm(dn_dxi_mtx,el_coo,jac_mtx_int)
    if ( present(jac_mtx) ) jac_mtx = jac_mtx_int
    ! Jacobian determinant
    associate( j => jac_mtx_int )
    det_jac_int = - j(1,3)*j(2,2)*j(3,1) + j(1,2)*j(2,3)*j(3,1) &
    &             + j(1,3)*j(2,1)*j(3,2) - j(1,1)*j(2,3)*j(3,2) &
    &             - j(1,2)*j(2,1)*j(3,3) + j(1,1)*j(2,2)*j(3,3)
    end associate
    if ( debug ) then
      if ( det_jac_int <= 0.0_rk ) then
        esta = -1
        emsg = 'Jacobian: determinant is zero'
        return
      end if
    end if
    if ( present(det_jac) ) det_jac = det_jac_int
    ! Jacobian inverse
    associate( j => jac_mtx_int )
    inv_jac_mtx_int = ( 1.0_rk / det_jac_int ) * reshape( [ &
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
    end associate
    if ( present(inv_jac_mtx) ) inv_jac_mtx = inv_jac_mtx_int
    !
    if ( present(b_mtx) ) call gemm(inv_jac_mtx_int,dn_dxi_mtx,b_mtx)
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine gradient_c3d10
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
! Get connectivity
!*****************************************************************************80
  pure subroutine get_connectivity_c3d10(self,connectivity,esta,emsg)
    class(c3d10_t), intent(in) :: self
    integer(ik), intent(out) :: connectivity(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    if ( debug ) then
      if ( size(connectivity) /= nelnod ) then
        esta = -1
        emsg = 'Connectivity vector wrong size'
      end if
    end if
    !
    connectivity = self%connectivity
    !
  end subroutine get_connectivity_c3d10
end module fe_c3d10
