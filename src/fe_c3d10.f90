module fe_c3d10
!*****************************************************************************80
  use types
  use ll
  use point, only: dom, point_3d_t
  use general_routines, only: exclude, real_interval, i2str, r2str
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
    procedure, nopass :: n_mtx_sca => n_mtx_sca_c3d10
    procedure, nopass :: dn_dxi_mtx_sca => dn_dxi_mtx_sca_c3d10
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
! N matrix c3d10
!*****************************************************************************80
  pure subroutine n_mtx_sca_c3d10(gp_num,xi_coo_pnt,n_mtx_sca,esta,emsg)
    integer(ik), optional, intent(in) :: gp_num
    type(point_3d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), intent(out) :: n_mtx_sca(:)
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
      if ( .not. size(n_mtx_sca)==nelnod ) then
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
    n_mtx_sca = [ &
    & lambda * ( 2.0_rk * lambda - 1.0_rk ), &
    & xi(1) * ( 2.0_rk * xi(1) - 1.0_rk ), &
    & xi(2) * ( 2.0_rk * xi(2) - 1.0_rk ), &
    & xi(3) * ( 2.0_rk * xi(3) - 1.0_rk ), &
    & 4.0_rk * xi(1) * lambda, &
    & 4.0_rk * xi(1) * xi(2),  &
    & 4.0_rk * xi(2) * lambda, &
    & 4.0_rk * xi(3) * lambda, &
    & 4.0_rk * xi(1) * xi(3),  &
    & 4.0_rk * xi(2) * xi(3) ]
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine n_mtx_sca_c3d10
!*****************************************************************************80
! dN dXi matrix c3d10
!*****************************************************************************80
  pure subroutine dn_dxi_mtx_sca_c3d10(gp_num,xi_coo_pnt,dn_dxi_mtx_sca,esta,&
  &emsg)
    integer(ik), optional, intent(in) :: gp_num
    type(point_3d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), intent(out) :: dn_dxi_mtx_sca(:,:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    real(rk) :: xi(dom)
    logical(lk) :: gp_num_p, xi_coo_p
    !
    gp_num_p = present(gp_num); xi_coo_p = present(xi_coo_pnt)
    if ( debug ) then
      if ( .not. exclude(gp_num_p,xi_coo_p) ) then
        esta = -1
        emsg ='dN dXi matrix scalar tet: Gauss points and xi coordinates must&
        & be excluded'
        return
      end if
      if ( (.not. size(dn_dxi_mtx_sca,dim=1)==dom) .and. &
      & (.not. size(dn_dxi_mtx_sca,dim=2)==nelnod) ) then
        esta = -1
        emsg ='dN dXi matrix scalar tet: matrix size incorrect'
        return
      end if
    end if
    ! Set xi value
    if( gp_num_p ) then
      if ( debug ) then
        if( gp_num > ngp ) then
          esta = -1
          emsg ='dN dXi matrix scalar tet: max five Gauss points'
          return
        end if
      end if
      xi = gp(gp_num,1:dom)
    else if( xi_coo_p ) then
      if ( debug ) then
        do i = 1, dom
          if ( .not. real_interval(0.0e0_rk,1.0e0_rk,xi_coo_pnt%x(i))) then
            esta = -1
            emsg ='dN dXi matrix scalar tet: xi coo vector component number '&
            & // trim(i2str(i)) // ' out of bounds'
            return
          end if
        end do
      end if
      xi = xi_coo_pnt%x
    end if
    !
    dn_dxi_mtx_sca = 0.0_rk
    !
    dn_dxi_mtx_sca(1,1) = 1.0_rk-4.0_rk*(1.0_rk-xi(1)-xi(2)-xi(3))
    dn_dxi_mtx_sca(2,1) = 1.0_rk-4.0_rk*(1.0_rk-xi(1)-xi(2)-xi(3))
    dn_dxi_mtx_sca(3,1) = 1.0_rk-4.0_rk*(1.0_rk-xi(1)-xi(2)-xi(3))
    !
    dn_dxi_mtx_sca(1,2) = -1.0_rk+4.0_rk*xi(1)
    !
    dn_dxi_mtx_sca(2,3) = -1.0_rk+4.0_rk*xi(2)
    !
    dn_dxi_mtx_sca(3,4) = -1.0_rk+4.0_rk*xi(3)
    !
    dn_dxi_mtx_sca(1,5) = 4.0_rk*xi(1)+4.0_rk*(1.0_rk-xi(1)-xi(2)-xi(3))
    dn_dxi_mtx_sca(2,5) = 4.0_rk*xi(1)
    dn_dxi_mtx_sca(3,5) = 4.0_rk*xi(1)
    !
    dn_dxi_mtx_sca(1,6) = 4.0_rk*xi(2)
    dn_dxi_mtx_sca(2,6) = 4.0_rk*xi(1)
    !
    dn_dxi_mtx_sca(1,7) = -4.0_rk*xi(2)
    dn_dxi_mtx_sca(2,7) = -4.0_rk*xi(2)+4.0_rk*(1.0_rk-xi(1)-xi(2)-xi(3))
    dn_dxi_mtx_sca(3,7) = -4.0_rk*xi(2)
    !
    dn_dxi_mtx_sca(1,8) = -4.0_rk*xi(2)
    dn_dxi_mtx_sca(2,8) = -4.0_rk*xi(2)
    dn_dxi_mtx_sca(3,8) = -4.0_rk*(1.0_rk-xi(1)-xi(2)-xi(3))-4.0_rk*xi(3)
    !
    dn_dxi_mtx_sca(1,9) = 4.0_rk*xi(3)
    dn_dxi_mtx_sca(3,9) = 4.0_rk*xi(1)
    !
    dn_dxi_mtx_sca(2,10) = 4.0_rk*xi(3)
    dn_dxi_mtx_sca(3,10) = 4.0_rk*xi(2)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine dn_dxi_mtx_sca_c3d10
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
