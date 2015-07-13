module types
  use iso_fortran_env, only: &
  & int8, &
  & ik => int32, &
  & int64, &
#ifndef DOUBLE
  & rk => real32, &
#else
  & rk => real64, &
#endif
  & stdout => output_unit, &
  & stdin => input_unit, &
  & stderr => error_unit
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: cl = 256, lk = 4
  real(rk), parameter :: eps = epsilon(1.0_rk)
!*****************************************************************************80
#ifndef DOUBLE
  character(len=*), parameter :: es = 'es14.6e3'
#else
  character(len=*), parameter :: es = 'es23.15e3'
#endif
!*****************************************************************************80
#ifndef DEBUG
  logical(lk), parameter :: debug = .false.
#else
  logical(lk), parameter :: debug = .true.
#endif
!*****************************************************************************80
  character(len=1), parameter :: nl = achar(10)
!*****************************************************************************80
  public :: ik, &
  & int8, &
  & int64, &
  & rk, &
  & cl, &
  & lk, &
  & eps, &
  & es, &
  & stdout, &
  & stdin, &
  & stderr, &
  & debug, &
  & nl
!*****************************************************************************80
end module types

module point
  use types, only: ik, rk
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: dom = 2
!*****************************************************************************80
  type :: point_2d_t
    real(rk) :: x(dom) = 0.0_rk
  end type point_2d_t
!*****************************************************************************80
  public :: dom, point_2d_t
contains

!*****************************************************************************80
end module point
module cp2d4
!*****************************************************************************80
  use blas95, only: gemm
  use types, only: ik, rk, lk, cl
  use point, only : dom, point_2d_t
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: nelnod = 4
  integer(ik), parameter :: ngp = 4
!*****************************************************************************80
  real(rk), parameter :: gp(ngp,dom) = reshape( [ &
  & - 0.57735027_rk, - 0.57735027_rk, &
  & + 0.57735027_rk, - 0.57735027_rk, &
  & + 0.57735027_rk, + 0.57735027_rk, &
  & - 0.57735027_rk, + 0.57735027_rk ], shape = [ngp,dom], order = [2,1] )
  real(rk), parameter :: w(ngp) = [ 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk ]
!*****************************************************************************80
  type :: cp2d4_t
    type(point_2d_t) :: nodes(nelnod)
    integer(ik) :: connectivity(nelnod) = 0
  contains
    procedure, nopass :: n_mtx => n_mtx_cp2d4
    procedure, nopass :: dn_dxi_mtx => dn_dxi_mtx_cp2d4
    procedure :: gradient => gradient_cp2d4
  end type cp2d4_t
!*****************************************************************************80
  public :: cp2d4_t, ngp, nelnod
!*****************************************************************************80
contains
!*****************************************************************************80
  pure subroutine n_mtx_cp2d4(gp_num,xi_coo_pnt,n_mtx,esta,emsg)
    integer(ik), optional, intent(in) :: gp_num
    type(point_2d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), intent(out) :: n_mtx(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    real(rk) :: xi(dom)
    logical(lk) :: gp_num_p, xi_coo_p
    !
    gp_num_p = present(gp_num); xi_coo_p = present(xi_coo_pnt)
    if( gp_num_p ) then
      xi = gp(gp_num,1:dom)
    else if( xi_coo_p ) then
      xi = xi_coo_pnt%x
    end if
    !
    n_mtx(1) = 2.5e-1_rk * ( 1.0e0_rk - xi(1) ) * ( 1.0e0_rk - xi(2) )
    n_mtx(2) = 2.5e-1_rk * ( 1.0e0_rk + xi(1) ) * ( 1.0e0_rk - xi(2) )
    n_mtx(3) = 2.5e-1_rk * ( 1.0e0_rk + xi(1) ) * ( 1.0e0_rk + xi(2) )
    n_mtx(4) = 2.5e-1_rk * ( 1.0e0_rk - xi(1) ) * ( 1.0e0_rk + xi(2) )
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine n_mtx_cp2d4
!*****************************************************************************80
  pure subroutine dn_dxi_mtx_cp2d4(gp_num,xi_coo_pnt,dn_dxi_mtx,esta,&
  &emsg)
    integer(ik), optional, intent(in) :: gp_num
    type(point_2d_t), optional, intent(in) :: xi_coo_pnt
    real(rk), intent(out) :: dn_dxi_mtx(:,:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    real(rk) :: xi(dom)
    logical(lk) :: gp_num_p, xi_coo_p
    !
    gp_num_p = present(gp_num); xi_coo_p = present(xi_coo_pnt)
    if( gp_num_p ) then
      xi = gp(gp_num,1:dom)
    else if( xi_coo_p ) then
      xi = xi_coo_pnt%x
    end if
    !
    dn_dxi_mtx(1,1) = 2.5e-1_rk * ( xi(2) - 1.0_rk )
    dn_dxi_mtx(2,1) = 2.5e-1_rk * ( xi(1) - 1.0_rk )
    dn_dxi_mtx(1,2) = 2.5e-1_rk * ( 1.0_rk - xi(2) )
    dn_dxi_mtx(2,2) = 2.5e-1_rk * ( - xi(1) - 1.0_rk )
    dn_dxi_mtx(1,3) = 2.5e-1_rk * ( 1.0_rk + xi(2) )
    dn_dxi_mtx(2,3) = 2.5e-1_rk * ( 1.0_rk + xi(1) )
    dn_dxi_mtx(1,4) = 2.5e-1_rk * ( - xi(2) - 1.0_rk )
    dn_dxi_mtx(2,4) = 2.5e-1_rk * ( 1.0_rk - xi(1) )
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine dn_dxi_mtx_cp2d4
!*****************************************************************************80
  pure subroutine gradient_cp2d4(self,gp_num,xi_coo_pnt,jac_mtx,det_jac &
  &,inv_jac_mtx,b_mtx,esta,emsg)
    class(cp2d4_t), intent(in) :: self
    integer(ik), optional, intent(in) :: gp_num
    type(point_2d_t), optional, intent(in) :: xi_coo_pnt
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
    det_jac_int = j(1,1) * j(2,2) - j(1,2) * j(2,1)
    end associate
    if ( det_jac_int <= 0.0_rk ) then
      esta = -1
      emsg = 'Jacobian: determinant is zero'
      return
    end if
    if ( present(det_jac) ) det_jac = det_jac_int
    ! Jacobian inverse
    associate( j => jac_mtx_int )
    inv_jac_mtx_int = ( 1.0_rk / det_jac ) * reshape( &
    & [ j(2,2), -j(1,2),         &
    &  -j(2,1),  j(1,1) ],       &
    & shape = [2,2], order = [2,1] )
    end associate
    if ( present(inv_jac_mtx) ) inv_jac_mtx = inv_jac_mtx_int
    !
    if ( present(b_mtx) ) call gemm(inv_jac_mtx_int,dn_dxi_mtx,b_mtx)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine gradient_cp2d4
end module cp2d4
!*****************************************************************************80
program main
  use blas95, only: gemm
  use types, only: ik, rk, cl
  use point, only: dom, point_2d_t
  use cp2d4, only: cp2d4_t, ngp, nelnod
!*****************************************************************************80
  implicit none
!*****************************************************************************80
  integer(ik) :: i, j, p
  integer(ik) :: connectivity(4)
  real(rk) :: det_jac, vol, h
  real(rk) :: b_mtx(dom,nelnod), b_mtx_int(dom,nelnod)
  real(rk) :: l_max_analy
  real(rk) :: rtmp1(nelnod,nelnod), l_mtx_int(nelnod,nelnod)
  type(point_2d_t) :: coordinates(4)
  type(cp2d4_t) :: element
!*****************************************************************************80
  integer(ik) :: esta
  character(len=cl) :: emsg
!*****************************************************************************80
  integer(ik) :: fpm(128)
  character(len=1), parameter :: UPLO = 'F'
  real  :: Emin, Emax
  real  :: epsout
  integer :: loop
  integer :: M0, M, info
  real :: E(nelnod)
  real :: X(nelnod,nelnod)
  real :: res(nelnod)
!*****************************************************************************80
  h = 0.5
  coordinates = [ point_2d_t([0.,0.]), point_2d_t([h,0.]),  &
  &point_2d_t([h+1,h+1]), point_2d_t([0.,h]) ]
  connectivity = [ 1,2,3,4 ]
  element = cp2d4_t(coordinates,connectivity)
!*****************************************************************************80
  vol = 0.0_rk
  b_mtx_int = 0.0_rk
  l_mtx_int = 0.0_rk
  do p = 1, ngp
    call element%gradient(gp_num=p,det_jac=det_jac,b_mtx=b_mtx,esta=esta,&
    &emsg=emsg)
    if (esta/=0) print*, trim(emsg)
    vol = vol + det_jac
    b_mtx_int = b_mtx_int + b_mtx * det_jac
    call gemm(transpose(b_mtx),b_mtx,rtmp1)
    l_mtx_int = l_mtx_int + rtmp1 * det_jac
  end do
  print*, 'vol=',vol
!  print*, '--b mtx avg--'
!  b_mtx_int = 1.0_rk / vol * b_mtx_int
!  do i = 1, dom
!    print*, (b_mtx_int(i,j),j=1,nelnod)
!  end do
!  l_max_analy = 0.0_rk
!  do i = 1, dom
!    do j=1,nelnod
!      l_max_analy = l_max_analy + b_mtx_int(i,j) * b_mtx_int(i,j)
!    end do
!  end do
!  l_max_analy = sqrt(l_max_analy)
!  print*, 'Maximal eigenvalue analyitcal upper bound = ', l_max_analy
!  print*, 'Characheristic element dimension = ', vol / l_max_analy
  print*, '--l mtx avg--'
  l_mtx_int = 1.0_rk / vol * l_mtx_int
  do i = 1, nelnod
    print*, (l_mtx_int(i,j),j=1,nelnod)
  end do
  l_max_analy = 0.0_rk
  do i = 1, nelnod
    do j=1,nelnod
      l_max_analy = l_max_analy + l_mtx_int(i,j) * l_mtx_int(i,j)
    end do
  end do
  l_max_analy = sqrt(l_max_analy)
  print*, 'Maximal eigenvalue analyitcal upper bound = ', l_max_analy
  print*, 'Characheristic element dimension = ', 1.0 / l_max_analy
!*****************************************************************************80
! Initialize
  call feastinit(fpm)
  fpm(1)=1
  Emin=tiny(0.0)
  Emax=10.0
  M0 = 4
  M = 4
!*****************************************************************************80
  call sfeast_syev(UPLO,nelnod,l_mtx_int,nelnod,fpm,epsout,loop, &
  & Emin,Emax,M0,E,X,M,res,info)
  print  *,' FEAST OUTPUT INFO ',info
  if(info.ne.0) stop 1
  print*, 'Eigenvalues of l mtx:'
  do i = 1, nelnod
    print*, e(i)
  end do
  print*, 'Characheristic element dimension = ', 1.0 /  maxval(e)
!*****************************************************************************80
end program main




















