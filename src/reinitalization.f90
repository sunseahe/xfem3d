module reinitalzation
!*****************************************************************************80
  use blas95, only: dot, gemv, gemm
  use types
  use general_routines, only: size_mtx
  use point, only: dom
  use fe_c3d10, only: nelnod, ngp, c3d10_t, w
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  type :: reinit_t
    private
    logical(lk) :: configured = .false.
    !type(spar_mtx) :: c_mtx
    !type(spar_lin_sys) :: lss
  contains
    procedure :: set
    procedure :: is_configured
    !procedure, private :: cmtx_setup
    !procedure :: calculate
  end type reinit_t
  type(reinit_t), save :: reinit
!*****************************************************************************80
  integer(ik) :: num_reinit = 200 ! Number of reinitalization equations
  logical(lk) :: par_status = .false. ! Parameters set
  real(rk) :: alpha = 0.5e0_rk ! Correction to time step
  real(rk) :: d_t = 0.0e0_rk ! Time step
  real(rk) :: er = 0.0e0_rk
  real(rk) :: c = 1.0e-1_rk ! Diffusion coeficient
  real(rk) :: rho = 1.0e4_rk ! Enforce Dirichlet boundary
  real(rk) :: sign_dist_tol = 1.0e-3_rk ! Tolerance for convergence
  logical(lk) :: write_par = .false. ! Write parameters to log file
!*****************************************************************************80
  public :: reinit
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine set(self,alpha_in,c_in,rho_in,num_reinit_in,sign_dist_tol_in)
    class(reinit_t), intent(out) :: self
    integer(ik), optional, intent(in) :: num_reinit_in
    real(rk), optional, intent(in) :: alpha_in
    real(rk), optional, intent(in) :: c_in
    real(rk), optional, intent(in) :: rho_in
    real(rk), optional, intent(in) :: sign_dist_tol_in
    ! set parameters
    if (present(num_reinit_in)) num_reinit = num_reinit_in
    if (present(alpha_in)) alpha = alpha_in
    if (present(c_in)) c = c_in
    if (present(rho_in)) rho = rho_in
    if (present(sign_dist_tol_in)) sign_dist_tol = sign_dist_tol_in
    self%configured = .true.
    !
  end subroutine set
  pure logical(lk) function is_configured(self)
    class(reinit_t), intent(in) :: self
    is_configured = self%configured
  end function is_configured
!*****************************************************************************80
  pure subroutine cal_el_cmtx(c3d10,g_mtx,el_cmtx,esta,emsg)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(in) :: g_mtx(:,:)
    real(rk), intent(out) :: el_cmtx(:,:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: p
    real(rk) :: rtmp1(nelnod,nelnod), rtmp2(nelnod,nelnod)
    real(rk) :: det_jac
    real(rk) :: n(nelnod), b(dom,nelnod)
    ! Checks
    if ( debug ) then
      if ( .not.size_mtx(g_mtx,nelnod,nelnod) ) then
        esta = -1
        emsg ='Cal el cmtx: G mtx size incorrect'
        return
      end if
      !
      if ( .not.size_mtx(el_cmtx,nelnod,nelnod) ) then
        esta = -1
        emsg ='Cal el cmtx: El cmtx size incorrect'
        return
      end if
    end if
    !
    el_cmtx = 0.0_rk
    do p = 1, ngp
      ! Set main values
      call c3d10%n_mtx(gp_num=p,n_mtx=n,esta=esta,emsg=emsg)
      call c3d10%gradient(gp_num=p,b_mtx=b,det_jac=det_jac,&
      &esta=esta,emsg=emsg)
      !call domain_fe%b_mtx_sca(el_coo=el_coo,gp_num=p,b_mtx_sca=b)
      !w = domain_fe%get_w(gp_num=p)
      !det_j = domain_fe%det_j(el_coo=el_coo,gp_num=p)
      !
      !call outer(n,n,rtmp1)
      !call gemm(transpose(b),b,rtmp2)
      ! integrate
      el_cmtx = el_cmtx + ( rtmp1 + d_t * er * rtmp2 +  &
      & rho * g_mtx )* w(p) * det_jac
    end do
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine cal_el_cmtx
!*****************************************************************************80
end module reinitalzation
