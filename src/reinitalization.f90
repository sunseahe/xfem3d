module reinitalzation
!*****************************************************************************80
  use blas95, only: dot, gemv, gemm
  use types, only: ik, int64, rk, lk, debug, stdout
  use general_routines, only: size_mtx, outer, pnorm
  use memory_storage, only: size_in_bytes, write_size_of_storage
  use point, only: dom
  use fe_c3d10, only: nelnod, ngp, c3d10_t, w
  use mesh_data, only: nnod, nfe, char_fe_dim, finite_elements
  use scalar_field, only: scalar_field_t
  use sparse, only: sparse_square_matrix_t, sparse_linear_system_t
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  logical(lk), protected :: configured = .false.
!*****************************************************************************80
  integer(ik) :: num_reinit = 200 ! Number of reinitalization equations
  real(rk) :: alpha = 0.5e0_rk ! Correction to time step
  real(rk) :: d_t = 0.0e0_rk ! Time step
  real(rk) :: er = 0.0e0_rk
  real(rk) :: c = 1.0e-1_rk ! Diffusion coeficient
  real(rk) :: rho = 1.0e1_rk ! Enforce Dirichlet boundary
  real(rk) :: sign_dist_tol = 1.0e-3_rk ! Tolerance for convergence
  logical(lk) :: write_par = .false. ! Write parameters to log file
!*****************************************************************************80
  type(scalar_field_t), save :: lsf_0
  type(scalar_field_t), pointer :: sdf => null()
  type(scalar_field_t), save :: r_vec
!*****************************************************************************80
  type(sparse_square_matrix_t), save :: c_mtx
  type(sparse_linear_system_t), save :: linear_system
!*****************************************************************************80
  integer(int64) :: mem_fac_c_mtx = 0
!*****************************************************************************80
  public :: set, configured, reinitalization_statistics
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine set(alpha_in,c_in,rho_in,num_reinit_in,sign_dist_tol_in)
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
    configured = .true.
    !
  end subroutine set
!*****************************************************************************80
! C matrix fe
!*****************************************************************************80
  pure subroutine calc_fe_c_mtx(c3d10,c_mtx,esta,emsg)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: c_mtx(:,:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: p
    real(rk) :: det_jac, lsf_0_gp, lsf_0_nv(nelnod)
    real(rk) :: rtmp1(nelnod,nelnod), rtmp2(nelnod,nelnod)
    real(rk) :: m_gam_mtx(nelnod,nelnod)
    real(rk) :: n(nelnod), b(dom,nelnod)
    ! Checks
    if ( debug ) then
      if ( .not.size_mtx(c_mtx,nelnod,nelnod) ) then
        esta = -1
        emsg ='Cal fe c mtx: El cmtx size incorrect'
        return
      end if
    end if
    !
    c_mtx = 0.0_rk
    do p = 1, ngp
      ! Set main values
      call c3d10%n_mtx(gp_num=p,n_mtx=n,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      call c3d10%gradient(gp_num=p,b_mtx=b,det_jac=det_jac,&
      &esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      !
      call outer(n,n,rtmp1)
      call gemm(transpose(b),b,rtmp2)
      ! M gamma
      call lsf_0%get_element_nodal_values(c3d10,lsf_0_nv,esta,emsg)
      lsf_0_gp = dot(n,lsf_0_nv)
      m_gam_mtx = rtmp1 * dirac_delta(lsf_0_gp)
      ! integrate
      c_mtx = c_mtx + ( rtmp1 + d_t * er * rtmp2 +  &
      & rho * m_gam_mtx )* w(p) * det_jac
    end do
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine calc_fe_c_mtx
!*****************************************************************************80
! Calculate c mtx
!*****************************************************************************80
  subroutine c_mtx_setup(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e, i, j, indx
    integer(ik) :: nnz
    integer(ik) :: connectivity(nelnod)
    integer(ik), allocatable :: crow(:), ccol(:)
    real(rk) :: fe_c_mtx(nelnod,nelnod)
    real(rk), allocatable :: cx(:)

    ! Allocate c mtx coordinate format
    nnz = nelnod * ( nelnod + 1 ) / 2 * nfe
    allocate(crow(nnz),ccol(nnz),cx(nnz),stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    indx = 1
    !$omp parallel do schedule(static,1) &
    !$omp private(e,fe_c_mtx) &
    !$omp shared(finite_elements,esta,emsg)
    do e = 1, nfe
      if ( esta == 0 ) then
        call calc_fe_c_mtx(finite_elements(e),fe_c_mtx,esta,emsg)
        !$omp critical
        ! Add to sparse matrix
        call finite_elements(e)%get_connectivity(connectivity,esta,emsg)
        do i = 1, nelnod
          do j = 1, nelnod
            ! symmetric matrix only upper part
            if ( connectivity(i) <= connectivity(j) ) then
              crow(indx) = connectivity(i)
              ccol(indx) = connectivity(j)
              cx(indx) = fe_c_mtx(i,j)
              indx = indx + 1
            end if
          end do
        end do
        !$omp end critical
      end if
    end do
    !$omp end parallel do
    if ( esta /= 0 ) return
    ! Allocate sparse matrix
    call c_mtx%set(nnod,crow,ccol,cx,esta,emsg)
    if ( esta /= 0 ) return
    ! Factorize
    call linear_system%solve(job=1,a=c_mtx,mem_used=mem_fac_c_mtx,&
    &esta=esta,emsg=emsg)
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine c_mtx_setup
!*****************************************************************************80
! R vector fe
!*****************************************************************************80
  pure subroutine calc_fe_r_vec(c3d10,r_vec,esta,emsg)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: r_vec(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: p
    real(rk) :: beta1
    real(rk) :: det_jac
    real(rk) :: n(nelnod), b(dom,nelnod), inv_jac_mtx(dom,dom)
    real(rk) :: lsf_0_nv(nelnod), sdf_nv(nelnod)
    real(rk) :: lsf_0_gp, sdf_gp
    real(rk) :: grad_sdf(dom), norm_sdf(dom)
    real(rk) :: s_0, v(dom), v_tilde(nelnod)
    real(rk) :: rtmp1(dom), rtmp2(nelnod)
    real(rk) :: r1(nelnod), r2(nelnod), r3(nelnod)
    !
    if ( debug ) then
      if ( .not.size(r_vec)==nelnod ) then
        esta = -1
        emsg ='Cal fe r vec: r vec size incorrect'
        return
      end if
    end if
    !
    r1 = 0.0e0_rk; r2 = 0.0e0_rk; r3 = 0.0e0_rk
    do p = 1, ngp
      ! Set main values
      call c3d10%n_mtx(gp_num=p,n_mtx=n,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      call c3d10%gradient(gp_num=p,b_mtx=b,det_jac=det_jac,&
      &inv_jac_mtx=inv_jac_mtx,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      !
      call lsf_0%get_element_nodal_values(c3d10,lsf_0_nv,esta,emsg)
      call sdf%get_element_nodal_values(c3d10,sdf_nv,esta,emsg)
      lsf_0_gp = dot(n,lsf_0_nv)
      sdf_gp = dot(n,sdf_nv)
      call gemv(b,sdf_nv,grad_sdf)
      norm_sdf = + grad_sdf / reg_pnorm(2,grad_sdf)
      s_0 = smooth_sign(lsf_0_gp)
      v = s_0 * norm_sdf
      !
      call gemv(inv_jac_mtx,v,rtmp1)
      beta1 = 1.0e0_rk / ( 2.0e0_rk * sqrt( d_t**(-2) + &
      & reg_pnorm(2,rtmp1) ) )
      call gemv(transpose(b),v,rtmp2)
      v_tilde = n + beta1 * rtmp2
      !
      r1 = r1 + s_0 * v_tilde * w(p) * det_jac
      r2 = r2 + dot(v,grad_sdf) * v_tilde * w(p) * det_jac
      r3 = r3 + sdf_gp * n * w(p) * det_jac
      !
    end do
    r_vec = d_t * (r1 - r2) + r3
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine calc_fe_r_vec
!*****************************************************************************80
! Setup r vector
!*****************************************************************************80
  subroutine r_vec_setup(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e
    real(rk) :: fe_r_vec(nelnod)
    !
    call r_vec%set(esta=esta,emsg=emsg)
    if ( esta /= 0 ) return
    !
    !$omp parallel do schedule(static,1) &
    !$omp private(e,fe_r_vec) &
    !$omp shared(finite_elements,esta,emsg)
    do e = 1, nfe
      if ( esta == 0 ) then
        call calc_fe_r_vec(finite_elements(e),fe_r_vec,esta,emsg)
        !$omp critical
        call r_vec%assemble_element_nodal_values(finite_elements(e),&
        &fe_r_vec,esta,emsg)
        !$omp end critical
      end if
    end do
    !$omp end parallel do
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine r_vec_setup
!*****************************************************************************80
! Dirac delta function
!*****************************************************************************80
  pure function dirac_delta(x) result(res)
    real(rk), intent(in) :: x
    real(rk) :: res
    associate( delta => 2.0_rk * char_fe_dim )
    if( abs(x) <= delta ) then
      res = 3.0_rk / 4.0_rk * delta * &
      & ( 1.0_rk - x**2 / delta**2 )
    else
      res = 0.0_rk
    end if
    end associate
  end function dirac_delta
!*****************************************************************************80
  pure function smooth_sign(lsf) result(res)
    real(rk), intent(in) :: lsf
    real(rk) :: res
    res = lsf / sqrt( lsf**2 + char_fe_dim**2 )
  end function smooth_sign
!*****************************************************************************80
  pure function reg_pnorm(p,x) result(res)
    integer(ik), intent(in) :: p
    real(rk), intent(in) :: x(:)
    real(rk) :: res
    res = sqrt( pnorm(p,x)**2 + char_fe_dim**2 )
  end function reg_pnorm
!*****************************************************************************80
! Reinitalization statistics
!*****************************************************************************80
  subroutine reinitalization_statistics()
    !
    integer(int64) :: storage
    !
    write(stdout,'(a)') '*** Reinitalization statistics ***'
    storage = mem_fac_c_mtx
    write(stdout,'(a,a)') 'Data allocated: ', trim(write_size_of_storage( &
    &storage))
    !
  end subroutine reinitalization_statistics
!*****************************************************************************80
end module reinitalzation
