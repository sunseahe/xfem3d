module reinitalzation
!*****************************************************************************80
  use blas95, only: dot, gemv, gemm
  use types, only: ik, rk, lk, es, debug, stdout, nl, log_file, cl
  use general_routines, only: size_mtx, outer, pnorm, time, &
  &write_dense_mtx_real
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
  integer(ik) :: num_conv_iter = 3 ! Number of convergence iterations
  real(rk) :: alpha = 0.5e0_rk ! Correction to time step
  real(rk) :: d_t = 0.0e0_rk ! Time step - to be calculated
  real(rk) :: er = 0.0e0_rk ! Diffusion related - to be calculated
  logical(lk) :: status_par = .false.
  real(rk) :: c = 1.0e-1_rk ! Diffusion coeficient
  real(rk) :: rho = 1.0e0_rk ! Enforce Dirichlet boundary
  real(rk) :: sign_dist_tol = 1.0e-3_rk ! Tolerance for convergence
  logical(lk) :: write_par = .false. ! Write parameters to log file
  logical(lk) :: write_iter_t = .false. ! Reinitalization iteration time
  ! written
!*****************************************************************************80
  logical(lk) :: iter_sol = .false. ! Iterative solver
  integer(ik) :: iter_niter = 50 ! Number of iterations for iterative solver
  real(rk) :: iter_tol = 1.0e-8_rk ! Tolerance
!*****************************************************************************80
  type(scalar_field_t), save :: sdf_0
  type(scalar_field_t), pointer :: sdf => null()
  type(scalar_field_t), save :: r_vec
!*****************************************************************************80
  type(sparse_square_matrix_t), save :: c_mtx
  type(sparse_linear_system_t), save :: linear_system
!*****************************************************************************80
  public :: set_reinitalization, configured, calculate_reinitalization, &
  & reg_pnorm
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine set_reinitalization(alpha_in,c_in,rho_in,num_reinit_in,&
    &num_conv_iter_in,sign_dist_tol_in,iter_sol_in)
    integer(ik), optional, intent(in) :: num_reinit_in
    integer(ik), optional, intent(in) :: num_conv_iter_in
    real(rk), optional, intent(in) :: alpha_in
    real(rk), optional, intent(in) :: c_in
    real(rk), optional, intent(in) :: rho_in
    real(rk), optional, intent(in) :: sign_dist_tol_in
    logical(lk), optional, intent(in) :: iter_sol_in
    ! set parameters
    if (present(num_reinit_in)) num_reinit = num_reinit_in
    if (present(num_conv_iter_in)) num_conv_iter = num_conv_iter_in
    if (present(alpha_in)) alpha = alpha_in
    if (present(c_in)) c = c_in
    if (present(rho_in)) rho = rho_in
    if (present(sign_dist_tol_in)) sign_dist_tol = sign_dist_tol_in
    if (present(iter_sol_in)) iter_sol = .true.
    configured = .true.
    !
  end subroutine set_reinitalization
!*****************************************************************************80
! C matrix fe
!*****************************************************************************80
  pure subroutine calc_fe_c_mtx(c3d10,c_mtx)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: c_mtx(:,:)
    !
    integer(ik) :: p
    real(rk) :: det_jac, sdf_0_gp, sdf_0_nv(nelnod)
    real(rk) :: rtmp1(nelnod,nelnod), rtmp2(nelnod,nelnod)
    real(rk) :: m_gam_mtx(nelnod,nelnod)
    real(rk) :: n(nelnod), b(dom,nelnod)
    !
    c_mtx = 0.0_rk
    do p = 1, ngp
      ! Get main values
      call c3d10%main_values(gp_num=p,n_mtx=n,b_mtx=b,det_jac=det_jac)
      !
      call outer(n,n,rtmp1)
      call gemm(transpose(b),b,rtmp2)
      ! M gamma
      call sdf_0%get_element_nodal_values(c3d10,sdf_0_nv)
      sdf_0_gp = dot(n,sdf_0_nv)
      m_gam_mtx = rtmp1 * dirac_delta(sdf_0_gp)
      ! integrate
      c_mtx = c_mtx + ( rtmp1 + d_t * er * rtmp2 +  &
      & rho * m_gam_mtx )* w(p) * det_jac
    end do
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
    !$omp shared(finite_elements)
    do e = 1, nfe
      call calc_fe_c_mtx(finite_elements(e),fe_c_mtx)
      !$omp critical
      ! Add to sparse matrix
      do i = 1, nelnod
        do j = 1, nelnod
          ! symmetric matrix only upper part
          if ( finite_elements(e)%connectivity(i) <= &
          & finite_elements(e)%connectivity(j) ) then
            crow(indx) = finite_elements(e)%connectivity(i)
            ccol(indx) = finite_elements(e)%connectivity(j)
            cx(indx) = fe_c_mtx(i,j)
            indx = indx + 1
          end if
        end do
      end do
      !$omp end critical
    end do
    !$omp end parallel do
    ! Allocate sparse matrix
    call c_mtx%set(nnod,crow,ccol,cx,esta,emsg)
    if ( esta /= 0 ) return
    ! Factorize if direct solution
    if ( .not. iter_sol ) then
      call linear_system%solve_dir(job=1,a=c_mtx,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
    end if
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine c_mtx_setup
!*****************************************************************************80
! R vector fe
!*****************************************************************************80
  pure subroutine calc_fe_r_vec(c3d10,r_vec)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: r_vec(:)
    !
    integer(ik) :: p
    real(rk) :: beta1
    real(rk) :: det_jac
    real(rk) :: n(nelnod), b(dom,nelnod), inv_jac_mtx(dom,dom)
    real(rk) :: sdf_0_nv(nelnod), sdf_nv(nelnod)
    real(rk) :: sdf_0_gp, sdf_gp
    real(rk) :: grad_sdf(dom), norm_sdf(dom)
    real(rk) :: s_0, v(dom), v_tilde(nelnod)
    real(rk) :: rtmp1(dom), rtmp2(nelnod), rtmp3
    real(rk) :: r1(nelnod), r2(nelnod), r3(nelnod)
    !
    r1 = 0.0e0_rk; r2 = 0.0e0_rk; r3 = 0.0e0_rk
    do p = 1, ngp
      ! Get main values
      call c3d10%main_values(gp_num=p,n_mtx=n,b_mtx=b,det_jac=det_jac,&
      &inv_jac_mtx=inv_jac_mtx)
      call sdf_0%get_element_nodal_values(c3d10,sdf_0_nv)
      call sdf%get_element_nodal_values(c3d10,sdf_nv)
      !
      sdf_0_gp = dot(n,sdf_0_nv)
      sdf_gp = dot(n,sdf_nv)
      call gemv(b,sdf_nv,grad_sdf)
      norm_sdf = + grad_sdf / reg_pnorm(2,grad_sdf)
      s_0 = smooth_sign(sdf_0_gp)
      v = s_0 * norm_sdf
      !
      call gemv(inv_jac_mtx,v,rtmp1)
      call gemv(transpose(b),v,rtmp2)
      rtmp3 = reg_pnorm(2,rtmp1) !sqrt( pnorm(2,rtmp1)**2 + char_fe_dim**2 ) ! Bug
      beta1 = 1.0e0_rk / ( 2.0e0_rk * sqrt( d_t**(-2) + rtmp3 ) )
      v_tilde = n + beta1 * rtmp2
      !
      r1 = r1 + s_0 * v_tilde * w(p) * det_jac
      r2 = r2 + dot(v,grad_sdf) * v_tilde * w(p) * det_jac
      r3 = r3 + sdf_gp * n * w(p) * det_jac
      !
    end do
    r_vec = d_t * (r1 - r2) + r3
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
    type(time) :: t
    ! Reset the right hand side vector
    call t%start_timer()
    call r_vec%set(esta=esta,emsg=emsg)
    if ( esta /= 0 ) return
    !
    !print*, 'here'
    !$omp parallel do & !schedule(static,1)
    !$omp private(e,fe_r_vec) &
    !$omp shared(finite_elements)
    do e = 1, nfe
      call calc_fe_r_vec(finite_elements(e),fe_r_vec)
      !$omp critical
      call r_vec%assemble_element_nodal_values(finite_elements(e),&
      &fe_r_vec)
      !$omp end critical
    end do
    !$omp end parallel do
    !call t%write_elapsed_time(stdout)
    !stop 'test'
    !
  end subroutine r_vec_setup
!*****************************************************************************80
! Calculate signed distance tolerance
!*****************************************************************************80
  pure subroutine calc_fe_sdf_tol(c3d10,sdf_tol)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: sdf_tol
    !
    integer(ik) :: p
    real(rk) :: det_jac
    real(rk) :: sdf_nv(nelnod), grad_sdf(dom)
    real(rk) :: b(dom,nelnod)
    !
    sdf_tol = 0.0e0_rk
    do p = 1, ngp
      call c3d10%main_values(gp_num=p,b_mtx=b,det_jac=det_jac)
      call sdf%get_element_nodal_values(c3d10,sdf_nv)
      call gemv(b,sdf_nv,grad_sdf)
      sdf_tol = sdf_tol + (reg_pnorm(2,grad_sdf) - 1.0e0_rk)**2 &
      & * w(p) * det_jac
    end do
    !
  end subroutine calc_fe_sdf_tol
!*****************************************************************************80
  subroutine calc_sdf_tol(sdf_tol)
    real(rk), intent(out) :: sdf_tol
    !
    integer(ik) :: e
    real(rk) :: fe_sdf_tol
    !
    sdf_tol = 0.0e0_rk
    !$omp parallel do schedule(static,1)     &
    !$omp private(e,fe_sdf_tol) &
    !$omp shared(finite_elements) &
    !$omp reduction(+:sdf_tol)
    do e = 1, nfe
      call calc_fe_sdf_tol(finite_elements(e),fe_sdf_tol)
      sdf_tol = sdf_tol + fe_sdf_tol
    end do
    !$omp end parallel do
    sdf_tol = sqrt( sdf_tol )
    !
  end subroutine calc_sdf_tol
!*****************************************************************************80
! Calculate reinitalization
!*****************************************************************************80
  subroutine calculate_reinitalization(sdf_inout,esta,emsg)
    type(scalar_field_t), target, intent(inout) :: sdf_inout
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    integer(ik) :: conv_iter
    real(rk) :: sd_tol_previous, sd_tol_current, rel_tol
    logical(lk) :: converged
    character(len=cl) :: info
    type(time) :: t_complete, t_iter
    ! Check
    if ( .not. configured ) then
      esta = -1
      emsg = 'Reinitalization procedure not configured'
      return
    end if
    ! Copy scalar fields
    sdf => sdf_inout
    call sdf_0%copy(sdf_inout,esta,emsg); if ( esta /= 0 ) return
    ! Setup parameters
    if ( .not. status_par ) then
      d_t = alpha * char_fe_dim
      er = c * char_fe_dim**2 / d_t
      status_par = .true.
    end if
    ! Parameters
    if ( .not.write_par ) then
      write(log_file,'(a)')         'Solution parameters '
      write(log_file,'(a,i0)')      ' - num of steps is: ', num_reinit
      write(log_file,'(a,'//es//')') ' - stable time inc is: ', d_t
      write(log_file,'(a,'//es//')') ' - step corr par is: ', alpha
      write(log_file,'(a,'//es//')') ' - diffusion par is: ', c
      write(log_file,'(a,'//es//')') ' - enforce dirichlet: ', rho
      write(log_file,'(a,'//es//')') ' - solution tolerance: ', sign_dist_tol
      write_par = .true.
    end if
    ! Solve
    if ( .not. iter_sol ) then
      write(log_file,'(a)') 'Solving the reinitalization &
      &equation with direct solver ...'
    else
      write(log_file,'(a)') 'Solving the reinitalization &
      &equation with iterative solver ...'
    end if
    call t_complete%start_timer()
    ! Calculate c matrix
    call c_mtx_setup(esta,emsg); if ( esta /= 0 ) return
    ! Set up right hand side vector
    call r_vec%set(esta=esta,emsg=emsg); if ( esta /= 0 ) return
    !
    sd_tol_previous = huge(1.0_rk)
    converged = .false.
    conv_iter = 0
    do i = 1, num_reinit
      if ( i == 1 .and. .not. write_iter_t ) call t_iter%start_timer()
      call r_vec_setup(esta,emsg); if ( esta /= 0 ) return
      ! Backsubstitution or iterative solve
      if ( .not. iter_sol ) then
        call linear_system%solve_dir(job=2,a=c_mtx,b=r_vec%values,&
        &x=sdf%values,esta=esta,emsg=emsg)
        if ( esta /= 0 ) return
      else
        call linear_system%solve_iter(a=c_mtx,b=r_vec%values,&
        &x=sdf%values,niter=iter_niter,tol=iter_tol,info=info,&
        & esta=esta,emsg=emsg); if ( esta /= 0 ) return
        !if ( i == 1 .and. .not. write_iter_t ) then
          write(log_file,'(a)') trim(info)
        !end if
      end if
      ! Check for convergence
      call calc_sdf_tol(sd_tol_current)
      rel_tol = abs(sd_tol_current-sd_tol_previous) / sd_tol_current
      if( rel_tol <= sign_dist_tol ) then
        conv_iter = conv_iter + 1
        if( conv_iter >= num_conv_iter ) then
          converged = .true.
          exit
        end if
      else
        conv_iter = 0
      end if
      sd_tol_previous = sd_tol_current
      if ( i == 1 .and. .not. write_iter_t ) then
        call t_complete%write_elapsed_time(log_file,'One reinitalzation &
        &iteration time')
        write_iter_t = .true.
      end if
    end do
    ! Write status
    if ( .not. converged ) then
      write(log_file,'(a,'//es//',a,i0,a,a,a,'//es//',a)') &
      & '**Warning: solution relative tolerance (', rel_tol, &
      & ' ) after ', num_reinit , ' iterations, ', nl, 'is larger than &
      & specified (', sign_dist_tol ,' ).'
      write(log_file,'(a,'//es//',a)') 'Acheved a value of ', sd_tol_current, &
      ' for signed distance norm approximation.'
    else
      write(log_file,'(a,'//es//',a,a,a,'//es//',a,i0,a)') &
      & 'Acheved a value of ', sd_tol_current,' for signed distance norm &
      &approximation and', nl, 'relative tolerance of ', rel_tol, ' in ', i, &
      & ' iterations.'
    end if
    call t_complete%write_elapsed_time(log_file,'Reinitalization finished, &
    &elapsed time')
    ! Clean
    if ( .not. iter_sol ) then
      call linear_system%solve_dir(job=3,a=c_mtx,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
    end if
    call c_mtx%delete(esta,emsg); if ( esta /= 0 ) return
    call r_vec%delete(esta,emsg); if ( esta /= 0 ) return
    sdf => null()
    call sdf_0%delete(esta,emsg); if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine calculate_reinitalization
!*****************************************************************************80
! Dirac delta function
!*****************************************************************************80
  pure function dirac_delta(x) result(res)
    real(rk), intent(in) :: x
    real(rk) :: res
    associate( delta => char_fe_dim )
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
end module reinitalzation
