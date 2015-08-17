module hamilton_jacobi
!*****************************************************************************80
  use blas95, only: dot, gemv, gemm
  use types, only: ik, rk, lk, cl, debug, stdout, log_file
  use reinitalzation, only: reg_pnorm
  use sparse, only: sparse_square_matrix_t, sparse_linear_system_t
  use fe_c3d10, only: nelnod, ngp, c3d10_t, w
  use point, only: dom
  use general_routines, only: size_mtx, outer, time
  use mesh_data, only: nnod, nfe, finite_elements, char_fe_dim
  use scalar_field, only: scalar_field_t
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  logical(lk), protected :: configured = .false.
!*****************************************************************************80
  integer(ik) :: max_num_advect = 50 ! Maximal number of advection steps
  real(rk) :: max_alpha_t = 1.0_rk ! Maximal alpha correction
  real(rk) :: d_t = 0.0e0_rk ! Time step - to be calculated
  logical(lk) :: write_iter_t = .false. ! Advection iteration time
  ! written
!*****************************************************************************80
  logical(lk) :: iter_sol = .false. ! Iterative solver
  integer(ik) :: iter_niter = 50 ! Number of iterations for iterative solver
  real(rk) :: iter_tol = 1.0e-16_rk ! Tolerance
!*****************************************************************************80
  type(sparse_square_matrix_t), save :: m_mtx
  type(sparse_linear_system_t), save :: linear_system
!*****************************************************************************80
  type(scalar_field_t), pointer :: lsf
  type(scalar_field_t), save :: v
  type(scalar_field_t), save :: r_vec
!*****************************************************************************80
  public :: configured, calculate_advection
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine set_hamilton_jacobi(max_num_advect_in,max_alpha_t_in,&
    &iter_sol_in)
    integer(ik), optional, intent(in) :: max_num_advect_in
    real(rk), optional, intent(in) :: max_alpha_t_in
    logical(lk), optional, intent(in) :: iter_sol_in
    !
    if (present(max_num_advect_in)) max_num_advect = max_num_advect_in
    if (present(max_alpha_t_in)) max_alpha_t = max_alpha_t_in
    if (present(iter_sol_in)) iter_sol = .true.
    configured = .true.
    !
  end subroutine set_hamilton_jacobi
!*****************************************************************************80
! M matrix fe
!*****************************************************************************80
  pure subroutine calc_fe_m_mtx(c3d10,m_mtx)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: m_mtx(:,:)
    !
    integer(ik) :: p
    real(rk) :: det_jac
    real(rk) :: rtmp1(nelnod,nelnod)
    real(rk) :: n(nelnod)
    !
    m_mtx = 0.0_rk
    do p = 1, ngp
      ! Set main values
      call c3d10%main_values(gp_num=p,n_mtx=n,det_jac=det_jac)
      call outer(n,n,rtmp1)
      ! integrate
      m_mtx = m_mtx + rtmp1 * w(p) * det_jac
    end do
    !
  end subroutine calc_fe_m_mtx
!*****************************************************************************80
! Calculate m mtx
!*****************************************************************************80
  subroutine m_mtx_setup(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e, i, j, indx
    integer(ik) :: nnz
    integer(ik), allocatable :: crow(:), ccol(:)
    real(rk) :: fe_m_mtx(nelnod,nelnod)
    real(rk), allocatable :: cx(:)
    ! Allocate c mtx coordinate format
    nnz = nelnod * ( nelnod + 1 ) / 2 * nfe
    allocate(crow(nnz),ccol(nnz),cx(nnz),stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    indx = 1
    !$omp parallel do schedule(static,1) &
    !$omp private(e,fe_m_mtx) &
    !$omp shared(finite_elements)
    do e = 1, nfe
      call calc_fe_m_mtx(finite_elements(e),fe_m_mtx)
      !$omp critical
      ! Add to sparse matrix
      do i = 1, nelnod
        do j = 1, nelnod
          ! symmetric matrix only upper part
          if ( finite_elements(e)%connectivity(i) <= &
          & finite_elements(e)%connectivity(j) ) then
            crow(indx) = finite_elements(e)%connectivity(i)
            ccol(indx) = finite_elements(e)%connectivity(j)
            cx(indx) = fe_m_mtx(i,j)
            indx = indx + 1
          end if
        end do
      end do
      !$omp end critical
    end do
    !$omp end parallel do
    ! Allocate sparse matrix
    call m_mtx%set(nnod,crow,ccol,cx,esta,emsg)
    if ( esta /= 0 ) return
    ! Factorize
    if ( .not. iter_sol ) then
      call linear_system%solve_dir(job=1,a=m_mtx,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
    end if
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine m_mtx_setup
!*****************************************************************************80
! R vector fe
!*****************************************************************************80
  pure subroutine calc_fe_r_vec(c3d10,r_vec)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: r_vec(:)
    !
    integer(ik) :: p
    real(rk) :: beta1, det_jac
    real(rk) :: n(nelnod), b(dom,nelnod), inv_jac_mtx(dom,dom)
    real(rk) :: v_n(dom)
    real(rk) :: lsf_nv(nelnod), v_nv(nelnod)
    real(rk) :: grad_lsf(dom), norm_lsf(dom)
    real(rk) :: rtmp1(dom), rtmp2(nelnod)
    !
    r_vec = 0.0_rk
    do p = 1, ngp
      ! Set main values
      call c3d10%main_values(gp_num=p,n_mtx=n,b_mtx=b,det_jac=det_jac,&
      &inv_jac_mtx=inv_jac_mtx)
      call lsf%get_element_nodal_values(c3d10,lsf_nv)
      call v%get_element_nodal_values(c3d10,v_nv)
      !
      call gemv(b,lsf_nv,grad_lsf)
      norm_lsf = - grad_lsf / reg_pnorm(2,grad_lsf)
      v_n = dot(n,v_nv) * norm_lsf
      call gemv(inv_jac_mtx,v_n,rtmp1)
      call gemv(transpose(b),v_n,rtmp2)
      beta1 = 1.0e0_rk / ( 2.0e0_rk * sqrt( d_t**(-2) + &
      & reg_pnorm(2,rtmp1) ) )
      rtmp2 = ( n + beta1 * rtmp2 ) * dot(v_n,grad_lsf)
      ! Integrate
      r_vec = r_vec - d_t * rtmp2 * w(p) * det_jac
    end do
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
    ! Reset the right hand side vector
    call r_vec%set(esta=esta,emsg=emsg)
    if ( esta /= 0 ) return
    !
    !$omp parallel do schedule(static,1) &
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
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine r_vec_setup
!*****************************************************************************80
! Calculate advection
!*****************************************************************************80
  subroutine calculate_advection(lsf_inout,velocity,num_advect,alpha_t,&
    &esta,emsg)
    type(scalar_field_t), target, intent(inout) :: lsf_inout
    integer(ik), intent(in) :: num_advect
    real(rk), intent(in) :: alpha_t
    interface
      subroutine velocity(lsf,v,esta,emsg)
        import :: ik, scalar_field_t
        type(scalar_field_t), intent(in) :: lsf
        type(scalar_field_t), intent(out) :: v
        integer(ik), intent(out) :: esta
      character(len=*), intent(out) :: emsg
      end subroutine velocity
    end interface
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    character(len=cl) :: info
    type(time) :: t_complete, t_iter
    ! Check
    if ( .not. configured ) then
      esta = -1
      emsg = 'Advection procedure not configured'
      return
    end if
    !
    lsf => lsf_inout
    call r_vec%set(esta=esta,emsg=emsg); if ( esta /= 0 ) return
    ! Advect
    write(stdout,'(a)') 'Solving the Hamilton Jacobi equation ...'
    call t_complete%start_timer()
    ! Calculate stable time increment
    d_t = alpha_t * char_fe_dim / v%max_value()
    !
    do i = 1, min(max_num_advect,num_advect)
      if ( i == 1 .and. .not. write_iter_t ) call t_iter%start_timer()
      ! Get advection velocity
      call velocity(lsf,v,esta,emsg); if ( esta /= 0 ) return
      call r_vec_setup(esta,emsg); if ( esta /= 0 ) return
      ! Backsubstitution or iterative solve
      if ( .not. iter_sol ) then
        call linear_system%solve_dir(job=2,a=m_mtx,b=r_vec%values,&
        &x=lsf%values,esta=esta,emsg=emsg)
        if ( esta /= 0 ) return
      else
        call linear_system%solve_iter(a=m_mtx,b=r_vec%values,&
        &x=lsf%values,niter=iter_niter,tol=iter_tol,info=info,&
        & esta=esta,emsg=emsg); if ( esta /= 0 ) return
        !if ( i == 1 .and. .not. write_iter_t ) then
          write(log_file,'(a)') trim(info)
        !end if
      end if
      !
      if ( i == 1 .and. .not. write_iter_t ) then
        call t_complete%write_elapsed_time(log_file,'One advection &
        &iteration time')
        write_iter_t = .true.
      end if
    end do
    !
    call t_complete%write_elapsed_time(stdout,'Advection finished, &
    &elapsed time')
    ! Clean
    if ( .not. iter_sol ) then
      call linear_system%solve_dir(job=3,a=m_mtx,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
    end if
    lsf => null()
    call m_mtx%delete(esta,emsg); if ( esta /= 0 ) return
    call v%delete(esta,emsg); if ( esta /= 0 ) return
    call r_vec%delete(esta,emsg); if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine calculate_advection
!*****************************************************************************80
end module hamilton_jacobi
