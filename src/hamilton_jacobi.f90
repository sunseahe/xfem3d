module hamilton_jacobi
!*****************************************************************************80
  use blas95, only: dot, gemv, gemm
  use types, only: ik, rk, lk, int64, debug, stdout
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
  logical(lk) :: iter_sol = .false. ! Iterative solver
!*****************************************************************************80
  type(sparse_square_matrix_t), save :: m_mtx
  type(sparse_linear_system_t), save :: linear_system
!*****************************************************************************80
  type(scalar_field_t), pointer :: lsf
  type(scalar_field_t), pointer :: v
  type(scalar_field_t), save :: r_vec
!*****************************************************************************80
  integer(int64) :: mem_fac_c_mtx = 0
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
  pure subroutine calc_fe_m_mtx(c3d10,m_mtx,esta,emsg)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: m_mtx(:,:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: p, p_tmp
    real(rk) :: det_jac
    real(rk) :: rtmp1(nelnod,nelnod)
    real(rk) :: n(nelnod)
    ! Checks
    if ( debug ) then
      if ( .not.size_mtx(m_mtx,nelnod,nelnod) ) then
        esta = -1
        emsg ='Cal fe m mtx: El mmtx size incorrect'
        return
      end if
    end if
    !
    m_mtx = 0.0_rk
    do p = 1, ngp
      p_tmp = p ! Gfortran bug
      ! Set main values
      call c3d10%n_mtx(gp_num=p_tmp,n_mtx=n,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      call c3d10%gradient(gp_num=p_tmp,det_jac=det_jac,&
      &esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      !
      call outer(n,n,rtmp1)
      ! integrate
      m_mtx = m_mtx + rtmp1 * w(p) * det_jac
    end do
    ! Sucess
    esta = 0
    emsg = ''
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
    !$omp shared(finite_elements,esta,emsg)
    do e = 1, nfe
      if ( esta == 0 ) then
        call calc_fe_m_mtx(finite_elements(e),fe_m_mtx,esta,emsg)
        !$omp critical
        ! Add to sparse matrix
        !associate( c => finite_elements(e)%connectivity )
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
        !end associate
        !$omp end critical
      end if
    end do
    !$omp end parallel do
    if ( esta /= 0 ) return
    ! Allocate sparse matrix
    call m_mtx%set(nnod,crow,ccol,cx,esta,emsg)
    if ( esta /= 0 ) return
    ! Factorize
    if ( .not. iter_sol ) then
      call linear_system%solve_dir(job=1,a=m_mtx,mem_used=mem_fac_c_mtx,&
      &esta=esta,emsg=emsg)
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
  pure subroutine calc_fe_r_vec(c3d10,r_vec,esta,emsg)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: r_vec(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: p, p_tmp
    real(rk) :: beta1, det_jac
    real(rk) :: n(nelnod), b(dom,nelnod), inv_jac_mtx(dom,dom)
    real(rk) :: v_n(dom)
    real(rk) :: lsf_nv(nelnod), v_nv(nelnod)
    real(rk) :: grad_lsf(dom), norm_lsf(dom)
    real(rk) :: rtmp1(dom), rtmp2(nelnod)
    !
    if ( debug ) then
      if ( .not.size(r_vec)==nelnod ) then
        esta = -1
        emsg ='Cal fe r vec: r vec size incorrect'
        return
      end if
    end if
    !
    r_vec = 0.0_rk
    do p = 1, ngp
      p_tmp = p ! Gfortran bug
      ! Set main values
      call c3d10%n_mtx(gp_num=p_tmp,n_mtx=n,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      call c3d10%gradient(gp_num=p_tmp,b_mtx=b,det_jac=det_jac,&
      &inv_jac_mtx=inv_jac_mtx,esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      call lsf%get_element_nodal_values(c3d10,lsf_nv,esta,emsg)
      if ( esta /= 0 ) return
      call v%get_element_nodal_values(c3d10,v_nv,esta,emsg)
      if ( esta /= 0 ) return
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
    ! Reset the right hand side vector
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
! Calculate advection
!*****************************************************************************80
  subroutine caluclate_advection(lsf_inout,v_in,num_advect,alpha_t,&
    &esta,emsg)
    type(scalar_field_t), target, intent(inout) :: lsf_inout
    type(scalar_field_t), target, intent(in) :: v_in
    integer(ik), intent(in) :: num_advect
    real(rk), intent(in) :: alpha_t
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    type(time) :: t_complete
    ! Check
    if ( .not. configured ) then
      esta = -1
      emsg = 'Advection procedure not configured'
      return
    end if
    !
    lsf => lsf_inout
    v   => v_in
    call r_vec%set(esta=esta,emsg=emsg); if ( esta /= 0 ) return
    ! Advect
    write(stdout,'(a)') 'Solving the Hamilton Jacobi equation ...'
    call t_complete%start_timer()
    ! Calculate stable time increment
    d_t = alpha_t * char_fe_dim / v%max_value()
    !
    do i = 1, min(max_num_advect,num_advect)
      call r_vec_setup(esta,emsg); if ( esta /= 0 ) return
    end do
    !
    call t_complete%write_elapsed_time(stdout,'Advection finished, &
    &elapsed time')
    !
    lsf => null()
    v   => null()
    call r_vec%delete(esta,emsg); if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine caluclate_advection
!*****************************************************************************80
end module hamilton_jacobi
