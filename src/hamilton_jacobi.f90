module hamilton_jacobi
!*****************************************************************************80
  use types, only: ik, rk
  use reinitalzation, only: reg_pnorm
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  logical(lk), protected :: configured = .false.
!*****************************************************************************80
  integer(ik) :: num_max_advect = 50 ! Maximal number of advection steps
!*****************************************************************************80
  type(sparse_square_matrix_t), save :: m_mtx
  type(sparse_linear_system_t), save :: linear_system
!*****************************************************************************80
  integer(int64) :: mem_fac_c_mtx = 0
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine set_hamilton_jacobi(num_max_advect_in)
    integer(ik), optional, intent(in) :: num_max_advect_in
    if (present(num_max_advect_in)) num_max_advect = num_max_advect_in
    configured = .true.
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
    real(rk) :: m_mtx(nelnod,nelnod), rtmp_1(nelnod,nelnod)
    real(rk) :: n(nelnod)
    ! Checks
    if ( debug ) then
      if ( .not.size_mtx(c_mtx,nelnod,nelnod) ) then
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
      call outer(n,n,rtmp_1)
      ! integrate
      m_mtx = m_mtx + rtmp_1 * w(p) * det_jac
    end do
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine calc_fe_c_mtx
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
    !$omp private(e,fe_c_mtx) &
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
              cx(indx) = fe_c_mtx(i,j)
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
    call linear_system%solve(job=1,a=c_mtx,mem_used=mem_fac_c_mtx,&
    &esta=esta,emsg=emsg)
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine c_mtx_setup
end module hamilton_jacobi
