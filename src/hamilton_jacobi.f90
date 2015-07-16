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
    real(rk) :: m_mtx(nelnod,nelnod)
    real(rk) :: n(nelnod)
    ! Checks
    if ( debug ) then
      if ( .not.size_mtx(c_mtx,nelnod,nelnod) ) then
        esta = -1
        emsg ='Cal fe c mtx: El cmtx size incorrect'
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
      call c3d10%gradient(gp_num=p_tmp,b_mtx=b,det_jac=det_jac,&
      &esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      !
      call outer(n,n,rtmp1)
      call gemm(transpose(b),b,rtmp2)
      ! M gamma
      call sdf_0%get_element_nodal_values(c3d10,sdf_0_nv,esta,emsg)
      sdf_0_gp = dot(n,sdf_0_nv)
      m_gam_mtx = rtmp1 * dirac_delta(sdf_0_gp)
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
end module hamilton_jacobi
