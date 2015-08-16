module mesh_data
!*****************************************************************************80
  use types, only: ik, rk, stdout, es, log_file
  use point, only: dom, zero_pnt, point_3d_t
  use fe_c3d10, only: nelnod, c3d10_t, ngp, w
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  character(len=*), parameter :: fe_type = 'c3d10'
  integer(ik), protected :: nnod = 0, nfe = 0
  real(rk), protected :: char_fe_dim = 0.0_rk
  type(point_3d_t), allocatable, protected :: nodes(:)
  type(c3d10_t), allocatable, protected :: finite_elements(:)
!*****************************************************************************80
  public :: fe_type, nnod, nfe, nodes, set_nodes, finite_elements, &
  & set_finite_elements, mesh_data_finish, mesh_data_statistics, char_fe_dim
!*****************************************************************************80
  contains
!*****************************************************************************80
! Set nodes
!*****************************************************************************80
  subroutine set_nodes(nodes_in,esta,emsg)
    type(point_3d_t), allocatable, intent(inout) :: nodes_in(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    if ( .not. allocated(nodes_in) ) then
      esta = -1
      emsg = 'Set nodes: nodes_in not allocated'
      return
    end if
    nnod = size(nodes_in)
    call move_alloc(nodes_in,nodes)
    !
  end subroutine set_nodes
!*****************************************************************************80
! Set elements
!*****************************************************************************80
  subroutine set_finite_elements(finite_elements_in,esta,emsg)
    type(c3d10_t), allocatable, intent(inout) :: finite_elements_in(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i, j
    integer(ik) :: el_conn(nelnod)
    !
    if ( .not. allocated(finite_elements_in) ) then
      esta = -1
      emsg = 'Set finite elements: finite_elements_in not allocated'
      return
    end if
    nfe = size(finite_elements_in)
    call move_alloc(finite_elements_in,finite_elements)
    ! Copy nodes to finite elements
    do i = 1, nfe
      el_conn = finite_elements(i)%connectivity
      do j = 1, nelnod
        finite_elements(i)%nodes(j) = nodes(el_conn(j))
      end do
    end do
    ! Check if jacobian determinant is positive for all elements
    call jac_det_check(esta,emsg)
    if ( esta /= 0 ) return
    ! Calculate minimal characteristic dimension
    call min_char_fe_dim()
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine set_finite_elements
!*****************************************************************************80
! Jacobian determinant check
!*****************************************************************************80
  subroutine jac_det_check(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e
    real(rk), parameter :: ot = 1.0_rk/3.0_rk
    real(rk) :: det_jac
    type(point_3d_t), parameter :: centroid = point_3d_t([ ot, ot, ot ])
    !
    check: do e = 1, nfe
      call finite_elements(e)%main_values(xi_coo_pnt=centroid,&
      &det_jac=det_jac)
      if ( det_jac <= 0.0_rk ) then
        esta = -1
        write(emsg,'(a,i0)') 'Zero or negative determinant for finite &
        &element no:', e
        exit check
      end if
    end do check
    esta = 0
    emsg = ''
    !
  end subroutine jac_det_check
!*****************************************************************************80
! Characteristic finite element dimension
!*****************************************************************************80
  pure subroutine calc_char_fe_dim(c3d10,le)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: le
    !
    integer(ik) :: p
    real(rk) :: det_jac, vol
    !
    vol = 0.0_rk
    !
    do p = 1, ngp
      call c3d10%main_values(gp_num=p,det_jac=det_jac)
      vol = vol + w(p) * det_jac
    end do
    !print*, vol
    !le = (12.0_rk * vol)**(1.0_rk/3.0_rk) / 2.0_rk
    le = vol**(1.0_rk/3.0_rk)
    !
  end subroutine calc_char_fe_dim
  subroutine min_char_fe_dim()
    !
    integer(ik) :: e
    real(rk) :: le, le_sum
    !
    le_sum = 0.0_rk
    !$omp parallel do schedule(static,1) &
    !$omp private(e,le) &
    !$omp shared(finite_elements,le_sum)
    do e = 1, nfe
      call calc_char_fe_dim(finite_elements(e),le)
      !$omp critical
      le_sum = le_sum + le
      !$omp end critical
    end do
    !$omp end parallel do
    char_fe_dim = le_sum / nfe
    !
  end subroutine min_char_fe_dim
!*****************************************************************************80
! Read input statistics
!*****************************************************************************80
  subroutine mesh_data_statistics()
    !
    write(log_file,'(a)') '*** Mesh data statistics ***'
    write(log_file,'(a,i0)') 'Number of nodes is: ', nnod
    write(log_file,'(a,i0)') 'Number of finite elements is: ', nfe
    write(log_file,'(a,'//es//')') 'Characteristic finite element dimension &
    &is: ', char_fe_dim
    !
  end subroutine mesh_data_statistics
!*****************************************************************************80
! Finish
!*****************************************************************************80
  subroutine mesh_data_finish(esta,emsg)
    !
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    if ( allocated(nodes) ) deallocate(nodes,stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    if ( allocated(finite_elements) ) deallocate(finite_elements,&
    &stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    esta = 0
    emsg = ''
    !
  end subroutine mesh_data_finish
!*****************************************************************************80
end module mesh_data
