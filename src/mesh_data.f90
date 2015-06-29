module mesh_data
!*****************************************************************************80
  use types, only: ik, rk, stdout, es, int64
  use memory_storage, only: size_in_bytes, write_size_of_storage
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
    ! Calculate minimal characteristic dimension
    call min_char_fe_dim(esta,emsg)
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine set_finite_elements
!*****************************************************************************80
! Characteristic finite element dimension
!*****************************************************************************80
  pure subroutine calc_char_fe_dim(c3d10,le,esta,emsg)
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: le
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: a, i, p, p_tmp
    real(rk) :: det_jac, vol_fe
    real(rk) :: b(dom,nelnod), usf(dom,nelnod)
    !
    vol_fe = 0.0_rk
    usf = 0.0_rk
    le = 0.0_rk
    !
    do p = 1, ngp
      p_tmp = p ! Gfortran bug
      call c3d10%gradient(gp_num=p_tmp,b_mtx=b,det_jac=det_jac,&
      &esta=esta,emsg=emsg)
      if ( esta /= 0 ) return
      !
      vol_fe = vol_fe + w(p) * det_jac
      usf = usf + b * w(p) * det_jac
    end do
    usf = 1.0_rk / vol_fe * usf
    do i = 1, dom
      do a = 1, nelnod
        le = le + usf(i,a) * usf(i,a)
      end do
    end do
    le = vol_fe / sqrt(le)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine calc_char_fe_dim
  subroutine min_char_fe_dim(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e
    real(rk) :: le_all(nfe)
    !
    le_all = 0.0_rk
    !$omp parallel do schedule(static,1) &
    !$omp private(e) &
    !$omp shared(finite_elements,le_all,esta,emsg)
    do e = 1, nfe
      if ( esta == 0 ) then
        call calc_char_fe_dim(finite_elements(e),le_all(e),esta,emsg)
      end if
    end do
    !$omp end parallel do
    if ( esta /= 0 ) return
    char_fe_dim = minval(le_all)
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine min_char_fe_dim
!*****************************************************************************80
! Read input statistics
!*****************************************************************************80
  subroutine mesh_data_statistics()
    !
    integer(int64) :: storage
    !
    write(stdout,'(a)') '*** Mesh data statistics ***'
    write(stdout,'(a,i0)') 'Number of nodes is: ', nnod
    write(stdout,'(a,i0)') 'Number of elements is: ', nfe
    write(stdout,'(a,'//es//')') 'Characteristic finite element dimension &
    &is: ', char_fe_dim
    storage = size_in_bytes(nodes) + size_in_bytes(finite_elements)
    write(stdout,'(a,a)') 'Data allocated: ', trim(write_size_of_storage( &
    &storage))
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
