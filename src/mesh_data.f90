module mesh_data
!*****************************************************************************80
  use types, only: ik, rk, stdout, int64
  use memory_storage, only: size_in_bytes, write_size_of_storage
  use point, only: dom, zero_pnt, point_3d_t, point_3d_t_ll
  use fe_c3d10, only: nelnod, c3d10_t, c3d10_t_ll
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  character(len=*), parameter :: fe_type = 'c3d10'
  integer(ik), protected :: nnod = 0, nfe = 0
  real(rk), protected :: char_length = 0.0_rk
  type(point_3d_t), allocatable, protected :: nodes(:)
  type(c3d10_t), allocatable, protected :: finite_elements(:)
!*****************************************************************************80
  public :: fe_type, nnod, nfe, nodes, set_nodes, finite_elements, &
  & set_finite_elements, mesh_data_finish, mesh_data_statistics
!*****************************************************************************80
  contains
!*****************************************************************************80
! Set nodes
!*****************************************************************************80
  subroutine set_nodes(nodes_ll,esta,emsg)
    type(point_3d_t_ll), intent(inout) :: nodes_ll
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    call nodes_ll%fill_array(nodes,esta,emsg); if ( esta /= 0 ) return
    nnod = size(nodes)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine set_nodes
!*****************************************************************************80
! Set elements
!*****************************************************************************80
  subroutine set_finite_elements(fe_ll,esta,emsg)
    type(c3d10_t_ll), intent(inout) :: fe_ll
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i, j
    integer(ik) :: el_conn(nelnod)
    !
    call fe_ll%fill_array(finite_elements,esta,emsg)
    if ( esta /= 0 ) return
    nfe = size(finite_elements)
    ! Copy nodes to finite elements
    do i = 1, nfe
      el_conn = finite_elements(i)%connectivity
      do j = 1, nelnod
        finite_elements(i)%nodes(j) = nodes(el_conn(j))
      end do
    end do
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine set_finite_elements
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
