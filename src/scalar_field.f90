module scalar_field
  use types, only: ik, rk, lk, cl, debug
  use mesh_data, only : nnod
  use fe_c3d10, only: nelnod, c3d10_t
!*****************************************************************************80
  implicit none
!*****************************************************************************80
  type :: scalar_field_t
    real(rk), allocatable :: values(:)
  contains
    procedure :: set
    procedure :: get_element_nodal_values
    procedure :: copy
    procedure :: assemble_element_nodal_values
  end type  scalar_field_t
!*****************************************************************************80
  public :: scalar_field_t
!*****************************************************************************80
  contains
!*****************************************************************************80
  pure subroutine set(self,values,esta,emsg)
    class(scalar_field_t), intent(inout) :: self
    real(rk), optional, intent(in) :: values(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    if ( debug ) then
      if ( present(values) ) then
        if ( .not. size(values)==nnod ) then
          esta = -1
          emsg ='Set scalar field: vector size incorrect'
          return
        end if
      end if
    end if
    !
    if ( .not. allocated (self%values) ) then
      allocate(self%values(nnod),stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
    end if
    if ( present(values) ) then
      self%values = values
    else
      self%values = 0.0_rk
    end if
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine set
!*****************************************************************************80
  pure subroutine copy(self,other,esta,emsg)
    class(scalar_field_t), intent(inout) :: self
    type(scalar_field_t), intent(in) :: other
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: n
    !
    if ( .not. allocated (self%values) ) then
      n = size(other%values)
      allocate(self%values(n),source=other%values,stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
    else
      self%values = other%values
    end if
    !
  end subroutine copy
!*****************************************************************************80
  pure subroutine get_element_nodal_values(self,c3d10,nodal_values,esta,emsg)
    class(scalar_field_t), intent(in) :: self
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: nodal_values(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: connectivity(nelnod)
    ! Checks
    if ( debug ) then
      if ( .not. size(nodal_values)==nelnod ) then
        esta = -1
        emsg ='Get element nodal values: vector size incorrect'
        return
      end if
      if ( .not.allocated (self%values) ) then
        esta = -1
        emsg ='Get element nodal values: field values not allocated'
        return
      end if
    end if
    !
    call c3d10%get_connectivity(connectivity,esta,emsg)
    if ( debug ) then
      if ( esta /= 0 ) return
    end if
    nodal_values = self%values(connectivity)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine  get_element_nodal_values
!*****************************************************************************80
  pure subroutine assemble_element_nodal_values(self,c3d10,nodal_values,&
  &esta,emsg)
    class(scalar_field_t), intent(inout) :: self
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(in) :: nodal_values(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: connectivity(nelnod)
    !
    ! Checks
    if ( debug ) then
      if ( .not. size(nodal_values)==nelnod ) then
        esta = -1
        emsg ='Get element nodal values: vector size incorrect'
        return
      end if
      if ( .not.allocated (self%values) ) then
        esta = -1
        emsg ='Get element nodal values: field values not allocated'
        return
      end if
    end if
    call c3d10%get_connectivity(connectivity,esta,emsg)
    if ( debug ) then
      if ( esta /= 0 ) return
    end if
    ! Assemble
    self%values(connectivity) = self%values(connectivity) + nodal_values
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine assemble_element_nodal_values
!*****************************************************************************80
end module scalar_field
