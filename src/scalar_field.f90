module scalar_field
  use types, only: ik, rk, lk, cl, debug
  use mesh_data, only : nnod
  use fe_c3d10, only: nelnod, c3d10_t
!*****************************************************************************80
  implicit none
!*****************************************************************************80
  type :: scalar_field_t
    private
    real(rk), allocatable :: values(:)
  contains
    procedure :: set
    procedure :: get_element_nodal_values
    procedure :: copy
    generic :: assignment(=) => copy
  end type  scalar_field_t
!*****************************************************************************80
  public :: scalar_field_t
!*****************************************************************************80
  contains
!*****************************************************************************80
  subroutine set(self,values,esta,emsg)
    class(scalar_field_t), intent(out) :: self
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
  subroutine copy(self,other)
    class(scalar_field_t), intent(out) :: self
    type(scalar_field_t), intent(in) :: other
    allocate(self%values,source=other%values)
  end subroutine copy
!*****************************************************************************80
  pure subroutine get_element_nodal_values(self,c3d10,nodal_values,esta,emsg)
    class(scalar_field_t), intent(in) :: self
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(out) :: nodal_values(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
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
    nodal_values = self%values(c3d10%connectivity)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine  get_element_nodal_values
!*****************************************************************************80
end module scalar_field
