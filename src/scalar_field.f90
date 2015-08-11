module scalar_field
  use types, only: ik, rk, lk, cl, debug, es
  use mesh_data, only : nnod
  use fe_c3d10, only: nelnod, c3d10_t
!*****************************************************************************80
  implicit none
!*****************************************************************************80
  type :: scalar_field_t
    real(rk), allocatable :: values(:)
  contains
    procedure :: set
    procedure :: set_nodal_value
    procedure :: get_element_nodal_values
    procedure :: copy
    procedure :: assemble_element_nodal_values
    procedure :: delete
    procedure :: write_field
    procedure :: max_value
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
  pure subroutine set_nodal_value(self,n,nodal_value,esta,emsg)
    class(scalar_field_t), intent(inout) :: self
    integer(ik), intent(in) :: n
    real(rk), intent(in) :: nodal_value
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    if ( debug ) then
      if ( .not.allocated (self%values) ) then
        esta = -1
        emsg ='Set nodal value: field values not allocated'
        return
      end if
      if ( n > size(self%values) ) then
        esta = -1
        emsg ='Set nodal value: node not found'
        return
      end if
    end if
    !
    self%values(n) = nodal_value
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine set_nodal_value
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
  pure subroutine assemble_element_nodal_values(self,c3d10,nodal_values,&
  &esta,emsg)
    class(scalar_field_t), intent(inout) :: self
    type(c3d10_t), intent(in) :: c3d10
    real(rk), intent(in) :: nodal_values(:)
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
    ! Assemble
    associate( c => c3d10%connectivity )
    self%values(c) = self%values(c) + nodal_values
    end associate
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine assemble_element_nodal_values
!*****************************************************************************80
  pure subroutine delete(self,esta,emsg)
    class(scalar_field_t), intent(inout) :: self
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    if ( debug ) then
      if ( .not.allocated (self%values) ) then
        esta = -1
        emsg ='Delete scalar field: field values not allocated'
        return
      end if
    end if
    deallocate(self%values,stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine delete
!*****************************************************************************80
  subroutine write_field(self,write_unit,esta,emsg)
    class(scalar_field_t), intent(in) :: self
    integer(ik), intent(in) :: write_unit
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    !
    if ( debug ) then
      if ( .not.allocated (self%values) ) then
        esta = -1
        emsg ='Write scalar field: field values not allocated'
        return
      end if
    end if
    !
    do i = 1, size(self%values)
      write(write_unit,'(i0,1x,'//es//')') i, self%values(i)
    end do
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine write_field
!*****************************************************************************80
  pure real(rk) function max_value(self)
    class(scalar_field_t), intent(in) :: self
    if ( .not. allocated(self%values) ) then
      max_value = 0.0_rk
      return
    else
      max_value = maxval(self%values)
    end if
  end function max_value
!*****************************************************************************80
end module scalar_field
