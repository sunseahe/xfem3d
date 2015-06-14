module scalar_field
  use types, only: ik, rk, debug
  use mesh_data, only : nnod
!*****************************************************************************80
  implicit none
!*****************************************************************************80
  type :: scalar_field_t
    real(rk), allocatable :: values(:)
  contains
    procedure :: get_nodal_values
  end type  scalar_field_t
  interface scalar_field_t
    module procedure constructor
  end interface scalar_field_t
!*****************************************************************************80
  public :: scalar_field_t
!*****************************************************************************80
  contains
!*****************************************************************************80
  pure type(scalar_field_t) function constructor()
    scalar_field_t%values
  end function constructor
!*****************************************************************************80
  pure subroutine get_nodal_values(self,e,nodal_values,esta,emsg)
    class(scalar_field_t), intent(in) :: self
    integer(ik), intent(in) :: e
    real(rk), intent(out) :: nodal_values(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    if ( debug ) then
      if ( .not. size(nodal_values)==nelnod ) then
        esta = -1
        emsg ='Get nodal values: Vector size incorrect'
        return
      end if
    end if
    !
  end subroutine get_nodal_values
!*****************************************************************************80
end module scalar_field
