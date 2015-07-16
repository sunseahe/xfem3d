module ll
  use types, only: ik, rk, lk, cl
  implicit none
  private
!*****************************************************************************80
  type :: link
    class(*), allocatable :: any_item
    type(link), pointer :: next => null()
  end type link
!*****************************************************************************80
  type, abstract :: list
    private
    integer(ik) :: nitem = 0
    type(link), pointer :: head => null()
    type(link), pointer :: current => null()
    type(link), pointer :: tail => null()
    contains
      procedure, non_overridable :: add
      procedure, non_overridable :: get_nitem
      procedure, non_overridable :: clean
      procedure, non_overridable :: is_empty
  end type list
!*****************************************************************************80
  contains
!*****************************************************************************80
  pure subroutine add(self,any_item,esta,emsg)
    class(list), intent(inout) :: self
    class(*), intent(in) :: any_item
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    if ( .not. associated(self%head) ) then
      allocate(self%head,stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      self%tail => self%head
      allocate(self%tail%any_item,source=any_item,stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      self%nitem = self%nitem + 1
    else
      allocate(self%tail%next,stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      self%tail => self%tail%next
      self%tail%next => null()
      allocate(self%tail%any_item,source=any_item,stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      self%nitem = self%nitem + 1
    end if
    esta = 0
    emsg = ''
    !
  end subroutine add
  !
  pure integer(ik) function get_nitem(self)
    class(list), intent(in) :: self
    get_nitem = self%nitem
  end function get_nitem
  !
  pure logical(lk) function is_empty(self)
    class(list), intent(in) :: self
    is_empty = .true.
    if ( associated(self%head) ) is_empty = .false.
  end function is_empty
  !
  pure subroutine fill_array(self,arr,esta,emsg)
    class(list), intent(inout) :: self
    class(*), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    type(link), pointer :: curr
    !
    if ( self%is_empty() ) then
      esta = -1
      emsg = 'Linked list: empty list'
      return
    end if
    allocate(arr(self%nitem),stat=esta,errmsg=emsg)
    if( esta /= 0 ) return
    !
    do i = 1, self%nitem
      allocate(curr,source=self%get_current())
      select type(curr); type is( il )
        arr(i) = curr
      end select
      call self%set_next()
    end do
    !
  end subroutine fill_array
  !
  pure subroutine clean(self,esta,emsg)
    class(list), intent(inout) :: self
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    type(link), pointer :: current, next
    !
    if ( .not. associated(self%head) ) return
    current => self%head
    next => current%next
    do
      deallocate(current,stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      current => next
      if( .not. associated(current) ) exit
      next => current%next
    end do
    self%head => null()
    self%current => null()
    self%tail => null()
    self%nitem = 0
    !
    esta = 0
    emsg = ''
    !
  end subroutine clean

!*****************************************************************************80
end module ll
