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
      procedure, non_overridable :: add_list
      procedure, non_overridable :: get_nitem
      procedure, non_overridable :: clean
      procedure, non_overridable :: is_empty
      procedure, non_overridable :: get_first
      procedure, non_overridable :: get_current
      procedure, non_overridable :: reset
      procedure, non_overridable :: set_next
  end type list
!*****************************************************************************80
! Character linked list
!*****************************************************************************80
  type, extends(list) :: char_ll
    private
  contains
    procedure :: add => add_char_ll
    procedure :: fill_array => fill_array_char_ll
  end type char_ll
!*****************************************************************************80
! Integer linked list
!*****************************************************************************80
  type, extends(list) :: integer_ll
    private
  contains
    procedure :: add => add_integer_ll
    procedure :: fill_array => fill_array_integer_ll
  end type integer_ll
!*****************************************************************************80
! Real linked list
!*****************************************************************************80
  type, extends(list) :: real_ll
    private
  contains
    procedure :: add => add_real_ll
    procedure :: fill_array => fill_array_real_ll
  end type real_ll
!*****************************************************************************80
  public :: link, list, char_ll, integer_ll, real_ll
!*****************************************************************************80
contains
  !
  pure subroutine add_list(self,any_item,esta,emsg)
    class(list), intent(inout) :: self
    class(*), intent(in) :: any_item
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    if ( .not. associated(self%head) ) then
      allocate(self%head,stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      self%tail => self%head
      allocate(self%tail%any_item,source=any_item,&
      &stat=esta,errmsg=emsg); if( esta /= 0 ) return
      self%nitem = self%nitem + 1
    else
      allocate(self%tail%next,stat=esta,errmsg=emsg)
      if( esta /= 0 ) return
      self%tail => self%tail%next
      self%tail%next => null()
      allocate(self%tail%any_item,source=any_item,&
      & stat=esta,errmsg=emsg); if( esta /= 0 ) return
      self%nitem = self%nitem + 1
    end if
    esta = 0
    emsg = ''
    !
  end subroutine add_list
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
  pure subroutine get_first(self,first,esta,emsg)
    class(list), intent(in) :: self
    class(*), allocatable, intent(inout) :: first
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    if ( .not. associated(self%head) ) then
      esta = -1
      emsg = 'Linked list: no items'
    else
      if ( allocated(first) ) deallocate(first,&
      &stat=esta,errmsg=emsg); if( esta /= 0 ) return
      allocate(first,source=self%head%any_item,&
      &stat=esta,errmsg=emsg); if( esta /= 0 ) return
    end if
    !
  end subroutine get_first
  !
  pure subroutine get_current(self,current,esta,emsg)
    class(list), intent(in) :: self
    class(*), allocatable, intent(inout) :: current
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    if ( .not. associated(self%head) ) then
      esta = -1
      emsg = 'Linked list: no items'
    else
      if ( allocated(current) ) deallocate(current,&
      &stat=esta,errmsg=emsg); if( esta /= 0 ) return
      allocate(current,source=self%current%any_item,&
      &stat=esta,errmsg=emsg); if( esta /= 0 ) return
    end if
    !
  end subroutine get_current
  !
  pure subroutine reset(self)
    class(list), intent(inout) :: self
    if ( .not. associated(self%head) ) return
    self%current => self%head
  end subroutine reset
  !
  pure subroutine set_next(self)
    class(list), intent(inout) :: self
    if ( .not. associated(self%head) ) return
    self%current => self%current%next
  end subroutine set_next
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
  pure subroutine add_char_ll(self,arg,esta,emsg)
    class(char_ll), intent(inout) :: self
    character(len=*), intent(in) :: arg
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,esta,emsg)
    !
  end subroutine add_char_ll
!*****************************************************************************80
  pure subroutine fill_array_char_ll(self,arr,esta,emsg)
    class(char_ll), intent(inout) :: self
    character(len=cl), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      esta = -1
      emsg = 'List is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=esta,errmsg=emsg)
    if( esta /= 0 ) return
    call self%reset()
    do i = 1, self%get_nitem()
      call self%get_current(curr,esta,emsg)
      if( esta /= 0 ) return
      select type(curr)
      type is( character(*) )
        arr(i) = curr
      class default
        esta = -1
        emsg = 'Wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
    ! Clean
    call self%clean(esta,emsg)
    if( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine fill_array_char_ll
!*****************************************************************************80
  pure subroutine add_integer_ll(self,arg,esta,emsg)
    class(integer_ll), intent(inout) :: self
    integer(ik), intent(in) :: arg
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,esta,emsg)
    !
  end subroutine add_integer_ll
!*****************************************************************************80
  pure subroutine fill_array_integer_ll(self,arr,esta,emsg)
    class(integer_ll), intent(inout) :: self
    integer(ik), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      esta = -1
      emsg = 'Integer linked list: list is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=esta,errmsg=emsg)
    if( esta /= 0 ) return
    call self%reset()
    do i = 1, self%get_nitem()
      call self%get_current(curr,esta,emsg)
      if( esta /= 0 ) return
      select type(curr)
      type is( integer(ik) )
        arr(i) = curr
      class default
        esta = -1
        emsg = 'Integer linked list: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
    ! Clean
    call self%clean(esta,emsg)
    if( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine fill_array_integer_ll
!*****************************************************************************80
  pure subroutine add_real_ll(self,arg,esta,emsg)
    class(real_ll), intent(inout) :: self
    real(rk), intent(in) :: arg
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,esta,emsg)
    !
  end subroutine add_real_ll
!*****************************************************************************80
  pure subroutine fill_array_real_ll(self,arr,esta,emsg)
    class(real_ll), intent(inout) :: self
    real(rk), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      esta = -1
      emsg = 'Real linked list: list is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=esta,errmsg=emsg)
    if( esta /= 0 ) return
    !
    call self%reset()
    do i = 1, self%get_nitem()
      call self%get_current(curr,esta,emsg)
      if( esta /= 0 ) return
      select type(curr)
      type is( real(rk) )
        arr(i) = curr
      class default
        esta = -1
        emsg = 'Real linked list: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
     ! Clean
    call self%clean(esta,emsg)
    if( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
  end subroutine fill_array_real_ll
!*****************************************************************************80
end module ll
