module ll
  use types
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
  type :: real_arr_t
    real(rk), allocatable :: dat(:)
  end type real_arr_t
!*****************************************************************************80
  public :: link, list, char_ll, integer_ll, real_ll
!*****************************************************************************80
contains
  !
  pure subroutine add_list(self,any_item,istat,emsg)
    class(list), intent(inout) :: self
    class(*), intent(in) :: any_item
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    if ( .not. associated(self%head) ) then
      allocate(self%head,stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      self%tail => self%head
      allocate(self%tail%any_item,source=any_item,stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      self%nitem = self%nitem + 1
    else
      allocate(self%tail%next,stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      self%tail => self%tail%next
      self%tail%next => null()
      allocate(self%tail%any_item,source=any_item,stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      self%nitem = self%nitem + 1
    end if
    istat = 0
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
  pure function get_first(self)
    class(list), intent(in) :: self
    class(*), allocatable :: get_first
    if ( .not. associated(self%head) ) then
      allocate(get_first,source=0) ! Something for error
    else
      allocate(get_first,source=self%head%any_item)
    end if
  end function get_first
  !
  pure function get_current(self)
    class(list), intent(in) :: self
    class(*), allocatable :: get_current
    if ( .not. associated(self%head) ) then
      allocate(get_current,source=0) ! Something for error
    else
      allocate(get_current,source=self%current%any_item)
    end if
  end function get_current
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
  pure subroutine clean(self,istat,emsg)
    class(list), intent(inout) :: self
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    type(link), pointer :: current, next
    !
    if ( .not. associated(self%head) ) return
    current => self%head
    next => current%next
    do
      deallocate(current,stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      current => next
      if( .not. associated(current) ) exit
      next => current%next
    end do
    self%head => null()
    self%current => null()
    self%tail => null()
    self%nitem = 0
    !
    istat = 0
    emsg = ''
    !
  end subroutine clean
!*****************************************************************************80
  pure subroutine add_char_ll(self,arg,istat,emsg)
    class(char_ll), intent(inout) :: self
    character(len=*), intent(in) :: arg
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,istat,emsg)
    !
  end subroutine add_char_ll
!*****************************************************************************80
  pure subroutine fill_array_char_ll(self,arr,istat,emsg)
    class(char_ll), intent(inout) :: self
    character(len=cl), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      istat = -1
      emsg = 'List is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=istat,errmsg=emsg)
    if( istat /= 0 ) return
    call self%reset()
    do i = 1, self%get_nitem()
      if(allocated(curr)) then
        deallocate(curr,stat=istat,errmsg=emsg)
        if( istat /= 0 ) return
      end if
      allocate(curr,source=self%get_current(),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      select type(curr)
      type is( character(*) )
        arr(i) = curr
      class default
        istat = -1
        emsg = 'Wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
  end subroutine fill_array_char_ll
!*****************************************************************************80
  pure subroutine add_integer_ll(self,arg,arg_arr,istat,emsg)
    class(integer_ll), intent(inout) :: self
    integer(ik), optional, intent(in) :: arg
    integer(ik), optional, intent(in) :: arg_arr(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    if ( present(arg) ) then
      call self%add_list(arg,istat,emsg)
    else if( present(arg_arr) ) then
    else
      istat = -1
      emsg = 'Integer linked list: incorrect input arguments'
      return
    end if
    !
  end subroutine add_integer_ll
!*****************************************************************************80
  pure subroutine fill_array_integer_ll(self,arr,istat,emsg)
    class(integer_ll), intent(inout) :: self
    integer(ik), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      istat = -1
      emsg = 'List is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=istat,errmsg=emsg)
    if( istat /= 0 ) return
    call self%reset()
    do i = 1, self%get_nitem()
      if(allocated(curr)) then
        deallocate(curr,stat=istat,errmsg=emsg)
        if( istat /= 0 ) return
      end if
      allocate(curr,source=self%get_current(),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      select type(curr)
      type is( integer(ik) )
        arr(i) = curr
      class default
        istat = -1
        emsg = 'Wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
  end subroutine fill_array_integer_ll
!*****************************************************************************80
  pure subroutine add_real_ll(self,arg,arg_arr,istat,emsg)
    class(real_ll), intent(inout) :: self
    real(rk), optional, intent(in) :: arg
    real(rk), optional, intent(in) :: arg_arr(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    type(real_arr_t) :: real_arr
    !
    if ( present(arg) ) then
      call self%add_list(arg,istat,emsg)
    else if( present(arg_arr) ) then
      istat = -1
      emsg = 'Real linked list: not implemented'
      return
      allocate(real_arr%dat(size(arg_arr)),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      call self%add_list(real_arr,istat,emsg)
    else
      istat = -1
      emsg = 'Real linked list: incorrect input arguments'
      return
    end if
    !
  end subroutine add_real_ll
!*****************************************************************************80
  pure subroutine fill_array_real_ll(self,arr,istat,emsg)
    class(real_ll), intent(inout) :: self
    real(rk), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      istat = -1
      emsg = 'Real linked list: list is empty'
      return
    end if
    allocate(arr(self%get_nitem()),stat=istat,errmsg=emsg)
    if( istat /= 0 ) return
    !
    call self%reset()
    do i = 1, self%get_nitem()
      if(allocated(curr)) then
        deallocate(curr,stat=istat,errmsg=emsg)
        if( istat /= 0 ) return
      end if
      allocate(curr,source=self%get_current(),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      select type(curr)
      type is( real(rk) )
        arr(i) = curr
      class default
        istat = -1
        emsg = 'Real linked list: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
  end subroutine fill_array_real_ll
!*****************************************************************************80
end module ll

