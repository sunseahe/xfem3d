module tokenize_string
  use types, only: ik, cl
  implicit none
  private
!*****************************************************************************80
  character(len=1), parameter :: tab   = char(9)
  character(len=1), parameter :: space = char(32)
  character(len=1), parameter :: comma = char(44)
  character(len=1), parameter :: equal = char(61)
  character(len=1), parameter :: semicolon = char(59)
!*****************************************************************************80
  public :: tab, space, comma, equal, semicolon
  public :: tokenize
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine tokenize(str,delimiter,words,esta,emsg)
    character(len=*), intent(in) :: str
    character(len=1), intent(in) :: delimiter
    character(len=cl), allocatable, intent(out) :: words(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i, ns, pos1
    character(len=cl) :: current_word
    !type(char_ll) :: words_list
    !
    ns = len_trim(str); pos1 = 1
    current_word = ''
    !
    do i = 1, ns
      if( .not. ( str(i:i) == delimiter ) ) then
        current_word(pos1:pos1) = str(i:i)
        pos1 = pos1 + 1
      else
        if( pos1 > 1 ) then
          if( len_trim(current_word) /= 0 ) then
            current_word = adjustl(current_word)
            call add_char(words,current_word,esta,emsg)
            if ( esta /= 0 ) return
          end if
        end if
        current_word = ''
        pos1 = 1
      end if
    end do
    ! Set the last word if not empty string
    if( len_trim(current_word) /= 0 ) then
      current_word = adjustl(current_word)
      call add_char(words,current_word,esta,emsg)
      if ( esta /= 0 ) return
    end if
    ! Empty list
    if ( size(words) == 0 ) then
      esta = 1
      emsg = 'Tokenize: empty list'
    end if
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine tokenize
!*****************************************************************************80
  pure subroutine add_char(vector,element,esta,emsg)
    character(len=cl), allocatable, intent(inout) :: vector(:)
    character(len=cl), intent(in) :: element
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: n
    character(len=cl), allocatable :: tmp_vector(:)
    !
    if ( allocated(vector) ) then
      n = size(vector,dim=1)
      allocate(tmp_vector(n+1),stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
      tmp_vector(1:n) = vector
      tmp_vector(n+1) = element
      call move_alloc(tmp_vector,vector)
    else
      allocate(vector(1),source=element,stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
    end if
    !
  end subroutine add_char
!*****************************************************************************80
end module tokenize_string
