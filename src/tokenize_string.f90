module tokenize_string
  use types, only: ik, cl
  use ll, only: char_ll
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
    type(char_ll) :: words_list
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
            call words_list%add(current_word,esta,emsg)
            if( esta /= 0 ) return
          end if
        end if
        current_word = ''
        pos1 = 1
      end if
    end do
    ! Set the last word if not empty string
    if( len_trim(current_word) /= 0 ) then
      current_word = adjustl(current_word)
      call words_list%add(current_word,esta,emsg)
      if( esta /= 0 ) return
    end if
    ! Copy to array
    if( .not. words_list%is_empty() ) then
      call words_list%fill_array(words,esta,emsg)
      if( esta /= 0 ) return
      ! Empty list
      call words_list%clean(esta,emsg)
      if( esta /= 0 ) return
      ! Sucess
      esta = 0
      emsg = ''
    else
      ! Empty list
      esta = 1
      emsg = 'Tokenize: empty list'
    end if
    !
  end subroutine tokenize
!*****************************************************************************80
end module tokenize_string
