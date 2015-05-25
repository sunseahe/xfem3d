module types
  use iso_fortran_env, only: &
  & ik => int32, &
  & int64, &
#ifndef DOUBLE
  & rk => real32, &
#else
  & rk => real64, &
#endif
  & stdout => output_unit, &
  & stdin => input_unit, &
  & stderr => error_unit
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: cl = 256, lk = 4
!*****************************************************************************80
#ifndef DOUBLE
  character(len=*), parameter :: es = 'es14.6e3'
#else
  character(len=*), parameter :: es = 'es20.15e3'
#endif
!*****************************************************************************80
#ifndef DEBUG
  logical(lk), parameter :: debug = .false.
#else
  logical(lk), parameter :: debug = .true.
#endif
!*****************************************************************************80
  public :: ik, int64, rk, cl, lk, es, stdout, stdin, stderr, debug
!*****************************************************************************80
end module types

module general_routines
  use types
  implicit none
  private
  public :: to_lower, &
  & str2i, &
  & i2str, &
  & str2r, &
  & r2str, &
  & write_dense_mtx_real, &
  & print_error, &
  & f2c_char, &
  & exclude, &
  & real_interval
contains
!*****************************************************************************80
! Print error
!*****************************************************************************80
  subroutine print_error(esta,emsg)
    integer(ik), intent(in) :: esta
    character(len=*), intent(in) :: emsg
    !
    write(stderr,'(a,i0)') 'Error number: ', esta
    write(stderr,'(a,a)') 'Error reason: ', trim(emsg)
    stop -1
    !
  end subroutine print_error
!*****************************************************************************80
! To lower
!*****************************************************************************80
  pure subroutine to_lower(str)
    character(len=*), intent(inout) :: str
    integer(ik) :: i
    do i = 1, len( str )
      select case( str(i:i) )
      case( 'A':'Z' )
        str(i:i) = achar( iachar( str(i:i) ) + 32 )
      end select
    end do
  end subroutine to_lower
!*****************************************************************************80
! String to integer
!*****************************************************************************80
  pure subroutine str2i(str,iout,esta,emsg)
    character(len=*), intent(in) :: str
    integer(ik), intent(out) :: iout
    integer, intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    character(len=cl) :: tmp_str, f
    integer :: n
    !
    tmp_str = adjustl(str)
    n = len_trim(tmp_str)
    write(unit=f,fmt='("(i",i0,")")',iostat=esta, &
    &iomsg=emsg) n
    iout = 0_ik
    read(unit=tmp_str,fmt=f,iostat=esta,iomsg=emsg) iout
    !
  end subroutine str2i
  pure character(len=cl) function i2str(i)
    integer(ik), intent(in) :: i
    write(unit=i2str,fmt='(i0)') i
  end function i2str
!*****************************************************************************80
! String to real
!*****************************************************************************80
  pure subroutine str2r(str,rout,esta,emsg)
    character(len=*), intent(in) :: str
    real(rk), intent(out) :: rout
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    character(len=cl) :: tmp_str, f
    integer :: dot, n
    !
    tmp_str = adjustl(str)
    n = len_trim(tmp_str)
    dot = index(tmp_str, '.')
    if( dot < 1 ) then
      write(unit=f,fmt='("(f",i0,".0)")',iostat=esta,iomsg=emsg) n
    else
      write(unit=f,fmt='("(f",i0,".",i0,")")',iostat=esta, &
      &iomsg=emsg) n, n - dot
    end if
    rout = 0.0e0_rk
    read(unit=tmp_str,fmt=f,iostat=esta,iomsg=emsg) rout
    !
  end subroutine str2r
  pure character(len=cl) function r2str(r)
    real(rk), intent(in) :: r
    write(unit=r2str,fmt='('//es//')') r
  end function r2str
!*****************************************************************************80
  subroutine write_dense_mtx_real(unit,mtx,name)
    integer(ik), intent(in) :: unit
    real(rk), intent(in) :: mtx(:,:)
    character(len=*), intent(in) :: name
    !
    integer(ik) :: i, m
    character(len=*), parameter :: wrf = '(*('//es//',:,","))'
    !
    m = size(mtx,dim=1)
    !
    write(unit,'(a,a)') 'Matrix: ', trim(name)
    do i =1, m
      write(unit,wrf) mtx(i,:)
    end do
    !
  end subroutine write_dense_mtx_real
!*****************************************************************************80
! Exclude
!*****************************************************************************80
  pure logical(lk) function exclude(x,y)
    logical(lk), intent(in) :: x,y
    exclude = ( x .and. .not.y ) .or. ( .not.x .and. y )
  end function exclude
!*****************************************************************************80
! Real interval
!*****************************************************************************80
  pure logical(lk) function real_interval(a,b,var)
    real(rk), intent(in) :: a, b, var
    real_interval = ( a <= var ) .and. ( var <= b )
  end function real_interval
!*****************************************************************************80
! Convert to c characters
!*****************************************************************************80
  pure subroutine f2c_char(f_char_in,c_char_out,esta,emsg)
    use iso_c_binding, only: c_char, c_null_char
    character(len=*), intent(in) :: f_char_in
    character(len=1,kind=c_char), &
    & allocatable, intent(out) :: c_char_out(:)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: n, i
    !
    n = len_trim(f_char_in)
    allocate(c_char_out(n+1),stat=esta,errmsg=emsg)
    if( esta /= 0 ) return
    do i = 1, n
      c_char_out(i) = f_char_in(i:i)
    end do
    c_char_out(n+1) = c_null_char
    !
  end subroutine f2c_char
!*****************************************************************************80
end module general_routines

module memory_storage
  use types
  implicit none
  private
!*****************************************************************************80
  interface size_in_bits
    module procedure size_in_bits_scalar
    module procedure size_in_bits_vector
    module procedure size_in_bits_matrix
  end interface size_in_bits
!*****************************************************************************80
  public :: size_in_bits, write_size_of_storage
!*****************************************************************************80
contains
!*****************************************************************************80
  pure integer(int64) function size_in_bits_scalar(scalar)
    class(*), intent(in) :: scalar
    size_in_bits_scalar = storage_size(scalar,kind=int64)
  end function size_in_bits_scalar
!*****************************************************************************80
  pure integer(int64) function size_in_bits_vector(vector)
    class(*), intent(in) :: vector(:)
    size_in_bits_vector = size(vector) * storage_size(vector,kind=int64)
  end function size_in_bits_vector
!*****************************************************************************80
  pure integer(int64) function size_in_bits_matrix(matrix)
    class(*), intent(in) :: matrix(:,:)
    size_in_bits_matrix = size(matrix,dim=1) * &
    & size(matrix,dim=2) * storage_size(matrix,kind=int64)
  end function size_in_bits_matrix
!*****************************************************************************80
  pure character(len=cl) function write_size_of_storage(x)
    integer(int64), intent(in) :: x
    real :: storage
    storage = 1.25e-7 * x
    !
    if ( storage < 1.0e+3 ) then
      write(write_size_of_storage,'(f4.1,a)') storage, ' MB'
    else
      storage = 1.0e-3 * storage
      write(write_size_of_storage,'(f0.1,a)') storage, ' GB'
    end if
    !
  end function write_size_of_storage
!*****************************************************************************80
end module memory_storage


