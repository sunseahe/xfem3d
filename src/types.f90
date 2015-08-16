module types
  use iso_fortran_env, only: &
  & int8, &
  & ik => int32, &
  & int64, &
#ifndef _DOUBLE
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
  real(rk), parameter :: eps = epsilon(1.0_rk)
!*****************************************************************************80
#ifndef _DOUBLE
  character(len=*), parameter :: es = 'es14.6e3'
#else
  character(len=*), parameter :: es = 'es23.15e3'
#endif
!*****************************************************************************80
#ifndef _DEBUG
  logical(lk), parameter :: debug = .false.
#else
  logical(lk), parameter :: debug = .true.
#endif
!*****************************************************************************80
  character(len=1), parameter :: nl = achar(10) ! New line character
!*****************************************************************************80
  integer(ik), parameter :: log_file = 11
!*****************************************************************************80
  public :: ik, &
  & int8, &
  & int64, &
  & rk, &
  & cl, &
  & lk, &
  & eps, &
  & es, &
  & stdout, &
  & stdin, &
  & stderr, &
  & debug, &
  & nl, &
  & log_file
!*****************************************************************************80
end module types

module general_routines
!*****************************************************************************80
  use types, only: ik, rk, cl, lk, es, stderr
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  type :: time
    real(rk) :: saved_time
    contains
      procedure :: start_timer => start_timer_sub
      procedure :: elapsed_time => elapsed_time_fnk
      procedure, nopass, private :: ctime_fnk
      procedure, nopass :: print_time => print_time_fnk
      procedure :: write_elapsed_time
  end type time
!*****************************************************************************80
  interface resize_vec
    module procedure :: resize_ivec
    module procedure :: resize_rvec
  end interface resize_vec
!*****************************************************************************80
  public :: to_lower, &
  & str2i, &
  & i2str, &
  & str2r, &
  & r2str, &
  & write_dense_mtx_real, &
  & print_error, &
  & f2c_char, &
  & size_mtx, &
  & exclude, &
  & real_interval, &
  & outer, &
  & resize_vec, &
  & pnorm, &
  & strip_char, &
  & time, &
  & change_extension
!*****************************************************************************80
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
    integer(ik) :: i
    real(rk) :: row(size(mtx,dim=2))
    character(len=*), parameter :: wrf = '(*('//es//',:,","))'
    !
    write(unit,'(a,a)') 'Matrix: ', trim(name)
    do i =1, size(mtx,dim=1)
      row = mtx(i,:)
      write(unit,wrf) row
    end do
    !
  end subroutine write_dense_mtx_real
!*****************************************************************************80
! Exclude
!*****************************************************************************80
  pure logical(lk) function size_mtx(mtx,a,b)
    class(*), intent(in) :: mtx(:,:)
    integer(ik), intent(in) :: a, b
    size_mtx = (size(mtx,dim=1)==a) .or. (size(mtx,dim=2)==b)
  end function size_mtx
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
  pure subroutine outer(x,y,res)
    real(rk), intent(in)  :: x(:), y(:)
    real(rk), intent(out) :: res(:,:)
    integer :: m,n
    m = size(x)
    n = size(y)
    res = 0.0_rk
    res = spread(x,dim=2,ncopies=n) * spread(y,dim=1,ncopies=m)
  end subroutine outer
!*****************************************************************************80
! Resize vector
!*****************************************************************************80
  pure subroutine resize_ivec(n,vec,esta,emsg)
    integer, intent(in) :: n ! new size
    integer(ik), allocatable, intent(inout) :: vec(:)  ! input vector
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik), allocatable :: tmp(:)
    !
    allocate(tmp(1:n),source=vec(1:n),stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    call move_alloc(from=tmp,to=vec)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine resize_ivec
  pure subroutine resize_rvec(n,vec,esta,emsg)
    integer, intent(in) :: n ! new size
    real(rk), allocatable, intent(inout) :: vec(:)  ! input vector
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    real(rk), allocatable :: tmp(:)
    !
    allocate(tmp(1:n),source=vec(1:n),stat=esta,errmsg=emsg)
    if ( esta /= 0 ) return
    call move_alloc(from=tmp,to=vec)
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine resize_rvec
!*****************************************************************************80
! P norm
!*****************************************************************************80
  pure function pnorm(p,x) result(res)
    integer(ik), intent(in) :: p
    real(rk), intent(in) :: x(:)
    real(rk) :: res
    res = (sum(abs(x)**p))**(1.0_rk/p)
  end function pnorm
!*****************************************************************************80
! Strip characters
!*****************************************************************************80
  elemental pure subroutine strip_char(string,set)
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: set
    integer(ik) :: old, new, stride
    old = 1; new = 1
    do
      stride = scan( string( old : ), set )
      if ( stride > 0 ) then
        string( new : new+stride-2 ) = string( old : old+stride-2 )
        old = old+stride
        new = new+stride-1
      else
        string( new : ) = string( old : )
        return
      end if
    end do
  end subroutine strip_char
!*****************************************************************************80
! Time routines
!*****************************************************************************80
  subroutine start_timer_sub(self)
    class(time), intent(inout) :: self
    self%saved_time = self%ctime_fnk()
  end subroutine start_timer_sub
!
  character(len=cl) function elapsed_time_fnk(self)
    class(time), intent(in) :: self
    real(rk) :: et
    et = self%ctime_fnk() - self%saved_time
    write(elapsed_time_fnk,'(es9.1e3)') et
  end function elapsed_time_fnk
!
  real(rk) function ctime_fnk()
    integer :: datetime(8)
    call date_and_time(values = datetime)
    ctime_fnk = 86400.0_rk * datetime(3) + 3600.0_rk * datetime(5) &
    &           + 60.0_rk * datetime(6) + datetime(7) + 0.001_rk * datetime(8)
  end function ctime_fnk
!
  function print_time_fnk() result(print_time)
    character(len=8) :: print_time
    integer :: datetime(8)
    call date_and_time(values = datetime)
    write(print_time,'(i2.2,a,i2.2,a,i2.2)') datetime(5), ':',&
    & datetime(6), ':', datetime(7)
  end function print_time_fnk
!
  subroutine write_elapsed_time(self,write_unit,specific)
    class(time), intent(in) :: self
    integer(ik), intent(in) :: write_unit
    character(len=*), optional, intent(in) :: specific
    !
    real(rk) :: elapsed_time
    !
    elapsed_time = self%ctime_fnk() - self%saved_time
    if ( present(specific) ) then
      write(write_unit,'(a,es9.1e3,a)') specific, elapsed_time, ' s.'
    else
      write(write_unit,'(a,es9.1e3,a)') 'Elapsed time', elapsed_time, ' s.'
    end if
    !
  end subroutine write_elapsed_time
!*****************************************************************************80
! Change extension
!*****************************************************************************80
  pure character(len=cl) function change_extension(name,ext)
    character(len=*), intent(in) :: name
    character(len=3), intent(in) :: ext
    !
    integer(ik) :: slen
    !
    slen = len_trim(name)
    write(change_extension,'(a,a)') name(1:slen-3), ext
    !
  end function change_extension
!*****************************************************************************80
end module general_routines




