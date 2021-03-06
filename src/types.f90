module types
  use iso_fortran_env, only: &
  & ik => int32, &
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
  public :: ik, rk, cl, lk, es, stdout, stdin, stderr, debug
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
  subroutine print_error(istat,emsg)
    integer(ik), intent(in) :: istat
    character(len=*), intent(in) :: emsg
    !
    write(stderr,'(a,i0)') 'Error number: ', istat
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
  pure subroutine str2i(str,iout,istat,emsg)
    character(len=*), intent(in) :: str
    integer(ik), intent(out) :: iout
    integer, intent(out) :: istat
    character(len=cl), intent(out) :: emsg
    !
    character(len=cl) :: tmp_str, f
    integer :: n
    !
    tmp_str = adjustl(str)
    n = len_trim(tmp_str)
    write(unit=f,fmt='("(i",i0,")")',iostat=istat, &
    &iomsg=emsg) n
    iout = 0_ik
    read(unit=tmp_str,fmt=f,iostat=istat,iomsg=emsg) iout
    !
  end subroutine str2i
  pure character(len=cl) function i2str(i)
    integer(ik), intent(in) :: i
    write(unit=i2str,fmt='(i0)') i
  end function i2str
!*****************************************************************************80
! String to real
!*****************************************************************************80
  pure subroutine str2r(str,rout,istat,emsg)
    character(len=*), intent(in) :: str
    real(rk), intent(out) :: rout
    integer(ik), intent(out) :: istat
    character(len=cl), intent(out) :: emsg
    !
    character(len=cl) :: tmp_str, f
    integer :: dot, n
    !
    tmp_str = adjustl(str)
    n = len_trim(tmp_str)
    dot = index(tmp_str, '.')
    if( dot < 1 ) then
      write(unit=f,fmt='("(f",i0,".0)")',iostat=istat,iomsg=emsg) n
    else
      write(unit=f,fmt='("(f",i0,".",i0,")")',iostat=istat, &
      &iomsg=emsg) n, n - dot
    end if
    rout = 0.0e0_rk
    read(unit=tmp_str,fmt=f,iostat=istat,iomsg=emsg) rout
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
  pure subroutine f2c_char(f_char_in,c_char_out,istat,emsg)
    use iso_c_binding, only: c_char, c_null_char
    character(len=*), intent(in) :: f_char_in
    character(len=1,kind=c_char), &
    & allocatable, intent(out) :: c_char_out(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer :: n, i
    !
    n = len_trim(f_char_in)
    allocate(c_char_out(n+1),stat=istat,errmsg=emsg)
    if( istat /= 0 ) return
    do i = 1, n
      c_char_out(i) = f_char_in(i:i)
    end do
    c_char_out(n+1) = c_null_char
    !
  end subroutine f2c_char
!*****************************************************************************80
end module general_routines
