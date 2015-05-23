module lsf_test_functions
  use types
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: dim3 = 3
!*****************************************************************************80
  public :: s_1, s_2, s_3, s_4, s_5
!*****************************************************************************80
  contains
!*****************************************************************************80
  pure real(rk) function s_1(x)
    real(rk), intent(in) :: x(:)
    !
    s_1 = 0.0_rk
    if ( debug ) then
      if ( size(x) /= dim3 ) return
    end if
    s_1 = - ( -1._rk + 6.249999999999999_rk*(-0.5_rk + x(1))**2 + &
    & 24.999999999999996_rk*(-0.5_rk + x(2))**2 + &
    & 99.99999999999999_rk*(-0.5_rk + x(3))**2 )
    !
  end function s_1
!*****************************************************************************80
  pure real(rk) function s_2(x)
    real(rk), intent(in) :: x(:)
    !
    s_2 = 0.0_rk
    if ( debug ) then
      if ( size(x) /= dim3 ) return
    end if
    s_2 = - ( -0.01_rk + (-0.3_rk + sqrt((-0.5_rk + x(1))**2 + &
    & (-0.5_rk + x(2))**2))**2 + (-0.5_rk + x(3))**2 )
    !
  end function s_2
!*****************************************************************************80
  pure real(rk) function s_3(x)
    real(rk), intent(in) :: x(:)
    !
    s_3 = 0.0_rk
    if ( debug ) then
      if ( size(x) /= dim3 ) return
    end if
    s_3 = - ( -0.05_rk + ((-1.6667_rk + 3.3333333333333335_rk*x(1))**2*(1.2_rk - &
    & (-1.6667_rk + 3.3333333333333335_rk*x(1))**2) - &
    & (-1.6667_rk + 3.3333333333333335_rk*x(2))**2)**2 + &
    & (-1.6667_rk + 3.3333333333333335_rk*x(3))**2 )
    !
  end function s_3
!*****************************************************************************80
  pure real(rk) function s_4(x)
    real(rk), intent(in) :: x(:)
    !
    s_4 = 0.0_rk
    if ( debug ) then
      if ( size(x) /= dim3 ) return
    end if
    associate( xm => x(:) / 0.25_rk - 2.0_rk )
    s_4 = - ( ((-1.0_rk + xm(1)**2 + xm(2)**2)**2 + xm(3)**2)* &
    & (xm(2)**2 + (-1.0_rk + xm(1)**2 + xm(3)**2)**2)*(xm(1)**2 + &
    &(-1.0_rk + xm(2)**2 + xm(3)**2)**2) - &
    & 0.005625_rk*(1.0_rk + 3.0_rk*(xm(1)**2 + xm(2)**2 + xm(3)**2)) )
    end associate
    !
  end function s_4
!*****************************************************************************80
   pure real(rk) function s_5(x)
    real(rk), intent(in) :: x(:)
    !
    s_5 = 0.0_rk
    if ( debug ) then
      if ( size(x) /= dim3 ) return
    end if
    s_5 = 2.5e-1_rk - sqrt( (-0.5_rk + x(1))**2 + (-0.5_rk + x(2))**2 + &
    & (-0.5_rk + x(3))**2 )
    !
  end function s_5
!*****************************************************************************80
end module lsf_test_functions
