module determinant
  use types, only: ik, rk
  implicit none
  private
!*****************************************************************************80
  public :: det_sqr_mtx
!*****************************************************************************80
  contains
!*****************************************************************************80
  pure subroutine det_sqr_mtx(a_in,det,esta,emsg)
    !
    use lapack95, only: getrf
    !
    real(rk), intent(in) :: a_in(:,:)
    real(rk), intent(out) :: det
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    integer(ik) :: piv(size(a_in,1))
    real(rk) :: a(size(a_in,1),size(a_in,1))
    !
    if( size(a_in,1) /= size(a_in,2) ) then
      esta = -1
      emsg = 'Det sqr mtx: matrix must be square.'
      return
    end if
    ! LU decomposition
    a = a_in
    call getrf(a,piv,esta)
    ! Check sucess
    det = 0.0_rk
    if ( esta /= 0 ) then
      if( esta > 0 ) then
        emsg = 'Det sqr mtx: u matrix is singular.'
      else
        emsg = 'Det sqr mtx: piviot parameter has an illegal value.'
      end if
      return
    end if
    ! Calculate determinant
    det = 1.0_rk
    do i = 1 , size(piv)
      if ( piv(i) /= i ) then
        det = - det * a(i,i)
      else
        det = det * a(i,i)
      endif
    end do
    !
  end subroutine det_sqr_mtx
!*****************************************************************************80
end module determinant
