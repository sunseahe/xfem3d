module sparse
!*****************************************************************************80
  use types, only: ik, rk, cl, eps, debug
  use general_routines, only: resize_vec
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
  type :: sparse_square_matrix_t
    private
    integer(ik) :: n ! Size
    integer(ik) :: nnz ! Number of nonzero entries
    integer(ik), allocatable :: ai(:), aj(:) ! CSR format entries
    real(rk), allocatable :: ax(:) ! Values
  contains
    procedure :: set => set_ssm
    procedure, private :: sum_dup_rm_zero
    procedure :: get_val_from_indices
  end type sparse_square_matrix_t
!*****************************************************************************80
! Interface for matrix conversion
!*****************************************************************************80
  interface mkl_csrcoo
    subroutine mkl_scsrcoo(job,n,acsr,ajr,air,nnz,acoo,ir,jc,info)
      use iso_fortran_env, only: sp => real32
      integer  :: job(8)
      integer  :: n, nnz, info
      integer  :: ajr(*), air(n+1), ir(*), jc(*)
      real(sp) :: acsr(*), acoo(*)
    end subroutine mkl_scsrcoo
    subroutine mkl_dcsrcoo(job,n,acsr,ajr,air,nnz,acoo,ir,jc,info)
      use iso_fortran_env, only: dp => real64
      integer  :: job(8)
      integer  :: n, nnz, info
      integer  :: ajr(*), air(n+1), ir(*), jc(*)
      real(dp) :: acsr(*), acoo(*)
    end subroutine mkl_dcsrcoo
  end interface mkl_csrcoo
!*****************************************************************************80
contains
!*****************************************************************************80
! Set sparse matrix
!*****************************************************************************80
  subroutine set_ssm(self,n,ar,ac,ax,esta,emsg)
    class(sparse_square_matrix_t), intent(out) :: self
    integer(ik), intent(in) :: n
    integer(ik), intent(in) :: ar(:), ac(:)
    real(rk), intent(in) :: ax(:)
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: job(8)
    ! Check
    if ( debug ) then
      if ( .not.( size(ar,dim=1) == size(ac,dim=1) ) .or. &
      &.not.(size(ar,dim=1) == size(ax,dim=1)) ) then
       esta = -1
       emsg = 'Sparse square matrix: sizes do not match'
       return
      end if
    end if
    self%n = n
    self%nnz = size(ar,dim=1)
    ! Convert to csr format
    job(1) = 2 ! Conversion to csr
    job(2) = 1 ! One-based indexing for the matrix in CSR format is used
    job(3) = 1 ! One-based indexing for the matrix in COO format is used
    job(5) = self%nnz ! Number of entries
    job(6) = 0 ! All arrays are filled
    job(7) = 0 ! Reserved
    job(8) = 0 ! Reserved
    !
    if ( .not. allocated(self%ai) ) then
      allocate(self%ai(self%n+1),stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
    end if
    if ( .not. allocated(self%aj) ) then
      allocate(self%aj(self%nnz),stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
    end if
    if ( .not. allocated(self%ax) ) then
      allocate(self%ax(self%nnz),stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
    end if
    ! Convert
    CALL mkl_csrcoo(job,self%n,self%ax,self%aj,self%ai, &
    & self%nnz,ax,ar,ac,esta)
    if ( esta /= 0 ) then
      emsg = 'Sparse square matrix: conversion to csr format error'
      return
    end if
    ! Sum duplicate entries and remove zero entries
    call self%sum_dup_rm_zero(esta,emsg)
    if ( esta /= 0 ) return
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine set_ssm
!*****************************************************************************80
! Sum duplicate entries and remove zero entries
!*****************************************************************************80
  pure subroutine sum_dup_rm_zero(self,esta,emsg)
    class(sparse_square_matrix_t), intent(inout) :: self
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i, j, jj, k, row_end, r1
    real(rk)  :: x
    ! Sum duplicate
    k        = 1
    row_end  = 1
    do i = 1, self%n
      r1      = row_end
      row_end = self%ai(i+1)
      jj      = r1
      outter: do while ( jj < row_end )
        j  = self%aj(jj)
        x  = self%ax(jj)
        jj = jj + 1
        inner: do while ( jj < row_end  )
          if ( self%aj(jj) == j  ) then
            x  = x + self%ax(jj)
            jj = jj + 1
          else
            exit inner
          end if
        end do inner
        self%aj(k) = j
        self%ax(k)  = x
        k = k + 1
      end do outter
      self%ai(i + 1) = k
    end do
    ! delete zero values
    k       = 1
    row_end = 1
    do i = 1, self%n
      r1      = row_end
      row_end = self%ai(i+1)
      jj      = r1
      do while ( jj < row_end )
        j = self%aj(jj)
        x = self%ax(jj)
        if ( abs(x) >= eps ) then
          self%aj(k)  = j
          self%ax(k) = x
          k = k + 1
        end if
        jj = jj + 1
      end do
      self%ai(i+1) = k
    end do
    k = k - 1
    ! New sizes
    self%nnz = k
    call resize_vec(k,self%aj,esta,emsg); if ( esta /= 0 ) return
    call resize_vec(k,self%ax,esta,emsg); if ( esta /= 0 ) return
    !
  end subroutine sum_dup_rm_zero
!*****************************************************************************80
! Get sparse matrix value from indices
!*****************************************************************************80
  pure subroutine get_val_from_indices(self,i,j,val_out,indx)
    class(sparse_square_matrix_t), intent(in) :: self
    integer(ik), intent(in) :: i, j != indices
    real(rk), intent(out) :: val_out != output value
    integer(ik), intent(out) :: indx != position
    !
    integer(ik) :: ibeg, iend, imid
    !
    val_out = 0.0_rk
    indx = 0
    ibeg = self%ai(i)
    iend = self%ai(i+1)-1
    !
    search1: do
      imid = ( ibeg + iend ) / 2
      if( self%aj(imid) == j ) then
        indx = imid
        exit search1
      end if
      if( ibeg > iend ) exit search1
      if( self%aj(imid) > j ) then
        iend = imid - 1
      else
        ibeg = imid + 1
      end if
    end do search1
    if( indx /= 0 ) val_out = self%ax(indx)
    !
  end subroutine get_val_from_indices


end module sparse