include 'mkl_pardiso.f90'
module sparse
!*****************************************************************************80
  use blas95, only: dot
  use iso_fortran_env, only: sp => real32, dp => real64, int64
  use types, only: ik, rk, lk, cl, eps, debug, es, nl
  use general_routines, only: resize_vec, time, pnorm
  use mkl_pardiso, only: mkl_pardiso_handle, pardiso
!*****************************************************************************80
  implicit none
  private
!*****************************************************************************80
! Sparse square matrix
!*****************************************************************************80
  type :: sparse_square_matrix_t
    integer(ik) :: n ! Size
    integer(ik) :: nnz ! Number of nonzero entries
    integer(ik), allocatable :: ai(:), aj(:) ! CSR format entries
    real(rk), allocatable :: ax(:) ! Values
  contains
    procedure :: set => set_ssm
    procedure, private :: sum_dup_rm_zero
    procedure :: get_val_from_indices
    procedure :: delete => delete_ssm
    procedure :: write_matrix
  end type sparse_square_matrix_t
!*****************************************************************************80
! Interface for matrix conversion
!*****************************************************************************80
  interface mkl_csrcoo
    subroutine mkl_scsrcoo(job,n,acsr,ajr,air,nnz,acoo,ir,jc,info)
      import :: ik, sp
      integer(ik) :: job(8)
      integer(ik) :: n, nnz, info
      integer(ik) :: ajr(*), air(n+1), ir(*), jc(*)
      real(sp) :: acsr(*), acoo(*)
    end subroutine mkl_scsrcoo
    subroutine mkl_dcsrcoo(job,n,acsr,ajr,air,nnz,acoo,ir,jc,info)
      import :: ik, dp
      integer(ik) :: job(8)
      integer(ik) :: n, nnz, info
      integer(ik) :: ajr(*), air(n+1), ir(*), jc(*)
      real(dp) :: acsr(*), acoo(*)
    end subroutine mkl_dcsrcoo
  end interface mkl_csrcoo
!*****************************************************************************80
! Sparse linear system
!*****************************************************************************80
  type :: sparse_linear_system_t
    private
    logical(lk) :: configured_direct = .false.
    integer(ik) :: iparm(64) = 0
    type(mkl_pardiso_handle), pointer :: pt(:) => null()
  contains
    procedure :: solve_dir => solve_direct_sls
    procedure, nopass :: solve_iter => bicgstab_jac_prec
  end type sparse_linear_system_t
  ! Solution parameters
  integer(ik), parameter :: &
  &  msglvl = 0, & ! No printing of statistical information
  &  mtype  = 2, & ! Symmetric positive definite
  &  maxfct = 1, &
  &  mnum   = 1
  ! Paradiso storing routines
  integer(ik), parameter :: phs = 1, phr = 2, phd = 3
!*****************************************************************************80
  interface mkl_csrsymv
    subroutine mkl_scsrsymv(uplo,m,a,ia,ja,x,y)
      import :: ik, sp
      character(len=1) :: uplo
      integer(ik) :: m
      real(sp) :: a(*)
      integer(ik) :: ia(*), ja(*)
      real(sp) :: x(*), y(*)
    end subroutine mkl_scsrsymv
    subroutine mkl_dcsrsymv(uplo,m,a,ia,ja,x,y)
      import :: ik, dp
      character(len=1) :: uplo
      integer(ik) :: m
      real(dp) :: a(*)
      integer(ik) :: ia(*), ja(*)
      real(dp) :: x(*), y(*)
    end subroutine mkl_dcsrsymv
  end interface mkl_csrsymv
!*****************************************************************************80
  interface mkl_csrsv
    subroutine mkl_scsrsv( transa, m, alpha, matdescra, &
    &val, indx, pntrb, pntre,  x, y)
      import :: ik, sp
      character(len=1) :: transa
      character(len=1) :: matdescra(*)
      integer(ik) :: m
      real(sp) :: alpha
      integer(ik) :: indx(*), pntrb(*), pntre(*)
      real(sp) :: val(*)
      real(sp) :: y(*), x(*)
    end subroutine mkl_scsrsv
    subroutine mkl_dcsrsv( transa, m, alpha, matdescra, &
    &val, indx, pntrb, pntre,  x, y)
      import :: ik, dp
      character(len=1) :: transa
      character(len=1) :: matdescra(*)
      integer(ik) :: m
      real(dp) :: alpha
      integer(ik) :: indx(*), pntrb(*), pntre(*)
      real(dp) :: val(*)
      real(dp) :: y(*), x(*)
    end subroutine mkl_dcsrsv
  end interface mkl_csrsv
!*****************************************************************************80
  public :: sparse_square_matrix_t, sparse_linear_system_t
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
    call mkl_csrcoo(job,self%n,self%ax,self%aj,self%ai, &
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
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine sum_dup_rm_zero
!*****************************************************************************80
! Get sparse matrix value from indices
!*****************************************************************************80
  pure subroutine get_val_from_indices(self,i,j,val_out,indx,esta,emsg)
    class(sparse_square_matrix_t), intent(in) :: self
    integer(ik), intent(in) :: i, j != indices
    real(rk), intent(out) :: val_out != output value
    integer(ik), intent(out) :: indx != position
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: ibeg, iend, imid
    !
    if ( debug ) then
      if ( .not. allocated(self%ai) ) then
        esta = -1
        emsg = 'Sparse square matrix: not allocated'
        return
      end if
    end if
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
    ! Sucess
    esta = 0
    emsg = ''
    !
  end subroutine get_val_from_indices
!*****************************************************************************80
! Write sparse matrix
!*****************************************************************************80
  subroutine write_matrix(self,write_unit,esta,emsg)
    class(sparse_square_matrix_t), intent(in) :: self
    integer(ik), intent(in) :: write_unit
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    integer(ik) :: indx
    !
    if ( debug ) then
      if ( .not. allocated(self%ai) ) then
        esta = -1
        emsg = 'Sparse square matrix: not allocated'
        return
      end if
    end if
    ! First write sizes
    write(write_unit,'(i0,1x,i0)') self%n, self%nnz
    ! Write ai
    indx = 1
    row1: do
      do i = 1, 8
        if ( indx > self%n + 1 ) exit row1
        write(write_unit,'(i0,1x)',advance='no') self%ai(indx)
        indx = indx + 1
      end do
      write(write_unit,'(a)') ''
    end do row1
    write(write_unit,'(a)') ''
    ! Write aj
    indx = 1
    row2: do
      do i = 1, 8
        if ( indx > self%nnz ) exit row2
        write(write_unit,'(i0,1x)',advance='no') self%aj(indx)
        indx = indx + 1
      end do
      write(write_unit,'(a)') ''
    end do row2
    write(write_unit,'(a)') ''
    ! Write ax
    indx = 1
    row3: do
      do i = 1, 8
        if ( indx > self%nnz ) exit row3
        write(write_unit,'('//es//',1x)',advance='no') self%ax(indx)
        indx = indx + 1
      end do
      write(write_unit,'(a)') ''
    end do row3
    write(write_unit,'(a)') ''
    !
  end subroutine write_matrix
!*****************************************************************************80
! Delete routine
!*****************************************************************************80
  subroutine delete_ssm(self,esta,emsg)
    class(sparse_square_matrix_t), intent(inout) :: self
    integer(ik), intent(out) :: esta
    character(len=cl), intent(out) :: emsg
    !
    self%n = 0
    self%nnz = 0
    deallocate(self%ai,self%aj,self%ax,stat=esta,errmsg=emsg)
    !
  end subroutine delete_ssm
!*****************************************************************************80
! Direct solve routine
!*****************************************************************************80
  subroutine solve_direct_sls(self,job,a,b,x,mem_used,esta,emsg)
    class(sparse_linear_system_t), intent(inout) :: self
    integer(ik), intent(in) :: job
    type(sparse_square_matrix_t), intent(inout) :: a
    real(rk), optional, intent(inout) :: b(:)
    real(rk), optional, intent(out) :: x(:)
    integer(int64), optional, intent(out) :: mem_used
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik), parameter :: one = 1
    integer(ik) :: phase, idum(1)
    real(rk) :: rdum(1)
    ! Set solution parameters
    if ( .not. self%configured_direct ) then
      !
      allocate(self%pt(64),stat=esta,errmsg=emsg)
      if ( esta /= 0 ) return
      self%pt%dummy = 0
      !
      self%iparm = 0
      self%iparm(1) = 1 ! no solver default
      self%iparm(2) = 2 ! fill-in reordering from METIS
      self%iparm(4) = 0 ! no iterative-direct algorithm
      self%iparm(5) = 0 ! no user fill-in reducing permutation
      self%iparm(6) = 0 ! =0 solution on the first n compoments of x
      self%iparm(8) = 8 ! numbers of iterative refinement steps
      self%iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      self%iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      self%iparm(27) = 1 ! check input values
#ifndef _DOUBLE
      self%iparm(28) = 1 ! Single precision
#else
      self%iparm(28) = 0 ! Double precision
#endif
      !
      self%configured_direct = .true.
    end if
    idum = 1
    rdum = 1.0_rk
    ! Select job
    select case ( job )
!*****************************************************************************80
    case ( 1 ) ! Analyze and factorize
      ! Factorize
      phase = 12
      call pardiso(self%pt,maxfct,mnum,mtype,phase,a%n,a%ax,a%ai,a%aj,idum, &
      & one,self%iparm,msglvl,rdum,rdum,esta)
      if ( esta /= 0 ) then
        call pardiso_err(esta,emsg)
        emsg = trim(emsg) // ' - in job analyze and factorize'
        return
      end if
      if ( present(mem_used) ) mem_used = int(self%iparm(15)*1000,&
      & kind=int64)
!*****************************************************************************80
    case ( 2 ) ! Solve
      if ( .not. present(x) .or. .not. present(b) ) then
        call pardiso_err(-20,emsg)
        emsg = trim(emsg) // ' - in job solve'
        return
      end if
      ! Solve
      phase = 33
      CALL pardiso(self%pt,maxfct,mnum,mtype,phase,a%n,a%ax,a%ai,a%aj, &
      &idum,one,self%iparm,msglvl,b,x,esta)
      if ( esta /= 0 ) then
        call pardiso_err(esta,emsg)
        emsg = trim(emsg) // ' - in job solve'
        return
      end if
!*****************************************************************************80
    case ( 3 ) ! Clean
      phase = -1
      CALL pardiso(self%pt,maxfct,mnum,mtype,phase,one,rdum,idum,idum, &
      & idum,one,self%iparm,msglvl,rdum,rdum,esta)
      if ( esta /= 0 ) then
        call pardiso_err(esta,emsg)
        emsg = trim(emsg) // ' - in job clean'
        return
      end if
      deallocate(self%pt,stat=esta,errmsg=emsg); if ( esta /= 0 ) return
      self%configured_direct = .false.
    case default
      esta = -1
      emsg = 'Solve sls: unknown job'
    end select
    ! Sucess
    esta = 0
    emsg = ''
    !
    contains
!*****************************************************************************80
! Error handling
!*****************************************************************************80
    subroutine pardiso_err(esta,emsg)
      integer(ik), intent(in) :: esta
      character(len=*), intent(out) :: emsg
      !
      character(len=cl) :: adit_emsg
      !
      select case ( esta )
      case ( -1 )
        adit_emsg = 'Input inconsistent'
      case ( -2 )
        adit_emsg = 'Not enough memory'
      case ( -3 )
        adit_emsg = 'Reordering problem'
      case ( -4 )
        adit_emsg = 'Zero pivot, numerical factorization or iterative &
        &refinement problem'
      case ( -5 )
        adit_emsg = 'Unclassified (internal) error'
      case ( -6 )
        adit_emsg = 'Reordering failed (matrix types 11 and 13 only)'
      case ( -7 )
        adit_emsg = 'Diagonal matrix is singular'
      case ( -8 )
        adit_emsg = '32-bit integer overflow problem'
      case ( -9 )
        adit_emsg = 'Not enough memory for OOC'
      case ( -10 )
        adit_emsg = 'Problems with opening OOC temporary files'
      case ( -11 )
        adit_emsg = 'Read/write problems with the OOC data file'
      case ( -20 )
        adit_emsg = 'Rhs vector or solution vector not present'
      case default
        adit_emsg = 'Unknown error'
      end select
      emsg = 'Solve sls: pardiso error - '// trim(adit_emsg)
      !
    end subroutine pardiso_err
  end subroutine solve_direct_sls
!*****************************************************************************80
! Iterative solver
!*****************************************************************************80
  subroutine bicgstab_jac_prec(a,b,x,niter,tol,info,esta,emsg)
    type(sparse_square_matrix_t), intent(in) :: a
    real(rk), intent(in) :: b(:)
    real(rk), intent(out) :: x(:)
    integer(ik), intent(in) :: niter
    real(rk), intent(in) :: tol
    character(len=*), optional, intent(out) :: info
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i
    real(rk), parameter :: one = 1.0_rk
    real(rk) :: rho, rho_old, alpha, beta, omega
    real(rk), allocatable :: r(:), r_tilde(:), p(:), p_hat(:),&
    & v(:), s(:), s_hat(:), tmp(:)
    character(len=1), parameter :: transa = 'N'
    character(len=1), parameter :: matdescra(3) = ['D','U','N']
    character(len=1), parameter :: uplo = 'U'
    type(time) :: t
    !
    call t%start_timer()
    !
    allocate(r(a%n),r_tilde(a%n),p(a%n),p_hat(a%n),v(a%n),s(a%n),&
    &s_hat(a%n),tmp(a%n),stat=esta,errmsg=emsg); if ( esta /= 0 ) return
    !
    x = 0.0_rk
    call mkl_csrsymv(uplo,a%n,a%ax,a%ai,a%aj,x,tmp)
    r = b - tmp
    r_tilde = r
    omega = 1.0_rk
    !
    i = 1
    loop: do
      if ( i > niter ) then
        esta = -1
        emsg = 'biCGSTAB : not converged in the given number&
        & of iterations'
        exit loop
      end if
      rho = dot(r_tilde,r)
      if ( abs(rho) <= tiny(1.0_rk) ) then
        esta = -1
        emsg = 'biCGSTAB : breakdown - rho'
        exit loop
      end if
      !
      if ( i > 1 ) then
        beta = (rho / rho_old ) * ( alpha / omega )
        p = r + beta * ( p - omega * v )
      else
        p = r
      end if
      !
      call mkl_csrsv(transa,a%n,one,matdescra,a%ax,a%aj,a%ai(1:a%n),&
      &a%ai(2:a%n+1),p,p_hat)
      call mkl_csrsymv(uplo,a%n,a%ax,a%ai,a%aj,p_hat,v)
      alpha = rho / dot(r_tilde,v)
      s = r - alpha * v
      if ( pnorm(2,s) <= tol ) then
        x = x + alpha * p_hat
        esta = 0
        emsg = ''
        exit loop
      end if
      !
      call mkl_csrsv(transa,a%n,one,matdescra,a%ax,a%aj,a%ai(1:a%n),&
      &a%ai(2:a%n+1),s,s_hat)
      call mkl_csrsymv(uplo,a%n,a%ax,a%ai,a%aj,s_hat,tmp)
      omega = dot(tmp,s) / dot(tmp,tmp)
      if ( abs(omega) <= tiny(1.0_rk) ) then
        esta = -1
        emsg = 'biCGSTAB : breakdown - omega'
        exit loop
      end if
      !
      x = x + alpha * p_hat + omega * s_hat
      r = s - omega * tmp
      if ( pnorm(2,r) <= tol ) then
        esta = 0
        emsg = ''
        exit loop
      end if
      rho_old = rho
      i = i + 1
    end do loop
    !
    call mkl_csrsymv(uplo,a%n,a%ax,a%ai,a%aj,x,tmp)
    r = b - tmp
    !
    if ( present(info) ) then
      write(info,'(a,a,a,i0,a,a,'//es//',a,a,a,a)') &
      & 'biCGSTAB:', nl, &
      & ' - number of iterations: ', i, nl, &
      & ' - tolerance of the solution: ', pnorm(2,r), nl, &
      & ' - elapsed time: ', trim(t%elapsed_time()), ' s'
    end if
    !
  end subroutine bicgstab_jac_prec

!*****************************************************************************80
end module sparse
