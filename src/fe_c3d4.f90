module fe_c3d4
  use types
  use determinant
  use ll
  use general_routines, only: exclude, real_interval, i2str, r2str
  implicit none
  private
!*****************************************************************************80
  integer(ik), parameter :: dom = 3
  integer(ik), parameter :: nnod = 4
!*****************************************************************************80
! Point
!*****************************************************************************80
  type :: point_3d_t
    real(rk) :: x(dom) = 0.0_rk
  contains
    procedure :: subtract_pnt
    procedure :: dot_pnt
    procedure :: assign_pnt
    generic :: operator(-) => subtract_pnt
    generic :: operator(.dot.) => dot_pnt
    generic :: assignment(=) => assign_pnt
    procedure :: write => write_pnt
  end type point_3d_t
  type, extends(list) :: point_3d_t_ll
    private
  contains
    procedure :: add => add_point_3d_t_ll
    procedure :: fill_array => fill_array_point_3d_t_ll
  end type point_3d_t_ll
!*****************************************************************************80
  type(point_3d_t), parameter :: zero_pnt = point_3d_t([0.0_rk,0.0_rk, &
  &0.0_rk])
  integer(ik), parameter :: ngp_c3d4 = 1
  real(rk), parameter :: gp_c3d4 = 1.0_rk / 4.0_rk, w_c3d4 = 1.0_rk / 6.0_rk
!*****************************************************************************80
! c3d4
!*****************************************************************************80
  type :: c3d4_t
    type(point_3d_t) :: nodes(nnod) = zero_pnt
    integer(ik) :: connectivity(nnod) = 0
  contains
    procedure :: get_ver_coo => get_ver_coo_tet
    procedure, nopass :: n_mtx_sca => n_mtx_sca_tet
    procedure :: g_coo => g_coo_tet
    procedure :: volume => calculate_tetrahedral_volume
    procedure :: write => write_tet
  end type c3d4_t
  type, extends(list) :: c3d4_t_ll
    private
  contains
    procedure :: add => add_c3d4_t_ll
    procedure :: fill_array => fill_array_c3d4_t_ll
  end type c3d4_t_ll
!*****************************************************************************80
  public :: dom, nnod, zero_pnt
  public :: point_3d_t, point_3d_t_ll, c3d4_t, c3d4_t_ll, &
  & calc_tet_vol, deallocate_tet, calc_tet_strech
!*****************************************************************************80
  contains
!*****************************************************************************80
  pure subroutine deallocate_tet(tets, n)
    type(c3d4_t), intent(inout), allocatable :: tets(:)
    integer(ik), intent(in) :: n
    !
    type(c3d4_t), allocatable :: alt_tets(:)
    integer(ik) :: n_len, k, i
    !
    n_len=size(tets)
    allocate(alt_tets(n_len-1))
    k=0
    do i=1,n_len
      if (i.ne.n) then
        k=k+1
        alt_tets(k)=tets(i)
      end if
    end do
    deallocate(tets)
    if (n_len.gt.1) allocate(tets(n_len-1))
    !gfortran bug
    if(allocated(tets)) then
      do i = 1, size(tets)
        tets(i) = alt_tets(i)
      end do
    end if
  end subroutine deallocate_tet
!*****************************************************************************80
!Racunanje volumna tetraedra
!*****************************************************************************80
  pure subroutine calc_tet_vol(tet,volume,istat,emsg)
    type(c3d4_t), intent(in) :: tet
    real(rk), intent(out) :: volume
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i, j
    real(rk) :: det
    real(rk) :: det_mtx(dom,dom)
    type(point_3d_t), parameter :: e(dom) = [ &
    & point_3d_t([1.0_rk,0.0_rk,0.0_rk]), &
    & point_3d_t([0.0_rk,1.0_rk,0.0_rk]), &
    & point_3d_t([0.0_rk,0.0_rk,1.0_rk]) ]
    !
    volume = 0.0_rk
    do i = 1, dom
      do j = 1, dom
        det_mtx(i,j) = (tet%nodes(j+1) - tet%nodes(1)) .dot. e(i)
      end do
    end do
    ! Calculate determinant
    call det_sqr_mtx(det_mtx,det,istat,emsg)
    if( istat /= 0 ) return
    ! Volume
    volume = 0.16666667_rk * abs(det)
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine calc_tet_vol
!*****************************************************************************80
! Tetrahedral volume - type bound procedure
!*****************************************************************************80
  pure subroutine calculate_tetrahedral_volume(self,volume,istat,emsg)
    class(c3d4_t), intent(in) :: self
    real(rk), intent(out) :: volume
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i, j
    real(rk) :: det
    real(rk) :: det_mtx(dom,dom)
    type(point_3d_t), parameter :: e(dom) = [ &
    & point_3d_t([1.0_rk,0.0_rk,0.0_rk]), &
    & point_3d_t([0.0_rk,1.0_rk,0.0_rk]), &
    & point_3d_t([0.0_rk,0.0_rk,1.0_rk]) ]
    !
    volume = 0.0_rk
    do i = 1, dom
      do j = 1, dom
        det_mtx(i,j) = (self%nodes(j+1) - self%nodes(1)) .dot. e(i)
      end do
    end do
    ! Calculate determinant
    call det_sqr_mtx(det_mtx,det,istat,emsg)
    if( istat /= 0 ) return
    ! Volume
    volume = 0.16666667_rk * abs(det)
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine calculate_tetrahedral_volume
!*****************************************************************************80
!Racunanje deformiranosti tetraedra
  pure subroutine calc_tet_strech(tet,volume,strech,istat,emsg)
    type(c3d4_t), intent(in) :: tet
    real(rk), intent(in) :: volume
    real(rk), intent(out) :: strech
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: i, j
    integer(ik) :: pnts(6)
    real(rk) :: len_a(3), s, len_max, area
    type(point_3d_t) :: a(3)
    !
    area=0.0_rk
    pnts(1)=1; pnts(2)=2; pnts(3)=3; pnts(4)=4; pnts(5)=1; pnts(6)=2;
    do i = 1, 4
       len_a=0.0_rk; s=0.0_rk; len_max=0.0_rk
       a(1)=tet%nodes(pnts(i))-tet%nodes(pnts(i+1))
       a(2)=tet%nodes(pnts(i+1))-tet%nodes(pnts(i+2))
       a(3)=tet%nodes(pnts(i+2))-tet%nodes(pnts(i))
       do j=1,3
         len_a(j)=len_pnt(a(j))
         if (len_a(j).gt.len_max) len_max=len_a(j)
         s=s+len_a(j)/2.0_rk
       end do
     area=area+sqrt(s*(s-len_a(1))*(s-len_a(2))*(s-len_a(3)))
    end do
    strech=14.6969385_rk*volume/area/len_max
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine calc_tet_strech
!*****************************************************************************80
  pure function subtract_pnt(self,other) result(res)
    class(point_3d_t), intent(in) :: self
    type(point_3d_t), intent(in) :: other
    type(point_3d_t) :: res
    real(rk) :: sub(dom)
    sub = self%x - other%x
    res = point_3d_t(sub)
  end function subtract_pnt
!*****************************************************************************80
  pure function dot_pnt(self,other) result(res)
    class(point_3d_t), intent(in) :: self
    type(point_3d_t), intent(in) :: other
    real(rk) :: res
    res = sum(self%x * other%x)
  end function dot_pnt
!*****************************************************************************80
  pure function len_pnt(self) result(res)
    class(point_3d_t), intent(in) :: self
    real(rk) :: res
    res = sqrt(sum(self%x * self%x))
  end function len_pnt
!*****************************************************************************80
  pure subroutine assign_pnt(self,other)
    class(point_3d_t), intent(inout) :: self
    type(point_3d_t), intent(in) :: other
    self%x = other%x
  end subroutine assign_pnt
!*****************************************************************************80
  subroutine write_pnt(self)
    class(point_3d_t), intent(in) :: self
    character(len=*), parameter :: pw = '(a,3(' // &
    &es//',:,","))'
    write(stdout,pw,advance='no') 'Point coordinates: (', self%x
    write(stdout,'(a)') ' )'
  end subroutine write_pnt
!*****************************************************************************80
  pure subroutine add_point_3d_t_ll(self,arg,istat,emsg)
    class(point_3d_t_ll), intent(inout) :: self
    type(point_3d_t), intent(in) :: arg
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,istat,emsg)
    !
  end subroutine add_point_3d_t_ll
!*****************************************************************************80
  pure subroutine fill_array_point_3d_t_ll(self,arr,istat,emsg)
    class(point_3d_t_ll), intent(inout) :: self
    type(point_3d_t), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      istat = -1
      emsg = 'Fill array point 3d: list is empty'
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
      type is( point_3d_t )
        arr(i) = curr
      class default
        istat = -1
        emsg = 'Fill array point 3d: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
  end subroutine fill_array_point_3d_t_ll
!*****************************************************************************80
! Get vert coordinates
!*****************************************************************************80
  pure function get_ver_coo_tet(self,n) result(coo)
    class(c3d4_t), intent(in) :: self
    integer(ik), intent(in) :: n
    real(rk) :: coo(dom)
    !
    if (n > nnod ) then
      coo = 0.0_rk
      return
    end if
    coo = self%nodes(n)%x
    !
  end function get_ver_coo_tet
!*****************************************************************************80
  subroutine write_tet(self)
    class(c3d4_t), intent(in) :: self
    !
    integer(ik) :: i
    !
    write(stdout,'(a)') 'C3d4 coordinates:'
    do i = 1, nnod
      write(stdout,'(i0,a)',advance='no') i, '. '
      call self%nodes(i)%write()
    end do
    !
  end subroutine write_tet
!*****************************************************************************80
! N matrix c3d4
!*****************************************************************************80
  pure subroutine n_mtx_sca_tet(gp_num,xi_coo,n_mtx_sca,istat,emsg)
    integer(ik), optional, intent(in) :: gp_num
    real(rk), optional, intent(in) :: xi_coo(:)
    real(rk), intent(out) :: n_mtx_sca(:)
    integer(ik), intent(out) :: istat
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    real(rk) :: xi(dom)
    logical(lk) :: gp_num_p, xi_coo_p
    !
    gp_num_p = present(gp_num); xi_coo_p = present(xi_coo)
    if ( debug ) then
      if ( .not. exclude(gp_num_p,xi_coo_p) ) then
        istat = -1
        emsg ='N matrix scalar tet: Gauss points and xi coordinates must&
        & be excluded'
        return
      end if
      if ( .not. size(n_mtx_sca)==nnod ) then
        istat = -1
        emsg ='N matrix scalar tet: N matrix size incorrect'
        return
      end if
    end if
    ! Set xi value
    if( gp_num_p ) then
      if ( debug ) then
        if( gp_num > 1 ) then
          istat = -1
          emsg ='N matrix scalar tet: only one Gauss point'
          return
        end if
      end if
      xi = gp_c3d4
    else if( xi_coo_p ) then
      if ( debug ) then
        if(.not.size(xi_coo) == dom ) then
          istat = -1
          emsg ='N matrix scalar tet: xi coo matrix size incorrect'
          return
        end if
        do i = 1, size(xi_coo)
          if ( .not. real_interval(0.0e0_rk,1.0e0_rk,xi_coo(i))) then
            istat = -1
            emsg ='N matrix scalar tet: xi coo vector component number '&
            & // trim(i2str(i)) // ' out of bounds'
            return
          end if
        end do
      end if
      xi = xi_coo
    end if
    n_mtx_sca = [ 1.0e0_rk - xi(1) - xi(2) - xi(3), xi(1), xi(2), xi(3) ]
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine n_mtx_sca_tet
!*****************************************************************************80
! Global coordinates c3d4
!*****************************************************************************80
  pure subroutine g_coo_tet(self,gp_num,xi_coo_pnt,g_coo_pnt,istat,emsg)
    use blas95, only: gemv
    class(c3d4_t), intent(in) :: self
    integer(ik), optional, intent(in) :: gp_num
    type(point_3d_t), optional, intent(in) :: xi_coo_pnt
    type(point_3d_t), intent(out) :: g_coo_pnt
    integer(ik), intent(out) :: istat
    character(len=cl), intent(out) :: emsg
    !
    integer(ik) :: i
    real(rk) :: xi(dom), n_mtx_sca(nnod)
    real(rk) :: el_coo(nnod,dom), g_coo(dom), xi_coo(dom)
    logical(lk) :: gp_num_p, xi_coo_p
    !
    gp_num_p = present(gp_num); xi_coo_p = present(xi_coo_pnt)
    if ( debug ) then
      if ( .not. exclude(gp_num_p,xi_coo_p) ) then
        istat = -1
        emsg ='Global coordinates tet: Gauss points and xi coordinates must&
        & be excluded'
        return
      end if
    end if
    ! Set xi value
    if( gp_num_p ) then
      if ( debug ) then
        if( gp_num > 1 ) then
          istat = -1
          emsg ='Global coordinates tet: only one Gauss point'
          return
        end if
      end if
      xi = gp_c3d4
    else if( xi_coo_p ) then
      xi_coo = xi_coo_pnt%x
      if ( debug ) then
        do i = 1, size(xi_coo)
          if ( .not. real_interval(0.0e0_rk,1.0e0_rk,xi_coo(i))) then
            istat = -1
            emsg ='Global coordinates tet: xi coo vector component number ' &
            & // trim(i2str(i)) // ' out of bounds, with value >' //  &
            & trim(r2str(xi_coo(i))) // '<'
            return
          end if
        end do
      end if
      xi = xi_coo
    end if
    ! Calculate
    do i = 1, nnod
      el_coo(i,:) = self%get_ver_coo(i)
    end do
    call self%n_mtx_sca(xi_coo=xi,n_mtx_sca=n_mtx_sca,istat=istat,&
    &emsg=emsg)
    call gemv(transpose(el_coo),n_mtx_sca,g_coo)
    g_coo_pnt = point_3d_t(g_coo)
    ! Sucess
    istat = 0
    emsg = ''
    !
  end subroutine g_coo_tet
!*****************************************************************************80
  pure subroutine add_c3d4_t_ll(self,arg,istat,emsg)
    class(c3d4_t_ll), intent(inout) :: self
    type(c3d4_t), intent(in) :: arg
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    call self%add_list(arg,istat,emsg)
    !
  end subroutine add_c3d4_t_ll
!*****************************************************************************80
  pure subroutine fill_array_c3d4_t_ll(self,arr,istat,emsg)
    class(c3d4_t_ll), intent(inout) :: self
    type(c3d4_t), allocatable, intent(out) :: arr(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer :: i
    class(*), allocatable :: curr
    !
    if ( self%is_empty() ) then
      istat = -1
      emsg = 'Fill array tet dat: list is empty'
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
      type is( c3d4_t )
        arr(i) = curr
      class default
        istat = -1
        emsg = 'Fill array tet dat: wrong linked list item.'
        return
      end select
      call self%set_next()
    end do
  end subroutine fill_array_c3d4_t_ll
!*****************************************************************************80
end module fe_c3d4

