module volume_integral
  use types
  use tet_volume, only: dom, nver, point_3d_dat, tet_dat
  use read_input, only: nodes, finite_elements
  use xtet, only: set_sub_xtets
  use lsf_test_functions
  use write_odb, only: start_odb_api, finish_odb_api, write_model_data,  &
  & create_step, create_frame, write_scalar_field
  implicit none
!*****************************************************************************80
  public :: test_functions
!*****************************************************************************80
contains
!*****************************************************************************80
! Write functions to odb
!*****************************************************************************80
  subroutine write_test_functions(istat,emsg)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e, i
    integer(ik) :: active_node
    real(rk) :: x(dom)
    real(rk), allocatable :: lsf_field(:)
    character(len=cl), parameter :: step = 'Test step', field_name = 'LSF', &
    & field_description = 'Level set function values'
    !
    allocate(lsf_field(size(nodes)),stat=istat,errmsg=emsg)
    lsf_field = 0.0_rk
    ! Start odb api
    call start_odb_api()
    ! Model data
    call write_model_data(1,istat,emsg)
    if( istat /= 0 ) return
    ! Step
    call create_step(step,'LSF values',istat,emsg)
    ! First function
    call create_frame(step,'Ellipsoid',1,1.0_rk,&
    &istat,emsg)
    do e = 1, size(finite_elements)
      do i = 1, nver
        x = finite_elements(e)%get_ver_coo(i)
        active_node = finite_elements(e)%connectivity(i)
        lsf_field(active_node) = s_1(x)
      end do
    end do
    call write_scalar_field(step,1,field_name,&
    & field_description,lsf_field,.true.,istat,emsg)
    ! Second function
    call create_frame(step,'Torus',2,2.0_rk,&
    &istat,emsg)
    do e = 1, size(finite_elements)
      do i = 1, nver
        x = finite_elements(e)%get_ver_coo(i)
        active_node = finite_elements(e)%connectivity(i)
        lsf_field(active_node) = s_2(x)
      end do
    end do
    call write_scalar_field(step,2,field_name,&
    & field_description,lsf_field,.true.,istat,emsg)
    !
    ! Finish odb api
    call finish_odb_api()
    !
  end subroutine write_test_functions
!*****************************************************************************80
! Calculation of volume
!*****************************************************************************80
  subroutine calc_vol(test_fun,volume,istat,emsg)
    !
    interface
      pure real(rk) function test_fun(x)
        import :: rk
        real(rk), intent(in) :: x(:)
      end function test_fun
    end interface
    real(rk), intent(out) :: volume
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e, i, t
    real(rk) :: x(dom) ! Coordinates
    real(rk) :: lsf(nver) ! Level set values
    real(rk) :: volume_tet ! Volume of a single tetrahedral
    type(point_3d_dat) :: sub_tets_g_coo_pnts(nver)
    type(tet_dat), allocatable :: sub_tets(:)
    type(tet_dat) :: global_sub_tet
    !
    volume = 0.0_rk
    elements: do e = 1, size(finite_elements)
      ! Get element coordinates and calculate lsf
      do i = 1, nver
        x = finite_elements(e)%get_ver_coo(i)
        lsf(i) = test_fun(x)
      end do
      ! Find sub tets
      call set_sub_xtets(lsf,sub_tets,istat,emsg)
      if ( istat /= 0 ) return
      if ( .not. allocated(sub_tets) ) cycle elements
      !write(stdout,'(a,i0)') 'Element number: ', e
      !write(stdout,'(a,4('//es//',:,","))') 'Lsf values: ', lsf
      !write(stdout,'(a,i0)') 'Nuber of subtets: ', size(sub_tets)
      ! Calculate integral
      do t = 1, size(sub_tets)
        !
        !write(stdout,'(a,i0,a)') 'Sub tet ', t, ' local coordinates'
        !call sub_tets(t)%write()
        ! Global coordinates for
        do i = 1, nver
          call finite_elements(e)%g_coo(xi_coo_pnt=sub_tets(t)%vert(i),&
          &g_coo_pnt=sub_tets_g_coo_pnts(i),istat=istat,emsg=emsg)
          if ( istat /= 0 ) then
            write(stdout,'(a,i0)') 'Element number: ', e
            write(stdout,'(a,4('//es//',:,","))') 'Lsf values: ', lsf
            write(stdout,'(a,i0)') 'Nuber of subtets: ', size(sub_tets)
            write(stdout,'(a,i0,a)') 'Sub tet ', t, ' local coordinates'
            call sub_tets(t)%write()
            return
          end if
        end do
        ! Global_sub_tet
        global_sub_tet = tet_dat(vert=sub_tets_g_coo_pnts)
        !write(stdout,'(a,i0,a)') 'Sub tet ', t, ' global coordinates'
        !call global_sub_tet%write()
        call global_sub_tet%volume(volume_tet,istat,emsg)
        if ( istat /= 0 ) return
        !write(stdout,'(a,i0,a,'//es//')') 'Subtet number: ', t, ', volume: ', &
        !& volume_tet
        volume = volume + volume_tet
      end do
      ! Deallocate sub tets
      deallocate(sub_tets,stat=istat,errmsg=emsg)
      if ( istat /= 0 ) return
    end do elements
    ! Sucess
    istat = 0
    emsg = ''
  end subroutine calc_vol
!*****************************************************************************80
! Test functions
!*****************************************************************************80
  subroutine test_functions(istat,emsg)
    !
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    real(rk), parameter :: s_1_vol = 0.0334707_rk, s_2_vol = 0.0591926_rk
    real(rk) :: volume
    character(len=*), parameter :: w_vol = '(a,'//es//')'
    !
    write(stdout,'(a)') 'Writing test functions to odb ...'
    call write_test_functions(istat,emsg)
    if( istat /= 0 ) return
    !
    write(stdout,'(a)') 'Calculating volume for ellipsoid ...'
    call calc_vol(s_1,volume,istat,emsg)
    if ( istat /= 0 ) return
    write(stdout,w_vol) 'Relative error is: ', abs(volume- s_1_vol) /  s_1_vol
    !
    write(stdout,'(a)') 'Calculating volume for torus ...'
    call calc_vol(s_2,volume,istat,emsg)
    if ( istat /= 0 ) return
    write(stdout,w_vol) 'Relative error is: ', abs(volume- s_2_vol) /  s_2_vol


  end subroutine test_functions
!*****************************************************************************80
end module volume_integral
