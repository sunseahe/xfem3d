module volume_integral
  use types
  use point, only: dom, point_3d_t
  use fe_c3d10, only:  nelnod, c3d10_t, c3d10_t_ll
  use read_input, only: nodes, finite_elements
  use xtet, only: set_sub_xtets
  use lsf_test_functions
  use write_odb, only: start_odb_api, finish_odb_api, write_model_data,  &
  & create_step, create_frame, write_scalar_field, close_odb_file
  use read_input, only: input_file_name
  implicit none
!*****************************************************************************80
  ! For writing in input file and checking global sub tets
  integer(ik) :: gst_inp_file = 0
  logical(lk) :: wrt_gst_inpf = .true.
!*****************************************************************************80
  public :: test_functions
!*****************************************************************************80
contains
!*****************************************************************************80
! Write functions to odb
!*****************************************************************************80
  subroutine write_test_functions(esta,emsg)
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: e, i
    integer(ik) :: active_node
    real(rk) :: x(dom)
    real(rk), allocatable :: lsf_field(:)
    character(len=cl), parameter :: step = 'Test step', field_name = 'LSF', &
    & field_description = 'Level set function values'
    logical(lk), allocatable :: lsf_field_val_cal(:)
    !
    allocate(lsf_field(size(nodes)),lsf_field_val_cal(size(nodes)),&
    & stat=esta,errmsg=emsg)
    ! Start odb api
    call start_odb_api()
    ! Model data
    call write_model_data(1,esta,emsg)
    if( esta /= 0 ) return
    ! Step
    call create_step(step,'LSF values',esta,emsg)
    ! Functions
    call write_ind_fun('Ellipsoid',1,1.0_rk,s_1,esta,emsg)
    if ( esta /= 0 ) return
    call write_ind_fun('Torus',2,2.0_rk,s_2,esta,emsg)
    if ( esta /= 0 ) return
    call write_ind_fun('Genus two',3,3.0_rk,s_3,esta,emsg)
    if ( esta /= 0 ) return
    call write_ind_fun('Genus seven',4,4.0_rk,s_4,esta,emsg)
    if ( esta /= 0 ) return
    call write_ind_fun('Sphere',5,5.0_rk,s_5,esta,emsg)
    if ( esta /= 0 ) return
    ! Close odb file
    call close_odb_file(esta,emsg)
    if( esta /= 0 ) return
    ! Finish odb api
    call finish_odb_api()
    !
  contains
    subroutine write_ind_fun(name,frame_num,time,test_fun,esta,emsg)
      character(len=*), intent(in) :: name
      integer(ik), intent(in) :: frame_num
      real(rk), intent(in) :: time
      interface
        pure real(rk) function test_fun(x)
          import :: rk
          real(rk), intent(in) :: x(:)
        end function test_fun
      end interface
      integer(ik), intent(out) :: esta
      character(len=*), intent(out) :: emsg
      !
      call create_frame(step,name,frame_num,time,&
      &esta,emsg)
      if ( esta /= 0 ) return
      lsf_field = 0.0_rk
      lsf_field_val_cal = .false.
      !
      elements: do e = 1, size(finite_elements)
        element_nodes: do i = 1, nelnod
          x = finite_elements(e)%get_node_coo(i)
          active_node = finite_elements(e)%connectivity(i)
          if ( lsf_field_val_cal(active_node) ) cycle element_nodes
          lsf_field_val_cal(active_node) = .true.
          lsf_field(active_node) = test_fun(x)
        end do element_nodes
      end do elements
      !
      call write_scalar_field(step,frame_num,field_name,&
      & field_description,lsf_field,.true.,esta,emsg)
      if ( esta /= 0 ) return
      ! Sucess
      esta = 0
      emsg = ''
    end subroutine write_ind_fun
  end subroutine write_test_functions
!*****************************************************************************80
! Calculation of volume
!*****************************************************************************80
!  subroutine calc_vol(test_fun,volume,esta,emsg)
!    !
!    interface
!      pure real(rk) function test_fun(x)
!        import :: rk
!        real(rk), intent(in) :: x(:)
!      end function test_fun
!    end interface
!    real(rk), intent(out) :: volume
!    integer(ik), intent(out) :: esta
!    character(len=*), intent(out) :: emsg
!    !
!    integer(ik) :: e, i, t
!    integer(ik) :: connectivity_gst(nelnod)
!    integer(ik) :: count_gst ! Count global subtets
!    real(rk) :: x(dom) ! Coordinates
!    real(rk) :: lsf(nelnod) ! Level set values
!    real(rk) :: volume_tet ! Volume of a single tetrahedral
!    type(point_3d_t) :: sub_tets_g_coo_pnts(nelnod)
!    type(c3d4_t), allocatable :: sub_tets(:)
!    type(c3d4_t) :: global_sub_tet
!    ! For writing in input file and checking global sub tets
!    type(c3d4_t), allocatable :: global_sub_tets(:)
!    type(c3d4_t_ll) :: global_sub_tets_ll
!    !
!    volume = 0.0_rk
!    count_gst = 0
!    elements: do e = 1, size(finite_elements)
!      ! Get element coordinates and calculate lsf
!      do i = 1, nelnod
!        x = finite_elements(e)%get_node_coo(i)
!        lsf(i) = test_fun(x)
!      end do
!      ! Find sub tets
!      call set_sub_xtets(lsf,sub_tets,esta,emsg)
!      if ( esta /= 0 ) return
!      if ( .not. allocated(sub_tets) ) cycle elements
!      !write(stdout,'(a,i0)') 'Element number: ', e
!      !write(stdout,'(a,4('//es//',:,","))') 'Lsf values: ', lsf
!      !write(stdout,'(a,i0)') 'Nuber of subtets: ', size(sub_tets)
!      ! Calculate integral
!      do t = 1, size(sub_tets)
!
!        !
!        !write(stdout,'(a,i0,a)') 'Sub tet ', t, ' local coordinates'
!        !call sub_tets(t)%write()
!        ! Global coordinates for
!        do i = 1, nelnod
!          call finite_elements(e)%g_coo(xi_coo_pnt=sub_tets(t)%nodes(i),&
!          &g_coo_pnt=sub_tets_g_coo_pnts(i),esta=esta,emsg=emsg)
!          if ( esta /= 0 ) then
!            write(stdout,'(a,i0)') 'Element number: ', e
!            write(stdout,'(a,4('//es//',:,","))') 'Lsf values: ', lsf
!            write(stdout,'(a,i0)') 'Nuber of subtets: ', size(sub_tets)
!            write(stdout,'(a,i0,a)') 'Sub tet ', t, ' local coordinates'
!            call sub_tets(t)%write()
!            return
!          end if
!        end do
!        ! Global_sub_tet
!        do i = 1, nelnod
!          connectivity_gst(i) = nelnod * count_gst + i
!        end do
!        global_sub_tet = c3d4_t(nodes=sub_tets_g_coo_pnts,connectivity=&
!        &connectivity_gst)
!        !write(stdout,'(a,i0,a)') 'Sub tet ', t, ' global coordinates'
!        !call global_sub_tet%write()
!
!        if ( wrt_gst_inpf ) then
!          call global_sub_tets_ll%add(global_sub_tet,esta,emsg)
!          if ( esta /= 0 ) return
!        end if
!
!        call global_sub_tet%volume(volume_tet,esta,emsg)
!        if ( esta /= 0 ) return
!        !write(stdout,'(a,i0,a,'//es//')') 'Subtet number: ', t, ', volume: ', &
!        !& volume_tet
!        volume = volume + volume_tet
!        count_gst = count_gst + 1
!      end do
!      ! Deallocate sub tets
!      deallocate(sub_tets,stat=esta,errmsg=emsg)
!      if ( esta /= 0 ) return
!    end do elements
!!*****************************************************************************80
!    if ( wrt_gst_inpf ) then
!      call global_sub_tets_ll%fill_array(global_sub_tets,esta,emsg)
!      if ( esta /= 0 ) return
!      call global_sub_tets_ll%clean(esta,emsg)
!      if ( esta /= 0 ) return
!      ! Write nodes
!      write(gst_inp_file,fmt='(a)') '*node'
!      do e = 1, size(global_sub_tets)
!        do i = 1, nelnod
!          x = global_sub_tets(e)%get_node_coo(i)
!          write(gst_inp_file,fmt='(i0,",",3('//es//',:,","))') &
!          & global_sub_tets(e)%connectivity(i), x
!        end do
!      end do
!      ! Write elements
!      write(gst_inp_file,fmt='(a)') '*element,type=c3d4'
!      do e = 1, size(global_sub_tets)
!        write(gst_inp_file,fmt='(i0,",",4(i0,:,","))') e, &
!        & global_sub_tets(e)%connectivity
!      end do
!      close(gst_inp_file)
!    end if
!    ! Sucess
!    esta = 0
!    emsg = ''
!  end subroutine calc_vol
!*****************************************************************************80
! Test functions
!*****************************************************************************80
  subroutine test_functions(esta,emsg)
    !
    integer(ik), intent(out) :: esta
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: n
    real(rk), parameter :: s_1_vol = 0.0334707_rk, s_2_vol = 0.0591926_rk
    real(rk) :: volume
    character(len=cl) :: gst_inp_file_name
    !
    write(stdout,'(a)') 'Writing test functions to odb ...'
    call write_test_functions(esta,emsg)
    if( esta /= 0 ) return
    !
    !write(stdout,'(a)') 'Calculating volume for ellipsoid ...'
    ! Write to input file global subtets
!    if ( wrt_gst_inpf ) then
!      n = len_trim(input_file_name)
!      write(gst_inp_file_name,'(a,a,a)' ) input_file_name(1:n-4), &
!      &'_ellipsoid', '.inp'
!      open(newunit=gst_inp_file,file=gst_inp_file_name,iostat=esta,&
!      &iomsg=emsg,action='write',status='replace')
!      if ( esta /= 0 ) return
!    end if
!    call calc_vol(s_1,volume,esta,emsg)
!    if ( esta /= 0 ) return
!    write(stdout,'(a,'//es//')') 'Calculated volume is: ', volume
!    write(stdout,'(a,'//es//')') 'Relative error is: ', abs(volume - &
!    & s_1_vol) /  s_1_vol
    !
!    write(stdout,'(a)') 'Calculating volume for torus ...'
!    call calc_vol(s_2,volume,esta,emsg)
!    if ( esta /= 0 ) return
!    write(stdout,'(a,'//es//')') 'Relative error is: ', abs(volume- s_2_vol) /  s_2_vol


  end subroutine test_functions
!*****************************************************************************80
end module volume_integral