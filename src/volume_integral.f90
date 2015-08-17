module volume_integral
  use blas95, only: dot
  use types, only: ik, rk, lk, cl, stdout, log_file
  use point, only: dom, point_3d_t
  use fe_c3d10, only:  nelnod, c3d10_t, ngp, w
  use mesh_data, only: nnod, nodes, nfe, finite_elements, char_fe_dim
  use scalar_field, only: scalar_field_t
  !use xtet, only: set_sub_xtets
  use lsf_test_functions
  use write_odb, only: write_model_data, create_step, create_frame,  &
  &  write_scalar_field, close_odb_file
  use read_input, only: input_file_name
  use reinitalzation, only: calculate_reinitalization, r_conf => configured
  use hamilton_jacobi, only: calculate_advection, a_conf => configured
!*****************************************************************************80
  implicit none
  private
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
    type(scalar_field_t) :: lsf
    character(len=cl), parameter :: step = 'Test step', field_name = 'LSF', &
    & field_description = 'Level set function values'
    logical(lk), allocatable :: lsf_field_val_cal(:)
    ! Allocate already calculated values
    allocate(lsf_field_val_cal(nnod),stat=esta,errmsg=emsg)
    ! Model data
    call write_model_data(1,esta,emsg)
    if( esta /= 0 ) return
    ! Step
    call create_step(step,'LSF values',esta,emsg)
    ! Functions
    !call write_ind_fun('Sphere',1,1.0_rk,s_5,esta,emsg)
    !if ( esta /= 0 ) return
    call write_ind_fun('Ellipsoid',1,1.0_rk,s_1,esta,emsg)
    if ( esta /= 0 ) return
    !call write_ind_fun('Torus',3,3.0_rk,s_2,esta,emsg)
    !if ( esta /= 0 ) return
    !call write_ind_fun('Genus two',1,1.0_rk,s_3,esta,emsg)
    !if ( esta /= 0 ) return
    !call write_ind_fun('Genus seven',1,1.0_rk,s_4,esta,emsg)
    !if ( esta /= 0 ) return
    ! Close odb file
    call close_odb_file(esta,emsg)
    if( esta /= 0 ) return
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
      real(rk) :: volume
      !
      call create_frame(step,name,frame_num,time,&
      &esta,emsg)
      if ( esta /= 0 ) return
      ! Reset lsf field
      call lsf%set(esta=esta,emsg=emsg); if( esta /= 0 ) return
      lsf_field_val_cal = .false.
      !
      elements: do e = 1, nfe
        element_nodes: do i = 1, nelnod
          x = finite_elements(e)%get_node_coo(i)
          active_node = finite_elements(e)%connectivity(i)
          if ( lsf_field_val_cal(active_node) ) cycle element_nodes
          lsf_field_val_cal(active_node) = .true.
          call lsf%set_nodal_value(active_node,test_fun(x))
        end do element_nodes
      end do elements
      ! Calculate volume
      call calc_vol_heaviside(lsf,volume)
      print*, 'Volume before reinitalization =', volume
      ! Reinitalize function
      if ( r_conf ) then
        call calculate_reinitalization(lsf,esta,emsg)
        if ( esta /= 0 ) return
      end if
      ! Calculate volume
      call calc_vol_heaviside(lsf,volume)
      print*, 'Volume after reinitalization =', volume
      ! Write field
      call write_scalar_field(step,frame_num,field_name,&
      & field_description,lsf%values,.true.,esta,emsg)
      if ( esta /= 0 ) return
      ! Advect
      if ( a_conf ) then
        call advect_ind_fun(name,frame_num+1,time+1.0_rk,lsf,esta,emsg)
      end if
      ! Sucess
      esta = 0
      emsg = ''
    end subroutine write_ind_fun
!*****************************************************************************80
    subroutine advect_ind_fun(name,frame_num,time,lsf,esta,emsg)
      character(len=*), intent(in) :: name
      integer(ik), intent(in) :: frame_num
      real(rk), intent(in) :: time
      type(scalar_field_t), intent(inout) :: lsf
      integer(ik), intent(out) :: esta
      character(len=*), intent(out) :: emsg
      !
      call calculate_advection(lsf,velocity,5,1.0e-1_rk,&
      &esta,emsg); if ( esta /= 0 ) return
      !
      call write_scalar_field(step,frame_num,field_name,&
      & field_description,lsf%values,.true.,esta,emsg)
      if ( esta /= 0 ) return
      !
    end subroutine advect_ind_fun
!*****************************************************************************80
    subroutine velocity(lsf,v,esta,emsg)
      type(scalar_field_t), intent(in) :: lsf
      type(scalar_field_t), intent(out) :: v
      integer(ik), intent(out) :: esta
      character(len=*), intent(out) :: emsg
      call v%set(esta=esta,emsg=emsg)
      v%values = 1.0_rk
    end subroutine velocity
!*****************************************************************************80
  end subroutine write_test_functions

  pure function heaviside(x) result(res)
    real(rk), intent(in) :: x
    real(rk) :: res
    !
    real(rk), parameter :: alpha = 1.0e-16_rk
    !
    associate ( delta => char_fe_dim )
    if( x < -delta ) then
      res = alpha
    else if ( (-delta<= x).and.(x <= delta) ) then
      res = 3.0_rk*(1.0_rk-alpha)/4.0_rk*(x/delta-x**3 &
      & /(3.0_rk*delta**3)) + (1.0_rk+alpha)/2.0_rk
    else
      res = 1.0_rk
    end if
    end associate
    !
  end function heaviside

  pure subroutine calc_vol_heaviside_fe(c3d10,lsf,volume_fe)
    type(c3d10_t), intent(in) :: c3d10
    type(scalar_field_t), intent(in) :: lsf
    real(rk), intent(out) :: volume_fe
    !
    integer(ik) :: p
    real(rk) :: n(nelnod)
    real(rk) :: det_jac
    real(rk) :: lsf_nv(nelnod), lsf_gp
    !
    volume_fe = 0.0e0_rk
    do p = 1, ngp
      call c3d10%main_values(gp_num=p,n_mtx=n,det_jac=det_jac)
      call lsf%get_element_nodal_values(c3d10,lsf_nv)
      lsf_gp = dot(n,lsf_nv)
      volume_fe = volume_fe + heaviside(lsf_gp)  * w(p) * det_jac
    end do
    !
  end subroutine calc_vol_heaviside_fe

  subroutine calc_vol_heaviside(lsf,volume)
    type(scalar_field_t), intent(in) :: lsf
    real(rk), intent(out) :: volume
    !
    integer(ik) :: e
    real(rk) :: volume_fe
    !
    volume = 0.0_rk
    !$omp parallel do & !schedule(static,1)
    !$omp private(e,volume_fe) &
    !$omp shared(finite_elements,lsf)
    do e = 1, nfe
      call calc_vol_heaviside_fe(finite_elements(e),lsf,volume_fe)
      !$omp critical
      volume = volume + volume_fe
      !$omp end critical
    end do
    !$omp end parallel do
    !
  end subroutine calc_vol_heaviside





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
    !integer(ik) :: n
    !real(rk), parameter :: s_1_vol = 0.0334707_rk, s_2_vol = 0.0591926_rk
    !real(rk) :: volume
    !character(len=cl) :: gst_inp_file_name
    !
    write(log_file,'(a)') 'Writing test functions to odb ...'
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
