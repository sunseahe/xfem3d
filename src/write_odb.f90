module write_odb
  use iso_c_binding, only: c_int, c_float, c_char, &
  & c_null_char, c_ptr, c_loc, c_bool, c_null_ptr
  use types
  use general_routines, only: f2c_char
  use read_input, only: input_file_name, fe_type, nodes, finite_elements
  use fe_c3d4, only: dom, nnod
  implicit none
  private
!*****************************************************************************80
  character(len=cl) :: odb_file_name = ''
  type(c_ptr) :: odb = c_null_ptr
!*****************************************************************************80
! Interfaces to c++ routines
!*****************************************************************************80
  interface
    ! Start odb api
    subroutine start_odb_api() bind(c,name='odb_initializeAPI')
    end subroutine start_odb_api
    ! Finish odb api
    subroutine finish_odb_api() bind(c,name='odb_finalizeAPI')
    end subroutine finish_odb_api
    ! Model data
    function model_data(odb_name,fe_type,dom,nnod,&
      & coo_ptr,nel,nnel,conn_ptr) &
      & bind(c,name='model_data') result(odb_ptr)
      import :: c_int, c_char, c_ptr
      character(c_char), intent(in) :: odb_name(*)
      character(c_char), intent(in) :: fe_type(*)
      integer(c_int), value, intent(in) :: dom
      integer(c_int), value, intent(in) :: nnod
      type(c_ptr), value, intent(in) :: coo_ptr
      integer(c_int), value, intent(in) :: nel
      integer(c_int), value, intent(in) :: nnel
      type(c_ptr), value, intent(in) :: conn_ptr
      type(c_ptr) :: odb_ptr
    end function model_data
    ! Close odb
    subroutine close_odb(odb_ptr) bind(c,name='close_odb')
      import :: c_ptr
      type(c_ptr), value, intent(in) :: odb_ptr
    end subroutine close_odb
    ! Create step
    subroutine step(odb_ptr,step_name,&
      & step_desc) bind(c,name='step')
      import :: c_ptr, c_char
      type(c_ptr), value, intent(in) :: odb_ptr
      character(c_char), intent(in) :: step_name(*)
      character(c_char), intent(in) :: step_desc(*)
    end subroutine step
    ! Create frame
    subroutine frame(odb_ptr,step_name,&
      & frame_desc,fnum,analy_time) bind(c,name='frame')
      import :: c_int, c_float, c_char, c_ptr
      type(c_ptr), value, intent(in) :: odb_ptr
      character(c_char), intent(in) :: step_name(*)
      character(c_char), intent(in) :: frame_desc(*)
      integer(c_int), value, intent(in) :: fnum
      real(c_float), value, intent(in) :: analy_time
    end subroutine frame
    ! Write scalar field
    subroutine scalar_field(odb_ptr,step_name,&
      &fnum,field_name,field_desc,nnod,field_vals,df) &
      & bind(c,name='scalar_field')
      import :: c_int, c_char, c_ptr, c_bool
      type(c_ptr), value, intent(in) :: odb_ptr
      character(c_char), intent(in) :: step_name(*)
      integer(c_int), value, intent(in) :: fnum
      character(c_char), intent(in) :: field_name(*)
      character(c_char), intent(in) :: field_desc(*)
      integer(c_int), value, intent(in) :: nnod
      type(c_ptr), value, intent(in) :: field_vals
      logical(c_bool), value, intent(in) :: df
    end subroutine scalar_field
    !
  end interface
!*****************************************************************************80
  public :: start_odb_api, finish_odb_api, write_model_data, create_step, &
  & create_frame, write_scalar_field, close_odb_file
!*****************************************************************************80
  contains
!*****************************************************************************80
! Write model data
!*****************************************************************************80
  subroutine write_model_data(odb_num,istat,emsg)
    !
    integer(ik), intent(in) :: odb_num
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: n
    integer(ik) :: nnod, nel
    integer(c_int), allocatable, target :: connectivity(:,:)
    real(c_float), allocatable, target :: coordinates(:,:)
    character(len=1,kind=c_char), allocatable :: c_odb_file_name(:),&
    & c_fe_type(:)
    type(c_ptr) :: coordinates_ptr, connectivity_ptr
    ! Odb file name
    n = len_trim(input_file_name)
    write(odb_file_name,'(a,a,i3.3,a)' ) input_file_name(1:n-4), &
    &'_', odb_num, '.odb'
    call f2c_char(odb_file_name,c_odb_file_name,istat,emsg)
    ! Fe type
    call f2c_char(fe_type,c_fe_type,istat,emsg)
    ! Node coordinates and connectivity
    nnod = size(nodes)
    nel = size(finite_elements)
    allocate(coordinates(nnod,dom),stat=istat,errmsg=emsg)
    if ( istat /= 0 ) return
    allocate(connectivity(nel,nnod),stat=istat,errmsg=emsg)
    if ( istat /= 0 ) return
    ! Copy
    do n = 1, nnod
      coordinates(n,:) = real(nodes(n)%x,kind=c_float)
    end do
    do n = 1, nel
      connectivity(n,:) = int(finite_elements(n)%connectivity,kind=c_int)
    end do
    ! Set pointers
    coordinates_ptr = c_loc(coordinates(1,1))
    connectivity_ptr = c_loc(connectivity(1,1))
    ! Write model data
    odb = model_data(c_odb_file_name,c_fe_type,int(dom,kind=c_int), &
    &int(nnod,kind=c_int),coordinates_ptr,int(nel,kind=c_int), &
    & int(nnod,kind=c_int),connectivity_ptr)
    !
  end subroutine write_model_data
!*****************************************************************************80
  subroutine create_step(step_name,step_desc,istat,emsg)
    character(len=*),intent(in) :: step_name
    character(len=*),intent(in) :: step_desc
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    character(len=1,kind=c_char), allocatable :: c_step_name(:), &
    & c_step_desc(:)
    !
    call f2c_char(step_name,c_step_name,istat,emsg); if ( istat /= 0 ) return
    call f2c_char(step_desc,c_step_desc,istat,emsg); if ( istat /= 0 ) return
    call step(odb,c_step_name,c_step_desc)
    !
  end subroutine create_step
!*****************************************************************************80
! Create frame
!*****************************************************************************80
  subroutine create_frame(step_name,frame_desc,frame_num,time,&
    &istat,emsg)
    character(len=*),intent(in) :: step_name
    character(len=*),intent(in) :: frame_desc
    integer, intent(in) :: frame_num
    real(rk), intent(in) :: time
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    character(len=1,kind=c_char), allocatable :: c_step_name(:), &
    & c_frame_desc(:)
    !
    call f2c_char(step_name,c_step_name,istat,emsg); if ( istat /= 0 ) return
    call f2c_char(frame_desc,c_frame_desc,istat,emsg); if ( istat /= 0 ) return
    call frame(odb,c_step_name,c_frame_desc,&
    & int(frame_num,kind=c_int),real(time,kind=c_float))
    !
  end subroutine create_frame
!*****************************************************************************80
! Write scalar field
!*****************************************************************************80
  subroutine write_scalar_field(step_name,frame_num,field_name,&
    & field_description,field_values,field_default,istat,emsg)
    !
    character(len=*),intent(in) :: step_name
    integer, intent(in) :: frame_num
    character(len=*),intent(in) :: field_name
    character(len=*),intent(in) :: field_description
    real(rk), intent(in) :: field_values(:)
    logical, intent(in) :: field_default
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    integer(ik) :: nnod
    real(c_float), allocatable, target :: c_field_values(:)
    character(len=1,kind=c_char), allocatable :: c_step_name(:), &
    &c_field_name(:), c_field_description(:)
    type(c_ptr) :: field_values_ptr
    !
    call f2c_char(step_name,c_step_name,istat,emsg)
    if ( istat /= 0 ) return
    call f2c_char(field_name,c_field_name,istat,emsg)
    if ( istat /= 0 ) return
    call f2c_char(field_description,c_field_description,istat,emsg)
    if ( istat /= 0 ) return
    ! Field values
    nnod = size(nodes)
    allocate(c_field_values(size(field_values)),stat=istat,errmsg=emsg)
    if ( istat /= 0 ) return
    c_field_values = real(field_values,kind=c_float)
    field_values_ptr = c_loc(c_field_values(1))
    ! Write field
    call scalar_field(odb,c_step_name,&
    &int(frame_num,kind=c_int),c_field_name,c_field_description,&
    &int(nnod,kind=c_int),field_values_ptr, &
    &logical(field_default,kind=c_bool))
    !
  end subroutine write_scalar_field
!*****************************************************************************80
! Write scalar field
!*****************************************************************************80
  subroutine close_odb_file(istat,emsg)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    call close_odb(odb)
    odb = c_null_ptr
    istat = 0
    emsg = ''
  end subroutine close_odb_file
!*****************************************************************************80
end module write_odb
