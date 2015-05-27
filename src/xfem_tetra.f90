module xtet
  use types
  use tet_volume
  use determinant
  implicit none
  private
!*****************************************************************************80
  real(rk), parameter :: volume_trsh = 0.0_rk, strech_trsh=0.0_rk
!*****************************************************************************80
  public :: set_sub_xtets
!*****************************************************************************80
contains
!*****************************************************************************80
  subroutine set_sub_xtets(lsf,sub_tets,istat,emsg)
    !
    real(rk), intent(in) :: lsf(:)
    type(tet_dat), allocatable, intent(out) :: sub_tets(:)
    integer(ik), intent(out) :: istat
    character(len=*), intent(out) :: emsg
    !
    real(rk) :: volume, strech
    integer(ik) :: i, i1, eta, n_dealloc
    integer(ik) :: nod_in(4), nod_tet(4)

  ! Definiraj vozlisca v izoparamtericnih koordinatah
    real(rk), parameter :: xi(4,3) =           &
    &     reshape([ 0.0_rk, 0.0_rk, 0.0_rk,    &
    &               1.0_rk, 0.0_rk, 0.0_rk,    &
    &               0.0_rk, 1.0_rk, 0.0_rk,    &
    &               0.0_rk, 0.0_rk, 1.0_rk  ], &
    &     shape=[4,3],order=[2,1])

  ! Preveri stevilo pozitivnih vrednosti nivojske funkcije
    eta=0
    do i = 1, 4
      if(lsf(i) > 0.0_rk ) then
        eta = eta + 1
        nod_in(eta)=i   ! shrani vozlisce s pozitivnim lsf
      end if
    end do
    select case(eta)
!*****************************************************************************80
    case( 0 ) ! Prazen element
!*****************************************************************************80
    case( 1 ) ! Poln le vrh tetraedra
      allocate(sub_tets(1),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      do i = 1, 4
        if (i==nod_in(1)) then
          ! Ohrani vozlisce s pozitivno lsf vrednostjo
          sub_tets(1)%vert(i)%x = xi(i,:)
        else
          ! Delitev stranic
          sub_tets(1)%vert(i)%x = xi(i,:) - ((xi(i,:) -  &
          & xi(nod_in(1),:))/(lsf(i) - lsf(nod_in(1))))*lsf(i)
        endif
      end do
      call calc_tet_vol(sub_tets(1),volume,istat,emsg)
      call calc_tet_strech(sub_tets(1),volume,strech,istat,emsg)
      if ( (volume.lt.volume_trsh) .or. (strech.lt.strech_trsh) ) then
        deallocate(sub_tets,stat=istat,errmsg=emsg)
        if( istat /= 0 ) return
      endif
!*****************************************************************************80
    case( 2 ) != Polna polovica tetraedra
      allocate(sub_tets(3),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      n_dealloc=0

      if ((nod_in(1)==1).and.(nod_in(2)==2)) then
          nod_tet=[1,2,3,4]

      else if ((nod_in(1)==1).and.(nod_in(2)==3)) then
          nod_tet=[3,1,2,4]

      else if ((nod_in(1)==1).and.(nod_in(2)==4)) then
          nod_tet=[1,4,2,3]

      else if ((nod_in(1)==2).and.(nod_in(2)==3)) then
          nod_tet=[2,3,4,1]

      else if ((nod_in(1)==2).and.(nod_in(2)==4)) then
          nod_tet=[2,4,1,3]

      else if ((nod_in(1)==3).and.(nod_in(2)==4)) then
          nod_tet=[3,4,1,2]
      else
        istat = -1
        emsg = 'Set sub xtets: No variants'
        return
      end if
      ! Definiranje treh tetraedrov, ki delijo osnovnega
!*****************************************************************************80
      ! prvi tetraeder
      sub_tets(1)%vert(1)%x = xi(nod_tet(1),:)
      sub_tets(1)%vert(2)%x = xi(nod_tet(2),:)
      sub_tets(1)%vert(3)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) -  &
      & xi(nod_tet(3),:))/(lsf(nod_tet(1)) - lsf(nod_tet(3))))*lsf(nod_tet(1))
      sub_tets(1)%vert(4)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) -  &
      & xi(nod_tet(4),:))/(lsf(nod_tet(1)) - lsf(nod_tet(4))))*lsf(nod_tet(1))
      call calc_tet_vol(sub_tets(1),volume,istat,emsg)
      call calc_tet_strech(sub_tets(1),volume,strech,istat,emsg)
      if ((volume.lt.volume_trsh).or.(strech.lt.strech_trsh)) then
        call deallocate_tet(sub_tets, 1)
        n_dealloc=n_dealloc+1
      end if
!*****************************************************************************80
      ! drugi tetraeder
      sub_tets(2-n_dealloc)%vert(1)%x = xi(nod_tet(2),:)
      sub_tets(2-n_dealloc)%vert(2)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) &
      & - xi(nod_tet(3),:))/(lsf(nod_tet(1)) - lsf(nod_tet(3))))*lsf(nod_tet(1))
      sub_tets(2-n_dealloc)%vert(4)%x = xi(nod_tet(2),:) - ((xi(nod_tet(2),:) &
      & - xi(nod_tet(3),:))/(lsf(nod_tet(2)) - lsf(nod_tet(3))))*lsf(nod_tet(2))
      sub_tets(2-n_dealloc)%vert(3)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) &
      & - xi(nod_tet(4),:))/(lsf(nod_tet(1)) - lsf(nod_tet(4))))*lsf(nod_tet(1))
      call calc_tet_vol(sub_tets(2-n_dealloc),volume,istat,emsg)
      call calc_tet_strech(sub_tets(2-n_dealloc),volume,strech,istat,emsg)
      if ((volume.lt.volume_trsh).or.(strech.lt.strech_trsh)) then
        call deallocate_tet(sub_tets, 2-n_dealloc)
        n_dealloc=n_dealloc+1
      end if
!*****************************************************************************80
      ! tretji tetraeder
      sub_tets(3-n_dealloc)%vert(1)%x = xi(nod_tet(2),:)
      sub_tets(3-n_dealloc)%vert(2)%x = xi(nod_tet(2),:) - ((xi(nod_tet(2),:) &
      & - xi(nod_tet(3),:))/(lsf(nod_tet(2)) - lsf(nod_tet(3))))*lsf(nod_tet(2))
      sub_tets(3-n_dealloc)%vert(3)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) &
      & - xi(nod_tet(4),:))/(lsf(nod_tet(1)) - lsf(nod_tet(4))))*lsf(nod_tet(1))
      sub_tets(3-n_dealloc)%vert(4)%x = xi(nod_tet(2),:) - ((xi(nod_tet(2),:) &
      & - xi(nod_tet(4),:))/(lsf(nod_tet(1)) - lsf(nod_tet(4))))*lsf(nod_tet(1))
      call calc_tet_vol(sub_tets(3-n_dealloc),volume,istat,emsg)
      call calc_tet_strech(sub_tets(3-n_dealloc),volume,strech,istat,emsg)
      if ((volume.lt.volume_trsh).or.(strech.lt.strech_trsh)) then
        call deallocate_tet(sub_tets, 3-n_dealloc)
        n_dealloc=n_dealloc+1
      endif
!*****************************************************************************80
    CASE( 3 ) != Odreze vrh tetraedra
      n_dealloc=0
      allocate(sub_tets(3),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return

      if ((nod_in(1)==1).and.(nod_in(2)==2).and.(nod_in(3)==3)) then
          nod_tet=[1,2,3,4]

      else if ((nod_in(1)==2).and.(nod_in(2)==3).and.(nod_in(3)==4)) then
          nod_tet=[2,3,4,1]

      else if ((nod_in(1)==1).and.(nod_in(2)==2).and.(nod_in(3)==4)) then
          nod_tet=[1,2,4,3]

      else if ((nod_in(1)==1).and.(nod_in(2)==3).and.(nod_in(3)==4)) then
          nod_tet=[1,3,4,2]

      else
        istat = -1
        emsg = 'Set sub xtets: No variants'
        return
      endif
      ! Definiranje treh tetraedrov, ki delijo osnovnega
!*****************************************************************************80
      ! prvi tetraeder
      sub_tets(1)%vert(1)%x = xi(nod_tet(1),:)
      sub_tets(1)%vert(2)%x = xi(nod_tet(2),:)
      sub_tets(1)%vert(3)%x = xi(nod_tet(3),:)
      sub_tets(1)%vert(4)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) - &
      &xi(nod_tet(4),:))/(lsf(nod_tet(1)) - lsf(nod_tet(4))))*lsf(nod_tet(1))
!*****************************************************************************80
      ! drugi tetraeder
      call calc_tet_vol(sub_tets(1),volume,istat,emsg)
      call calc_tet_strech(sub_tets(1),volume,strech,istat,emsg)
      if ((volume.lt.volume_trsh).or.(strech.lt.strech_trsh)) then
        call deallocate_tet(sub_tets, 1-n_dealloc)
        n_dealloc=n_dealloc+1
      end if
      sub_tets(2-n_dealloc)%vert(1)%x = xi(nod_tet(2),:)
      sub_tets(2-n_dealloc)%vert(2)%x = xi(nod_tet(3),:)
      sub_tets(2-n_dealloc)%vert(3)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) -&
      & xi(nod_tet(4),:))/(lsf(nod_tet(1)) - lsf(nod_tet(4))))*lsf(nod_tet(1))
      sub_tets(2-n_dealloc)%vert(4)%x = xi(nod_tet(3),:) - ((xi(nod_tet(3),:) -&
      & xi(nod_tet(4),:))/(lsf(nod_tet(3)) - lsf(nod_tet(4))))*lsf(nod_tet(3))
      call calc_tet_vol(sub_tets(2-n_dealloc),volume,istat,emsg)
      call calc_tet_strech(sub_tets(2-n_dealloc),volume,strech,istat,emsg)
      if ((volume.lt.volume_trsh).or.(strech.lt.strech_trsh)) then
        call deallocate_tet(sub_tets, 2-n_dealloc)
        n_dealloc=n_dealloc+1
      end if
!*****************************************************************************80
      ! tretji tetraeder
      sub_tets(3-n_dealloc)%vert(1)%x = xi(nod_tet(2),:)
      sub_tets(3-n_dealloc)%vert(2)%x = xi(nod_tet(1),:) - ((xi(nod_tet(1),:) -&
      & xi(nod_tet(4),:))/(lsf(nod_tet(1)) - lsf(nod_tet(4))))*lsf(nod_tet(1))
      sub_tets(3-n_dealloc)%vert(3)%x = xi(nod_tet(2),:) - ((xi(nod_tet(2),:) -&
      & xi(nod_tet(4),:))/(lsf(nod_tet(2)) - lsf(nod_tet(4))))*lsf(nod_tet(2))
      sub_tets(3-n_dealloc)%vert(4)%x = xi(nod_tet(3),:) - ((xi(nod_tet(3),:) -&
      xi(nod_tet(4),:))/(lsf(nod_tet(3)) - lsf(nod_tet(4))))*lsf(nod_tet(3))
      call calc_tet_vol(sub_tets(3-n_dealloc),volume,istat,emsg)
      call calc_tet_strech(sub_tets(3-n_dealloc),volume,strech,istat,emsg)
      if ((volume.lt.volume_trsh).or.(strech.lt.strech_trsh)) then
        call deallocate_tet(sub_tets, 3-n_dealloc)
        n_dealloc=n_dealloc+1
      endif
!*****************************************************************************80
    case( 4 ) != poln tetraeder
      allocate(sub_tets(1),stat=istat,errmsg=emsg)
      if( istat /= 0 ) return
      do i1=1,4
        sub_tets(1)%vert(i1)%x = xi(i1,:)
      end do
    end select
!*****************************************************************************80
    istat = 0
    emsg = ''
  end subroutine
!*****************************************************************************80
end module xtet

