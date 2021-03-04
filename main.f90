module global_variables

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical parameter
  real(8),parameter :: fs=0.024189d0


! physical system
  integer :: nsite, nelec
  real(8), allocatable :: ham0(:,:), xx(:)
  real(8) :: delta_gap, t_hop, delta_gap0
  real(8) :: mass, lattice_a
  complex(8),allocatable :: zpsi(:,:)
  real(8),allocatable :: phi_gs(:,:),sp_energy(:)
  real(8),allocatable :: phi_gs_transpose(:,:)

! time propagation
  integer :: nt
  real(8) :: Tprop, dt

! laser fields
  real(8) :: omega0, Efield0, Tpulse0
  real(8),allocatable :: Efield_t(:)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call preparation

!  stop
  call time_propagation

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

! system parameters
! CdS: PRB 39, 10935 (1989)
  lattice_a = 5.82d0/0.529d0
  mass = 1d0/(1d0/0.18d0+1d0/0.53d0)
  delta_gap0= 2.5d0/27.2114d0
  delta_gap = delta_gap0

  write(*,*)"lattice_a=",lattice_a
  write(*,*)"mass     =",mass
  write(*,*)"delta_gap=",delta_gap

  t_hop     = 0.5d0/lattice_a*sqrt(delta_gap0/mass)
  write(*,*)"t_hop    =",t_hop

  nelec = 1
  nsite = 2*nelec+1

! laser fields
  omega0 = 0.35424d0/27.2114d0 ! 3.5 mum
  Efield0 = 20d6 *(0.529d-8/27.2114d0) ! MV/cm
!  Tpulse0 = 10d0*2d0*pi/omega0
  Tpulse0 = 0.5d0*pi*(80d0/fs)/acos((0.5d0)**(1d0/8d0))
  write(*,"(A,2x,e26.16e3)")"Tpulse0 [fs]=",Tpulse0*fs

! time-propagation
  Tprop = Tpulse0
  dt = 0.2d0
  nt = aint(Tprop/dt)+1
  dt = Tprop/nt
  write(*,"(A,2x,e26.16e3)")"refined dt=",dt
  write(*,"(A,2x,I7)")"nt        =",nt

end subroutine input
!-------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer :: isite, jsite, istate
  real(8) :: x0
!==LAPACK==
  real(8),allocatable :: amat_t(:,:)
  real(8),allocatable :: work(:)
  integer :: lwork, info

  lwork = 4*nsite**2 + 16*nsite
  allocate(work(lwork))
  allocate(amat_t(0:nsite-1,0:nsite-1))
!==LAPACK==

  allocate(ham0(0:nsite-1,0:nsite-1), xx(0:nsite-1))
  allocate(zpsi(0:nsite-1,0:nelec-1))
  allocate(phi_gs(0:nsite-1,0:nsite-1),sp_energy(0:nsite-1))
  allocate(phi_gs_transpose(0:nsite-1,0:nsite-1))

  ham0 = 0d0
  do isite = 0, nsite-1
    ham0(isite, isite) = 0.5d0*delta_gap*(-1)**isite

    jsite = isite + 1
    if(jsite <=  nsite-1)then
      ham0(isite, jsite) = -t_hop
      ham0(jsite, isite) = -t_hop
    end if
  end do
  
  amat_t = ham0
  call dsyev('V','U',nsite,amat_t,nsite,sp_energy,work,lwork,info)
  phi_gs = amat_t
  phi_gs_transpose = transpose(phi_gs)

  write(*,"(A,2x,1e26.16e3)")"Diameterp=",nsite*lattice_a*0.5d0
  write(*,"(A,2x,1e26.16e3)")"E_gap=",sp_energy(nelec)-sp_energy(nelec-1)
  write(*,"(A,2x,1e26.16e3)")"delta E_gap=",sp_energy(nelec)-sp_energy(nelec-1)-delta_gap0
  open(20,file='eigenvalue.out')
  do isite = 0, nsite-1
    write(20,"(I7,2x,999e26.16e3)")isite,sp_energy(isite),sum(phi_gs(:,isite)**2)
  end do
  close(20)

  zpsi(:,0:nelec-1) = phi_gs(:,0:nelec-1)

  do isite = 0, nsite-1
    xx(isite)= isite
  end do

  x0 = 0d0
  do istate = 0, nelec-1
    x0 = x0 + sum(xx*abs(zpsi(:,istate))**2)
  end do
  x0 = x0/nelec
  xx = xx -x0

end subroutine preparation
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: dip, ngs, nex

  call init_laser_field

  open(20,file='Et_Dt.out')
  open(21,file='nex_t.out')
  do it = 0,nt

    call calc_dipole(dip)
    write(20,"(999e26.16e3)")it*dt,Efield_t(it),dip

    if(mod(it,100)==0 .or. it == nt)then
      call calc_nex(ngs, nex)
      write(21,"(999e26.16e3)")it*dt,ngs,nex
    end if

    call dt_evolve(it)

  end do
  close(20)
  close(21)


end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine init_laser_field
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt, ss

  allocate(Efield_t(-1:nt+1))
  Efield_t = 0d0

  do it = 0, nt
    tt = dt*it
    ss = (tt - 0.5d0*Tpulse0)
    if(abs(ss)<= 0.5d0*Tpulse0)then
      Efield_t(it) = Efield0*sin(omega0*ss)*cos(pi*ss/Tpulse0)**4
    end if
        
  end do


end subroutine init_laser_field
!-------------------------------------------------------------------------------
subroutine calc_dipole(dip)
  use global_variables
  implicit none
  real(8),intent(out) :: dip
  integer :: istate

  dip = 0d0
!$omp parallel do reduction(+:dip)
  do istate = 0, nelec-1
    dip = dip + sum(xx*abs(zpsi(:,istate))**2)
  end do
  dip = dip/nelec

end subroutine calc_dipole
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: istate, jstate
  complex(8) :: zexp_ziVdt(0:nsite-1),zc(0:nsite-1, 0:nelec-1)

  
  zexp_ziVdt = exp(-zi*0.5d0*dt*Efield_t(it)*xx)
!$openmp parallel do
  do istate = 0, nelec-1
    zpsi(:,istate) = zexp_ziVdt(:)*zpsi(:,istate)
  end do

!$openmp parallel do private(istate,jstate) collapse(2)
  do istate = 0, nelec-1
    do jstate = 0, nsite-1
      zc(jstate,istate) = sum(phi_gs(:,jstate)*zpsi(:,istate))
    end do
  end do

!$openmp parallel do
  do jstate = 0, nsite-1
    zc(jstate,:) = zc(jstate,:)*exp(-zi*dt*sp_energy(jstate))
  end do

  zpsi(:,:) =  0d0
!$openmp parallel do private(istate,jstate)
  do istate = 0, nelec-1
    do jstate = 0, nsite-1
      zpsi(:,istate) = zpsi(:,istate) + phi_gs(:,jstate)*zc(jstate,istate)
    end do
  end do

  zexp_ziVdt = exp(-zi*0.5d0*dt*Efield_t(it+1)*xx)
!$openmp parallel do
  do istate = 0, nelec-1
    zpsi(:,istate) = zexp_ziVdt(:)*zpsi(:,istate)
  end do

end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine calc_nex(ngs, nex)
  use global_variables
  implicit none
  real(8),intent(out) :: ngs, nex
  real(8) :: occ_t(0:nelec-1)
  integer :: ib1, ib2

  occ_t = 0d0

!$openmp parallel do private(ib1, ib2)
  do ib1 =0 , nelec-1
    do ib2 = 0, nelec-1
      occ_t(ib1) = occ_t(ib1) + abs(sum(phi_gs(:,ib1)*zpsi(:,ib2)))**2
    end do
  end do

  ngs = sum(occ_t)/nelec
  nex = 1d0 -ngs
  
  
end subroutine calc_nex
!-------------------------------------------------------------------------------
