module global_variables

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical system
  integer :: nsite, nelec
  real(8), allocatable :: ham0(:,:), xx(:)
  real(8) :: delta_gap, t_hop
  complex(8),allocatable :: zpsi(:,:)
  real(8),allocatable :: phi_gs(:,:),sp_energy(:)
  real(8),allocatable :: phi_gs_transpose(:,:)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call preparation

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  delta_gap = 1d0
  t_hop     = 4d0*delta_gap

  nelec = 4
  nsite = 2*nelec+1

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
  allocate(zpsi(0:nsite-1,nelec))
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

  write(*,"(A,2x,1e26.16e3)")"E_gap=",sp_energy(nelec)-sp_energy(nelec-1)
  open(20,file='eivenvalue.out')
  do isite = 0, nsite-1
    write(20,"(I7,2x,999e26.16e3)")isite,sp_energy(isite)
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
  xx = xx -x0
  

end subroutine preparation
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
