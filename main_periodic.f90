module global_variables

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical parameter
  real(8),parameter :: fs=0.024189d0


! physical system
  integer :: nk
  real(8) :: delta_gap, t_hop
  real(8) :: mass, lattice_a
  real(8),allocatable :: ham_kt(:,:,:)
  complex(8),allocatable :: zpsi(:,:)
  real(8),allocatable :: phi_gs(:,:,:),sp_energy(:,:)
  real(8),allocatable :: kn(:)


! time propagation
  integer :: nt
  real(8) :: Tprop, dt

! laser fields
  real(8) :: omega0, Efield0, Tpulse0
  real(8),allocatable :: Efield_t(:), Afield_t(:)



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
  delta_gap = 9d0/27.2114d0

  write(*,*)"lattice_a=",lattice_a
  write(*,*)"mass     =",mass
  write(*,*)"delta_gap=",delta_gap

  t_hop     = 0.5d0/lattice_a*sqrt(delta_gap/mass)
  write(*,*)"t_hop    =",t_hop

  nk = 2

! laser fields
  omega0 = 1.55d0/27.2114d0 ! 3.5 mum
  Efield0 = 1d4 *(0.529d-8/27.2114d0) ! MV/cm
!  Efield0 = 0d0
!  Tpulse0 = 10d0*2d0*pi/omega0
  Tpulse0 = 0.5d0*pi*(10d0/fs)/acos((0.5d0)**(1d0/8d0))
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
  integer :: ik
  real(8) :: ham(2,2)

  allocate(ham_kt(2,2,0:nk-1))
  allocate(zpsi(2,0:nk-1), phi_gs(2,2,0:nk-1),sp_energy(2,0:nk-1))
  allocate(kn(0:nk-1))

  do ik = 0, nk-1
    kn(ik) = (2d0/lattice_a)*(pi/nk)*ik
  end do


  do ik = 0, nk-1
    ham(1,1) = -0.5d0*delta_gap
    ham(1,2) = -2d0*t_hop*cos(0.5d0*lattice_a*kn(ik))
    ham(2,1) = ham(1,2)
    ham(2,2) = 0.5d0*delta_gap

    call diag_2x2(ham, phi_gs(:,:,ik), sp_energy(:,ik))
    zpsi(:,ik) = phi_gs(:,2,ik)
    write(*,*)phi_gs(:,2,ik) ! debug
  end do

end subroutine preparation
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: current
  real(8) :: nex_bloch,nex_houston,nex_pol_houston

  call init_laser_field

  open(20,file='current.out')
  open(21,file='nex.out')
  do it = 0,nt

    call calc_current(current,it)
    write(20,"(999e26.16e3)")it*dt,Afield_t(it),Efield_t(it),current

    call calc_nex(nex_bloch,nex_houston,nex_pol_houston,it)
    write(21,"(999e26.16e3)")it*dt,nex_bloch,nex_houston,nex_pol_houston

    call dt_evolve(it)

  end do
  close(20)
  close(21)

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ik
  real(8) :: ham(2,2), vec(2,2), lambda(2)
  complex(8) :: zvec(2)
  real(8) :: tt, kt




  do ik = 0, nk-1

! propagation from dt*it to dt*(it+0.5)
    tt = dt*it

    kt = kn(ik) + Afield_t(it)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)
    
    zvec(1) = vec(1,1)*zpsi(1,ik)+ vec(2,1)*zpsi(2,ik)
    zvec(2) = vec(1,2)*zpsi(1,ik)+ vec(2,2)*zpsi(2,ik)
    zvec(1) = exp(-zi*0.5*dt*lambda(1))*zvec(1)
    zvec(2) = exp(-zi*0.5*dt*lambda(2))*zvec(2)

    zpsi(1,ik) = vec(1,1)*zvec(1) + vec(1,2)*zvec(2)
    zpsi(2,ik) = vec(2,1)*zvec(1) + vec(2,2)*zvec(2)


! propagation from dt*(it+0.5) to dt*(it+1)
    tt = dt*it

    kt = kn(ik) + Afield_t(it+1)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)
    
    zvec(1) = vec(1,1)*zpsi(1,ik)+ vec(2,1)*zpsi(2,ik)
    zvec(2) = vec(1,2)*zpsi(1,ik)+ vec(2,2)*zpsi(2,ik)
    zvec(1) = exp(-zi*0.5*dt*lambda(1))*zvec(1)
    zvec(2) = exp(-zi*0.5*dt*lambda(2))*zvec(2)

    zpsi(1,ik) = vec(1,1)*zvec(1) + vec(1,2)*zvec(2)
    zpsi(2,ik) = vec(2,1)*zvec(1) + vec(2,2)*zvec(2)

  end do


end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine calc_current(jt_t,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: jt_t
  integer :: ik
  real(8) :: pmat
  real(8) :: tt, kt

  tt = dt*it


  jt_t = 0d0
  do ik = 0, nk-1

    kt = kn(ik) + Afield_t(it)
    pmat = 2d0*t_hop*sin(0.5d0*lattice_a*kt)*0.5d0*lattice_a

    jt_t = jt_t + conjg(zpsi(1,ik))*pmat*zpsi(2,ik) + conjg(zpsi(2,ik))*pmat*zpsi(1,ik)
  end do

  jt_t = jt_t/nk

end subroutine calc_current
!-------------------------------------------------------------------------------
subroutine calc_nex(nex_bloch,nex_houston,nex_pol_houston,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: nex_bloch, nex_houston,nex_pol_houston
  integer :: ik
  real(8) :: kt, dkdt, phi, xx, yy, factor
  real(8) :: ham(2,2), vec(2,2), lambda(2)
  real(8) :: duc_dk(2), uv(2)

! Bloch projection
  nex_bloch = 0d0
  do ik = 0, nk-1
    nex_bloch = nex_bloch &
        + abs(phi_gs(1,1,ik)*zpsi(1,ik)+phi_gs(2,1,ik)*zpsi(2,ik))**2
  end do


! Houston projection
  nex_houston = 0d0
  do ik = 0, nk-1

    kt = kn(ik) + Afield_t(it)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)

    nex_houston = nex_houston + abs(vec(1,1)*zpsi(1,ik)+vec(2,1)*zpsi(2,ik))**2
  end do

! polarized Houston projection
  nex_pol_houston = 0d0
  do ik = 0, nk-1

    kt = kn(ik) + Afield_t(it)
    dkdt = 0.5d0*(Afield_t(it+1)-Afield_t(it-1))/dt
    phi = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    xx =  phi/(0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))
    yy = -phi/(0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))

    duc_dk(1)=-xx/(sqrt(1d0+xx**2))**3*xx + 1d0/sqrt(1d0+xx**2)**3
    duc_dk(2)=-xx/(sqrt(1d0+xx**2))**3 
    factor = 1d0/(0.5d0*delta_gap + sqrt(delta_gap**2/4d0 + phi**2))
    factor = factor - phi**2/( &
        (0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))**2 &
        *sqrt(delta_gap**2/4d0+phi**2) &
        )

    factor = factor *lattice_a*t_hop*sin(0.5d0*lattice_a*kt)

    duc_dk = duc_dk*factor
    

  end do

end subroutine calc_nex
!-------------------------------------------------------------------------------
subroutine init_laser_field
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt, ss

  allocate(Efield_t(-1:nt+1),Afield_t(-1:nt+1))
  Efield_t = 0d0
  Afield_t = 0d0

  do it = 0, nt+1
    tt = dt*it
    ss = (tt - 0.5d0*Tpulse0)
    if(abs(ss)<= 0.5d0*Tpulse0)then
      Afield_t(it) = -(Efield0/omega0)*cos(omega0*ss)*cos(pi*ss/Tpulse0)**4
    end if
        
  end do

  do it = 0, nt
    Efield_t(it) = 0.5d0*(Afield_t(it+1)-Afield_t(it-1))/dt
  end do

end subroutine init_laser_field
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine diag_2x2(mat, vec, lambda)
  implicit none
  real(8),intent(in) :: mat(2,2)
  real(8),intent(out) :: vec(2,2)
  real(8),intent(out) :: lambda(2)
  real(8) :: a, b, c
  real(8) :: ss

  vec = 0d0
  lambda = 0d0

  a = mat(1,1)
  b = mat(1,2)
  c = mat(2,2)
  

  lambda(1) = 0.5d0*((a+c) + sqrt((a-c)**2 + 4d0*b**2)) 
  lambda(2) = 0.5d0*((a+c) - sqrt((a-c)**2 + 4d0*b**2)) 


  if( abs(lambda(1) - a) > abs(lambda(1) - c)  ) then
    vec(2,1) = 1d0
    vec(1,1) = b/(lambda(1)-a)

    vec(1,2) = 1d0
    vec(2,2) = b/(lambda(2)-c)
  else
    vec(1,1) = 1d0
    vec(2,1) = b/(lambda(1)-c)

    vec(2,2) = 1d0
    vec(1,2) = b/(lambda(2)-a)
  end if

  ss = sum(abs(vec(:,1))**2)
  vec(:,1) = vec(:,1)/sqrt(ss)

  ss = sum(abs(vec(:,2))**2)
  vec(:,2) = vec(:,2)/sqrt(ss)

end subroutine diag_2x2
!-------------------------------------------------------------------------------
subroutine diag_2x2_complex(zmat, zvec, lambda)
  implicit none
  complex(8),intent(in) :: zmat(2,2)
  complex(8),intent(out) :: zvec(2,2)
  real(8),intent(out) :: lambda(2)
  real(8) :: a, c
  complex(8) :: zb
  real(8) :: ss

  zvec = 0d0
  lambda = 0d0

  a  = zmat(1,1)
  c  = zmat(2,2)
  zb = zmat(1,2)

  lambda(1) = 0.5d0*((a+c) + sqrt((a-c)**2 + 4d0*abs(zb)**2)) 
  lambda(2) = 0.5d0*((a+c) - sqrt((a-c)**2 + 4d0*abs(zb)**2)) 


  if( abs(lambda(1) - a) > abs(lambda(1) - c)  ) then
    zvec(2,1) = 1d0
    zvec(1,1) = zb/(lambda(1)-a)

    zvec(1,2) = 1d0
    zvec(2,2) = conjg(zb)/(lambda(2)-c)
  else
    zvec(1,1) = 1d0
    zvec(2,1) = conjg(zb)/(lambda(1)-c)

    zvec(2,2) = 1d0
    zvec(1,2) = zb/(lambda(2)-a)
  end if

  
  ss = sum(abs(zvec(:,1))**2)
  zvec(:,1) = zvec(:,1)/sqrt(ss)

  ss = sum(abs(zvec(:,2))**2)
  zvec(:,2) = zvec(:,2)/sqrt(ss)

end subroutine diag_2x2_complex
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
