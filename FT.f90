program main
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8),parameter :: fs=0.024189d0
  integer,parameter :: n_harmonics = 25, nw_harmonic = 100
!  integer,parameter :: nt = 48266, nw = (n_harmonics+1)*nw_harmonic
!  integer,parameter :: nt = 126619, nw = (n_harmonics+1)*nw_harmonic
  integer,parameter :: nt = 63310, nw = (n_harmonics+1)*nw_harmonic

  real(8), parameter ::   omega0 = 0.35424d0/27.2114d0 ! 3.5 mum
  real(8), parameter :: wi = 0d0, wf = (n_harmonics+1)*omega0
  real(8), parameter :: dw = (wf- wi)/nw
  real(8) :: tt(0:nt), Et(0:nt), Dip(0:nt), window(0:nt)
  real(8) :: tpulse0
  integer :: it, iw
  real(8) :: ww, ss, dt
  complex(8) :: zDw
  real(8) :: signal(0:nw), signal_harmonics(n_harmonics)
  integer :: ih

  do it = 0, nt
    read(*,*)tt(it),Et(it), Dip(it)
  end do
  dt = tt(1)-tt(0)

!  tpulse0 = 10d0*2d0*pi/omega0
  Tpulse0 = 0.5d0*pi*(80d0/fs)/acos((0.5d0)**(1d0/8d0))
  window = 0d0
  do it = 0,nt
    ss = tt(it) -0.5d0*tpulse0
    if(abs(ss)<=0.5d0*tpulse0)then
      window(it)= cos(pi*ss/tpulse0)**4
    end if
  end do

  do iw = 0, nw
    ww = wi + dw*iw
    zDw = 0d0
    do it = 0,nt
      zDw = zDw + Dip(it)*exp(zi*ww*tt(it))*window(it)
    end do
    zDw = zDw*dt
    write(*,"(999e26.16e3)")ww,ww**2*abs(zDw)**2
    signal(iw) = ww**2*abs(zDw)**2

  end do

  do ih = 1, n_harmonics, 2
    signal_harmonics(ih)= 0d0
    do iw = (ih-1)*nw_harmonic,(ih+1)*nw_harmonic
      signal_harmonics(ih)= signal_harmonics(ih) + signal(iw)
    end do
    signal_harmonics(ih)= signal_harmonics(ih)*dw
  end do

  open(31,file='hhg_intensity.out')
  write(31,"(999e26.16e3)")signal_harmonics(1) &
                          ,signal_harmonics(3) &
                          ,signal_harmonics(5) &
                          ,signal_harmonics(7) &
                          ,signal_harmonics(9) &
                          ,signal_harmonics(11) &
                          ,signal_harmonics(13)
  close(31)

end program main
