program main
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)
  integer,parameter :: nt = 96531, nw = 256
  real(8), parameter :: wi = 0d0, wf = 6d0/27.2114
  real(8), parameter :: dw = (wf- wi)/nw
  real(8) :: tt(0:nt), Et(0:nt), Dip(0:nt), window(0:nt)
  real(8) :: tpulse0, omega0
  integer :: it, iw
  real(8) :: ww, ss, dt
  complex(8) :: zDw

  do it = 0, nt
    read(*,*)tt(it),Et(it), Dip(it)
  end do
  dt = tt(1)-tt(0)

  omega0 = 0.35424d0/27.2114d0 
  tpulse0 = 10d0*2d0*pi/omega0
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

  end do


end program main
