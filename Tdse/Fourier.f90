! Fourier Transformation to calculate the time dependent schrodinger wave equation.

program Fourier
implicit none

! constants 
real :: alpha = 20.0, x0 = -0.5 , p0 = 20.0 
real :: m = 14500 , xmin = -2.0 , dx = 0.02, dt = 0.1 
real :: k0 = 20.0 , h = 1,L
integer :: time_steps = 5000 , N = 256

real,parameter :: pi = acos(-1.0)
complex,parameter :: iota=cmplx(0,1.0)
integer :: i,t,j,mj
complex :: hsi,ddif
character(10):: filename

! the x-grid and the psi
real,dimension(:),allocatable :: x,V
complex,dimension(:,:),allocatable :: psi,phi 
real,dimension(:),allocatable :: km

! allocate the grids 
allocate(x(N),V(N),km(N),psi(N,0:time_steps-1),phi(N,0:time_steps-1))

L = N*dx

! the x-grid
do i = 1,N
x(i) = xmin + dx*(i-1)
end do  

! the km-grid
open(unit=1,file="kgrid")
do i = 1,N
 if(i .le. N/2.0) km(i) =2*pi*(i-1)/L
 if(i .gt. N/2.0) km(i) = 2*pi*(i-1-N)/L
 write(1,*) x(i),km(i)
end do
 close(1)

! potential grid 
do i = 1,N
if(x(i) .lt. 0) V(i) = 0.0
if(x(i) .ge. 0) V(i) = 0.1
end do 

! initial wave packet t = 0
do i=1,N 
 psi(i,0) =  sqrt(sqrt(2*alpha/pi))*exp(iota*k0*(x(i)-x0))*exp(-alpha*((x(i)-x0)**2))
end do 

 !wavepacket at time t = 1 using euler method
do i=1,N
 ! Fourier Transformation
  do mj = 1,N
    phi(mj,0) = 0
    do j=1,N
        phi(mj,0) = phi(mj,0) + 1./sqrt(N*1.0)*psi(j,0)*exp(-iota*km(mj)*x(j))
    end do 
  end do
 ! Inverse Fourier Transformation
 ddif=cmplx(0,0)
 do mj=1,N
   ddif = ddif - 1./sqrt(N*1.0)*(km(mj)**2)*phi(mj,0)*exp(iota*km(mj)*x(i))
 end do  
  
 psi(i,1) = psi(i,0) + iota*dt*(h/(2*m)*ddif - V(i)*psi(i,0))
end do


 !wave packet propagation 
do t=2,time_steps-1
  
  ! Fourier Transformation
  do mj = 1 , N
    phi(mj,t-1) = 0
    do j=1,N
        phi(mj,t-1) = phi(mj,t-1) + 1./sqrt(N*1.0)*psi(j,t-1)*exp(-iota*km(mj)*x(j))
    end do 
  end do
 

  do i=1,N
    !Inverse Fourier Transformation
    ddif=cmplx(0,0)
    do mj=1,N
    ddif = ddif - 1./sqrt(N*1.0)*(km(mj)**2)*phi(mj,t-1)*exp(iota*km(mj)*x(i))
    end do 
    hsi = -1.0/(2*m)*ddif + V(i)*psi(i,t-1)
    psi(i,t) = psi(i,t-2) - 2*iota*dt*hsi
  end do 
end do

 ! Writing the output to separate files ... 
t=0
do 
  if(t .ge. time_steps) exit
  open(unit=1,file='util')
  write(1,*) t/100
  close(1)
  open(unit=1,file='util')
  read(1,*) filename
  filename="psift"//filename
  print*,filename
  close(1)
  open(unit=t,file=filename)   
  do i=1,N
	write(t,*) x(i),V(i),real(psi(i,t))**2 + aimag(psi(i,t))**2
  end do
  close(t)
  t = t + 100
end do 
 
end program 

