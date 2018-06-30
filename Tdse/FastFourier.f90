! Fourier Transformation to calculate the time dependent schrodinger wave equation.

program Fourier
implicit none

! constants 
real*8 :: alpha = 20.0, x0 = -0.5 , p0 = 50.0 
real*8 :: m = 14500 , xmin = -2.0 , dx = 0.02, dt = 0.1 
real*8 :: k0 = 50.0 , h = 1,L
integer :: time_steps = 5000 , N = 256

real*8,parameter :: pi = acos(-1.0)
complex*16,parameter :: iota=cmplx(0,1.0)
integer :: i,t,j,mj
complex*8 :: hsi,ddif
character(10):: filename
complex*8,dimension(256) :: dummy

! the x-grid and the psi
real*8,dimension(:),allocatable :: x,V
complex*8,dimension(:,:),allocatable :: psi,phi 
real*8,dimension(:),allocatable :: km

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
if(x(i) .ge. 0 .and. x(i) .le. 0.5) V(i) = 0.1
if(x(i) .gt. 0.5) V(i) = 0.0
end do 

! initial wave packet t = 0
do i=1,N 
 psi(i,0) =  sqrt(sqrt(2*alpha/pi))*exp(iota*k0*(x(i)-x0))*exp(-alpha*((x(i)-x0)**2))
end do 



 !wave packet propagation 
do t=1,time_steps-1
  
  ! Fourier Transformation
  ! do mj = 1 , N
    ! phi(mj,t-1) = 0
      do j=1,N
      !  phi(mj,t-1) = phi(mj,t-1) + 1./sqrt(N*1.0)*exp(-iota*V(j)/2.0*dt)*psi(j,t-1)*exp(-iota*km(mj)*x(j))
	dummy(j) = psi(j,t-1)*exp(-iota*V(j)/2.0*dt)            
      end do 
        call fft(dummy,N,1)
	
  ! end do
 

  ! do i=1,N
    !Inverse Fourier Transformation
    !ddif=cmplx(0,0)
    do mj=1,N
    dummy(mj)=dummy(mj)*exp(-iota*km(mj)*km(mj)/(2*m)*dt)
    end do
    call fft(dummy,N,-1)
    dummy=dummy/N 
   do i=1,N
    psi(i,t) = exp(-iota*V(i)/2.0*dt)*dummy(i)
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
  filename="psiso"//filename
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

