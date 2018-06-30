program Euler
implicit none 

real :: m=1.0,k=1.0,p0=1.0,x0=1.0
integer :: i
real :: t0=0,T,dt,w
real,parameter :: PI = acos(-1.0)
real :: p,x,E

w=sqrt(k/m)
T=2*PI/w
dt = 0.02*T 

open(unit=1,file="Euler")

! Euler method 
do i=1,200

p=p0 + dt*(-1.0*x0)
x=x0 + dt*(p0)
E = (p0**2 + x**2)/2.0
write(1,*) t0/T,x0,p0,2*E

p0=p
x0=x

t0=t0+dt
end do 
 
end program 
