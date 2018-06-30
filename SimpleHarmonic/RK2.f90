program RK2
implicit none 

real :: m=1.0,k=1.0,p0=0,x0=1.0
integer :: i
real :: t0=0,T,dt,w,sleft,sright,xnext,pnext
real :: kleft,kright
real,parameter :: PI = acos(-1.0)
real :: p,x,E

w=sqrt(k/m)
T=2*PI/w
dt = 0.02*T 

open(unit=1,file="rk2")

! sleft -- p
! kleft -- x

! Euler method 
do i=1,200

sleft = fp(x0)
kleft = fx(p0)
xnext = x0 + dt*fx(p0)
pnext = p0 + dt*fp(x0)
sright = fp(xnext)
kright = fx(pnext)
p=p0 + dt*(sleft+sright)/2.0
x=x0 + dt*(kleft+kright)/2.0


E = (p0**2 + x0**2)/2.0
write(1,*) t0/T,x0,p0,2*E

p0=p
x0=x

t0=t0+dt
end do 
contains
real function fp(x)
real,intent(in) :: x
fp = -1.0*x
end function

real function fx(p)
real,intent(in) :: p
fx = p
end function

end program 
