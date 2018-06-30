program Rk4
implicit none 


real :: m=1.0,k=1.0,p0=0,x0=1.0
integer :: i
real :: t0=0,T,dt,w,s1,s2,s3,s4
real :: k1,k2,k3,k4
real,parameter :: PI = acos(-1.0)
real :: p,x,E

w=sqrt(k/m)
T=2*PI/w
dt = 0.02*T 

open(unit=1,file="rk4")

! s -- p 
! k -- x

do i=1,200

s1 = fp(x0)
k1 = fx(p0)

s2 = fp(x0+k1/2.0)
k2 = fx(p0+s1/2.0)

s3 = fp(x0+k2/2.0)
k3 = fx(p0+s2/2.0)

s4 = fp(x0+k3)
k4 = fx(p0+s3)

p=p0 + dt*(s1+2*s2+2*s3+s4)/6.0

x=x0 + dt*(k1+2*k2+2*k3+k4)/6.0


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


