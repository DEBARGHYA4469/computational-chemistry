program Radial

implicit none 

real :: p0,rx0,rx,rmax,p,E0,R,R0,dr
integer :: np = 10000
real :: s1,s2,s3,s4
real :: k1,k2,k3,k4
real :: e = 0.1,dE = 0.01,conv = 0.1,Eval
integer :: flag = 0,maxConv=0

E0 = -0.6


open(unit=1,file="RadialDistribution")

do
if(E0 > -0.4) exit
flag = 0
conv = e

! all the initial values
p0 = -1000.0
rx0 = 0.0005
R0 = 0.000001
rx=rx0
p=p0
R=R0
rmax = 5.0
dr = (rmax - rx0)/np
! all the initial values

do 
if(rx > rmax) exit 

! RK4 to calculate the p

! k --- r
! p --- s

s1 = dr*fp(p0,rx,R0,E0)  
k1 = dr*fR(p0)  
s2 = dr*fp(p0+s1/2.0,rx+dr/2.0,R0+k1/2.0,E0)
k2 = dr*fR(p0+s1/2.0)
s3 = dr*fp(p0+s2/2.0,rx+dr/2.0,R0+k2/2.0,E0)
k3 = dr*fR(p0+s2/2.0)
s4 = dr*fp(p0+s3,rx+dr,R0+k3,E0)
k4 = dr*fR(p0+s3)

p = p0 + (s1+2*s2+2*s3+s4)/6.0
R = R0 + (k1+2*k2+2*k3+k4)/6.0
 
if(abs(R0 - R0*R0*rx) < e .and. abs(R0-R0*R0*rx) < conv) then
  flag = flag + 1 
  conv = abs(R0-R0*R0*rx) 
else  
flag = 0
conv = 10
end if

if(abs(E0+0.5) < 0.0001)  write(1,*) rx,R0,rx*R0*R0

p0 = p
R0 = R
rx = rx + dr
end do  

if(flag > maxConv) then 
maxConv=flag 
Eval=E0
end if
E0=E0 + dE
end do 
print*,"The correct value of E is : ",Eval
contains 

real function fR(p)
real,intent(in) :: p
fR = p
end function

real function fp(p,rx,R,E0)
real,intent(in) :: p,rx,R,E0
fp = (-2.0*p)/rx - (2.0*(E0 + 1.0/rx)*R)
end function

end program 
