program schrodinger
implicit none

real :: minx = 0,maxx = 10.0 
integer :: i,j,gp = 1000 
real,parameter :: PI = 3.1412
real :: h= 6.5821 * (10 **(-16.0))
real :: m= 9.1 * (10 ** (-31.0))
real :: delx,tp=0.0,psisq=0.0,minC,maxC
real*8 :: t 
real :: E = 1.0,val,Pavg
complex,allocatable,dimension(:) :: psi,k
real,allocatable,dimension(:) :: x
real,allocatable,dimension(:) :: V
complex :: a,iota
real,dimension(250) :: tr_probability
real :: delE,scalingfactor = 0.0

delE = 0.1
E = 1.0

iota=cmplx(0.0,1.0)

allocate(psi(gp),x(gp),V(gp),k(gp))
open(unit=30,file = "transmission_probability")

do j= 1,250

E = 1.0 + (j-1)*delE

open(unit=10,file="potential")

delx = (maxx-minx)/gp
t = (h**2)/(2*m*delx*delx) ! A constant

! Creating the potential barrier
do i=1,gp
   x(i) = minx + (i-1)*delx 
   if(x(i) .ge. 4.0 .and. x(i) .le. 5.0) V(i) = 9.0
   if(x(i) .lt. 4.0 .or.  x(i) .gt. 5.0) V(i) = 0.0
   write(10,*) x(i),V(i)
end do

! Store the values of the k
do i=1,gp
   val=2*m*(E-V(i))/(h*h)
   if(V(i) .eq. 0) k(i)=cmplx(sqrt(val),0)
   if(V(i) .ge. E) k(i)=cmplx(0,sqrt(val*(-1)))
end do

psi(1)=exp(-iota*k(1)*x(1))
psi(2)=exp(-iota*k(2)*x(2))


! A example plotting
print*,E
if( abs(E-5.0) .eq. 0) then 
 open(unit=50,file="waves")
end if
do i=3,gp 
  psi(i) = (1/t)*(V(i)+ 2*t - E)*psi(i-1) - psi(i-2)
  if(E .eq. 5.0) write(50,*) x(i),real(psi(i))   
end do


! Calculating the maximum and minimum interference...
minC=(real(psi(1))**2 + aimag(psi(1))**2)
maxC=(real(psi(1))**2 + aimag(psi(1))**2)
do i=600,900
  psisq = (real(psi(i))**2 + aimag(psi(i))**2)
  if(psisq > maxC) maxC = psisq 
  if(psisq < minC) minC = psisq
end do 

! Calculating the Pavg
Pavg = (minC + maxC) / 2 
tp = 2.0 / (1.0+Pavg) 
write(30,*) E,tp
 close(10) 
end do 




end program schrodinger
