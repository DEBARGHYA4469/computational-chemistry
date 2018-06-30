! code for chebychev interpolation
! fx = Pnx = sum(0,N) ckTkx 
! ck = 2/n+1 sum(0,N) fj*Tk xj
! calculate chebychev's nodes 
! Tj x = 2*x*Tj-1  - Tj-2  1,x

program Chebychev
implicit none 

real,parameter :: PI = acos(-1.0)
real,dimension(0:3) :: x,f 
integer :: i,j,k,N=3
real :: lb = 0 , ub = 10 ,dx = 0.01
real :: s

open(unit=1,file ="chebychev_plot.dat")

! chebychev nodes 
do i=0,N
x(i) = cos((2*i+1)*PI/(2*N + 2))
f(i) = exp(x(i))
print*,x(i),f(i)
end do 

do
if(lb > ub) exit
s = 0.0

do k=0,N
s = s + C(k,x,f,N)*T(k,x,lb,N) ! summation Ck*Tk(x) 
end do
write(1,*) lb,s,exp(lb),T(0,x,lb,N),T(1,x,lb,N),T(2,x,lb,N),T(3,x,lb,N)
lb =lb + dx
end do 


contains 

real function T(k,x,xval,N)
integer,intent(in) :: k,N
real,dimension(0:) :: x
real,intent(in) :: xval
real :: lb
real,allocatable,dimension(:) :: Tval
integer :: i,j
allocate(Tval(0:k))
Tval(0) = 1 
Tval(1) = xval
do i=2,k
Tval(i) = 2*xval*Tval(i-1) - Tval(i-2)
end do  
T = Tval(k)
end function

real function C(j,x,f,N)
integer,intent(in) :: j,N
real,dimension(0:),intent(in) :: x,f
integer :: k
real :: Cval

 Cval = 0.0
do k = 0,N
 Cval = Cval  + f(k)*T(j,x,x(k),N)
end do
if(j .eq. 0) Cval = Cval * (1.0/(N+1.0))  
if(j .ne. 0) Cval = Cval * (2.0/(N+1.0)) 
 C =Cval 
end function
end program
