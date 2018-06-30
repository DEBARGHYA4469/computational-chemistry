! Use Newton interpolation to interpolate the function 
! f(x) = sqrt(x)

program Newton
implicit none 
real :: lb = 0 , ub = 10 ,dx = 0.01
real,dimension(5) :: x,f,A
real,dimension(5,5) :: D
integer :: i,j,N = 5
real :: p = 1.0,s=0

Data x/1,2,3,4,5/
Data f/1,1.4142,1.732,2,2.2361/

open(unit=1,file="Newton_plot.dat")


! Newton Interpolation by divided difference
do i=1,N
D(i,1) = f(i)
end do 

do j=2,N
do i=1,N-j+1
D(i,j) = (D(i+1,j-1)-D(i,j-1))/(x(i+j-1) - x(i))
end do
end do

do i=1,N
A(i)=D(1,i)
end do

do 
if(lb > ub) exit
p=1.0
s=0
do i=1,5
s = s + p*A(i)
p=p*(lb - x(i))
end do
write(1,*) lb,s,sqrt(lb)
lb = lb + dx
end do 

end program 
