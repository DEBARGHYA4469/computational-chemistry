! code to implement lagranges interpolation method to interpolate a set of data points for approximating the curve square root of n
program Lagrange 
implicit none

real,dimension(5) :: x
real,dimension(5) :: f
integer :: i
real :: s = 0 , lb = 0 , ub = 10 , dx = 0.01

DATA  x/1,2,3,4,5/
DATA  f/1,1.4142,1.7321,2,2.2361/

open(unit = 1 , file="lagrange_plot.dat")

! calculate P4
do
if(lb > ub) exit
s = 0 
do i=1,5
 s = s + calc_li(i,lb,x,5)*f(i)
end do 
write(1,*) lb,s,sqrt(lb) 
lb = lb + dx 
end do 
contains 
real function calc_li(i,xval,x,N)
real,dimension(5),intent(in) :: x 
integer,intent(in) :: i,N
real,intent(in) :: xval
integer :: j
real :: p
p=1.0
do j=1,N
if(i .eq. j) cycle  
p = p * (xval-x(j))/(x(i)-x(j))
end do
calc_li = p
end function
end Program
