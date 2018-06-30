program False
implicit none 

real :: x1 = 2 ,x2 = 10 ,x3 
real :: e = 0.0001

do 
x3 = x1 - (x2-x1)*f(x1)/(f(x2)-f(x1))
if(f(x3)*f(x1) > 0) x1=x3
if(f(x3)*f(x2) > 0) x2=x3
if(abs(f(x3))< e) exit
if(abs((x2-x1)/x1) < e)  exit
end do 

print*,x3

contains 
real function f(x)
real,intent(in) :: x
f = x**2 - 8*x + 7
end function
end program 
