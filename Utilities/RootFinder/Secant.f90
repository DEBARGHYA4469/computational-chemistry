program Secant
implicit none 

real :: x1 = 12,x2,x3 
real :: e = 0.0001

do 
x3 = x2 - f(x2)*(x2 - x1)/(f(x2)-f(x1)) 
x1=x2
x2=x3
if(abs(f(x3))< e) exit
end do 

print*,x3

contains 
real function f(x)
real,intent(in) :: x
f = x**2 - 8*x + 7
end function
end program 
