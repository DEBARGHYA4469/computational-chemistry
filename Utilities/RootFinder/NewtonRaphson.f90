program Raphson
implicit none 

real :: x0 = 12,x1 
real :: e = 0.0001

do 
x1 = x0 - f(x0)/fderv(x0)
x0 = x1 
if(abs(f(x1))< e) exit
end do 

print*,x1

contains 
real function f(x)
real,intent(in) :: x
f = x**2 - 8*x + 7
end function
real function fderv(x)
real,intent(in) :: x 
real :: h = 0.0001
real :: fd
fd = (f(x+h)-f(x-h))/(2.0*h) 
fderv= fd
end function
end program 
