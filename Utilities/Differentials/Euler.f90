program Euler
implicit none 

real :: h,y0 = 1 
real :: t0 = 0,tf = 1

integer :: i,np = 4 

h = (tf - t0)/np

do i=1,np
y0 = y0 + h*f(t0,y0)
t0 = t0 + h 
end do 

print*,y0

contains 
real function f(t,y)
real,intent(in) :: t,y
f=y
end function
end program
