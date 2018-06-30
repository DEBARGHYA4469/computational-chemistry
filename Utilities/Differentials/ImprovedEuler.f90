! RK2
program ImprovedEuler
implicit none 
real :: h,t0 = 0,t1 = 1, y0 = 1 
integer :: i,np = 4
real :: sleft,sright,yeuler

h = (t1 - t0)/np 

do i=1,np
sleft  = f(t0,y0)
t0 = t0 + h
yeuler = y0 + h*f(t0-h,y0)
sright = f(t0,yeuler)
y0 = y0 + h*(sleft+sright)/2.0
end do 

print*,y0

contains
real function f(t,y)
real,intent(in) :: t,y
f = y
end function
end program 
