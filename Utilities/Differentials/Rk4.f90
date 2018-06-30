program RK4
implicit none 

real :: t0 = 0,t1 = 0.4,y0 = 0
real :: h,s1,s2,s3,s4 
integer :: i,np = 4

h=0.2
  
do i = 1,2
s1 = h*f(t0,y0) 
s2 = h*f(t0+h/2.0,y0+s1/2.0)
s3 = h*f(t0+h/2.0,y0+s2/2.0)
s4 = h*f(t0+h,y0+s3)
y0 = y0 + (s1+2*s2+2*s3+s4)/6.0
t0 = t0 + h
end do
print*,y0

contains 
real function f(t,y)
real,intent(in) :: t,y 
f=t**2 + y**2
end function 
end program 
