! find the area under the curve exp(x) from 0 to 1
! Actual Area = 1.71828

program Trapezoidal
implicit none 

real :: lb = 0,ub = 1
real :: h,s=0.0
integer :: i,np=1000

h=(ub-lb)/np
s=exp(lb)+exp(ub)
lb = lb + h
do i=2,np
s=s+2*exp(lb)
lb = lb + h
end do
s=s*h/2.0
print*,s
end program 
