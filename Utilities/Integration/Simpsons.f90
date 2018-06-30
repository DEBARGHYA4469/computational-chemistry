program Simpsons
implicit none

real :: lb = 0, ub = 1 , h
integer :: i,np = 1000
real :: s = 0.0

s = exp(lb) +  exp(ub)
h= (ub-lb)/np

lb = lb + h
do i=2,np
if(MOD(i,2) .eq. 0) s = s + 4*exp(lb)
if(MOD(i,2) .eq. 1) s = s + 2*exp(lb)
lb=lb+h
end do 

s = h*s/3.0
print*,s 


end program 

