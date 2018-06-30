program Bisection
implicit none 

real :: lb = 2 ,ub = 10 ,mid,root 
real :: e = 0.0001

do 
mid = (lb + ub)/2
if(f(mid)*f(lb) > 0) lb=mid
if(f(mid)*f(ub) > 0) ub=mid
if(abs(f(mid))< e) exit
if(abs((ub-lb)/lb) < e)  exit
end do 

print*,mid

contains 
real function f(x)
real,intent(in) :: x
f = x**2 - 8*x + 7
end function
end program 
