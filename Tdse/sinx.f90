program sinx
implicit none

complex*8,dimension(256) ::  x,y
integer n,i,isign

open(unit=10,file="sinxplot")

do i=1,256
x(i) = (i-1)*0.01
end do 

y=x 

call fft(y,256,1)
call fft(y,256,-1)
y=y/256


print*,x(70),y(70)
 



end program 

