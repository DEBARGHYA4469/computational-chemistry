program Kinetics
implicit none 

real :: A0,B0,C0,D0
real :: t=0.0,tmax = 60.0 , dt = 0.0001
real :: A,B,C,D
real :: s1,s2,s3,s4

 A0=5.0
 B0=0.0
 C0=0.0
 D0=0.0
open(unit=1,file="concentrations")
do 
if(t > tmax) exit


! concentration of A
s1=dt*fA(A0,B0) 
s2=dt*fA(A0+s1/2.0,B0+s1/2.0)
s3=dt*fA(A0+s2/2.0,B0+s2/2.0)
s4=dt*fA(A0+s3,B0+s3)

 A = A0 + (s1+2*s2+2*s3+s4)/6.0

! concentration of B

s1=dt*fB(A0,B0,C0) 
s2=dt*fB(A0+s1/2.0,B0+s1/2.0,C0+s1/2.0)
s3=dt*fB(A0+s2/2.0,B0+s2/2.0,C0+s2/2.0)
s4=dt*fB(A0+s3,B0+s3,C0+s3)

 B = B0 + (s1+2*s2+2*s3+s4)/6.0

! concentration of C
s1=dt*fC(B0,C0) 
s2=dt*fC(B0+s1/2.0,C0+s1/2.0)
s3=dt*fC(B0+s2/2.0,C0+s2/2.0)
s4=dt*fC(B0+s3,C0+s3)

 C = C0 + (s1+2*s2+2*s3+s4)/6.0

! concentration of D

s1=dt*fD(C0) 
s2=dt*fD(C0+s1/2.0)
s3=dt*fD(C0+s2/2.0)
s4=dt*fD(C0+s3)

 D = D0 + (s1+2*s2+2*s3+s4)/6.0

write(1,*) t,A0,B0,C0,D0

 A0 = A
 B0 = B
 C0 = C
 D0 = D
t = t + dt 
end do 

contains 
real function fA(A,B)
real :: A,B
fA = -1.0*A + 3.0*B
end function

real function fB(A,B,C)
real :: A,B,C
fB = A - 3.0*B - 4.2*B + 7.3*C*C
end function

real function fC(B,C)
real :: B,C
fC = 2*4.2*B - 2*7.3*C*C - 0.4*C
end function

real function fD(C)
real :: C
fD = 0.4*C
end function
end program 
