program deepwell
implicit none

integer n,i,j,k,iplusj,nrot,ifail,lwork

real*8 smatrix, hmatrix, diag, x, res,pi, eigvec
allocatable::smatrix(:,:),hmatrix(:,:),diag(:), eigvec(:,:)

real*8 :: amat,a1mat, hmat,work,kn,km,Vmat,V0=1.0,a=1.0,L=2.0
allocatable:: amat(:,:),a1mat(:,:),hmat(:,:),work(:)

write (*,*) 'Give nr. of basis states'
read (*, *) n

allocate(smatrix(0:2*n-1,0:2*n-1),hmatrix(0:2*n-1,0:2*n-1),diag(0:2*n-1),eigvec(0:2*n-1,0:2*n-1))
allocate(amat(0:2*n-1,0:2*n-1),a1mat(0:2*n-1,0:2*n-1),hmat(0:2*n-1,0:2*n-1),work(64*n*2))

pi = 4.0*datan(1.d0)

do i=0,2*n-1
  do j=0,2*n-1
         kn = (i+1-n)*pi/L
         km = (j+1-n)*pi/L
         smatrix(i,j) = delta(i,j)
         if(i .eq. j) Vmat = -V0*a/L
         if(i .ne. j) Vmat = -(V0/L)*sin((km-kn)*a)/(km-kn)
         hmatrix(i,j) =  -1.0*kn*kn*delta(i,j) + Vmat         
  ENDDO
ENDDO

!diagonalize S-mat first
lwork=64*n*2
call dsyev('v','u',n,smatrix,n,diag,work,lwork,ifail)

do i = 0, 2*n-1
 do j = 0, 2*n-1
  amat(i,j) = smatrix(i,j)/sqrt(diag(j))
 enddo
enddo

hmat=matmul(transpose(amat),matmul(hmatrix,amat))

!diagonalize h-mat 
call dsyev('v','u',n,hmat,n,diag,work,lwork,ifail)

!carry out Av=c
a1mat=matmul(amat,hmat)

! Output the variational eigenvalues to the screen  together with the exact ones
! exact = n^2*pi^2/length^2

write (6,*) 'Variational     Exact'
do i=0, 2*n-1
  write (6,'(5F12.4)') diag(i), ((i+1-n)*pi/L)**2 - V0
enddo

contains 
integer function delta(m,n)
integer,intent(in) :: m,n
integer :: k
if(m .eq. n) k=1 
if(m .ne. n) k=0
delta = k
end function

end

