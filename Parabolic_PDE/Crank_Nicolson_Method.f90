PROGRAM CN
integer, parameter :: n=20, E=1, Re=5000, t=360 ! the number of time steps
integer N1,T1, i, j
real dy, dt
real,dimension(n+1,t+1) :: u
real,dimension(n-1,1) :: A,B,C
real,dimension(n-1,t+1) :: K
real,dimension(n-1) :: x
real,dimension(1:7,1:n+1) :: uu

N1=n+1
T1=t+1
dy=1/n
dt=E*Re*(dy**2)

! Initial Condition
do i=1,n
	u(i,1)=0
enddo
u(N1,1)=1

! Coefficients
do i=1,n-1
    A(i,1)=(-0.5)*E                                              !   9.37a
    B(i,1)=1.0+E                                                 !   9.37b  
    C(i,1)=(-0.5)*E                                                
    K(i,1)=(1.0-E)*u(i+1,1)+(0.5*E)*(u(i+2,1)+u(i,1))            !   9.37c
    if(i==n-1) then
        K(i,1)=K(i,1)-C(i,1)
    endif
enddo

do j=2,T1
	u(1,j)=0
	u(N1,j)=1
	call thomas(A,B,C,K(:,j-1),x,n-1)
	do i=2,n
		u(i,j)=x(i-1)
	enddo
	do i=1,n-1
        K(i,j)=(1.0-E)*u(i,j)+(0.5*E)*(u(i+2,j)+u(i,j))
        if (i==n-1) then
            K(i,j)=K(i,j)-C(i,1)
        endif
    enddo
enddo

!   Print
    do i=1,N1
        uu(1,i)=i
        uu(2,i)=u(i,13)
        uu(3,i)=u(i,37)
        uu(4,i)=u(i,61)
        uu(5,i)=u(i,121)
        uu(6,i)=u(i,191)
        uu(7,i)=u(i,361)
    enddo

    write(*,*) uu(:,1:N1)

OPEN (20, FILE="u.txt")
write (20,'(21f10.4)', ADVANCE='NO') u

contains
subroutine thomas(a,b,c,r,x,n)
implicit none

!    a - sub-diagonal (means it is the diagonal below the main diagonal)
!    b - the main diagonal
!    c - sub-diagonal (means it is the diagonal above the main diagonal)
!    r - right part
!    x - the answer
!    n - number of equations

integer,intent(in) :: n
real,dimension(n),intent(in) :: a,b,c,r
real,dimension(n),intent(out) :: x
real,dimension(n) :: cp,rp
real :: tmp
integer i

! initialize c-prime and r-prime
        cp(1) = c(1)/b(1)
        rp(1) = r(1)/b(1)

        

! solve for vectors c-prime and r-prime
         do i = 2,n
           tmp = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/tmp
           rp(i) = (r(i)-rp(i-1)*a(i))/tmp
         enddo

         

! initialize x
         x(n) = rp(n)

! solve for x from the vectors c-prime and r-prime
        do i = n-1, 1, -1
          x(i) = rp(i)-cp(i)*x(i+1)
        end do

end subroutine thomas



END PROGRAM CN
