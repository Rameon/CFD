PROGRAM c

implicit none
!-------------------------------------------------------------------------------!
!---------------------Hoffmann CFD Chapter 3.5 Application----------------------!
!-----------Hoffmann CFD APPENDIX B : Tridiagonal System of Equations-----------!
!-------------------Parabolic Partial Differential Equations--------------------! 
!-------------------------The Laasonen implicit method--------------------------!
!-------------------------------------------------------------------------------!
real, parameter :: vis = 0.000217              ! [m^2/s] a kinematic viscosity
real, parameter :: h = 0.04                    ! [m] the spacing between plates
real, parameter :: U0 = 40                     ! [m/s] the velocity of the lower wall
real, parameter :: dy = 0.001                  ! [m] y-direction steps
real, parameter :: dt = 0.002                  ! [s] time steps
integer :: n, j
integer, parameter :: NM=541, JM=41
real, dimension(NM,JM) :: u
real :: a, b, c, d
real, dimension(JM-2) :: aa, ab, ac, ad
real,dimension(JM-2) :: X

! Initial Conditions
do j=1,JM
	u(1,j)=0        ! u=0 (0<y<=h)
enddo
u(1,1)=U0           ! u=U0 (y=0)

! Boundary Conditions
do n=2,NM
	u(n,1)=U0
	u(n,JM)=0
enddo

! Coefficients
d=vis*dt/(dy*dy)
a=d
b=-(2*d+1)
c=d
do j=1,JM-2
	aa(j)=d
	ab(j)=b
	ac(j)=d
	aD(j)=-u(1,j+1)
	aD(1)=-u(1,2)-a*u(1,1)
	aD(JM-2)=-u(1,JM-1)-c*u(1,JM)
enddo

! Calculation
do n=2,NM
	call thomas(aa,ab,ac,ad,x,JM-2)
	do j=2,JM-1
		u(n,j)=x(j-1)
	enddo
	do j=1,JM-2
		aD(j)=-u(n,j+1)
		aD(1)=-u(n,2)-aa(1)*u(n,1)
		aD(JM-2)=-u(n,JM-1)-ac(j)*u(n,JM)
	enddo
enddo

! Datafile
OPEN(20, FILE="c.txt")
write(20,'(41f10.4)') transpose(u)

contains
subroutine thomas(a,b,c,r,x,nn)
implicit none

!    a - sub-diagonal (means it is the diagonal below the main diagonal)
!    b - the main diagonal
!    c - sub-diagonal (means it is the diagonal above the main diagonal)
!    r - right part
!    x - the answer
!    n - number of equations

integer,intent(in) :: nn
real,dimension(nn),intent(in) :: a,b,c,r
real,dimension(nn),intent(out) :: x
real,dimension(nn) :: cp,rp
real :: tmp
integer i

! initialize c-prime and r-prime
        cp(1) = c(1)/b(1)
        rp(1) = r(1)/b(1)

! solve for vectors c-prime and r-prime
         do i = 2,nn
           tmp = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/tmp
           rp(i) = (r(i)-rp(i-1)*a(i))/tmp
         enddo

! initialize x
         x(nn) = rp(nn)

! solve for x from the vectors c-prime and r-prime
        do i = nn-1, 1, -1
          x(i) = rp(i)-cp(i)*x(i+1)
        end do

end subroutine thomas

END PROGRAM c
