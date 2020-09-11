PROGRAM b
implicit none
!-------------------------------------------------------------------------------!
!---------------------Hoffmann CFD Chapter 3.5 Application----------------------!
!-------------------Parabolic Partial Differential Equations--------------------! 
!----------------------The DuFort-Frankel explicit method-----------------------!
!-------------------------------------------------------------------------------!
real, parameter :: vis = 0.000217              ! [m^2/s] a kinematic viscosity
real, parameter :: h = 0.04                    ! [m] the spacing between plates
real, parameter :: U0 = 40                     ! [m/s] the velocity of the lower wall
real, parameter :: dy = 0.001                  ! [m] y-direction steps
real, parameter :: dt = 0.002                  ! [s] time steps
integer :: n, j
integer, parameter :: NM=541, JM=41
real, dimension(NM,JM) :: u

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

! Calculation
do n=2,NM-1
	do j=2,JM-1
		u(n+1,j) = ((1-2*vis*dt/(dy*dy))*u(n-1,j) &
					+ 2*vis*dt/(dy*dy)*(u(n,j+1)+u(n,j-1))) / (1+2*vis*dt/(dy*dy)) ! (3-12)
	enddo
enddo

! Datafile
OPEN(20, FILE="b.txt")
write(20,'(41f10.4)') transpose(u)

END PROGRAM b