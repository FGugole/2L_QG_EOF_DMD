! SUBROUTINE that computes the projection of a field w.r.t. another.
! 
! INPUT:
!	- psi : variable w.r.t. which we want to project;
!	- A   : variable that should be projected;
!	- NX  : number of grid points used in x-direction;
!	- NY  : number of grid points used in y-direction.
!
! OUTPUT:
!	- PA  : the resulting projected variable.

SUBROUTINE Projection(PA,psi,A,NX,NY)
	IMPLICIT NONE
	INTEGER :: NX, NY
	REAL(KIND=8), DIMENSION(NX,NY) :: PA, psi, A
	
	INTEGER :: i, j
	REAL(KIND=8) :: res, norm_psi2
	
	! Operator projecting A w.r.t. psi
	res = 0.0d0
	norm_psi2 = 0.0d0
	DO j=1, NY
		DO i=1, NX
			res = res + psi(i,j)*A(i,j)
		END DO
	END DO
	norm_psi2 = sum(psi**2) 
	
	PA = A - (res/norm_psi2)*psi
	
END SUBROUTINE Projection
