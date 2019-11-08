! Computes Dynamic Mode Decomposition (DMD) of X1, X2
!	FORTRAN code based on the MATLAB code written by N. Kutz at al as complementary material 
!	to the book ' Dynamic Mode Decomposition: Data-Driven Modelling of Complex System ', SIAM (2016)
!	These codes are publicly available on the following website 
!		https://faculty.washington.edu/kutz/page26/
!	Last accessed: November 6th, 2019

! --------------------------------------------------------------------------------------------------- !

! INPUTS: 
! 	- X1 	: X, data matrix
! 	- X2 	: X', shifted data matrix
! 		Columns of X1 and X2 are state snapshots 
! 	- m 	: amount of snapshots contained in X1 and in X2
! 	- r 	: target rank of SVD
! 	- dt 	: time step advancing X1 to X2 (X to X')
! 	- NX*NY : dimension of snapshots when rearranged as an array (instead of as a matrix)
! 	!- Dt 	: time step advancing X1(:,n) to X1(:,n+1) 
! 		Note that in principle dt!=Dt 
!
! OUTPUTS:
! 	- Phi 	: the DMD modes
! 	- omega : the continuous-time DMD eigenvalues
! 	- lambda: the discrete-time DMD eigenvalues
! 	!- b		: a vector of magnitudes of modes Phi
! 	!- Xdmd	: the data matrix reconstrcted by Phi, omega, b


! Remember that DMD likes tall and skinny matrices !
! This subroutine makes use of the LAPACK library; be sure to correctly link the code to 
! the LAPACK library when compiling

SUBROUTINE DMD(Phi,omega,lambda,X1,X2,m,r,dt,NX,NY)
	IMPLICIT NONE
	INTEGER, INTENT(IN)			:: m, r, NX, NY
	REAL(KIND=8), INTENT(IN)	:: dt
	REAL(KIND=8), DIMENSION(NX*NY,m), INTENT(IN)	:: X1, X2
	COMPLEX(KIND=8), DIMENSION(r)					:: omega, lambda
	COMPLEX(KIND=8), DIMENSION(NX*NY,r)				:: Phi
	
	REAL(KIND=8), DIMENSION(m)			:: S
	REAL(KIND=8), DIMENSION(r) 			:: S_r, DI, DR
	REAL(KIND=8), DIMENSION(NX*NY,m)	:: U
	REAL(KIND=8), DIMENSION(NX*NY,r) 	:: U_r
	REAL(KIND=8), DIMENSION(m,m)		:: V
	REAL(KIND=8), DIMENSION(m,r)		:: V_r
	REAL(KIND=8), DIMENSION(NX*NY,r)	:: X2VS_r
	REAL(KIND=8), DIMENSION(r,r)		:: Atilde, VR_r
	COMPLEX(KIND=8), DIMENSION(r,r)		:: WR_r
	REAL(KIND=8), DIMENSION(2*NX*NY+m)	:: WORK
	REAL(KIND=8), DIMENSION(8*r)		:: WORK_eig
	COMPLEX(KIND=8), DIMENSION(1)		:: WORK_b
	INTEGER, DIMENSION(NX*NY+3*m)		:: IWORK
	REAL(KIND=8)						:: WL
	INTEGER :: ok, i
	

	! Compute rank-r SVD
	CALL DGEJSV('F','U','V','R','N','N',NX*NY,m,X1,NX*NY,S,U,NX*NY,V,m,WORK,2*NX*NY+m,IWORK,ok)
	WRITE(*,*) 'Reduced rank SVD. Success? ', ok
	U_r = U(:,1:r)
	S_r = S(1:r)*( WORK(2)/WORK(1) )
	V_r = V(:,1:r)
		
	X2VS_r 	= matmul(X2,V_r)
	DO i=1, r
		X2VS_r(:,i) = X2VS_r(:,i)/S_r(i)
	END DO
	Atilde 	= matmul(transpose(U_r),X2VS_r)

	! Compute eigenvalues and eigenvectors of Atilde
	CALL DGEEV('N','V',r,Atilde,r,DR,DI,WL,1,VR_r,r,WORK_eig,8*r,ok)
	WRITE(*,*) 'Eigenvalues and eigenvectors of Atilde. Success? ', ok
	! Form complex eigenvectors
	WR_r = 0.0d0
	DO i=1, r
		IF(abs(WR_r(1,i)) < 1e-16) THEN
			! not yet defined complex eigenvector
			IF(DR(i)-DR(i+1) < 1e-16) THEN
				! Complex conjugate pair of eigenvalues
				WR_r(:,i)	= cmplx(VR_r(:,i),+VR_r(:,i+1))
				WR_r(:,i+1)	= cmplx(VR_r(:,i),-VR_r(:,i+1))		
			ELSE 
				! No complex conjugate pair
				WR_r(:,i)	= cmplx(VR_r(:,i),0.0d0)
			END IF
		END IF
	END DO
	! Form complex eigenvalues
	lambda	= cmplx(DR,DI)
	omega	= log(lambda)/dt

	! 
	Phi		= matmul(X2VS_r,WR_r) 
	DO i=1, r
		! Multiply by the inverse of \Lambda 
		Phi(:,i) = Phi(:,i)/lambda(i) ! DMD modes
	END DO
	
END SUBROUTINE 
