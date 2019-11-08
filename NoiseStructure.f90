! SUBROUTINE that builds the noise structure using the eigenvectors computed with EOF
! 
! INPUT FILES: these files contain the eigenvectors and the eigenvalues 
! 	The name of the file should be: 
!		- eigenvec_psiT_x.dat
!		- eigenval_psiT_x.dat
!	Where x is a number of three digits corresponding on the grid point used in x-direction. For instance, 
!	if NX=64 the file name for the eigenvectors will be eigenvec_psiT_064.dat . 
!	Possible values for NX are: 128 
! 
! INPUTS:
!	- NX		: number of grid points in x-direction
!	- NY		: number of grid points in y-direction
!	- nstart	: number of the first eigenvector to be considered
!	- neigen	: total number of eigenvectors to be used
! 	Note that the code has been written having in mind that NX==NY. Two distinct variables are used to read more easily the code and 
!	to allow a easier extension to the case NX!=NY .
! 
! OUTPUT: 
!	- CSigma : spatio-temporal structure of the noise.

SUBROUTINE EOFStructure(Csigma,NX,NY,nstart,neigen)
	USE mtmod
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)									:: NX, NY, nstart, neigen
	REAL(KIND=8), DIMENSION(NX,NY), INTENT(OUT) 		:: Csigma 					! Resulting combination of EOFs in matrix form
	INTEGER												:: i, j, counter
	CHARACTER(LEN=200)									:: InputFile 				! Variable that will contain the name of the input files	
	REAL(KIND=8), DIMENSION(nstart+2*neigen-1)			:: alpha 					! Weights of the convex combination of EOFs. Note that their sum must be 1 !
	REAL(KIND=8), DIMENSION(NX*NY)						:: comb_EOF 				! Convex combination of EOFs. Note: it is an array !
	REAL(KIND=8), DIMENSION(NX*NY,nstart+2*neigen-1)	:: EOFs 					! Matrix containing all EOFs. Note that they are sorted by column and each column has dimension NX*NY
	REAL(KIND=8) 										:: tmp1, tmp2
	
	WRITE(InputFile,'(a,i3.3,a)') 'eigenvec_psiT_', NX, '.dat'

	! Read EOFs
	EOFs = 0.0d0
	OPEN(UNIT=444, FILE=InputFile, RECL=8000)
		READ(444,*) 
		DO i=1, NX*NY
			READ(444,*) EOFs(i,1:nstart+2*neigen-1)
		END DO
	CLOSE(444)
	
	! Set weights equal to EOFs eigenvalues
	WRITE(InputFile,'(a,i3.3,a)') 'eigenval_psiT_', NX, '.dat'
	alpha = 0.0d0
	OPEN(UNIT=444, FILE=InputFile, RECL=8000)
		READ(444,*) 
		DO i=1, nstart+neigen-1
			READ(444,*) tmp1, alpha(i), tmp2
		END DO
	CLOSE(444)
	
	! Compute the noise structure 
	comb_EOF = 0.0d0
	DO j=nstart, nstart+neigen-1 
		comb_EOF = comb_EOF + sqrt(alpha(j))*EOFs(:,j)
	END DO
	! Now remap the resulting array (NX*NY,1) to be a matrix (NX,NY)
	! Remap according to how the data (on which the EOFs have been computed) had been written originaly (see IO_restart.f90)
	counter = 0
	DO j=1, NY
		DO i=1, NX
			counter = counter + 1
			Csigma(i,j) = comb_EOF(counter)
		END DO
	END DO
	WRITE(*,*) 'CSigma computed with EOF.'
	
END SUBROUTINE EOFStructure




! SUBROUTINE that builds the noise structure using the eigenmodes computed with DMD
! 
! INPUTS:
!	- NX		: number of grid points in x-direction
!	- NY		: number of grid points in y-direction
!	- nstart	: number of the first eigenvector to be considered
!	- neigen	: total number of eigenvectors to be used
!	- m			: number of snapshots to be used for DMD
!	- r			: raduced rank desired
!	- dt		: time interval between X1(:,n) and X2(:,n)
!	- X1, X2	: X and X' respectively when following the notation in the DMD book of Kutz et al. 
! 	Note that the code has been written having in mind that NX==NY. Two distinct variables are used to read more easily the code and 
!	to allow a easier extension to the case NX!=NY .
! 
! OUTPUT: 
!	- CSigma : spatio-temporal structure of the noise.

SUBROUTINE DMDStructure(Csigma,NX,NY,nstart,neigen,m,r,dt,X1,X2)
	USE mtmod
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)							:: NX, NY, nstart, neigen, m, r
	REAL(KIND=8), INTENT(IN)					:: dt
	REAL(KIND=8), DIMENSION(NX*NY,m)			:: X1, X2
	REAL(KIND=8), DIMENSION(NX,NY), INTENT(OUT) :: Csigma 						! Resulting combination of DMDs in matrix form
	COMPLEX(KIND=8), DIMENSION(r)				:: omega, lambda, lambda_copy
	COMPLEX(KIND=8), DIMENSION(NX*NY,r)			:: Phi
	INTEGER, DIMENSION(r)						:: lambda_order
	INTEGER										:: i, j, counter, idx
	CHARACTER(LEN=200)							:: InputFile 					! Variable that will contain the name of the input files	
	REAL(KIND=8), DIMENSION(NX*NY,1)			:: comb_DMD 					! Convex combination of DMDs. Note: it is an array !
	REAL(KIND=8), DIMENSION(NX*NY,NX*NY)		:: WWT							! Covariance matrix given by comb_DMD before being scaled to have the same strenght as EOFs
	REAL(KIND=8), DIMENSION(NX*NY,r)			:: RDMDs, IDMDs					! Matrix containing all DMDs. Note that they are sorted by column and each column has dimension NX*NY
	REAL(KIND=8)								:: tmp1, tmp2, cum_var, trWWT
	
	! Compute DMDs
	RDMDs = 0.0d0
	IDMDs = 0.0d0
	CALL DMD(Phi,omega,lambda,X1,X2,m,r,dt,NX,NY)
	!!! Ensure to consider only eigenvalues on the unit circle and not the spurious ones outside
	DO i=1, r
		IF(abs(abs(lambda(i))-1)<5e-2) THEN
			lambda_copy(i) = lambda(i)
		ELSE
			lambda_copy(i) = -90
		END IF
	END DO
	!!! Reorder the DMD modes such that the first column corresponds to the mean mode, and then smaller and smaller scales
	DO i=1, r
		idx					= minloc(abs(lambda_copy-cmplx(1,0)),DIM=1)
		lambda_order(i)		= idx
		lambda_copy(idx)	= -99
		RDMDs(:,i)			= REAL(Phi(:,idx),KIND=8)
		IDMDs(:,i)			= DIMAG(Phi(:,idx))
	END DO
	
		
	! According to which EOFs have been used store their cumulated variance 
	! in order to set the DMD covariance matrix to have the same strenght as for the EOF case
	WRITE(InputFile,'(a,i3.3,a)') 'eigenval_psiT_', NX, '.dat'
	!WRITE(*,'(a,a,a,a)' ) ' Input file ', TRIM(InputFile)
	OPEN(UNIT=444, FILE=InputFile, RECL=8000)
		READ(444,*) ! first line = header
		READ(444,*) 
		READ(444,*) 
		READ(444,*) 
		READ(444,*)
		READ(444,*) cum_var, tmp1, tmp2
	CLOSE(444)
	
	! Compute the noise structure 
	comb_DMD = 0.0d0
	DO j=nstart, nstart+2*neigen-1, 2 
		comb_DMD(:,1) = comb_DMD(:,1) + ( REAL(lambda(lambda_order(j)))*RDMDs(:,j) - DIMAG(lambda(lambda_order(j)))*IDMDs(:,j) )
	END DO
	! Normalize it such that it has the same strength as with EOFs
	WWT = matmul(comb_DMD,transpose(comb_DMD))
	trWWT = 0.0d0
	DO i=1, NX*NY
		trWWT = trWWT + WWT(i,i)
	END DO
	comb_DMD = comb_DMD*sqrt(cum_var/trWWT)
	! Check that now it has the desired amount of variance
	!WWT = matmul(comb_DMD,transpose(comb_DMD))
	!trWWT = 0.0d0
	!DO i=1, NX*NY
	!	trWWT = trWWT + WWT(i,i)
	!END DO
	!WRITE(*,*) 'Trace of covariance matrix = ', trWWT
	
	! Now remap the resulting array (NX*NY,1) to be a matrix (NX,NY)
	! Remap according to how the data (on which the EOFs have been computed) had been written originaly (see IO_restart.f90)
	counter = 0
	DO j=1, NY
		DO i=1, NX
			counter = counter + 1
			Csigma(i,j) = comb_DMD(counter,1)
		END DO
	END DO
	WRITE(*,*) 'CSigma computed with DMD.'
	
END SUBROUTINE DMDStructure



! SUBROUTINE that computes and re-orders the DMD modes
! 
! INPUTS:
!	- NX			: number of grid points in x-direction
!	- NY			: number of grid points in y-direction
!	- m				: number of snapshots to be used for DMD
!	- r				: raduced rank desired
!	- dt			: time interval between X1(:,n) and X2(:,n)
!	- X1, X2		: X and X' respectively when following the notation in the DMD book of Kutz et al. 
! 	Note that the code has been written having in mind that NX==NY. Two distinct variables are used to read more easily the code and 
!	to allow a easier extension to the case NX!=NY .
! 
! OUTPUT: 
!	- RDMDs,IDMDs	: Real and imaginary parts of the DMD modes
!	- lambda		: DMD eigenvalues
!	- omega			: omega=log(lambda)/dt
!	- lambda_order	: order of the eigenmodes (i.e. the first corresponds to the mean and then according to how far the eigenvalues is w.r.t. 1+0i)

SUBROUTINE DMDModes(RDMDs,IDMDs,omega,lambda,lambda_order,NX,NY,m,r,dt,X1,X2)
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)									:: NX, NY, m, r
	REAL(KIND=8), INTENT(IN)							:: dt
	REAL(KIND=8), DIMENSION(NX*NY,m)					:: X1, X2
	COMPLEX(KIND=8), DIMENSION(r), INTENT(OUT)			:: omega, lambda
	COMPLEX(KIND=8), DIMENSION(r)						:: lambda_copy
	COMPLEX(KIND=8), DIMENSION(NX*NY,r)					:: Phi
	INTEGER, DIMENSION(r), INTENT(OUT)					:: lambda_order
	INTEGER												:: i, j, counter, idx
	REAL(KIND=8), DIMENSION(NX*NY,r), INTENT(OUT)		:: RDMDs, IDMDs					! Matrix containing all DMDs. Note that they are sorted by column and each column has dimension NX*NY
	
	! Compute DMDs
	RDMDs = 0.0d0
	IDMDs = 0.0d0
	CALL DMD(Phi,omega,lambda,X1,X2,m,r,dt,NX,NY)
	!!! Ensure to consider only eigenvalues on the unit circle and not the spurious ones outside
	!lambda_copy = lambda
	DO i=1, r
		IF(abs(abs(lambda(i))-1)<5e-2) THEN
			lambda_copy(i) = lambda(i)
		ELSE
			lambda_copy(i) = -90
		END IF
	END DO
	!!! Reorder them so that the first column corresponds to the mean and then smaller and smaller scales
	DO i=1, r
		idx					= minloc(abs(lambda_copy-cmplx(1,0)),DIM=1)
		lambda_order(i)		= idx
		lambda_copy(idx)	= -99
		RDMDs(:,i)			= REAL(Phi(:,idx),KIND=8)
		IDMDs(:,i)			= DIMAG(Phi(:,idx))
	END DO
	WRITE(*,*) 'DMD modes updated.'
		
END SUBROUTINE DMDModes




! SUBROUTINE that builds the noise structure and updates the DMD modes
! 
! INPUTS:
!	- RDMDs, IDMDs	: Real and imaginary parts of the DMD modes
!	- omega			: DMD eigenvalues; omega=log(lambda)/dt
!	- lambda_order	: order of the eigenmodes (i.e. the first corresponds to the mean and then according to how far the eigenvalues is w.r.t. 1+0i)
!	- nstart		: number of the first eigenvector to be considered
!	- neigen		: total number of eigenvectors to be used
!	- dt			: time interval between X1(:,n) and X2(:,n)
!	- ll			: for how many time steps the DMD modes need to be propagated
!	- cum_var		: corresponding amplitude of the noise when using EOFs
!	- NX			: number of grid points in x-direction
!	- NY			: number of grid points in y-direction
!	- r				: raduced rank desired
! 	Note that the code has been written having in mind that NX==NY. Two distinct variables are used to read more easily the code and 
!	to allow a easier extension to the case NX!=NY .
! 
! OUTPUT: 
!	- CSigma : spatio-temporal structure of the noise.

SUBROUTINE DMDCovariance(Csigma,RDMDs,IDMDs,omega,lambda_order,nstart,neigen,dt,ll,cum_var,NX,NY,r)
	IMPLICIT NONE
	
	INTEGER, INTENT(IN)								:: nstart, neigen, ll, NX, NY, r
	REAL(KIND=8), INTENT(IN)						:: dt
	REAL(KIND=8), DIMENSION(NX,NY), INTENT(OUT) 	:: Csigma 						! Resulting combination of DMDs in matrix form
	COMPLEX(KIND=8), DIMENSION(r), INTENT(IN)		:: omega
	INTEGER, DIMENSION(r), INTENT(IN)				:: lambda_order
	INTEGER											:: i, j, counter, idx
	REAL(KIND=8), DIMENSION(NX*NY,1)				:: comb_DMD 					! Convex combination of DMDs. Note: it is an array !
	REAL(KIND=8), DIMENSION(NX*NY,NX*NY)			:: WWT							! Covariance matrix given by comb_DMD before being scaled to have the same strenght as EOFs
	REAL(KIND=8), DIMENSION(NX*NY,r), INTENT(IN)	:: RDMDs, IDMDs					! Matrix containing all DMDs. Note that they are sorted by column and each column has dimension NX*NY
	REAL(KIND=8), INTENT(IN)						:: cum_var
	REAL(KIND=8)									:: trWWT	
		
	! Compute the noise structure 
	comb_DMD = 0.0d0
	DO j=nstart, nstart+2*neigen-1, 2 
		comb_DMD(:,1) = comb_DMD(:,1) + ( REAL(exp(omega(lambda_order(j))*dt*ll))*RDMDs(:,j) - DIMAG(exp(omega(lambda_order(j))*dt*ll))*IDMDs(:,j) )
	END DO
	! Normalize it such that it has the same strength as with EOFs
	WWT = matmul(comb_DMD,transpose(comb_DMD))
	trWWT = 0.0d0
	DO i=1, NX*NY
		trWWT = trWWT + WWT(i,i)
	END DO
	comb_DMD = comb_DMD*sqrt(cum_var/trWWT)
	! Check that now it has the desired amount of variance
	!WWT = matmul(comb_DMD,transpose(comb_DMD))
	!trWWT = 0.0d0
	!DO i=1, NX*NY
	!	trWWT = trWWT + WWT(i,i)
	!END DO
	!WRITE(*,*) 'Trace of covariance matrix = ', trWWT
	
	! Now remap the resulting array (NX*NY,1) to be a matrix (NX,NY)
	! Remap according to how the data (on which the EOFs have been computed) had been written originaly (see IO_restart.f90)
	counter = 0
	DO j=1, NY
		DO i=1, NX
			counter = counter + 1
			Csigma(i,j) = comb_DMD(counter,1)
		END DO
	END DO
	WRITE(*,*) 'DMD modes and Csigma propagated.'
	
END SUBROUTINE DMDCovariance
