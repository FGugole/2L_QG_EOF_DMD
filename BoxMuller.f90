! Box-Muller method to generate Gaussian random numbers with a certain mean and standard deviation.
! INPUT:
!	- N  		: number of normally distributed values that will be created;
!	- mu 		: mean of the distribution from which the values will be sampled;
!	- sd 		: standard deviation of the distribution from which the values will be sampled;
!	- X,Y 		: arrays of uniformly distributed values between 0 and 1 created with Mersenne-Twister (reference).
!					Note that X and Y should contain no zeros !
! OUTPUT : 
!	- Output 	: array of normally distributed random values sampled from a distribution N(mu,sd^2)


SUBROUTINE BM_MT(Output,N,mu,sd,X,Y)
	IMPLICIT NONE
	INTEGER :: N
	REAL(KIND=8) :: mu, sd
	REAL(KIND=8), DIMENSION(N) :: Output
	
	REAL(KIND=8), PARAMETER :: pi=3.14159265359d0
	REAL(KIND=8), DIMENSION(N+1) :: X, Y, N1, N2
	INTEGER :: i
	
	DO i=1, N
		N1(i) = sqrt(-2*log(X(i)))*cos(2*pi*Y(i))
		N2(i) = sqrt(-2*log(X(i)))*sin(2*pi*Y(i))
	END DO
	DO i=1, N
		Output(i) = mu + N1(i)*sd
	END DO
	
END SUBROUTINE BM_MT


