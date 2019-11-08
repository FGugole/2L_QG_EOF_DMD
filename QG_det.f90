! This is the main code. Be sure to link all modules/subroutines when compiling: 
!	AJ4.f90		EnergyOutput.f90	IO_restart.f90		Spectra.f90
! It also uses the FFTW module for the computation of the Fast Fourier Transform.


PROGRAM QG_FFT
	IMPLICIT NONE
	INCLUDE 'fftw3.f'
	
	!!! We solve the non-dimensional 2-layer QG equations reported in Vallis, 'Atmospheric and Oceanic Fluid Dynamics. Fundamentals and large-scale circulation' section 5.6 
	!!! We also include bottom friction, hyperviscosity and forcing 

	! NX, NY 	: number of elements in X and Y direction; NH is for the FFT
	! NMAX 		: maximum number of time iterations (to be sure that the program ends)
	! NPLOT 	: for memory reasons I do not save the value of energy & co at every time step but every 10 time steps
	! nTec 		: the number of iterations every which I want to save the field variables and the spectra
	INTEGER, PARAMETER :: NX=512, NY=512, NH=NX/2+1, NMAX=1e9, NPLOT=1e7, nTec=1e4

	! xL, xR, Lx, dx : domain extension and grid spacing in x-direction
	! yL, yR, Ly, dy : domain extension and grid spacing in y-direction
	! tend 			 : ending time of the simulation
	! f0, beta 		 : parameters for the beta-plane approximation
	! N 			 : Brunt-Vaisala frequency
	! InitEn 		 : initial energy of the system
	! InitEns		 : initial enstrophy of the system
	REAL(KIND=8), PARAMETER :: pi=3.14159265359d0, xL=-pi, xR=pi, yL=-pi, yR=pi
	REAL(KIND=8), PARAMETER :: Lx=xR-xL, Ly=yR-yL, dx=Lx/Nx, dy=Ly/Ny
	REAL(KIND=8), PARAMETER :: tend=5e3, f0=1.0d0, N=1.2e-2, H=200.0d0
	REAL(KIND=8), PARAMETER :: beta=0.509d0, kd=sqrt(2.0)*(2*f0/(N*H))
	REAL(KIND=8) :: time, dt, InitEn, Energy, t1, t2, InitEns, Enstrophy
	REAL(KIND=8) :: nu, tauf, cLeith, U, Energy_B, Energy_T, y0
	LOGICAL :: restart
		
	! Runge-Kutta parameters 
	REAL(KIND=8), DIMENSION(4)   :: b, c
	REAL(KIND=8), DIMENSION(4,4) :: a

	!!!!!!! Variables required for computations
	! To be known:
	! qB, qT 		: potential vorticity of barotropic and baroclinic mode
	! psiB, psiT 	: stream function of barotropic and baroclinic mode
	! LBx, LTx		: Laplacian of psiBx, psiTx respectively (x stands for stage number x of RK)
	! LqBx, LqTx	: Laplacian of qBx, qTx respectively (x stands for stage number x of RK) 
	! LLPx, LLMx	: Second Laplacian of plusx, minus x respectively (x stands for stage number x of RK)
	! 					plus and minus refers to the viscosity coefficient (psiB+psiT or psiB-psiT; see eq.s)
	! EnergyErr		: array containing the difference between the current energy of the system and the initial energy
	! SysEnergy 	: array containing the total energy of the system as time develops
	! EnstrophyErr	: array containing the difference between the current enstrophy of the system and the initial enstrophy
	! SysEnstrophy 	: array containing the total enstrophy of the system as time develops
	! kx, ky 		: wavenumbers in x and y directions respectively
	REAL(KIND=8), DIMENSION(NX) :: x, kx 
	REAL(KIND=8), DIMENSION(NY) :: y, ky
	REAL(KIND=8), DIMENSION(NPLOT) :: EnergyErr, SysEnergy, BEnergy, TEnergy 
	REAL(KIND=8), DIMENSION(NPLOT) :: SysEnstrophy, EnstrophyErr, tplot
	REAL(KIND=8), DIMENSION(NX,NY) :: qB, psiB, rhsB, psiB2, csiB1, csiB2, psiBnew, qBnew    
	REAL(KIND=8), DIMENSION(NX,NY) :: qT, psiT, rhsT, psiT2, csiT1, csiT2, psiTnew, qTnew  
	REAL(KIND=8), DIMENSION(NX,NY) :: psiBU, psiTU, psiBU2, psiTU2, psiBU3, psiTU3, psiBU4, psiTU4
	REAL(KIND=8), DIMENSION(NX,NY) :: DerBx, DerBy, DerBx2, DerTx, DerTy, DerTx2, randT, randB   
	REAL(KIND=8), DIMENSION(NX,NY) :: JBB1, JTT1, JTB1, JBT1, JBB2, JTT2, JTB2, JBT2
	REAL(KIND=8), DIMENSION(NX,NY) :: csiB3, csiT3, csiB4, csiT4, DerBx3, DerTx3, DerBx4, DerTx4
	REAL(KIND=8), DIMENSION(NX,NY) :: psiB3, psiB4, psiT3, psiT4, JBB3, JBB4, JTT3, JTT4
	REAL(KIND=8), DIMENSION(NX,NY) :: JTB3, JTB4, JBT3, JBT4, EnergyGrid 
	REAL(KIND=8), DIMENSION(NX,NY) :: LB1, LB2, LB3, LB4, LT1, LT2, LT3, LT4, LBplusT, LBminusT
	REAL(KIND=8), DIMENSION(NX,NY) :: LqB1, LqB2, LqB3, LqB4, LqT1, LqT2, LqT3, LqT4
	REAL(KIND=8), DIMENSION(NX,NY) :: LLP1, LLP2, LLP3, LLP4, LLM1, LLM2, LLM3, LLM4
	REAL(KIND=8), DIMENSION(2*NH)  :: TotSpectra, BKinSpectra, TKinSpectra
	REAL(KIND=8), DIMENSION(NH)    :: EnergyGridY 

	!!!!!!! Variables necessary for FFT
	COMPLEX(KIND=8), DIMENSION(NH,NY) :: fuB, fuT, frhsB, frhsT, fEnergy, fDBx, fDBy, fDTx, fDTy  
	   
	INTEGER(KIND=8) :: plan_backward
	INTEGER(KIND=8) :: plan_forward
	
	INTEGER :: i, j, nit, k, l, kk, niter

	PRINT *, '--------------------------------------------------------------'
	PRINT *, '- Welcome to QG deterministic model with forcing and damping -'
	PRINT *, '--------------------------------------------------------------'

    !!!!!!! QG deterministic model on a cartesian grid (Finite-Differences) with double-period BC
    ! 1 - Use Arakawa (4th order) discretization of the Jacobian operator to get qB and qT
    !       using the relation J(psiB,qB) = J(psiB,lap(psiB)) + beta*d(psiB)/dx
    ! 		where lap(psiB) stands for the Laplacian of psiB
    ! 2 - Solve with FFT for psiB and for psiT
    ! 3 - Use explicit Runge-Kutta 4 for the time integration
    ! 4 - Compute the energy/enstrophy and compare it with the initial energy/enstrophy
    ! 5 - Compute energy spectra

	! The forcing is introduced through a vertically sheared, baroclinically unstable background flow 

	!! Runge-Kutta parameters
	a 	 = 0.0d0
	a(2,1) = dble(1)/dble(2) 	
	a(3,1) = 0.0d0 				
	a(3,2) = dble(1)/dble(2) 	
	a(4,1) = 0.0d0 				
	a(4,2) = 0.0d0 				
	a(4,3) = 1.0d0 				
	c(1) = 0.0d0 				
	c(2) = dble(1)/dble(2) 		
	c(3) = dble(1)/dble(2) 		
	c(4) = 1.0d0 				
	b(1) = dble(1)/dble(6) 		
	b(2) = dble(1)/dble(3) 		
	b(3) = dble(1)/dble(3) 		
	b(4) = dble(1)/dble(6) 		

	DO i=1, NX
		x(i) = xL + dx/2 + (i-1)*dx
	END DO
	DO j=1, NY
		y(j) = yL + dy/2 + (j-1)*dy
	END DO
		
	!!!!!!! Define IC and time step
	dt 		  = 1e-4 	! Time step
	restart   = .TRUE. ! Set this variable equal to TRUE if you want to resume the simulation from previous (possibly interrupted) simulations. Check that the input file is properly written!
	IF(restart == .TRUE.) THEN
		OPEN(UNIT=333, FILE='QG-0000004000.dat', RECL=8000)
		READ(333,*)
		READ(333,*) time, InitEn, InitEns
		READ(333,*)
		READ(333,*)
		DO j=1, NY
			DO i=1, NX
				READ(333,*) x(i), y(j), qB(i,j), psiB(i,j), qT(i,j), psiT(i,j)
			END DO
		END DO
		CLOSE(333)
		niter = ceiling(time/dt)
		PRINT *, time, niter, InitEn
	ELSE
		DO i=1, NX
			x(i) = xL + dx/2 + (i-1)*dx
		END DO
		DO j=1, NY
			y(j) = yL + dy/2 + (j-1)*dy
		END DO
		time = 0.0d0
		niter = 0
		OPEN(UNIT=333, FILE='IC_512.dat', RECL=8000)
		READ(333,*)
		READ(333,*)
		READ(333,*)
		DO j=1, NY
			DO i=1, NX
				READ(333,*) qB(i,j), psiB(i,j), qT(i,j), psiT(i,j)
			END DO
		END DO
		CLOSE(333)
	END IF
	
	!!!! Define hyperviscosity \nabla^2(nu\nabla^4\psi) and friction 
	! Define tauf directly as the inverse of tau so I can set tauf=0 without damaging anything
	tauf = 1.0d0/10.0d0 
	cLeith = 0.005d0
	! In nu I put only the scalar factors appearing in the definition of nu used in Jansen & Held (2014)
	nu = cLeith*dx**6
	! Background-flow zonal velocity 
	U = 0.6d0
	
	EnergyErr 		= 0.0d0
	SysEnergy 		= 0.0d0
	BEnergy  		= 0.0d0
	TEnergy  		= 0.0d0
	EnstrophyErr 	= 0.0d0
	SysEnstrophy 	= 0.0d0
	
	! Write a file with IC
	CALL RestartOutput(x,y,qB,psiB,qT,psiT,Nx,Ny,niter/nTec,time,InitEn,InitEns) 
	
    CALL CPU_TIME(t1)
	
	!!!!!!! Time loop
	DO nit=niter+1, NMAX
		! Compute initial energy/enstrophy
		IF(nit==1) THEN
			InitEn  = 0.0d0
			InitEns = 0.0d0
			DO j=1, NY
				DO i=1, NX
					! x-direction
					IF(i==1) THEN
						DerBx(i,j) = (psiB(i+1,j)-psiB(NX,j))/(2*dx)
						DerTx(i,j) = (psiT(i+1,j)-psiT(NX,j))/(2*dx)
					ELSEIF(i==NX) THEN
						DerBx(i,j) = (psiB(1,j)-psiB(i-1,j))/(2*dx)
						DerTx(i,j) = (psiT(1,j)-psiT(i-1,j))/(2*dx)
					ELSE
						DerBx(i,j) = (psiB(i+1,j)-psiB(i-1,j))/(2*dx)
						DerTx(i,j) = (psiT(i+1,j)-psiT(i-1,j))/(2*dx)
					END IF
					! y-direction
					IF(j==1) THEN
						DerBy(i,j) = (psiB(i,j+1)-psiB(i,NY))/(2*dy)
						DerTy(i,j) = (psiT(i,j+1)-psiT(i,NY))/(2*dy)
					ELSEIF(j==NY) THEN
						DerBy(i,j) = (psiB(i,1)-psiB(i,j-1))/(2*dy)
						DerTy(i,j) = (psiT(i,1)-psiT(i,j-1))/(2*dy)
					ELSE
						DerBy(i,j) = (psiB(i,j+1)-psiB(i,j-1))/(2*dy)
						DerTy(i,j) = (psiT(i,j+1)-psiT(i,j-1))/(2*dy)
					END IF
					InitEn  = InitEn  + 0.5*dx*dy*( (DerBx(i,j)**2+DerBy(i,j)**2) + (DerTx(i,j)**2+DerTy(i,j)**2) + kd**2*psiT(i,j)**2 ) 
					InitEns = InitEns + 0.5*dx*dy*( (qB(i,j))**2 + (qT(i,j))**2 )
				END DO
			END DO

			PRINT *, ' Initial Energy = ', InitEn
			PRINT *, ' Initial Enstrophy = ', InitEns
		
		END IF
		
		IF(time<3e3) THEN
			y0 = 0.0d0
		ELSE
			y0 = pi/2.0d0
		END IF
		
		IF(time+dt>tend) THEN
			dt = tend - time
		END IF
		IF(time>=tend) THEN
			! Write a file with the energy values against time
			CALL WriteEnergyErr(EnergyErr,SysEnergy,BEnergy,TEnergy,tplot,(nit-1)/10,NMAX)
			! Write a file with the enstrophy values against time
			CALL WriteEnstrophyErr(EnstrophyErr,SysEnstrophy,tplot,(nit-1)/10,NMAX)
			EXIT
		END IF
		
		!! Solve for qB and qT with RK4 
		! First stage
		csiB1 = qB
		csiT1 = qT
		
		! Second stage
		! Get the laplacian of psiB from qB
		IF(nit==1) THEN
			DO j=1, NY
				DO i=1, NX
					LB1(i,j) = qB(i,j) - beta*y(j) ! Actually this is rhsB used to get psiBi i={1,2,3,4} from the laplacian, so after the first iteration I do not need to re-compute it.
				END DO
			END DO
		ELSE
			LB1 = rhsB
		END IF
		! Compute friction and hyperviscosity terms 
		LqB1 = 0.0d0
		LqT1 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				LT1(i,j) = qT(i,j) + kd**2*psiT(i,j)
				! x-derivatives
				IF(i==1) THEN
					LqB1(i,j) = LqB1(i,j) + (qB(i+1,j)-2*qB(i,j)+qB(NX,j))/(dx**2)
					LqT1(i,j) = LqT1(i,j) + (qT(i+1,j)-2*qT(i,j)+qT(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LqB1(i,j) = LqB1(i,j) + (qB(1,j)-2*qB(i,j)+qB(i-1,j))/(dx**2)
					LqT1(i,j) = LqT1(i,j) + (qT(1,j)-2*qT(i,j)+qT(i-1,j))/(dx**2)
				ELSE
					LqB1(i,j) = LqB1(i,j) + (qB(i+1,j)-2*qB(i,j)+qB(i-1,j))/(dx**2)
					LqT1(i,j) = LqT1(i,j) + (qT(i+1,j)-2*qT(i,j)+qT(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LqB1(i,j) = LqB1(i,j) + (qB(i,j+1)-2*qB(i,j)+qB(i,NY))/(dy**2)
					LqT1(i,j) = LqT1(i,j) + (qT(i,j+1)-2*qT(i,j)+qT(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LqB1(i,j) = LqB1(i,j) + (qB(i,1)-2*qB(i,j)+qB(i,j-1))/(dy**2)
					LqT1(i,j) = LqT1(i,j) + (qT(i,1)-2*qT(i,j)+qT(i,j-1))/(dy**2)
				ELSE
					LqB1(i,j) = LqB1(i,j) + (qB(i,j+1)-2*qB(i,j)+qB(i,j-1))/(dy**2)
					LqT1(i,j) = LqT1(i,j) + (qT(i,j+1)-2*qT(i,j)+qT(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add correction to acually get \nabla^4\psi_T
		LqT1 = LqT1 + kd**2*LT1
		! Multiply by viscosity coefficient
		!LqB1 = nu*abs(LqB1)*LqB1
		!LqT1 = nu*abs(LqT1)*LqT1
		LBplusT  = nu*abs(LqB1+LqT1)*(LqB1+LqT1)
		LBminusT = nu*abs(LqB1-LqT1)*(LqB1-LqT1)
		! Compute once more the Laplacian
		LLP1 = 0.0d0
		LLM1 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				! x-derivatives
				IF(i==1) THEN
					LLP1(i,j) = LLP1(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(NX,j))/(dx**2)
					LLM1(i,j) = LLM1(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LLP1(i,j) = LLP1(i,j) + (LBplusT(1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM1(i,j) = LLM1(i,j) + (LBminusT(1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				ELSE
					LLP1(i,j) = LLP1(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM1(i,j) = LLM1(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LLP1(i,j) = LLP1(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,NY))/(dy**2)
					LLM1(i,j) = LLM1(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LLP1(i,j) = LLP1(i,j) + (LBplusT(i,1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM1(i,j) = LLM1(i,j) + (LBminusT(i,1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				ELSE
					LLP1(i,j) = LLP1(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM1(i,j) = LLM1(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add background-flow zonal velocities (Note that it is only on layer 1! -> I have to add the factor 0.5 for barotropic and baroclinic modes)
		DO j=1, NY
			psiBU(:,j) = psiB(:,j) - 0.5*U*(y(j)-y0)
			psiTU(:,j) = psiT(:,j) - 0.5*U*(y(j)-y0)
		END DO
		! Compute the Jacobians
		CALL AJ4(JBB1,psiBU,LB1,NX,NY,dx,dy)
		JBB1 = JBB1 + beta*DerBx
		CALL AJ4(JTT1,psiTU,csiT1,NX,NY,dx,dy)
		CALL AJ4(JTB1,psiTU,LB1,NX,NY,dx,dy)
		JTB1 = JTB1 + beta*DerTx
		CALL AJ4(JBT1,psiBU,csiT1,NX,NY,dx,dy)
		! Second stage
		csiB2 = qB - a(2,1)*dt*(JBB1+JTT1) - a(2,1)*dt*tauf*(LB1-LT1)/2 - a(2,1)*dt*( LLP1+LLM1 )/2
		csiT2 = qT - a(2,1)*dt*(JTB1+JBT1) + a(2,1)*dt*tauf*(LB1-LT1)/2 - a(2,1)*dt*( LLP1-LLM1 )/2

		!! Now, for the second stage of qB and qT, get the corresponding psiB and psiT with FFT
		! Let's start with psiB2
		! Remember: qB = \nabla**2 psiB + beta*y
		DO j=1, NY
			DO i=1, NX
				rhsB(i,j) = csiB2(i,j) - beta*y(j)
			END DO
		END DO
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsB,frhsB,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuB
		! Remember that in position 1 and 1 there are the zero modes!!!
		DO j=1, NY
			DO i=1, NH
				IF(i==1 .AND. j==1) THEN
					fuB(i,j) = frhsB(i,j)
				ELSE
					fuB(i,j) = frhsB(i,j)/(-kx(i)**2-ky(j)**2)
				END IF
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuB,psiB2,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiB2 = psiB2/(NX*NY)
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)
		
		! Now get psiT2
		! Remember: qT = \nabla**2 psiT -kd**2*psiT
		rhsT = csiT2
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsT,frhsT,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuT
		DO j=1, NY
			DO i=1, NH
				fuT(i,j) = frhsT(i,j)/(-kd**2-kx(i)**2-ky(j)**2)
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuT,psiT2,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiT2 = psiT2/(NX*NY)
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)

		! Now evaluate the Jacobians for the new stages
		! First compute the derivatives in x of psiB2 and psiT2
		DO j=1, NY
			DO i=1, NX
				IF(i==1) THEN
					DerBx2(i,j) = (psiB2(i+1,j)-psiB2(NX,j))/(2*dx)
					DerTx2(i,j) = (psiT2(i+1,j)-psiT2(NX,j))/(2*dx)
				ELSEIF(i==NX) THEN
					DerBx2(i,j) = (psiB2(1,j)-psiB2(i-1,j))/(2*dx)
					DerTx2(i,j) = (psiT2(1,j)-psiT2(i-1,j))/(2*dx)
				ELSE
					DerBx2(i,j) = (psiB2(i+1,j)-psiB2(i-1,j))/(2*dx)
					DerTx2(i,j) = (psiT2(i+1,j)-psiT2(i-1,j))/(2*dx)
				END IF
			END DO
		END DO
		! Get the laplacian of psiB2
		! Actually this is rhsB used to get psiBi i={1,2,3,4} from the laplacian
		LB2 = rhsB
		! Compute friction and hyperviscosity terms 
		LqB2 = 0.0d0
		LqT2 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				LT2(i,j) = csiT2(i,j) + kd**2*psiT2(i,j)
				! x-derivatives
				IF(i==1) THEN
					LqB2(i,j) = LqB2(i,j) + (csiB2(i+1,j)-2*csiB2(i,j)+csiB2(NX,j))/(dx**2)
					LqT2(i,j) = LqT2(i,j) + (csiT2(i+1,j)-2*csiT2(i,j)+csiT2(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LqB2(i,j) = LqB2(i,j) + (csiB2(1,j)-2*csiB2(i,j)+csiB2(i-1,j))/(dx**2)
					LqT2(i,j) = LqT2(i,j) + (csiT2(1,j)-2*csiT2(i,j)+csiT2(i-1,j))/(dx**2)
				ELSE
					LqB2(i,j) = LqB2(i,j) + (csiB2(i+1,j)-2*csiB2(i,j)+csiB2(i-1,j))/(dx**2)
					LqT2(i,j) = LqT2(i,j) + (csiT2(i+1,j)-2*csiT2(i,j)+csiT2(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LqB2(i,j) = LqB2(i,j) + (csiB2(i,j+1)-2*csiB2(i,j)+csiB2(i,NY))/(dy**2)
					LqT2(i,j) = LqT2(i,j) + (csiT2(i,j+1)-2*csiT2(i,j)+csiT2(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LqB2(i,j) = LqB2(i,j) + (csiB2(i,1)-2*csiB2(i,j)+csiB2(i,j-1))/(dy**2)
					LqT2(i,j) = LqT2(i,j) + (csiT2(i,1)-2*csiT2(i,j)+csiT2(i,j-1))/(dy**2)
				ELSE
					LqB2(i,j) = LqB2(i,j) + (csiB2(i,j+1)-2*csiB2(i,j)+csiB2(i,j-1))/(dy**2)
					LqT2(i,j) = LqT2(i,j) + (csiT2(i,j+1)-2*csiT2(i,j)+csiT2(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add correction to acually get \nabla^4\psi_T
		LqT2 = LqT2 + kd**2*LT2
		! Multiply by viscosity coefficient
		!LqB2 = nu*abs(LqB2)*LqB2
		!LqT2 = nu*abs(LqT2)*LqT2
		LBplusT  = nu*abs(LqB2+LqT2)*(LqB2+LqT2)
		LBminusT = nu*abs(LqB2-LqT2)*(LqB2-LqT2)
		! Compute once more the Laplacian
		LLP2 = 0.0d0
		LLM2 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				! x-derivatives
				IF(i==1) THEN
					LLP2(i,j) = LLP2(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(NX,j))/(dx**2)
					LLM2(i,j) = LLM2(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LLP2(i,j) = LLP2(i,j) + (LBplusT(1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM2(i,j) = LLM2(i,j) + (LBminusT(1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				ELSE
					LLP2(i,j) = LLP2(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM2(i,j) = LLM2(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LLP2(i,j) = LLP2(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,NY))/(dy**2)
					LLM2(i,j) = LLM2(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LLP2(i,j) = LLP2(i,j) + (LBplusT(i,1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM2(i,j) = LLM2(i,j) + (LBminusT(i,1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				ELSE
					LLP2(i,j) = LLP2(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM2(i,j) = LLM2(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add background-flow zonal velocities (Note that it is only on layer 1! -> I have to add the factor 0.5 for barotropic and baroclinic modes)
		DO j=1, NY
			psiBU2(:,j) = psiB2(:,j) - 0.5*U*(y(j)-y0)
			psiTU2(:,j) = psiT2(:,j) - 0.5*U*(y(j)-y0)
		END DO
		!
		CALL AJ4(JBB2,psiBU2,LB2,NX,NY,dx,dy)
		JBB2 = JBB2 + beta*DerBx2
		CALL AJ4(JTT2,psiTU2,csiT2,NX,NY,dx,dy)
		CALL AJ4(JTB2,psiTU2,LB2,NX,NY,dx,dy)
		JTB2 = JTB2 + beta*DerTx2
		CALL AJ4(JBT2,psiBU2,csiT2,NX,NY,dx,dy)
		! Third stage
		csiB3 = qB - dt*( a(3,1)*(JBB1+JTT1) + a(3,2)*(JBB2+JTT2) ) & 
				   - dt*tauf*( a(3,1)*(LB1-LT1) + a(3,2)*(LB2-LT2) )/2 & 	! friction
				   - dt*( a(3,1)*(LLP1+LLM1) + a(3,2)*(LLP2+LLM2) )/2 		! hyperviscosity 
		csiT3 = qT - dt*( a(3,1)*(JTB1+JBT1) + a(3,2)*(JTB2+JBT2) ) &
				   + dt*tauf*( a(3,1)*(LB1-LT1) + a(3,2)*(LB2-LT2) )/2 & 	! friction
				   - dt*( a(3,1)*(LLP1-LLM1) + a(3,2)*(LLP2-LLM2) )/2		! hyperviscosity
		!! Now, for the third stage of qB and qT, get the corresponding psiB and psiT with FFT
		! Let's start with psiB3
		! Remember: qB = \nabla**2 psiB + beta*y
		DO j=1, NY
			DO i=1, NX
				rhsB(i,j) = csiB3(i,j) - beta*y(j)
			END DO
		END DO
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsB,frhsB,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuB
		! Remember that in position 1 and 1 there are the zero modes!!!
		DO j=1, NY
			DO i=1, NH
				IF(i==1 .AND. j==1) THEN
					fuB(i,j) = frhsB(i,j)
				ELSE
					fuB(i,j) = frhsB(i,j)/(-kx(i)**2-ky(j)**2)
				END IF
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuB,psiB3,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiB3 = psiB3/(NX*NY)
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)
		
		! Now get psiT3
		! Remember: qT = \nabla**2 psiT -kd**2*psiT
		rhsT = csiT3
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsT,frhsT,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuT
		DO j=1, NY
			DO i=1, NH
				fuT(i,j) = frhsT(i,j)/(-kd**2-kx(i)**2-ky(j)**2)
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuT,psiT3,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiT3 = psiT3/(NX*NY)
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)

		! Now evaluate the Jacobians for the new stages
		! First compute the derivatives in x of psiB3 and psiT3
		DO j=1, NY
			DO i=1, NX
				IF(i==1) THEN
					DerBx3(i,j) = (psiB3(i+1,j)-psiB3(NX,j))/(2*dx)
					DerTx3(i,j) = (psiT3(i+1,j)-psiT3(NX,j))/(2*dx)
				ELSEIF(i==NX) THEN
					DerBx3(i,j) = (psiB3(1,j)-psiB3(i-1,j))/(2*dx)
					DerTx3(i,j) = (psiT3(1,j)-psiT3(i-1,j))/(2*dx)
				ELSE
					DerBx3(i,j) = (psiB3(i+1,j)-psiB3(i-1,j))/(2*dx)
					DerTx3(i,j) = (psiT3(i+1,j)-psiT3(i-1,j))/(2*dx)
				END IF
			END DO
		END DO
		! Get the laplacian of psiB3
		! Actually this is rhsB used to get psiBi i={1,2,3,4} from the laplacian
		LB3 = rhsB
		! Compute friction and hyperviscosity terms 
		LqB3 = 0.0d0
		LqT3 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				LT3(i,j) = csiT3(i,j) + kd**2*psiT3(i,j)
				! x-derivatives
				IF(i==1) THEN
					LqB3(i,j) = LqB3(i,j) + (csiB3(i+1,j)-2*csiB3(i,j)+csiB3(NX,j))/(dx**2)
					LqT3(i,j) = LqT3(i,j) + (csiT3(i+1,j)-2*csiT3(i,j)+csiT3(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LqB3(i,j) = LqB3(i,j) + (csiB3(1,j)-2*csiB3(i,j)+csiB3(i-1,j))/(dx**2)
					LqT3(i,j) = LqT3(i,j) + (csiT3(1,j)-2*csiT3(i,j)+csiT3(i-1,j))/(dx**2)
				ELSE
					LqB3(i,j) = LqB3(i,j) + (csiB3(i+1,j)-2*csiB3(i,j)+csiB3(i-1,j))/(dx**2)
					LqT3(i,j) = LqT3(i,j) + (csiT3(i+1,j)-2*csiT3(i,j)+csiT3(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LqB3(i,j) = LqB3(i,j) + (csiB3(i,j+1)-2*csiB3(i,j)+csiB3(i,NY))/(dy**2)
					LqT3(i,j) = LqT3(i,j) + (csiT3(i,j+1)-2*csiT3(i,j)+csiT3(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LqB3(i,j) = LqB3(i,j) + (csiB3(i,1)-2*csiB3(i,j)+csiB3(i,j-1))/(dy**2)
					LqT3(i,j) = LqT3(i,j) + (csiT3(i,1)-2*csiT3(i,j)+csiT3(i,j-1))/(dy**2)
				ELSE
					LqB3(i,j) = LqB3(i,j) + (csiB3(i,j+1)-2*csiB3(i,j)+csiB3(i,j-1))/(dy**2)
					LqT3(i,j) = LqT3(i,j) + (csiT3(i,j+1)-2*csiT3(i,j)+csiT3(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add correction to acually get \nabla^4\psi_T
		LqT3 = LqT3 + kd**2*LT3
		! Multiply by viscosity coefficient
		!LqB3 = nu*abs(LqB3)*LqB3
		!LqT3 = nu*abs(LqT3)*LqT3
		LBplusT  = nu*abs(LqB3+LqT3)*(LqB3+LqT3)
		LBminusT = nu*abs(LqB3-LqT3)*(LqB3-LqT3)
		! Compute once more the Laplacian
		LLP3 = 0.0d0
		LLM3 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				! x-derivatives
				IF(i==1) THEN
					LLP3(i,j) = LLP3(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(NX,j))/(dx**2)
					LLM3(i,j) = LLM3(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LLP3(i,j) = LLP3(i,j) + (LBplusT(1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM3(i,j) = LLM3(i,j) + (LBminusT(1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				ELSE
					LLP3(i,j) = LLP3(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM3(i,j) = LLM3(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LLP3(i,j) = LLP3(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,NY))/(dy**2)
					LLM3(i,j) = LLM3(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LLP3(i,j) = LLP3(i,j) + (LBplusT(i,1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM3(i,j) = LLM3(i,j) + (LBminusT(i,1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				ELSE
					LLP3(i,j) = LLP3(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM3(i,j) = LLM3(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add background-flow zonal velocities (Note that it is only on layer 1! -> I have to add the factor 0.5 for barotropic and baroclinic modes)
		DO j=1, NY
			psiBU3(:,j) = psiB3(:,j) - 0.5*U*(y(j)-y0)
			psiTU3(:,j) = psiT3(:,j) - 0.5*U*(y(j)-y0)
		END DO
		!
		CALL AJ4(JBB3,psiBU3,LB3,NX,NY,dx,dy)
		JBB3 = JBB3 + beta*DerBx3
		CALL AJ4(JTT3,psiTU3,csiT3,NX,NY,dx,dy)
		CALL AJ4(JTB3,psiTU3,LB3,NX,NY,dx,dy)
		JTB3 = JTB3 + beta*DerTx3
		CALL AJ4(JBT3,psiBU3,csiT3,NX,NY,dx,dy)
		! Fourth stage
		csiB4 = qB - dt*( a(4,1)*(JBB1+JTT1) + a(4,2)*(JBB2+JTT2) + a(4,3)*(JBB3+JTT3) ) & 
				   - dt*tauf*( a(4,1)*(LB1-LT1) + a(4,2)*(LB2-LT2) + a(4,3)*(LB3-LT3))/2 & 	! friction
				   - dt*( a(4,1)*(LLP1+LLM1) + a(4,2)*(LLP2+LLM2) + a(4,3)*(LLP3+LLM3) )/2	! hyperviscosity
		csiT4 = qT - dt*( a(4,1)*(JTB1+JBT1) + a(4,2)*(JTB2+JBT2) + a(4,3)*(JTB3+JBT3) ) &
				   + dt*tauf*( a(4,1)*(LB1-LT1) + a(4,2)*(LB2-LT2) + a(4,3)*(LB3-LT3))/2 & 	! friction
				   - dt*( a(4,1)*(LLP1-LLM1) + a(4,2)*(LLP2-LLM2) + a(4,3)*(LLP3-LLM3) )/2 	! hyperviscosity
		!! Now, for the fourth stage of qB and qT, get the corresponding psiB and psiT with FFT
		! Let's start with psiB4
		! Remember: qB = \nabla**2 psiB + beta*y
		DO j=1, NY
			DO i=1, NX
				rhsB(i,j) = csiB4(i,j) - beta*y(j)
			END DO
		END DO
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsB,frhsB,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuB
		! Remember that in position 1 and 1 there are the zero modes!!!
		DO j=1, NY
			DO i=1, NH
				IF(i==1 .AND. j==1) THEN
					fuB(i,j) = frhsB(i,j)
				ELSE
					fuB(i,j) = frhsB(i,j)/(-kx(i)**2-ky(j)**2)
				END IF
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuB,psiB4,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiB4 = psiB4/(NX*NY)
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)
		
		! Now get psiT3
		! Remember: qT = \nabla**2 psiT -kd**2*psiT
		rhsT = csiT4
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsT,frhsT,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuT
		DO j=1, NY
			DO i=1, NH
				fuT(i,j) = frhsT(i,j)/(-kd**2-kx(i)**2-ky(j)**2)
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuT,psiT4,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiT4 = psiT4/(NX*NY)
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)

		! Now evaluate the Jacobians for the new stages
		! First compute the derivatives in x of psiB3 and psiT3
		DO j=1, NY
			DO i=1, NX
				IF(i==1) THEN
					DerBx4(i,j) = (psiB4(i+1,j)-psiB4(NX,j))/(2*dx)
					DerTx4(i,j) = (psiT4(i+1,j)-psiT4(NX,j))/(2*dx)
				ELSEIF(i==NX) THEN
					DerBx4(i,j) = (psiB4(1,j)-psiB4(i-1,j))/(2*dx)
					DerTx4(i,j) = (psiT4(1,j)-psiT4(i-1,j))/(2*dx)
				ELSE
					DerBx4(i,j) = (psiB4(i+1,j)-psiB4(i-1,j))/(2*dx)
					DerTx4(i,j) = (psiT4(i+1,j)-psiT4(i-1,j))/(2*dx)
				END IF
			END DO
		END DO
		! Get the laplacian of psiB4
		! Actually this is rhsB used to get psiBi i={1,2,3,4} from the laplacian
		LB4 = rhsB
		! Compute friction and hyperviscosity terms 
		LqB4 = 0.0d0
		LqT4 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				LT4(i,j) = csiT4(i,j) + kd**2*psiT4(i,j)
				! x-derivatives
				IF(i==1) THEN
					LqB4(i,j) = LqB4(i,j) + (csiB4(i+1,j)-2*csiB4(i,j)+csiB4(NX,j))/(dx**2)
					LqT4(i,j) = LqT4(i,j) + (csiT4(i+1,j)-2*csiT4(i,j)+csiT4(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LqB4(i,j) = LqB4(i,j) + (csiB4(1,j)-2*csiB4(i,j)+csiB4(i-1,j))/(dx**2)
					LqT4(i,j) = LqT4(i,j) + (csiT4(1,j)-2*csiT4(i,j)+csiT4(i-1,j))/(dx**2)
				ELSE
					LqB4(i,j) = LqB4(i,j) + (csiB4(i+1,j)-2*csiB4(i,j)+csiB4(i-1,j))/(dx**2)
					LqT4(i,j) = LqT4(i,j) + (csiT4(i+1,j)-2*csiT4(i,j)+csiT4(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LqB4(i,j) = LqB4(i,j) + (csiB4(i,j+1)-2*csiB4(i,j)+csiB4(i,NY))/(dy**2)
					LqT4(i,j) = LqT4(i,j) + (csiT4(i,j+1)-2*csiT4(i,j)+csiT4(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LqB4(i,j) = LqB4(i,j) + (csiB4(i,1)-2*csiB4(i,j)+csiB4(i,j-1))/(dy**2)
					LqT4(i,j) = LqT4(i,j) + (csiT4(i,1)-2*csiT4(i,j)+csiT4(i,j-1))/(dy**2)
				ELSE
					LqB4(i,j) = LqB4(i,j) + (csiB4(i,j+1)-2*csiB4(i,j)+csiB4(i,j-1))/(dy**2)
					LqT4(i,j) = LqT4(i,j) + (csiT4(i,j+1)-2*csiT4(i,j)+csiT4(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add correction to acually get \nabla^4\psi_T
		LqT4 = LqT4 + kd**2*LT4
		! Multiply by viscosity coefficient
		!LqB4 = nu*abs(LqB4)*LqB4
		!LqT4 = nu*abs(LqT4)*LqT4
		LBplusT  = nu*abs(LqB4+LqT4)*(LqB4+LqT4)
		LBminusT = nu*abs(LqB4-LqT4)*(LqB4-LqT4)
		! Compute once more the Laplacian
		LLP4 = 0.0d0
		LLM4 = 0.0d0
		DO j=1, NY
			DO i=1, NX
				! x-derivatives
				IF(i==1) THEN
					LLP4(i,j) = LLP4(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(NX,j))/(dx**2)
					LLM4(i,j) = LLM4(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(NX,j))/(dx**2)
				ELSEIF(i==NX) THEN
					LLP4(i,j) = LLP4(i,j) + (LBplusT(1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM4(i,j) = LLM4(i,j) + (LBminusT(1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				ELSE
					LLP4(i,j) = LLP4(i,j) + (LBplusT(i+1,j) -2*LBplusT(i,j) +LBplusT(i-1,j))/(dx**2)
					LLM4(i,j) = LLM4(i,j) + (LBminusT(i+1,j)-2*LBminusT(i,j)+LBminusT(i-1,j))/(dx**2)
				END IF
				! y-derivatives
				IF(j==1) THEN
					LLP4(i,j) = LLP4(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,NY))/(dy**2)
					LLM4(i,j) = LLM4(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,NY))/(dy**2)
				ELSEIF(j==NY) THEN
					LLP4(i,j) = LLP4(i,j) + (LBplusT(i,1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM4(i,j) = LLM4(i,j) + (LBminusT(i,1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				ELSE
					LLP4(i,j) = LLP4(i,j) + (LBplusT(i,j+1) -2*LBplusT(i,j) +LBplusT(i,j-1))/(dy**2)
					LLM4(i,j) = LLM4(i,j) + (LBminusT(i,j+1)-2*LBminusT(i,j)+LBminusT(i,j-1))/(dy**2)
				END IF
			END DO
		END DO
		! Add background-flow zonal velocities (Note that it is only on layer 1! -> I have to add the factor 0.5 for barotropic and baroclinic modes)
		DO j=1, NY
			psiBU4(:,j) = psiB4(:,j) - 0.5*U*(y(j)-y0)
			psiTU4(:,j) = psiT4(:,j) - 0.5*U*(y(j)-y0)
		END DO
		!
		CALL AJ4(JBB4,psiBU4,LB4,NX,NY,dx,dy)
		JBB4 = JBB4 + beta*DerBx4
		CALL AJ4(JTT4,psiTU4,csiT4,NX,NY,dx,dy)
		CALL AJ4(JTB4,psiTU4,LB4,NX,NY,dx,dy)
		JTB4 = JTB4 + beta*DerTx4
		CALL AJ4(JBT4,psiBU4,csiT4,NX,NY,dx,dy)

		!! qB and qT at n+1
		qBnew = qB - dt*( b(1)*(JBB1+JTT1) + b(2)*(JBB2+JTT2) + b(3)*(JBB3+JTT3) + b(4)*(JBB4+JTT4) ) & 
				   - dt*tauf*( b(1)*(LB1-LT1) + b(2)*(LB2-LT2) + b(3)*(LB3-LT3) + b(4)*(LB4-LT4))/2 & 	! friction
				   - dt*( b(1)*(LLP1+LLM1) + b(2)*(LLP2+LLM2) + b(3)*(LLP3+LLM3) + b(4)*(LLP4+LLM4) )/2 ! hyperviscosity 
		qTnew = qT - dt*( b(1)*(JTB1+JBT1) + b(2)*(JTB2+JBT2) + b(3)*(JTB3+JBT3) + b(4)*(JTB4+JBT4) ) &
				   + dt*tauf*( b(1)*(LB1-LT1) + b(2)*(LB2-LT2) + b(3)*(LB3-LT3) + b(4)*(LB4-LT4))/2 & 	! friction
				   - dt*( b(1)*(LLP1-LLM1) + b(2)*(LLP2-LLM2) + b(3)*(LLP3-LLM3) + b(4)*(LLP4-LLM4) )/2 ! hyperviscosity
		qB = qBnew
		qT = qTnew
		
		!! Get psiB at n+1
		! Remember: qB = \nabla**2 psiB + beta*y
		DO j=1, NY
			DO i=1, NX
				rhsB(i,j) = qB(i,j) - beta*y(j)
			END DO
		END DO
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsB,frhsB,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuB
		! Remember that in position 1 and 1 there are the zero modes!!!
		DO j=1, NY
			DO i=1, NH
				IF(i==1 .AND. j==1) THEN
					fuB(i,j) = frhsB(i,j)
				ELSE
					fuB(i,j) = frhsB(i,j)/(-kx(i)**2-ky(j)**2)
				END IF
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuB,psiBnew,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiBnew = psiBnew/(NX*NY)
		psiB 	= psiBnew
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)
		
		!! Get psiT at n+1
		! Remember: qT = \nabla**2 psiT -kd**2*psiT
		rhsT = qT
		! Define wave numbers kx and ky
		! x-direction
		DO i=0, NX/2 !-1
			kx(i+1) = (2*pi/Lx)*i
		END DO
		j = NX/2+1
		DO i=-NX/2+1, -1
			j = j+1
			kx(j) = (2*pi/Lx)*i
		END DO
		! y-direction
		DO j=0, NY/2 !-1
			ky(j+1) = (2*pi/Ly)*j
		END DO
		i = NY/2+1
		DO j=-NY/2+1, -1
			i = i+1
			ky(i) = (2*pi/Ly)*j
		END DO
		! Transform forward
		CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,rhsT,frhsT,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_forward)
		! Get fuT
		! Remember that in position 1 and 1 there are the zero modes!!!
		DO j=1, NY
			DO i=1, NH
				fuT(i,j) = frhsT(i,j)/(-kd**2-kx(i)**2-ky(j)**2)
			END DO
		END DO
		! Transform backward
		CALL dfftw_plan_dft_c2r_2d(plan_backward,NX,NY,fuT,psiTnew,FFTW_ESTIMATE)
		CALL dfftw_execute(plan_backward)		
		! Remember to normalize!!!
		psiTnew = psiTnew/(NX*NY)
		psiT 	= psiTnew
		CALL dfftw_destroy_plan(plan_forward)
		CALL dfftw_destroy_plan(plan_backward)

		! Compute the energy of the system at the new time step
		Energy		= 0.0d0
		Energy_B	= 0.0d0
		Energy_T	= 0.0d0
		Enstrophy	= 0.0d0
		DO j=1, NY
			DO i=1, NX
				! x-direction
				IF(i==1) THEN
					DerBx(i,j) = (psiB(i+1,j)-psiB(NX,j))/(2*dx)
					DerTx(i,j) = (psiT(i+1,j)-psiT(NX,j))/(2*dx)
				ELSEIF(i==NX) THEN
					DerBx(i,j) = (psiB(1,j)-psiB(i-1,j))/(2*dx)
					DerTx(i,j) = (psiT(1,j)-psiT(i-1,j))/(2*dx)
				ELSE
					DerBx(i,j) = (psiB(i+1,j)-psiB(i-1,j))/(2*dx)
					DerTx(i,j) = (psiT(i+1,j)-psiT(i-1,j))/(2*dx)
				END IF
				! y-direction
				IF(j==1) THEN
					DerBy(i,j) = (psiB(i,j+1)-psiB(i,NY))/(2*dy)
					DerTy(i,j) = (psiT(i,j+1)-psiT(i,NY))/(2*dy)
				ELSEIF(j==NY) THEN
					DerBy(i,j) = (psiB(i,1)-psiB(i,j-1))/(2*dy)
					DerTy(i,j) = (psiT(i,1)-psiT(i,j-1))/(2*dy)
				ELSE
					DerBy(i,j) = (psiB(i,j+1)-psiB(i,j-1))/(2*dy)
					DerTy(i,j) = (psiT(i,j+1)-psiT(i,j-1))/(2*dy)
				END IF
				Energy		= Energy 	+ 0.5*dx*dy*( (DerBx(i,j)**2+DerBy(i,j)**2) + (DerTx(i,j)**2+DerTy(i,j)**2) + kd**2*psiT(i,j)**2 ) 
				Energy_B	= Energy_B	+ 0.5*dx*dy*( (DerBx(i,j)**2+DerBy(i,j)**2) )
				Energy_T	= Energy_T	+ 0.5*dx*dy*( (DerTx(i,j)**2+DerTy(i,j)**2) + kd**2*psiT(i,j)**2 )
				Enstrophy	= Enstrophy + 0.5*dx*dy*( (qB(i,j))**2 + (qT(i,j))**2 )
			END DO
		END DO
		! Update time
		time = time + dt
		! Save the energy/enstrophy values every 10 time steps (together with the time)
		IF(mod(nit-niter,10)==0) THEN
			SysEnergy((nit-niter)/10) 		= Energy
			EnergyErr((nit-niter)/10) 		= InitEn - Energy
			BEnergy((nit-niter)/10)			= Energy_B
			TEnergy((nit-niter)/10)			= Energy_T
			SysEnstrophy((nit-niter)/10) 	= Enstrophy
			EnstrophyErr((nit-niter)/10) 	= InitEns - Enstrophy
			tplot((nit-niter)/10) 			= time
		END IF
		
		PRINT *, '- Energy of the system = ', Energy
		PRINT *, '- Energy error = ', InitEn - Energy
		PRINT *, '- Relative energy error = ', (InitEn - Energy)/InitEn
        PRINT *, ' ------------------------------------------------------------ ' 
		PRINT *, '- Enstrophy of the system = ', Enstrophy
		PRINT *, '- Enstrophy error = ', InitEns - Enstrophy
		PRINT *, '- Relative enstrophy error = ', (InitEns - Enstrophy)/InitEns
		
        PRINT *, ' ------------------------------------------------------------ ' 
        PRINT *, '  Current time = ', time         
        PRINT *, ' ------------------------------------------------------------ ' 

		IF(mod(nit,nTec)==0) THEN
			! Print fields values
			CALL RestartOutput(x,y,qB,psiB,qT,psiT,Nx,Ny,nit/nTec,time,InitEn,InitEns) 
			
!			!!!!!!! Compute energy spectra
!			! Compute the energy for every grid point
!			EnergyGrid 	= 0.0d0
!			DO i=1, NX
!				DO j=1, NY
!					EnergyGrid(i,j) = 0.5*dx*dy*( (DerBx(i,j)**2+DerBy(i,j)**2) + (DerTx(i,j)**2+DerTy(i,j)**2) + kd**2*psiT(i,j)**2 )
!				END DO
!			END DO
!			! Go to spectral space
!			! Total energy
!			CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,EnergyGrid,fEnergy,FFTW_ESTIMATE)
!			CALL dfftw_execute(plan_forward)
!
!			! Barotropic velocities
!			! Remember that u=-dpsi/dy and v=dpsi/dx
!			CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,-DerBy,fDBx,FFTW_ESTIMATE)
!			CALL dfftw_execute(plan_forward)
!			CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,DerBx,fDBy,FFTW_ESTIMATE)
!			CALL dfftw_execute(plan_forward)
!
!			! Baroclinic velocities
!			! Remember that u=-dpsi/dy and v=dpsi/dx
!			CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,-DerTy,fDTx,FFTW_ESTIMATE)
!			CALL dfftw_execute(plan_forward)
!			CALL dfftw_plan_dft_r2c_2d(plan_forward,NX,NY,DerTx,fDTy,FFTW_ESTIMATE)
!			CALL dfftw_execute(plan_forward)
!
!			! Compute spectra for the total wave number
!			TotSpectra 	= 0.0d0
!			BKinSpectra = 0.0d0
!			TKinSpectra = 0.0d0
!			DO k=0, NH-1
!				DO l=0, NH-1
!					IF (k==0 .AND. l==0) THEN
!						CYCLE
!					ELSE
!						kk = int( sqrt(dble(k)**2+dble(l)**2) )
!						TotSpectra(kk) 	= TotSpectra(kk)  + 0.5*2*pi*kk*abs( fEnergy(k+1,l+1) )
!						BKinSpectra(kk) = BKinSpectra(kk) + 0.5*4*pi*kk**2*abs( fDBx(k+1,l+1)**2 + fDBy(k+1,l+1)**2 )
!						TKinSpectra(kk) = TKinSpectra(kk) + 0.5*4*pi*kk**2*abs( fDTx(k+1,l+1)**2 + fDTy(k+1,l+1)**2 )
!					END IF
!				END DO
!				DO l=1, NH-2 ! negative wave numbers
!					kk = int( sqrt(dble(k)**2+dble(l)**2) )
!					TotSpectra(kk) 	= TotSpectra(kk)  + 0.5*2*pi*kk*abs( fEnergy(k+1,NY-l+1) )
!					BKinSpectra(kk) = BKinSpectra(kk) + 0.5*4*pi*kk**2*abs( fDBx(k+1,NY-l+1)**2 + fDBy(k+1,NY-l+1)**2 )
!					TKinSpectra(kk) = TKinSpectra(kk) + 0.5*4*pi*kk**2*abs( fDTx(k+1,NY-l+1)**2 + fDTy(k+1,NY-l+1)**2 )
!				END DO
!			END DO
!
!			CALL dfftw_destroy_plan(plan_forward)
!			CALL WriteEnergySpectra(TotSpectra,BKinSpectra,TKinSpectra,NH,nit)
!
		END IF
		
	END DO ! time loop
	
    CALL CPU_TIME(t2)

    PRINT *, ' ----------------------------------------------------------- '
    PRINT *, '  Program terminated correctly. Bye.  '
    PRINT *, '  CPU time = ', t2-t1
    PRINT *, ' ----------------------------------------------------------- '
	
END PROGRAM
