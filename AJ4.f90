! Subroutine computing the fourth order Arakawa discritization of the Jacobian. 
! See ' Computational Design for Long-Term Numerical Integration of the Equations of Fluid Motion: Two-Dimensional
! Incompressible Flow. Part I ', Akio Arakawa, Journ. Comp. Physics 1, 119-143 (1966)
!
! INPUT : 
!	- psi : first variable involved in the Jacobian;
!	- q   : second variable involved in the Jacobian;
!	- NX  : number of grid points used in x-direction;
!	- NY  : number of grid points used in y-direction;
!	- dx  : grid spacing in x-direction;
!	- dy  : grid spacing in y-direction;
! 
! OUTPUT:
!	- Jac : Jacobian of variables psi and q.

SUBROUTINE AJ4(Jac,psi,q,NX,NY,dx,dy)
    IMPLICIT NONE
    INTEGER :: NX, NY
    REAL(KIND=8) :: dx, dy
    REAL(KIND=8), DIMENSION(NX,NY) ::Jac, psi, q

    INTEGER :: i, j
    REAL(KIND=8), DIMENSION(NX+4,NY+4)    :: qBC, psiBC
    REAL(KIND=8), DIMENSION(NX,NY) :: J1, J2

    ! Function for the Arakawa discretization of the Jacobian
    ! Create a bigger matrix (containing boundary conditions) for q and for psi
    qBC(3:NX+2,3:NY+2) = q(:,:)
    ! Double periodic BC
    qBC(2,3:NY+2) = q(NX,:)
    qBC(1,3:NY+2) = q(NX-1,:)
    
    qBC(NX+3,3:NY+2) = q(1,:)
    qBC(NX+4,3:NY+2) = q(2,:)
    
    qBC(3:NX+2,2) = q(:,NY)
    qBC(3:NX+2,1) = q(:,NY-1)
    
    qBC(3:NX+2,NY+3) = q(:,1)
    qBC(3:NX+2,NY+4) = q(:,2)
    
    qBC(2,1:2) = q(NX,NY-1:NY)
    qBC(1,1:2) = q(NX-1,NY-1:NY)
    
    qBC(2,NY+3:NY+4) = q(NX,1:2)
    qBC(1,NY+3:NY+4) = q(NX-1,1:2)
    
    qBC(NX+3,1:2) = q(1,NY-1:NY)
    qBC(NX+4,1:2) = q(2,NY-1:NY)
    
    qBC(NX+3,NY+3:NY+4) = q(1,1:2)
    qBC(NX+4,NY+3:NY+4) = q(2,1:2)
    
    ! And now for psi
    psiBC(3:NX+2,3:NY+2) = psi(:,:)
    ! Double periodic BC
    psiBC(2,3:NY+2) = psi(NX,:)
    psiBC(1,3:NY+2) = psi(NX-1,:)
    
    psiBC(NX+3,3:NY+2) = psi(1,:)
    psiBC(NX+4,3:NY+2) = psi(2,:)
    
    psiBC(3:NX+2,2) = psi(:,NY)
    psiBC(3:NX+2,1) = psi(:,NY-1)
    
    psiBC(3:NX+2,NY+3) = psi(:,1)
    psiBC(3:NX+2,NY+4) = psi(:,2)
    
    psiBC(2,1:2) = psi(NX,NY-1:NY)
    psiBC(1,1:2) = psi(NX-1,NY-1:NY)
    
    psiBC(2,NY+3:NY+4) = psi(NX,1:2)
    psiBC(1,NY+3:NY+4) = psi(NX-1,1:2)
    
    psiBC(NX+3,1:2) = psi(1,NY-1:NY)
    psiBC(NX+4,1:2) = psi(2,NY-1:NY)
    
    psiBC(NX+3,NY+3:NY+4) = psi(1,1:2)
    psiBC(NX+4,NY+3:NY+4) = psi(2,1:2)
    
    !! Compute the Arakawa discretization of the Jacobian
    ! Compute J1
    DO i=3, NX+2
        DO j=3, NY+2
            J1(i-2,j-2) = -1/(12*dx*dy)*( (psiBC(i,j-1)+psiBC(i+1,j-1)-psiBC(i,j+1)-psiBC(i+1,j+1))*(qBC(i+1,j)-qBC(i,j)) + &
                                          (psiBC(i-1,j-1)+psiBC(i,j-1)-psiBC(i-1,j+1)-psiBC(i,j+1))*(qBC(i,j)-qBC(i-1,j)) + &
                                          (psiBC(i+1,j)+psiBC(i+1,j+1)-psiBC(i-1,j)-psiBC(i-1,j+1))*(qBC(i,j+1)-qBC(i,j)) + &
                                          (psiBC(i+1,j-1)+psiBC(i+1,j)-psiBC(i-1,j-1)-psiBC(i-1,j))*(qBC(i,j)-qBC(i,j-1)) + &
                                          (psiBC(i+1,j)-psiBC(i,j+1))*(qBC(i+1,j+1)-qBC(i,j)) + &
                                          (psiBC(i,j-1)-psiBC(i-1,j))*(qBC(i,j)-qBC(i-1,j-1)) + &
                                          (psiBC(i,j+1)-psiBC(i-1,j))*(qBC(i-1,j+1)-qBC(i,j)) + &
                                          (psiBC(i+1,j)-psiBC(i,j-1))*(qBC(i,j)-qBC(i+1,j-1)) )
        ENDDO
    ENDDO
    ! Compute J2
    DO i=3, NX+2
		DO j=3, NY+2
			J2(i-2,j-2) = 1/(24*dx*dy)*( (qBC(i+1,j+1)-qBC(i-1,j-1))*(psiBC(i-1,j+1)-psiBC(i+1,j-1)) - &
										 (qBC(i-1,j+1)-qBC(i+1,j-1))*(psiBC(i+1,j+1)-psiBC(i-1,j-1)) + &
										  qBC(i+1,j+1)*(psiBC(i,j+2)-psiBC(i+2,j)) - qBC(i-1,j-1)*(psiBC(i-2,j)-psiBC(i,j-2)) - &
										  qBC(i-1,j+1)*(psiBC(i,j+2)-psiBC(i-2,j)) + qBC(i+1,j-1)*(psiBC(i+2,j)-psiBC(i,j-2)) + &
										  qBC(i+2,j)*(psiBC(i+1,j+1)-psiBC(i+1,j-1)) - qBC(i-2,j)*(psiBC(i-1,j+1)-psiBC(i-1,j-1)) - &
										  qBC(i,j+2)*(psiBC(i+1,j+1)-psiBC(i-1,j+1)) + qBC(i,j-2)*(psiBC(i+1,j-1)-psiBC(i-1,j-1)) )
		END DO
    END DO

	Jac = - (2*J1 - J2) 

END SUBROUTINE AJ4
