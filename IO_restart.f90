! OUTPUT SUBROUTINE. It writes the field variables values in a file .dat whose name starts with QG- and contains the time  
! at which it has been written. Note that the field variables are rearranged to be arrays (NX*NY,1) and not matrixes (NX,NY).
! The time interval between each writing is set in the main code (variable nTec).
! 
! The four lines are as follows:
!	- 1st : header for the time at which the file has been written, the initial total energy and the initial enstrophy;
!	- 2nd : values of the previous header; these values are necessary in order to restart the computations from a desired point in time;
!	- 3rd : header for coordinates and field variables
!	- 4th : info to load the data in Tecplot (you may have to cancel the first two lines in order to actually load it with Tecplot) 
! 			(Tecplot reference : www.tecplot.com)
! From the 5th line on, the output is structured as follows:
!	- 1st column : x-coordinate;
!	- 2nd column : y-coordinate;
!	- 3rd column : Baroclinic potential vorticity (PV);
!	- 4th column : Baroclinic stream function;
!	- 5th column : Barotropic PV;
!	- 6th column : Barotropic stream function.

SUBROUTINE RestartOutput(x,y,qB,psiB,qT,psiT,NX,NY,nit,time,InitEn,InitEns)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NX, NY, nit
   REAL(KIND=8), INTENT(IN) :: time, InitEn, InitEns
   REAL(KIND=8), DIMENSION(NX), INTENT(IN) :: x
   REAL(KIND=8), DIMENSION(NY), INTENT(IN) :: y
   REAL(KIND=8), DIMENSION(NX,NY), INTENT(IN) :: qB, qT, psiB, psiT
   INTEGER              :: i,j 
   CHARACTER(LEN=200)   :: FileName 
   ! Generate filename, open unit and write header   
   WRITE(FileName,'(a,i10.10,a)') 'QG-', nit, '.dat' 
   WRITE(*,'(a,f12.5,a,a)' ) ' t = ', time, ' Output into file ', TRIM(FileName) 
   OPEN(FILE = TRIM(FileName), UNIT=444, RECL = 2000) 
   WRITE(444,*) ' VARIABLES = "time" "InitEn" "InitEns"'
   WRITE(444,'(8e20.10)') time, InitEn, InitEns 
   WRITE(444,*) ' VARIABLES = "x" "y" "qB" "psiB" "qT" "psiT" ' 
   WRITE(444,*) ' ZONE I = ', NX, ' J = ', NY, ' DATAPACKING=POINT ' 
   ! write output into file 
   DO j = 1, NY
      DO i = 1, NX
            WRITE(444,'(8e20.10)') x(i), y(j), qB(i,j), psiB(i,j), qT(i,j), psiT(i,j) 
      ENDDO 
   ENDDO 
   ! close the file 
   CLOSE(444) 
   ! 
END SUBROUTINE RestartOutput 
