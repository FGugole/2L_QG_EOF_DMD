! OUTPUT SUBROUTINE. It writes energy spectra in a file QG_spectra_*.dat . The file name contains the time ('*') at which 
! it has been written. The time interval between each writing is set in the main code (variable nTec). 
! The output file is structured as follows:
!	- 1st column : Total energy corresponding to the total wave number;
!	- 2nd column : Barotropic kinetic energy corresponding to the total wave number;
!	- 3rd column : Baroclinic kinetic energy corresponding to the total wave number.

SUBROUTINE WriteEnergySpectra(TotSpectra,BKinSpectra,TKinSpectra,NH,nit)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NH, nit
	REAL(KIND=8), DIMENSION(2*NH), INTENT(IN)	:: TotSpectra, BKinSpectra, TKinSpectra
	INTEGER :: k
	CHARACTER(LEN=200) :: Filename
   ! Generate filename, open unit and write header   
   WRITE(FileName,'(a,i10.10,a)') 'QG_spectra_', nit, '.dat' 
   WRITE(*,'(a,a)' ) ' Output into file ', TRIM(FileName) 
   OPEN(FILE = TRIM(FileName), UNIT=444, RECL = 2000) 
   WRITE(444,*) ' VARIABLES = "Total energy" "Barotropic kin-energy" "Baroclinic kin-energy"'
   WRITE(444,*) ' ZONE K = ', 2*NH, ' DATAPACKING=POINT ' 

   DO k=1, 2*NH
       WRITE(444,'(3e20.10)') TotSpectra(k), BKinSpectra(k), TKinSpectra(k)
   ENDDO

   ! close the file
   CLOSE(444)
END SUBROUTINE WriteEnergySpectra
