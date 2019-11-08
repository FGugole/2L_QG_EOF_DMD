! OUTPUT SUBROUTINE. The output file QG_energy.dat is structured as follows:
!	- 1st column : time;
!	- 2nd column : total energy error w.r.t. the total initial energy of the system (i.e. Energy(t=0) - Energy(t>0) );
!	- 3rd column : relative total energy error, i.e. it is equal to the second column diveded by the initial total energy
!	- 4th column : total energy of the system at the corresponding time.

SUBROUTINE WriteEnergyErr(EnergyErr,SysEnergy,BEnergy,TEnergy,tplot,nit,NMAX)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nit, NMAX
   REAL(KIND=8), DIMENSION(NMAX), INTENT(IN) :: EnergyErr, SysEnergy, tplot, BEnergy, TEnergy
   INTEGER :: i
   CHARACTER(LEN=200)   :: FileName 
   ! Generate filename, open unit and write header   
   WRITE(FileName,'(a,i10.10,a)') 'QG_energy_', nit, '.dat' 
   WRITE(*,'(a,a)' ) ' Output into file ', TRIM(FileName) 
   OPEN(FILE = TRIM(FileName), UNIT=444, RECL = 2000) 
   WRITE(444,*) ' VARIABLES = "Time" "Energy Error" "Energy of the system" "Barotropic energy" "Baroclinic Energy"'
   WRITE(444,*) ' ZONE I = ', nit, ' DATAPACKING=POINT ' 

   DO i=1, nit
       WRITE(444,'(5e20.10)') tplot(i), EnergyErr(i), SysEnergy(i), BEnergy(i), TEnergy(i)
   ENDDO

   ! close the file
   CLOSE(444)
END SUBROUTINE WriteEnergyErr



! OUTPUT SUBROUTINE. The output file QG_enstrophy.dat is structured as follows:
!	- 1st column : time;
!	- 2nd column : total enstrophy error w.r.t. the total initial enstrophy of the system (i.e. Enstrophy(t=0) - Enstrophy(t>0) );
!	- 3rd column : relative total enstrophy error, i.e. it is equal to the second column diveded by the initial total enstrophy
!	- 4th column : total enstrophy of the system at the corresponding time.

SUBROUTINE WriteEnstrophyErr(EnstrophyErr,SysEnstrophy,tplot,nit,NMAX)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nit, NMAX
   REAL(KIND=8), DIMENSION(NMAX), INTENT(IN) :: EnstrophyErr, SysEnstrophy, tplot
   INTEGER :: i
   CHARACTER(LEN=200)   :: FileName 
   ! Generate filename, open unit and write header   
   WRITE(FileName,'(a,i10.10,a)') 'QG_enstrophy_', nit, '.dat' 
   WRITE(*,'(a,a)' ) ' Output into file ', TRIM(FileName) 
   OPEN(FILE = TRIM(FileName), UNIT=444, RECL = 2000) 
   WRITE(444,*) ' VARIABLES = "Time" "Enstrophy Error" "Enstrophy of the system"'
   WRITE(444,*) ' ZONE I = ', nit, ' DATAPACKING=POINT ' 

   DO i=1, nit
       WRITE(444,'(3e20.10)') tplot(i), EnstrophyErr(i), SysEnstrophy(i)
   ENDDO

   ! close the file
   CLOSE(444)
END SUBROUTINE WriteEnstrophyErr
