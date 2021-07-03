!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! *************************************************************************
! *   Program to calculate structural properties of amorphous silicon
! *   The program reads a configuration of amorphous silicon model as 
! *   a single XYZ file and calculates Rair-Correlation Function (PCF)
! *   Structure Factor (S(k)), Bond-Angle Distribution (BAD), and   
! *   generate Coordination Number (CN) and Nearest Neighor Map (NMAP)
! *
! *   DIL LIMBU, USM
! *   APRIL 2018
! *
! *   TO COMPILE:: Use Makefile provided
! *
! *   USAGE:: ./structure.x INPUT.XYZ
! *
! *   OUTPUTS:: 
! *            gr.dat   <- Pair Correlation Function
! *            Sk.dat   <- Structure Factor
! *            bad.dat  <- Bond Angle Distribution
! *            Cn.dat   <- Coordination Number
! *            nmap.dat <- List of ALL Nearest Neighbors
! *
! * ************************************************************************* 

      PROGRAM StructureAnalysis
      USE DataVal
      IMPLICIT NONE

      integer  :: i,j, k
      integer  :: rbin, qbin
      integer  :: ntheta

      real(dp),allocatable,dimension(:) :: ri, gr
      real(dp) :: dr
      real(dp) :: theta,var     

      character (len=2)  :: atm
      character (len=50) :: infile
         
      CALL get_command_argument(1, infile)

      if(command_argument_count() /= 1)then
         write(0,*) ''
         write(0,*) '============ ERROR !!! ============'
         write(0,*) '= SOMETHING WRONG, PLEASE CHECK!! ='
         write(0,*) '== Usage: ./structure.x INPUT.XYZ ='
         write(0,*) '===== Exiting the program.... ====='
         write(0,*) ''
         stop
      endif
    
      write(0,*) ''
      write(0,*) '=================================================='       
      write(0,*) '=== PROGRAM TO CALCULATE STRUCTURAL PROPERTIES ==='
      write(0,*) '============== OF AMORPHOUS SILICON =============='
      write(0,*) '========== USAGE: ./structure.x INPUT.XYZ ========'
      write(0,*) '==================== OUTPUTS ====================='
      write(0,*) '===== gr.dat   <- Pair Correlation Function ======'
      write(0,*) '===== Sk.dat   <- Structure Factor ==============='
      write(0,*) '===== bad.dat  <- Bond Angle Distribution ========'
      write(0,*) '===== Cn.dat   <- Coordination Number ============'
      write(0,*) '===== nmap.dat <- List of ALL Nearest Neighbors =='
      write(0,*) '=================================================='

      dr = 0.05
      qbin = 1250

      open (1, file = infile, status = 'old')
      read(1,*) N
      read(1,*) L 

      rbin = int(.5*L/dr)
      allocate(x(N),y(N),z(N))
      allocate(nn(N), nmap(nmax,N))
      allocate(ri(rbin),gr(rbin))
               
      do i = 1, N
         read (1,*) atm, x(i), y(i), z(i)   
!         write(*,'(3F12.6)') x(i), y(i), z(i)
      end do
      close(1)

      CALL RDFcal(dr,rbin, ri, gr)

      CALL StructureFactorCal(rbin,2*rbin,qbin,ri,gr,rho)
      
      CALL BADcal(ntheta,theta,var)

      write(0,1000) 'Average theta:', theta, 'Variance: ', var
 1000 format(A18,F9.4,A12,F8.4) 

!!! NN-distribution
      open(1, file='Cn.dat', status='unknown')
      write(0,*) ' ========================'
      write(0,*) '    n   Cn   Cn_percent'
      write(1,'(A)') '#   n   Cn   Cn_percent'
      do i = 2,6
         write(0,'(2I6,f8.2)') i, count(nn==i), (count(nn==i)*100.0/N)
         write(1,'(2I6,f8.2)') i, count(nn==i), (count(nn==i)*100.0/N)
      end do
      close(1)

      write(0,*) ' ========================'
      write(0,*) 'Calculation Completed!!'
      write(0,*) '==== CHECK OUTPUTS ======'
      write(0,*) ''

      deallocate(x,y,z)
      deallocate(nn,nmap)
      deallocate(ri,gr)

      END PROGRAM
