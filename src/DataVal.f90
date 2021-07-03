!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ==================================================       
!     === PROGRAM TO CALCULATE STRUCTURAL PROPERTIES ===
!     ============== OF AMORPHOUS SILICON ==============
!     ========== USAGE: ./structure.x INPUT.XYZ ========     
!     April 04, 2018 :
!     DIL LIMBU, USM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Data Module !!!!!!!!!!!
      MODULE DataVal

      integer             :: N
      integer,parameter   :: dp=kind(1.d0)
      integer,parameter   :: nmax=10, nbond=180
      integer,allocatable :: nn(:), nmap(:,:)

      real(dp)            :: L, rho    
      real(dp),parameter  :: pi = 3.14159265
      real(dp),allocatable,dimension(:) :: x, y, z

      END MODULE DataVal
