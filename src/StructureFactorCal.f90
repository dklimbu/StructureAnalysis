!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ==================================================       
!     === PROGRAM TO CALCULATE STRUCTURAL PROPERTIES ===
!     ============== OF AMORPHOUS SILICON ==============
!     ========== USAGE: ./structure.x INPUT.XYZ ========     
!     April 04, 2018 :: Gaussian-Quadrature
!     DIL LIMBU, USM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Calculation of S(q) !!!!!!!!!!!
      SUBROUTINE StructureFactorCal(bin,bin2,qbin,ri,gr,den)
      IMPLICIT NONE

      integer, parameter :: dp=kind(1.d0)
      integer :: i,bin,bin2,qbin
      real(dp),dimension(bin)  :: ri, gr
      real(dp),dimension(qbin) :: qi, Sq
      real(dp) :: qmax, dq, q, intgr, den 
      real,parameter :: pi = 3.141593

      REAL     :: a,b,xi(bin2),wi(bin2)

!      finding gaussian-legendre nodes and weights
      a = minval(ri(1:bin))
      b = ri(bin)
      qmax = 25.0

      call gauleg(a,b,xi,wi,bin2)

      dq = (qmax)/qbin
      open(1, file='Sk.dat', status='unknown')
      do i = 1,qbin
         q = 0.5 + (i-1)*dq
         Sq(i) = 0.0
         call guassquadint(xi,wi,ri,gr,q,intgr,bin,bin2)
         Sq(i) = 1 + intgr*4.0*pi*den
         write(1,'(f8.4,f12.4)') q, Sq(i)
      end do
      close(1)

      END SUBROUTINE
