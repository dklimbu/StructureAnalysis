!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ==================================================       
!     === PROGRAM TO CALCULATE STRUCTURAL PROPERTIES ===
!     ============== OF AMORPHOUS SILICON ==============
!     ========== USAGE: ./structure.x INPUT.XYZ ========     
!     April 04, 2018 :
!     DIL LIMBU, USM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Gaussian-Quadrature Integral !!!!!!!!!!!
      SUBROUTINE guassquadint(xi,wi,r,gr,q,Sq,bin,bin2)
      IMPLICIT NONE

      integer,parameter :: dp=kind(1.d0)
      integer  :: i, j, ii,bin,bin2
      real(dp) :: r(bin),gr(bin),fr(bin),fx(bin2)
      real(dp) ::  q, Sq, xg
      real     :: xi(bin2),wi(bin2)
       
      do i = 1,bin
         fr(i) = r(i)*(gr(i)-1)*sin(q*r(i))/q
      end do

      ii = 1
      do i = 1,bin2
         xg = xi(i)
         do j = ii,bin-1
            if((xg > r(j)).and.(xg < r(j+1)))then
               fx(i) = fr(j) + (fr(j+1)-fr(j))/(r(j+1)-r(j))*(xg-r(j))
               ii = j
               exit
            end if
         end do
      end do

      Sq = 0.0
      do i = 1,bin2
         Sq = Sq + wi(i)*fx(i)
      end do

      return
      END SUBROUTINE
