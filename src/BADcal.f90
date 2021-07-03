!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ==================================================       
!     === PROGRAM TO CALCULATE STRUCTURAL PROPERTIES ===
!     ============== OF AMORPHOUS SILICON ==============
!     ========== USAGE: ./structure.x INPUT.XYZ ========     
!     April 04, 2018 :
!     DIL LIMBU, USM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Calculate Bond Angle !!!!!!!!!!!
      SUBROUTINE BADcal(ntheta,theta,var)
      USE DataVal
      IMPLICIT NONE

      integer    :: i, j, k
      integer    :: k1, k2, bin, ntheta

      real(dp)   :: rx1, ry1, rz1, rx2, ry2, rz2
      real(dp)   :: RtD, angle, r1, r2, rij, rik
      real(dp)   :: Adist(nbond)
      real(dp)   :: tot, tot2, theta, var
      real(dp)   :: R2D

      RtD = 180.0/pi     
        ntheta = 0
        tot = 0.0
        tot2 = 0.0
        Adist = 0
        do i = 1, N
           if(nn(i).ge.2) then
              do j = 1, nn(i)-1
                 k1 = nmap(j,i)
                 rx1 = x(k1) - x(i)
                 ry1 = y(k1) - y(i)
                 rz1 = z(k1) - z(i)
                 rx1 = rx1 - anint(rx1/L)*L
                 ry1 = ry1 - anint(ry1/L)*L
                 rz1 = rz1 - anint(rz1/L)*L
                 rij = sqrt(rx1*rx1 + ry1*ry1 + rz1*rz1) 
                 do k = j+1, nn(i)
                    k2 = nmap(k,i)
                    rx2 = x(k2) - x(i)
                    ry2 = y(k2) - y(i)
                    rz2 = z(k2) - z(i)
                    rx2 = rx2 - anint(rx2/L)*L
                    ry2 = ry2 - anint(ry2/L)*L
                    rz2 = rz2 - anint(rz2/L)*L
                    rik = sqrt(rx2*rx2 + ry2*ry2 + rz2*rz2)

                    angle = (rx1*rx2 + ry1*ry2 + rz1*rz2)/(rij*rik)
                    angle = acos(angle)*(RtD)
                    ntheta = ntheta + 1
                    bin = int(angle) + 1
                    Adist(bin) = Adist(bin) + 1
                    tot = tot + angle
                    tot2 = tot2 + angle*angle                          
                 end do
              end do
           end if
        end do
!
      theta = tot/ntheta
      var = sqrt((ntheta*tot2 - tot*tot)/ntheta/ntheta)

      open(1, file='bad.dat', status='unknown')
      write(1,'(A8,f9.4, A6,f8.4)') '#Theta: ', theta,'  Std:', var
      do i= 1,nbond
         write(1,'(f10.4,2f12.4)') (i-1)+.5, Adist(i), Adist(i)/N
      end do
      close(1)
  
      END SUBROUTINE BADCal
