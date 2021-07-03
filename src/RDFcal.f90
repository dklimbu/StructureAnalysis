!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ==================================================       
!     === PROGRAM TO CALCULATE STRUCTURAL PROPERTIES ===
!     ============== OF AMORPHOUS SILICON ==============
!     ========== USAGE: ./structure.x INPUT.XYZ ========     
!     April 04, 2018 :
!     DIL LIMBU, USM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Calculation of g(r) !!!!!!!!!!!
      SUBROUTINE RDFcal(dr, rbin, ri, gr)
      USE DataVal
      IMPLICIT NONE

      integer  :: i,j, bin, rbin

      real(dp) :: dx, dy, dz
      real(dp) :: dr, r2, rc2
      real(dp) :: l2
      real(dp) :: rl, ru, fact, const
      real(dp) :: ri(rbin), gr(rbin)

      rc2 = 2.8          
      nn = 0
      nmap = 0
      gr = 0
      l2 = 0.5*L

      do i = 1, N-1           
         do j = i+1, N
            dx = x(j) - x(i)
            dy = y(j) - y(i)
            dz = z(j) - z(i)
            dx = dx - anint(dx/L)*L
            dy = dy - anint(dy/L)*L
            dz = dz - anint(dz/L)*L   
            r2 = sqrt(dx*dx + dy*dy + dz*dz)

            if(r2 <= l2)then
               bin = int(r2/dr) + 1
               gr(bin) = gr(bin) + 2.0
            endif

            if(r2 <= rc2) then
               nn(i) = nn(i) + 1
               nn(j) = nn(j) + 1
               nmap(nn(i),i) = j
               nmap(nn(j),j) = i
            end if
         end do           
      end do

      open(1, file='nmap.dat', status='unknown')
      write(1,'(A)') '#  atm  nn    NeibhorList'
      do i = 1, N
         write(1,'(I6,I4,8I8)') i, nn(i), (nmap(j,i), j=1,nn(i))
      enddo
      close(1)
         
!  Normalization of gr_i
      rho = N/(L*L*L)
      const = 4.0/3.0*pi*N*rho

      open(1, file='gr.dat', status='unknown')
      write(1,*) '# ', rbin 
      do bin = 1, rbin
         rl = (bin-1.)*dr
         ru = rl + dr
         fact = const*((ru*ru*ru)-(rl*rl*rl))
         ri(bin) = rl+.5*dr
         gr(bin) = gr(bin)/fact
         write(1,'(f8.4,f12.4)') ri(bin),gr(bin)
      enddo
      close(1)

      END SUBROUTINE
