!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Program to Calculate Radial Distribution Function using MPI
!      Initially created on Spring 2016 for Serial version
!      March 21, 2018 :: modified by removing nbin loop
!      December 13, 2020 :: MPI implementation
!
!      DIL LIMBU, USM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program RDF2MPI
      use mpi
      implicit none
      integer, parameter :: nmax = 500000, gmax = 5000
      integer, parameter :: root = 0

      integer  :: i, j, N, nbin
      real     :: pos(3,nmax)
      real     :: gr(gmax), gr_local(gmax)
      real     :: rho, pi, L, r1, r2, dr
      real     :: rl, ru
      real     :: const, fact, tt0, tt1 
      character (len=2)  :: atm
      character (len=20) :: infil, argv

      integer  :: taskid, np, ierr
      integer  :: comm = MPI_COMM_WORLD

      pi = 3.141593

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(comm, np, ierr)
      call MPI_COMM_RANK(comm, taskid, ierr)

      tt0 = MPI_Wtime()

      if(taskid == root)then
         call get_command_argument(1, infil)
         call get_command_argument(2, argv)

         read(argv,*) dr

         write(0,*) '*********  PROGRAM TO CALCULATE g(r)  ************'
         write(0,*) '****************** USING MPI *********************'
         write(0,*) '*******************  USAGE:  *********************'
         write(0,*) '****mpirun -n 4 ./rdf2mpi file.xyz dr(~0.05) *****'
         write(0,*) '************* Output: gr_mpi.dat i****************'

         open (Unit = 1, file = infil, status = 'old')
         read(1,*) N
         read(1,*) L                            
         do i = 1, N
            read (1,*) atm, (pos(j,i), j=1,3)      !Read the coordinate of atoms
         end do
         close(1)

         rho = real(N)/(L*L*L)
         r1 = 0.0  
         r2 = real(L)/2.0
         nbin = int((r2 - r1)/dr)        
      endif

      call MPI_BCAST(N, 1, MPI_INT, root, comm,ierr)
      call MPI_BCAST(pos, 3*N, MPI_REAL, root, comm,ierr)

      call MPI_BCAST(nbin, 1, MPI_INT, root, comm,ierr)
      call MPI_BCAST(dr, 1, MPI_REAL, root, comm,ierr)
      call MPI_BCAST(L, 1, MPI_REAL, root, comm,ierr)

      call RDF2MPIBIN(pos, N, L, gr_local, nbin, dr)

      call MPI_BARRIER(comm, ierr)
      call MPI_REDUCE(gr_local,gr,nbin,MPI_REAL,MPI_SUM, root,comm,ierr)

      if(taskid == root)then 
         const = 4.0*pi*rho*N/3.0
         open(1,file='gr_mpi.dat', status='unknown')
         do i = 1,nbin
            rl = (i-1.0)*dr
            ru = rl + dr
            fact = const*((ru*ru*ru)-(rl*rl*rl))
            gr(i) = gr(i)/fact
            write(1,110) (rl+.5*dr), gr(i)
         end do
         close(1)
      endif
      110 format(1X, 4F12.4)

      tt1 = MPI_Wtime()

      if(taskid .eq. root) write(0,*) 'CPU Time:', (tt1-tt0),'seconds'

      call MPI_Finalize (ierr)

      end program RDF2MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE RDF2MPIBIN(pos, N, L, gr, nbin, dr)
      use mpi
      implicit none
     
      integer :: i, j, N, nbin, bin, istart
      integer :: start_atom, end_atom
      real    :: pos(3,N)
      real    :: gr(nbin)
      real    :: dis, dx, dy, dz
      real    :: L, L2, dr
      real    :: alpha    ! h*a + h*2a+h*3a+ ... = N  ==> a = 2/(p+1); 1.6/(p+1) for better speedup
      real    :: atom_per_proc1

      integer :: taskid, np, ierr
      integer  :: comm = MPI_COMM_WORLD

      call MPI_COMM_SIZE(comm, np, ierr)
      call MPI_COMM_RANK(comm, taskid, ierr)

      atom_per_proc1 = (real(N)/real(np))
      alpha = real(1.65/real(np+1))

      start_atom = 0 
      do i = 1, taskid
          start_atom = start_atom +  int(i*alpha*atom_per_proc1)
      enddo 
      start_atom = start_atom + 1
      end_atom = start_atom + int((taskid+1)*alpha*atom_per_proc1) - 1

      if(taskid == (np-1))then
         end_atom = N-1
      endif
!      start_atom = (taskid*atom_per_proc) + 1
!      end_atom = start_atom + atom_per_proc - 1
!      left_over = mod(N, np)

!      if(left_over .ne. 0)then
!         do i = 0, left_over - 1
!            if(taskid .eq. i) end_atom = end_atom + 1
!            if(i .lt. (np-1))then
!               do j = i+1,np-1
!                  if(taskid .eq. j)then
!                     start_atom = start_atom + 1
!                     end_atom = end_atom + 1
!                  endif
!                enddo
!            endif
!         enddo
!      endif

      gr = 0.0
      L2 = 0.5*L
      do i = start_atom, end_atom
         do j = i+1,N
            dx = pos(1,i) - pos(1,j)
            dy = pos(2,i) - pos(2,j)
            dz = pos(3,i) - pos(3,j)
            dx = dx - anint(dx/L)*L
            dy = dy - anint(dy/L)*L
            dz = dz - anint(dz/L)*L
            dis = sqrt(dx*dx + dy*dy + dz*dz)
            if(dis .lt. L2) then
               bin = int(dis/dr)+1
               gr(bin) = gr(bin) + 2.0                    
            end if
         end do
      end do

      end subroutine
