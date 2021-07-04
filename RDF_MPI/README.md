THIS IS THE PARALLEL IMPLEMENTATION USING MPI FOR
PAIR CORRELATION FUCNTION, g(r), CALCULATION


DIL LIMBU, USM

TO COMPILE:: mpiifort -O2 rdf2mpi.f90 -o rdf2mpi

TO RUN:: mpirun  -n  4  ./rdf2mpi  INPUT.XYZ rbin_width(~0.05)

OUTPUTS:: 
           gr_mpi.dat   <- Pair Correlation Function
