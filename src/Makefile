FC=ifort

FFLAGS=-traceback -O2 -fPIC 

FILES= 	Structure.f90 \
	DataVal.f90 \
	RDFcal.f90 \
	BADcal.f90 \
	StructureFactorCal.f90\
	guassquadint.f90\
	gauleg.f90 

DataVal.o: DataVal.f90
	${FC} -c $< -o $@
	${FC} ${FFLAGS} ${FILES} -o structure.x

clean:
	@echo "==> Cleaning ALL intermediate and executable files"	
	rm -f *.mod *.o *.x 
