F90_FLAG=-O1
#F90_FLAG=-O0 -check bounds -g -traceback
LIB_FLAG=-llammps -ljsonfortran
OBJ_FILE=lammps.o
SRC_LOC=~/softwares/lammps-29Aug2024/src
FC=mpiifx

delta_Fmn: ${OBJ_FILE} consts.o datatypes.o random.o spherical_harmonics.o nucleus.o mcmoves.o montecarlo.o delta_Fmn.o
	${FC} -L${SRC_LOC} ${F90_FLAG} ${LIB_FLAG} ${OBJ_FILE} consts.o datatypes.o random.o spherical_harmonics.o nucleus.o mcmoves.o montecarlo.o delta_Fmn.o -o delta_Fmn

delta_Fmn.o: consts.o datatypes.o montecarlo.o delta_Fmn.f90
	${FC} -fPIC ${F90_FLAG} -c delta_Fmn.f90
montecarlo.o: ${OBJ_FILE} consts.o nucleus.o random.o datatypes.o mcmoves.o  montecarlo.f90
	${FC} -fPIC ${F90_FLAG} -c montecarlo.f90
mcmoves.o: ${OBJ_FILE} consts.o nucleus.o random.o  datatypes.o mcmoves.f90
	${FC} -fPIC ${F90_FLAG} -c mcmoves.f90
nucleus.o: consts.o spherical_harmonics.o datatypes.o nucleus.f90
	${FC} -fPIC ${F90_FLAG} -c nucleus.f90
random.o: consts.o random.f90
	${FC} -fPIC ${F90_FLAG} -c random.f90
spherical_harmonics.o: consts.o spherical_harmonics.f90
	${FC} -fPIC ${F90_FLAG} -c spherical_harmonics.f90
datatypes.o: consts.o datatypes.f90
	${FC} -fPIC ${F90_FLAG} -c datatypes.f90
consts.o: consts.f90
	${FC} -fPIC ${F90_FLAG} -c consts.f90
lammps.o: lammps.f90
	${FC} -fPIC ${F90_FLAG} -c lammps.f90

clean:
	rm -f *.mod *.o delta_Fmn
