# ============================================================
# Machine/user-specific settings — edit these for a new setup
# ============================================================
FC=mpiifx                                          # Fortran/MPI compiler wrapper (e.g. mpiifx, mpiifort, mpif90)
LAMMPS_SRC=$(HOME)/softwares/lammps-29Aug2024/src   # LAMMPS source/build dir providing liblammps (path + version)
JSONFORTRAN_INC=$(HOME)/.local/include              # dir containing json-fortran's json_module.mod
JSONFORTRAN_LIB=$(HOME)/.local/lib64                # dir containing json-fortran's libjsonfortran.so

F90_FLAG=-O1
#F90_FLAG=-O0 -check bounds -g -traceback

# ============================================================
# Build internals — derived from the settings above, should not
# need to change between machines/users
# ============================================================
INC_FLAG=-I$(JSONFORTRAN_INC)
LIB_FLAG=-L$(LAMMPS_SRC) -L$(JSONFORTRAN_LIB) -Wl,-rpath,$(LAMMPS_SRC) -Wl,-rpath,$(JSONFORTRAN_LIB) -llammps -ljsonfortran
OBJ_FILE=lammps.o

delta_Fmn: ${OBJ_FILE} consts.o datatypes.o random.o spherical_harmonics.o nucleus.o mcmoves.o montecarlo.o delta_Fmn.o
	${FC} ${F90_FLAG} ${OBJ_FILE} consts.o datatypes.o random.o spherical_harmonics.o nucleus.o mcmoves.o montecarlo.o delta_Fmn.o ${LIB_FLAG} -o delta_Fmn

delta_Fmn.o: consts.o datatypes.o montecarlo.o delta_Fmn.f90
	${FC} -fPIC ${F90_FLAG} -c delta_Fmn.f90
montecarlo.o: ${OBJ_FILE} consts.o nucleus.o random.o datatypes.o mcmoves.o montecarlo.f90
	${FC} -fPIC ${F90_FLAG} ${INC_FLAG} -c montecarlo.f90
mcmoves.o: ${OBJ_FILE} consts.o nucleus.o random.o datatypes.o mcmoves.f90
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
