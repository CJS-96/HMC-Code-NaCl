# Parallel Code usage

This version of the code launches multiple simulations in parallel, each with its own set of input settings and configuration as described below:
1. input.json.${SimNo}: As many files as the number of simulations to be launched. Each file contains per simulation settings and is read by the code. See the input.json files in examples/ to know what settings are available.
2. data.lmp.${SimNo}: As many files as the number of simulations to be launched. Each contains the starting configuration for each simulation.
3. pair.lmp: This contains the force-field parameters for the NaCl-water system. Single file.
4. in.1.lmp and in.2.lmp: Contains LAMMPS settings, only 2 files.

The parallel simulation is launched in the following manner:

mpiexec -np ${TotalSimulations}*${ProcessorsPerSimulations} ./delta_Fmn1 ${TotalSimulations}
