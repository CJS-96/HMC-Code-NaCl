# HMC-Code-NaCl

This repository contains the hybrid monte-carlo codes used to study NaCl crystallisation from aqueous solutions.

## Working Principle

These codes invoke LAMMPS to run a short MD simulation for generating new configuration of the system and accept or reject it using the metropolis criterion satisfying the detailed balance. It implements the use of biasing potentials for the crystallinity (m) and solvation state (s) along with ensuring that the cluster size remains constant during the simulation.

## Code Versions

- `HMC_Parallel_Code/` — parallel version: This version runs multiple simulations at the same time, each with its own initial configuration.
- `HMC_Sequential_Code/` — sequential version: This version runs simulations one after the other, using the final configuration of the previous simulation as the starting point.

See the README.md files in the respective folders to understand how to use it.

## Compilation

### LAMMPS Pre-requisites:
1. LAMMPS should be built with -DBUILD_SHARED_LIBS=yes in a separate build/ directory using cmake.
2. Change the variables in Makefile according to the installation location.
3. The current code is tied to "lammps-29Aug2024" version. If there is a version mismatch, please use fortran/lammps.f90 of the LAMMPS version installed.

### JSON-Fortran:
1. Simple installation of JSON-Fortran is sufficient.
2. Appropriate paths should be set in Makefile such that the json module files and libraries are readily available.

### HMC-Code:
1. Each version of the code is self-contained and builds independently via its own `Makefile`.
2. To compile, just run 'make' in the directory itself. This generates the executable "delta_Fmn".
3. To clean up a previous complilation, run 'make clean'. This will delete all the module and object files as well as the executable.
