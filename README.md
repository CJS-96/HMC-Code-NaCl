# HMC-Code-NaCl

This repository contains the hybrid monte-carlo codes used to study NaCl crystallisation from aqueous solutions.

- `HMC_Parallel_Code/` — parallel version
- `HMC_Sequential_Code/` — sequential version

Each folder is self-contained and builds independently via its own `Makefile`.

# Working Principle

This code invokes LAMMPS to run a short MD simulation to generate new configuration of the system and accepts or rejects it using the metropolis criterion based on detailed balance.

# How to use:

The code is available 
