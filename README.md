# HMC-Code-NaCl

This repository contains the hybrid monte-carlo codes used to study NaCl crystallisation from aqueous solutions.

# Working Principle

These codes invoke LAMMPS to run a short MD simulation for generating new configuration of the system and accept or reject it using the metropolis criterion satisfying the detailed balance. It implements the use of biasing potentials for the crystallinity (m) and solvation state (s) along with ensuring that the cluster size remains constant during the simulation.

# Code Versions

- `HMC_Parallel_Code/` — parallel version: This version runs multiple simulations at the same time, each with its own initial configuration.
- `HMC_Sequential_Code/` — sequential version: This version runs simulations one after the other, using the final configuration of the previous simulation as the starting point.

See the README.md files in the respective folders to understand how to use it.
Note: Each folder is self-contained and builds independently via its own `Makefile`.

# Compilation
