# Optimization potential interaction parameters for ion mobility spectrometry

This repository contains a set of source code and input scripts that allows to optimization of Lennard-Jones parameters to perform collision cross section
(CCS) calculations. We developed optimization methods for force field parameters in order to improve ["MassCCS"](https://github.com/cces-cepid/massccs) accuracy. Our optimization method uses the differential evolution algorithm and ["OMPC"](https://ompcluster.gitlab.io/) to reduce computation time.

## The repository contents:
[`src`](src): This directory contains massccs-ompc source codes to perform a bunch of CCS calculations using MassCCS and OMPC.

[`N2_gas`](paper): This directory contains inputs, scripts and outputs for optimization Lennard-Jones parameters using N2 buffer gas.

[`CO2_gas`](paper): This directory contains inputs, scripts and outputs for optimization Lennard-Jones parameters using CO2 buffer gas.

## Installation 

Download the MassCCS-OMPC or clone the repository on your computer:
 
```bash
git clone https://github.com/cces-cepid/opt-massccs.git
```
## Required Software

MassCCS-OMPC depends on the following software:

* Singularity
* OpenMP Cluster container image 
* MPI compiler

On Ubuntu/Debian, you can simply run the following commands to dowlonad the last version of OMPC:
```bash
sudo apt-get install singularity-container  
singularity pull docker://ompcluster/runtime
```

## Installing

On your terminal, run the following commands to compile the source code:

```bash
$ cd massccs-ompc
$ singularity shell ./runtime_latest.sif # execute the OMPC container
Singularity> mkdir build # Create build directory 
Singularity> cd build
Singularity> export CC=clang CXX=clang++ # export Clang compiler 
Singularity> cmake .. # Generate Makefiles
Singularity> make  # Compile the program to create the massccs-ompc executable called massccs
```

Author & Contact:
--------------
Samuel Cajahuaringa - samuelcm@unicamp.br

