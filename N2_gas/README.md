# MassCSS-OMPC: A high performance Collision Cross-Section software for Cluster execution

MassCCS-OMPC is a extended version of ["MassCCS"](https://github.com/cces-cepid/massccs) to perform culations at the inter-node level using a novel programming model, called ["OpenMP Cluster (OMPC)"](https://ompcluster.gitlab.io/). OMPC is a task-parallel model that extends OpenMP for cluster programming, allowing applications to use the same programming model to exploit intra- and inter-node parallelism. That simplifies the development and software maintenance of HPC applications. 

This extension speedup even more the CCS calculations allowing expensive CCS calculations using nitrogen buffer for large systems such as human adenovirus with ~11 million atoms in just 4 min. 

If you use MassCCS-OMPC in your research please, cite the following papers:

["S. Cajahuaringa, L. N. Zanotto,S. Rigo, H. Yviquel, M. S. Skaf and G. Araujo, Ion-Molecule Collision Cross-Section Calculations Using Trajectory Parallelization in Distributed Systems, J.
Parallel Distr. Com. 191 (2024) 104902"](https://doi.org/10.1016/j.jpdc.2024.104902)

## The repository contents:
[`src`](src): This directory contains massccs source codes.

[`paper`](paper): This directory contains inputs, scripts and outputs of published paper.

## Installation 

Download the MassCCS-OMPC or clone the repository on your computer:
 
```bash
git clone https://github.com/cces-cepid/massccs-ompc.git
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
Singularity> make  # Compile the program to create the massccs-ompc executable
```

## Run

Finally, after compiled, the program can be launched using MPI's utility mpirun passing the number of processes as argument. Other flags of mpirun can be used normally (e.g. --hostfile).

In this case is used 3 physical nodes but only 2 nodes are used as worker nodes
```bash
cd .. # need to execute massccs from root dir
mpirun -np 3 singularity exec runtime_latest.sif ./build/bin/massccs-ompc input_Hellium.json > He_04_1ubq_charge_4e_min.log
```
where the input file input_Hellium.json define only 2 tasks for 2 workes nodes, thats means one task for node.

 ```bash
cd .. # need to execute massccs from root dir
mpirun -np 3 singularity exec runtime_latest.sif ./build/bin/massccs-ompc input_Nitrogen.json > N2_04_1ubq_charge_4e_min.log
```
where the input file input_Nitrogen.json define 4 tasks for 2 workes nodes available, that is a dynamics task schedule.




Author & Contact:
--------------
Samuel Cajahuaringa - samuelcm@unicamp.br

