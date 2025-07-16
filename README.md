# GCRYSTAL
This repository contains a first development release of the CRYSTAL code aimed at incorporating GPU support to improve efficiency on modern HPC architectures for ab initio solid-state materials simulations.

Two versions are provided. The first version, GCRYSTAL_1.0, directly integrates CUDA math libraries for accelerated linear algebra operations, building on the standard parallel version of the CRYSTAL code (PCRYSTAL). GCRYSTAL_1.0 has undergone extensive testing on accelerated HPC systems.

The second version, GCRYSTAL_2.0, is currently under development and involves a more extensive refactoring of the code. Its aim is to enhance the use of accelerated libraries by minimizing input/output operations and reducing frequent data transfers between the host and device. This version has to be regarded as an alpha release, as it has undergone limited testing, and not all features are fully implemented.

Both versions are available as executables, along with the source code for the newly developed code, and have limitations regarding the number of MPI processes and number of atoms in the unit cell. Additional details can be found in the relative sections below.

##July 2025 Update
The GCRYSTAL_2.0 version has receved an update (here GCRYSTAL_2.5), the code base is the same, but new features have been implemented: data memory storage options, construction of the density matrix on device, new k point distribution to processes. Also the code is now compatible with more advanced computational drivers, such as geometry optimization, molecular dynamics, frequencies calculaions, etc.

## Prerequisites
Install NVIDIA HPC SDK 23.9 for Linux x86_64 architecture
- [NVIDIA HPC SDK 23.9](https://developer.nvidia.com/nvidia-hpc-sdk-239-downloads)

## Using GCRYSTAL
To run the code the input file must be renamed “INPUT”, and the code must be executed in the directory of this file, as follows:
```
cp test_quartz.d12 INPUT
mpirun -np <number of processes between 1 and 2> /path/to/GCRYSTAL_X.0_executables
```

## License

This project is licensed under the MIT License
