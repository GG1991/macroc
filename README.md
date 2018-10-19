
Readme
------

MacroC is a FE code used to model composite material macrostructures. It works
jointly with MicroPP library and performs FE2 calculations. This is a procedure
where in each Gauss point of the macroscopic mesh a FE calculation is taken into
account to get the average properties that has the heterogeneous macrostructure.

Compilation
-----------

1. Download a PETSc Library from [www.mcs.anl.gov/petsc](www.mcs.anl.gov/petsc) and install.
2. Set the environmental variables `PETSC_DIR` and `PETSC_ARCH` with

```bash
   export PETSC_DIR=<path>
   export PETSC_ARCH=<architecture>
```

Build MacroC with CMake:
-----------------------

1. Clone the repository [GitHub](https://github.com/GG1991/macroc)
2. `cd <cloned directory>`
3. `mkdir build` (can be also `build + anything`)
4. `cd build`
5. `cmake ..`
6. `make`

To build the optimized version:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

and the debug version:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Execution
---------

The mesh size and the processors in each direction are given by command line as:

```
-da_grid_x <int> - number of grid points in x direction
-da_grid_y <int> - number of grid points in y direction
-da_grid_z <int> - number of grid points in z direction
-da_processors_x <int> number of processors in x direction
-da_processors_y <int> number of processors in y direction
-da_processors_z <int> number of processors in z direction
-ts <int> number of time steps
-dt <double> time step
-lx <double> length in x direction
-ly <double> length in y direction
-lz <double> length in z direction
```

Examples
--------

```bash
mpirun -np 1 ./macroc -da_grid_x 4 -da_grid_y 4 -da_grid_z 2
```
