
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
-da_grid_x <int> - number of grid points in x direction (10)
-da_grid_y <int> - number of grid points in y direction (10)
-da_grid_z <int> - number of grid points in z direction (10)
-da_processors_x <int> number of processors in x direction
-da_processors_y <int> number of processors in y direction
-da_processors_z <int> number of processors in z direction
-micro_n <int> number nodes in each direction of microscopic mesh
-ts <int> number of time steps (1)
-dt <double> time step (0.001)
-lx <double> length in x direction (10.0)
-ly <double> length in y direction (1.0)
-lz <double> length in z direction (1.0)
-vtu_freq <int> frequency of VTU plots (-1)
-new_its <int> maximum number of Newton Iterations
-new_tol <double> minimum tolerance of Newton residue
```

Examples
--------

```bash
mpirun -np 1 ./macroc -da_grid_x 4 -da_grid_y 4 -da_grid_z 2
```

Testing
---------

In order follow DevOps a testing methodology is used. **MacroC** uses CTest toolt to perform tests. The recommended use of it consists in going first to a stable version of the code in the master branch and execute:

```bash
ctest --VV -O output_stable.dat
```

Then in the file `output_stable.dat` the stable reference output would have be written for all the tests cases. Finally after doing modifications if it is necessary run the tests cases and compare the solution running the same command but changing the output file name, for example:

```bash
ctest -VV -O output_develop.dat
```

You can compare the solution using Linux tools like `diff` or `vimdiff`.
