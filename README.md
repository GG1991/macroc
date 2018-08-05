
Readme
------

MacroPP is a FE code used to model composite material macrostructures. It works
jointly with MicroPP library and performs FE2 calculations. This is procedure
where in each Gauss point of the macroscopic mesh a FE calculation is taken into
account to get the average properties that has the heterogeneous macrostructure.

Build steps with CMake:
-----------------------

1. Clone the repository 
2. cd cloned directory
3. mkdir build (can be also build+anything)
4. cd build
5. cmake .. (important the 2 points)
6. make
