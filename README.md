# NACA four-digit tool for collecting lift, drag, and moment coefficients
Compile code in basilisk with MPI
```
CC99='mpicc -std=c99' qcc -autolink -Wall -O2 -D_MPI=1 naca.c -o naca -lm -lfb_tiny
```
Run simulation with desired parameters with MPI
```
mpirun -np <# of cores> ./naca <NACA> <AoA> <Re> <LEVEL>
```
