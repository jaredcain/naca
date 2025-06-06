# NACA four-digit tool for aerodynamic performance data, for use in basilisk
Ensure you are using the bash file
```
source ~/.bashrc
```
Change to the proper directory, either the embed or gcm directory
```
cd ~/basilisk/src/examples/<embed or gcm>
```
Compile code in basilisk with MPI
```
CC99='mpicc -std=c99' qcc -autolink -disable-dimensions -Wall -O2 -D_MPI=1 naca.c -o naca -lm
```
Run simulation with desired parameters with MPI
```
mpirun -np <# of cores> ./naca <NACA> <AoA> <Re> <LEVEL>
```
For example, to simulate a NACA 0012 airfoil at a 2 degree angle of attack, at level 13, run
```
CC99='mpicc -std=c99' qcc -autolink -disable-dimensions -Wall -O2 -D_MPI=1 naca.c -o naca -lm
mpirun -np 4 ./naca 0012 2 10000 13
```
In a separate command window you can see live logs by running
```
tail -f <NACA>_<AoA>_<Re>_<LEVEL>.log
```
So the previous example, you would run
```
tail -f 0012_2_10000_13.log
```
