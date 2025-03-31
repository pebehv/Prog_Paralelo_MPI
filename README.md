# Prog_Paralelo_MPI
Ejercicio de practica Usando el lenguaje C , MPI y OpenMPI


camandos :


mpicc -fopenmp -o nbody 01-nbody.c -lm

mpirun -np 2 ./nbody
