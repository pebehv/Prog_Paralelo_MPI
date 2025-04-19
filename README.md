# Prog_Paralelo_MPI

Este fue el proyecto final de la universidad para la materia `Fundamentos de la Programación Paralela`. En este proyecto, se utilizó el lenguaje `C` junto con `MPI` y `OpenMP`, lo que permite la programación paralela mediante hilos o procesos.


## Instalacion 

Para descargar el codigo, debe ejecutar:

git clone https://github.com/pebehv/Prog_Paralelo_MPI.git

Una vez que hayas descargado el código, se debe copilar el programa con el siguiente comando `mpicc -fopenmp -o nbody 01-nbody.c -lm`  y luego  ejecutar el programa con el siguiente comando `mpirun -np 2 ./nbody` .  

## Observación
El comando proporcionado anteriormente para ejecutar el programa `mpirun -np 2 ./nbody`  utiliza 2 hilos. Puedes cambiar el número de hilos que desees utilizar para ejecutar el programa.

