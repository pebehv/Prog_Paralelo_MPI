#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include "files.h"
#include <mpi.h>
#include <omp.h>


#define SOFTENING 1e-9f
#define theread_  2


typedef struct { float x, y, z, vx, vy, vz; } Body;

Body *recvbuf;

void bodyForce(Body *p, float dt, int n, Body *recvbuf, int cantN, int my_id, int iter) {

 omp_set_num_threads(theread_); 
 #pragma omp parallel for 
  for (int i = 0; i < cantN; ++i) {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - recvbuf[i].x;
      float dy = p[j].y - recvbuf[i].y;
      float dz = p[j].z - recvbuf[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3;  Fz += dz * invDist3;
    }

    recvbuf[i].vx += dt*Fx; recvbuf[i].vy += dt*Fy; recvbuf[i].vz += dt*Fz;
    printf("\n\n iter: %d -- Proceso  %d, hilo %d ,i: %d - %4.2f, %4.2f, %4.2f\n", iter, my_id,  omp_get_thread_num(),i, recvbuf[i].x,recvbuf[i].y,recvbuf[i].z);

  }
}


int main( int argc,  char** argv) {

  int i;
	int my_id, nproc;
  Body *p;
  int bytes;
  double totalTime = 0.0;
  const float dt = 0.01f; // Time step
  const int nIters = 10;  // Simulation iterations
  const char * initialized_values;
  const char * solution_values;
  float *buf;
  int nBodies;
  
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  
  if(my_id == 0){  
       nBodies = 8;
       //nBodies = 2<<11;
       if (argc > 1) nBodies = 2<<atoi(argv[1]);
      bytes = nBodies * sizeof(Body);
      buf = (float *)malloc(bytes);
     
       if (nBodies == 2<<11) {
         initialized_values = "files/initialized_4096";
         solution_values = "files/solution_4096";
       } else { // nBodies == 2<<15
         initialized_values = "files/initialized_65536";
         solution_values = "files/solution_65536";
       }
     
     
       if (argc > 2) initialized_values = argv[2];
       if (argc > 3) solution_values = argv[3];
     
       p = (Body*)buf;

      read_values_from_file(initialized_values, buf, nBodies, sizeof(Body) );

      for (int i = 0 ; i < nBodies; i++) { // integrate position
        
       printf("i: %d - %4.2f, %4.2f, %4.2f",i, p[i].x,p[i].y,p[i].z);
      }
      printf("\n");
     

  }

  MPI_Bcast(&nBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int bytes2 = (nBodies/nproc) * sizeof(Body);
  float *buf2;
  buf2 = (float *)malloc(bytes2);

  if(my_id != 0){
    buf = (float *)malloc(bytes);
    p = (Body*)buf;
  }
  recvbuf = (Body*)buf2;
  for (int iter = 0; iter < 10; iter++) {
    
    StartTimer(); 
   
   
   MPI_Bcast(p, bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
   MPI_Scatter(p, bytes2, MPI_BYTE, recvbuf, bytes2, MPI_BYTE, 0 , MPI_COMM_WORLD);
   
   bodyForce(p, dt, nBodies, recvbuf, (nBodies/nproc), my_id, iter ); // compute interbody forces
   
   
   for (int i = 0 ; i < nBodies; i++) { // integrate position
    p[i].x += p[i].vx*dt;
    p[i].y += p[i].vy*dt;
    p[i].z += p[i].vz*dt;
  }
  printf("\n");
   MPI_Allgather(recvbuf, bytes2, MPI_BYTE, p, bytes2, MPI_BYTE, MPI_COMM_WORLD);
  
	   
   for (int i = 0 ; i < nBodies; i++) { // integrate position
    p[i].x += p[i].vx*dt;
    p[i].y += p[i].vy*dt;
    p[i].z += p[i].vz*dt;
   
    //printf("\n\n AllG-- iter: %d -- Proceso  %d, hilo %d ,i: %d - %4.2f, %4.2f, %4.2f\n", iter, my_id,  omp_get_thread_num(),i, p[i].x,p[i].y,p[i].z);
  }
  printf("\n");
    
  
      const double tElapsed = GetTimer() / 1000.0;
      totalTime += tElapsed;
  }

  if (my_id == 0) {

    double avgTime = totalTime / (double)(nIters);
    float billionsOfOpsPerSecond = 1e-9 * nBodies * nBodies / avgTime;
    write_values_to_file(solution_values, buf, nBodies);
  }

  MPI_Finalize(); 
  free(buf);
}
