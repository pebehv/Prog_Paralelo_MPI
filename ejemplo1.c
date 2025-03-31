#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include "files.h"
#include <mpi.h>
#include <omp.h>


#define SOFTENING 1e-9f
#define theread_  10

/*
* Each body contains x, y, and z coordinate positions,
* as well as velocities in the x, y, and z directions.
*/

typedef struct { float x, y, z, vx, vy, vz; } Body;

Body *recvbuf;
/*
 * Calculate the gravitational impact of all bodies in the system
 * on all others.
 */

void bodyForce(Body *p, float dt, int n) {
  //printf("bodyForce");
  for (int i = 0; i < n; ++i) {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3;  Fz += dz * invDist3;
    }

    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
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
       if (argc > 1) nBodies = 2<<atoi(argv[1]);
      bytes = nBodies * sizeof(Body);
      buf = (float *)malloc(bytes);
     // recvbuf = (Body*)buf;
       // The assessment will test against both 2<11 and 2<15.
       // Feel free to pass the command line argument 15 when you generate ./nbody report files
       //int nBodies = 2<<11;
     
       // The assessment will pass hidden initialized values to check for correctness.
       // You should not make changes to these files, or else the assessment will not work.
       /*const char * initialized_values;
       const char * solution_values;*/
     
       if (nBodies == 2<<11) {
         initialized_values = "files/initialized_4096";
         solution_values = "files/solution_4096";
       } else { // nBodies == 2<<15
         initialized_values = "files/initialized_65536";
         solution_values = "files/solution_65536";
       }
     
     
       if (argc > 2) initialized_values = argv[2];
       if (argc > 3) solution_values = argv[3];
     
       /*const float dt = 0.01f; // Time step
       const int nIters = 10;  // Simulation iterations*/
     
      // bytes = nBodies * sizeof(Body);
       //int bytes = nBodies * sizeof(Body);
       //float *buf;
     
       //buf = (float *)malloc(bytes);
     
       p = (Body*)buf;
      // Body *p = (Body*)buf;
     
       read_values_from_file(initialized_values, buf, nBodies);
       printf("Proceso Padre\n");

      for (int i = 0 ; i < nBodies; i++) { // integrate position
        p[i].x += p[i].vx*dt;
        p[i].y += p[i].vy*dt;
        p[i].z += p[i].vz*dt;
       printf("i: %d - %4.2f, %4.2f, %4.2f",i, p[i].x,p[i].y,p[i].z);
      }
      printf("\n");
      // read_values_from_file(initialized_values, buf, bytes);
     
       //printf("X,Y,Z,VX,VY,VZ\n");
        // for (int i=0;i<nBodies;i++)
         //printf("%10.2f, %10.2f, %10.2f\n",p[i].x,p[i].y,p[i].z);
         //printf("%10.2f, %10.2f, %10.2f, %10.2f, %10.2f, %10.2f \n",p[i].x,p[i].y,p[i].z,p[i].vx,p[i].vy,p[i].vz);
     
     //exit(0);
       //double totalTime = 0.0;
     
       /*
        * This simulation will run for 10 cycles of time, calculating gravitational
        * interaction amongst bodies, and adjusting their positions to reflect.
        */

  }

  MPI_Bcast(&nBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&buf, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  bytes = nBodies * sizeof(Body);
  buf = (float *)malloc(bytes);
  recvbuf = (Body*)buf;
  MPI_Scatter(p, nBodies/nproc, MPI_BYTE, recvbuf, nBodies/nproc, MPI_BYTE, 0 , MPI_COMM_WORLD);
	printf("Proceso  %d\n", my_id);
  omp_set_num_threads(theread_); 
  
#pragma omp parallel for shared(totalTime) schedule (static, 2  )

  for (int iter = 0; iter < nIters; iter++) {
    #pragma omp critical
    {

      StartTimer(); 
  
    /*
     * You will likely wish to refactor the work being done in `bodyForce`,
     * and potentially the work to integrate the positions.
     */
  
      bodyForce(recvbuf, dt, nBodies); // compute interbody forces
  
    /*
     * This position integration cannot occur until this round of `bodyForce` has completed.
     * Also, the next round of `bodyForce` cannot begin until the integration is complete.
     */
  
      for (int i = 0 ; i < nBodies/nproc; i++) { // integrate position
        recvbuf[i].x += recvbuf[i].vx*dt;
        recvbuf[i].y += recvbuf[i].vy*dt;
        recvbuf[i].z += recvbuf[i].vz*dt;
        if(i==0){
  
          printf("proceso %d, hilo %d :::%10.2f, %10.2f, %10.2f\n", my_id, omp_get_thread_num(),recvbuf[i].x,recvbuf[i].y,recvbuf[i].z);
        }
      }
  
      const double tElapsed = GetTimer() / 1000.0;
      totalTime += tElapsed;
    }
  }


  MPI_Gather(recvbuf, nBodies/nproc, MPI_BYTE, p, nBodies/nproc, MPI_BYTE,  0, MPI_COMM_WORLD);

  if (my_id == 0) {


    double avgTime = totalTime / (double)(nIters);
    float billionsOfOpsPerSecond = 1e-9 * nBodies * nBodies / avgTime;
    write_values_to_file(solution_values, buf, nBodies);
  }

  //read_values_from_file(solution_values, buf, nBodies);


 /* printf("X,Y,Z,VX,VY,VZ\n");
  for (int i=0;i<nBodies;i++)
    printf(":%10.2f, %10.2f, %10.2f\n",p[i].x,p[i].y,p[i].z);*/


  // You will likely enjoy watching this value grow as you accelerate the application,
  // but beware that a failure to correctly synchronize the device might result in
  // unrealistically high values.
//  printf("%0.3f Billion Interactions / second\n", billionsOfOpsPerSecond);

MPI_Finalize(); 
  free(buf);
}
