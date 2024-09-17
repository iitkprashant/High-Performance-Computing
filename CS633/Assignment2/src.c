#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <stdbool.h>

int main( int argc, char *argv[])
{
  int myrank, size, position, count, num_elem=2; 
  MPI_Status status;
  double sTime, eTime, time;

  // input arguments
  int Px = atoi (argv[1]);
  int N2 = atoi (argv[2]);
  int time_step = atoi (argv[3]);
  int seed = atoi(argv[4]);
  int stencil = atoi(argv[5]);
  
  int N = sqrt(N2);  // Number of grid points in one direction per processes

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  // initialize arrays  
  int Py = size/Px;
  int overlap = (stencil-1)/4;
  float array2D[N+overlap*2][N+overlap*2];
  float array2D_0[N][N];
  float temp_r[N][overlap];
  float temp_l[N][overlap];
  float temp_t[overlap][N];
  float temp_b[overlap][N];

  float buff_r[N*overlap];
  float buff_l[N*overlap];
  float buff_t[N*overlap];
  float buff_b[N*overlap];
	  
  float recvtemp_r[N][overlap];
  float recvtemp_l[N][overlap];
  float recvtemp_t[overlap][N];
  float recvtemp_b[overlap][N];

  float recvbuff_r[N*overlap];
  float recvbuff_l[N*overlap];
  float recvbuff_t[N*overlap];
  float recvbuff_b[N*overlap];
  
 // Filling of the array
  for (int i=0; i<N+2*overlap; i++){
   for (int j=0; j<N+2*overlap; j++){
       array2D[i][j] = 0.0;
   }
  }
  srand(seed*(myrank+10));
  for (int i=overlap; i<N+overlap; i++){
   for (int j=overlap; j<N+overlap; j++){
      array2D[i][j] = abs(rand()+(i*rand()+j*myrank))/100; 
   }
  }

  // define neighbour logics
  bool left = true;
  bool right = true;
  bool top = true;
  bool bottom = true;

  if (myrank%Px == Px-1){
   right = false;
  }
  if (myrank%Px == 0){
   left = false;
  }
  if (myrank >= size-Px){
   bottom = false;
  }
  if (myrank < Px){
   top = false;
  }

  sTime = MPI_Wtime();
 
 //time steps: Iteration begins
for (int i = 0; i < time_step; i++){

 // Packing
 position = 0;
 if (right){
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=N; i<N+overlap; i++) {
    temp_r[j-overlap][i-N] = array2D[j][i];
  }
  }
    MPI_Pack (&temp_r, N*overlap, MPI_FLOAT, &buff_r, N*overlap*sizeof(float), &position, MPI_COMM_WORLD);
 }

 position = 0;
 if (left){
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=overlap; i<2*overlap; i++) {
    temp_l[j-overlap][i-overlap] = array2D[j][i];
  }
  }
    MPI_Pack (&temp_l, N*overlap, MPI_FLOAT, &buff_l, N*overlap*sizeof(float), &position, MPI_COMM_WORLD);
 }
 
 position = 0;
 if (top){
  for (int i=overlap; i<2*overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    temp_t[i-overlap][j-overlap] = array2D[i][j];
  }
  }
    MPI_Pack (&temp_t, N*overlap, MPI_FLOAT, &buff_t, N*overlap*sizeof(float), &position, MPI_COMM_WORLD);
 }

 
 position = 0;
 if (bottom){
  for (int i=N; i<N+overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    temp_b[i-N][j-overlap] = array2D[i][j];
  }
  }
    MPI_Pack (&temp_b, N*overlap, MPI_FLOAT, &buff_b, N*overlap*sizeof(float), &position, MPI_COMM_WORLD);
 }


// 2D communication
  if (right)
  {
      // Send/recv right neighbor
      MPI_Send (buff_r, N*overlap*sizeof(float), MPI_PACKED, myrank+1, myrank+1, MPI_COMM_WORLD);
      MPI_Recv (recvbuff_r, N*overlap*sizeof(float), MPI_PACKED, myrank+1, myrank, MPI_COMM_WORLD, &status);
  }

  if (left)
  {
      // Send/recv left neighbor
      MPI_Recv (recvbuff_l, N*overlap*sizeof(float), MPI_PACKED, myrank-1, myrank, MPI_COMM_WORLD, &status);
      MPI_Send (buff_l, N*overlap*sizeof(float), MPI_PACKED, myrank-1, myrank-1, MPI_COMM_WORLD);
  }

  if (bottom)
  {
      // Send/recv bottom neighbor
      MPI_Send (buff_b, N*overlap*sizeof(float), MPI_PACKED, myrank+Px, myrank+Px, MPI_COMM_WORLD);
      MPI_Recv (recvbuff_b, N*overlap*sizeof(float), MPI_PACKED, myrank+Px, myrank, MPI_COMM_WORLD, &status);
  }

  if (top)
  {
      // Send/recv top neighbor
      MPI_Recv (recvbuff_t, N*overlap*sizeof(float), MPI_PACKED, myrank-Px, myrank, MPI_COMM_WORLD, &status);
      MPI_Send (buff_t, N*overlap*sizeof(float), MPI_PACKED, myrank-Px, myrank-Px, MPI_COMM_WORLD);
  }


// Unboxing
 position = 0;
 if (right){
  MPI_Unpack(recvbuff_r, N*overlap*sizeof(float), &position, recvtemp_r, N*overlap, MPI_FLOAT, MPI_COMM_WORLD);
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=N+overlap; i<N+2*overlap; i++) {
    array2D[j][i] = recvtemp_r[j-overlap][i-N-overlap];
  }
  }
 }

 position = 0;
 if (left){
  MPI_Unpack(recvbuff_l, N*overlap*sizeof(float), &position, recvtemp_l, N*overlap, MPI_FLOAT, MPI_COMM_WORLD);
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=0; i<overlap; i++) {
    array2D[j][i] = recvtemp_l[j-overlap][i];
  }
  }
 }

 position = 0;
 if (top){
  MPI_Unpack(recvbuff_t, N*overlap*sizeof(float), &position, recvtemp_t, N*overlap, MPI_FLOAT, MPI_COMM_WORLD);
  for (int i=0; i<overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    array2D[i][j] = recvtemp_t[i][j-overlap];
  }
  }
 }

 position = 0;
 if (bottom){
  MPI_Unpack(recvbuff_b, N*overlap*sizeof(float), &position, recvtemp_b, N*overlap, MPI_FLOAT, MPI_COMM_WORLD);
  for (int i=N+overlap; i<N+2*overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    array2D[i][j] = recvtemp_b[i-N-overlap][j-overlap];
  }
  }
 }

  MPI_Barrier(MPI_COMM_WORLD); 

// Computation
  for (int i=overlap; i<N+overlap; i++){
  for (int j=overlap; j<N+overlap; j++){
    float sum = 0.0;
    for (int k=1; k<overlap+1; k++){
      sum += array2D[i+k][j]+ array2D[i-k][j] + array2D[i][j+k] + array2D[i][j-k];
    }
    array2D_0[i-overlap][j-overlap] = (array2D[i][j]+sum)/stencil;
  }
  }
  // copy back
  for (int i=overlap; i<N+overlap; i++){
  for (int j=overlap; j<N+overlap; j++){
    array2D[i][j] = array2D_0[i-overlap][j-overlap];
  }
  }

}
eTime = MPI_Wtime();
  time = eTime - sTime;
  double finaltime;
  MPI_Reduce(&time, &finaltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!myrank) { 
          printf ("Without Leader: Max time : %lf\n", finaltime);
  }
  
  //***************************************************************************************************//
  //****************** WITH LEADER = 1 ****************************************************************//
  
  // Filling of the array
  for (int i=0; i<N+2*overlap; i++){
  for (int j=0; j<N+2*overlap; j++){
      array2D[i][j] = 0.0;
  }
  }
  srand(seed*(myrank+10));
  for (int i=overlap; i<N+overlap; i++){
  for (int j=overlap; j<N+overlap; j++){
      array2D[i][j] = abs(rand()+(i*rand()+j*myrank))/100; 
  }
  }
  
 
  sTime = MPI_Wtime();
 
//  time steps: Iteration begins
for (int i = 0; i < time_step; i++){

 // Packing
 position = 0;
 if (right){
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=N; i<N+overlap; i++) {
    temp_r[j-overlap][i-N] = array2D[j][i];
  }
  }
    MPI_Pack (&temp_r, N*overlap, MPI_FLOAT, &buff_r, N*overlap*sizeof(float), &position, MPI_COMM_WORLD);
 }

 position = 0;
 if (left){
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=overlap; i<2*overlap; i++) {
    temp_l[j-overlap][i-overlap] = array2D[j][i];
  }
  }
    MPI_Pack (&temp_l, N*overlap, MPI_FLOAT, &buff_l, N*overlap*sizeof(float), &position, MPI_COMM_WORLD);
 }
 

 if (top){
  for (int i=overlap; i<2*overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    temp_t[i-overlap][j-overlap] = array2D[i][j];
  }
  }
 }

 if (bottom){
  for (int i=N; i<N+overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    temp_b[i-N][j-overlap] = array2D[i][j];
  }
  }
 }
 
 // make subcommunicator
 int color, new_rank, new_size;
 MPI_Comm new_comm;
 // Determine the color for each process
 color = myrank / Px;

 // Split the global communicator into subcommunicators based on color
 MPI_Comm_split(MPI_COMM_WORLD, color, myrank, &new_comm);
    
 MPI_Comm_rank(new_comm, &new_rank);
 MPI_Comm_size(new_comm, &new_size);
 
 
 // GATHER 
 float ledtemp_t[overlap*new_size][N];
 float ledtemp_b[overlap*new_size][N];
 float ledrecvtemp_t[overlap*new_size][N];
 float ledrecvtemp_b[overlap*new_size][N];
 
 position = 0;
 if(bottom){
 MPI_Gather(temp_b, N*overlap, MPI_FLOAT, ledtemp_b, N*overlap, MPI_FLOAT, 0, new_comm);
 if(!new_rank){
 MPI_Pack (&ledtemp_b, N*overlap*new_size, MPI_FLOAT, &buff_b, N*overlap*new_size*sizeof(float), &position, MPI_COMM_WORLD);
 }
 }
 
  position = 0;
  if(top){
 MPI_Gather(temp_t, N*overlap, MPI_FLOAT, ledtemp_t, N*overlap, MPI_FLOAT, 0, new_comm);
 if(!new_rank){
 MPI_Pack(&ledtemp_t, N*overlap*new_size, MPI_FLOAT, &buff_t, N*overlap*new_size*sizeof(float), &position, MPI_COMM_WORLD);
 }
 }

// 2D communication
  if (right)
  {
      // Send/recv right neighbor
      MPI_Send (buff_r, N*overlap*sizeof(float), MPI_PACKED, myrank+1, myrank+1, MPI_COMM_WORLD);
      MPI_Recv (recvbuff_r, N*overlap*sizeof(float), MPI_PACKED, myrank+1, myrank, MPI_COMM_WORLD, &status);
  }

  if (left)
  {
      // Send/recv left neighbor
      MPI_Recv (recvbuff_l, N*overlap*sizeof(float), MPI_PACKED, myrank-1, myrank, MPI_COMM_WORLD, &status);
      MPI_Send (buff_l, N*overlap*sizeof(float), MPI_PACKED, myrank-1, myrank-1, MPI_COMM_WORLD);
  }
  
  if (!new_rank){
  if (bottom)
  {
      // Send/recv bottom neighbor
      MPI_Send (buff_b, N*overlap*new_size*sizeof(float), MPI_PACKED, myrank+Px, myrank+Px, MPI_COMM_WORLD);
      MPI_Recv (recvbuff_b, N*overlap*new_size*sizeof(float), MPI_PACKED, myrank+Px, myrank, MPI_COMM_WORLD, &status);
  }

  if (top )
  {
      // Send/recv top neighbor
      MPI_Recv (recvbuff_t, N*overlap*new_size*sizeof(float), MPI_PACKED, myrank-Px, myrank, MPI_COMM_WORLD, &status);
      MPI_Send (buff_t, N*overlap*new_size*sizeof(float), MPI_PACKED, myrank-Px, myrank-Px, MPI_COMM_WORLD);
  }
  }
  
//  // Scatter 

 
 position = 0;
 if(bottom){
 if(!new_rank){
 MPI_Unpack(recvbuff_b, N*overlap*new_size*sizeof(float), &position, ledrecvtemp_b, N*overlap*new_size, MPI_FLOAT, MPI_COMM_WORLD);
 }
 MPI_Scatter(ledrecvtemp_b, N*overlap, MPI_FLOAT, recvtemp_b, N*overlap, MPI_FLOAT, 0, new_comm);
 }
 
  position = 0;
  if(top){
  if(!new_rank){
  MPI_Unpack(recvbuff_t, N*overlap*new_size*sizeof(float), &position, ledrecvtemp_t, N*overlap*new_size, MPI_FLOAT, MPI_COMM_WORLD);  
  }
 MPI_Scatter(ledrecvtemp_t, N*overlap, MPI_FLOAT, recvtemp_t, N*overlap, MPI_FLOAT, 0, new_comm);
 }


// Unboxing
 position = 0;
 if (right){
  MPI_Unpack(recvbuff_r, N*overlap*sizeof(float), &position, recvtemp_r, N*overlap, MPI_FLOAT, MPI_COMM_WORLD);
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=N+overlap; i<N+2*overlap; i++) {
    array2D[j][i] = recvtemp_r[j-overlap][i-N-overlap];
  }
  }
 }

 position = 0;
 if (left){
  MPI_Unpack(recvbuff_l, N*overlap*sizeof(float), &position, recvtemp_l, N*overlap, MPI_FLOAT, MPI_COMM_WORLD);
  for (int j=overlap; j<N+overlap; j++) {
  for (int i=0; i<overlap; i++) {
    array2D[j][i] = recvtemp_l[j-overlap][i];
  }
  }
 }


 if (top){
  for (int i=0; i<overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    array2D[i][j] = recvtemp_t[i][j-overlap];
  }
  }
 }


 if (bottom){
  for (int i=N+overlap; i<N+2*overlap; i++) {
  for (int j=overlap; j<N+overlap; j++) {
    array2D[i][j] = recvtemp_b[i-N-overlap][j-overlap];
  }
  }
 }
 

  MPI_Barrier(MPI_COMM_WORLD); 

// Computation
  for (int i=overlap; i<N+overlap; i++){
  for (int j=overlap; j<N+overlap; j++){
    float sum = 0.0;
    for (int k=1; k<overlap+1; k++){
      sum += array2D[i+k][j]+ array2D[i-k][j] + array2D[i][j+k] + array2D[i][j-k];
    }
    array2D_0[i-overlap][j-overlap] = (array2D[i][j]+sum)/stencil;
  }
  }
  // copy back
  for (int i=overlap; i<N+overlap; i++){
  for (int j=overlap; j<N+overlap; j++){
    array2D[i][j] = array2D_0[i-overlap][j-overlap];
  }
  }

}
 
eTime = MPI_Wtime();
  time = eTime - sTime;
  double finaltime_led;
  MPI_Reduce(&time, &finaltime_led, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (!myrank) { 
          printf ("With Leader : Max time : %lf\n", finaltime_led);
            printf("Data: %f\n",array2D[overlap][overlap]);
  }

  MPI_Finalize();
  return 0;

}

