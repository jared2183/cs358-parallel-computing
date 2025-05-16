/* mm.cpp */

//
// Matrix multiplication implementation, computing C=A*B where A and B
// are NxN matrices. The resulting matrix C is therefore NxN.
// 
// MM is split across the main process and the N-1 worker processes.
//
#include <iostream>
#include <string>
#include <omp.h>
#include <mpi.h>

#include "alloc2D.h"

using namespace std;



#ifdef NEVER
//
// Matrix multiplication:
//
for (int i = 0; i < N; i++)
  for (int j = 0; j < N; j++)
    C[i][j] = 0.0;

for (int i = 0; i < N; i++)
{
  for (int j = 0; j < N; j++)
  {
    for (int k = 0; k < N; k++)
    {
      C[i][j] += (A[i][k] * B[k][j]);
    }
  }
}
#endif



//
// main MPI process:
//
double** main_process(double** A, double** B, int N, int numProcs)
{
  cout << "main starting..." << endl;
  cout.flush();

  // 
  // main process orchestrates the computation of C = A * B. Steps:
  //
  // 
  for (int w = 1; w < numProcs; w++) {
    int dest = w;
    int count = 1;  // number of data type we are receiving (in this case, ints)
    int tag = 0;
    
    // 1. send N --- matrix size --- to each worker
    MPI_Send(&N, count, MPI_INT, dest, tag, MPI_COMM_WORLD);

    // 2. send matrix A to each worker: hint A[0]
    count = N * N;
    // remember A is flattened in memory (rows are next to each other), 
    // so A[0] would point to the beginning of the flattened array
    // can NOT just directly send a pointer to array of pointers (double **) in MPI since pointers would not be valid on another machine
    // can only send direct chunks of memory, NO POINTERS
    MPI_Send(A[0], count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

    // 3. send matrix B to each worker: hint B[0]
    MPI_Send(B[0], count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  }

  //
  // 4. allocate memory for matrix C
  double** C = New2dMatrix<double>(N, N);
  // 
  // 5. instead of waiting, compute chunk size and compute first set of rows
  //
  int chunk_size = N / numProcs;
  int start_row = 0;
  int end_row = start_row + chunk_size;

  for (int i = start_row; i < end_row; i++)
    for (int j = 0; j < N; j++)
      C[i][j] = 0.0;

  for (int i = start_row; i < end_row; i++)
  {
    for (int j = 0; j < N; j++)
    {
      for (int k = 0; k < N; k++)
      {
        C[i][j] += (A[i][k] * B[k][j]);
      }
    }
  }
  // 6. receive results (chunk) from each worker: hint C[chunkSize*worker]
  //
  for (int w = 1; w < numProcs; w++) {
    int src = w;
    int count = chunk_size * N;
    int tag = 0;

    int row = chunk_size * w;

    MPI_Recv(C[row], count, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  // 7. return pointer to result matrix:
  // return C;
}


//
// worker MPI process(es):
//
void worker_process(int myRank, int numProcs)
{
  cout << "worker " << myRank << " starting..." << endl;
  cout.flush();
 
  // 
  // worker process performs a share of C = A * B. Steps:
  //
  // 1. receive N --- matrix size --- from main
  //
  int N;
  int src = 0;
  int count = 1;
  int tag = 0;

  MPI_Recv(&N, count, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
  // 2. allocate memory for matrices A and B
  double** A = New2dMatrix<double>(N, N);
  double** B = New2dMatrix<double>(N, N);
  // 
  // 3. receive matrix A from main: hint A[0]
  // 
  count = N * N;
  MPI_Recv(&A, count, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // 4. receive matrix B from main: hint B[0]
  MPI_Recv(&B, count, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //
  // 5. allocate memory for matrix C
  double** C = New2dMatrix<double>(N, N);
  // 
  // 6. compute chunk size, our start / end rows
  // 
  int chunk_size = N / numProcs;
  int start_row = myRank * chunk_size;
  int end_row = start_row + chunk_size;

  // 7. perform MM for our set of rows
  for (int i = start_row; i < end_row; i++)
    for (int j = 0; j < N; j++)
      C[i][j] = 0.0;
  
  for (int i = start_row; i < end_row; i++)
  {
    for (int j = 0; j < N; j++)
    {
      for (int k = 0; k < N; k++)
      {
        C[i][j] += (A[i][k] * B[k][j]);
      }
    }
  }
  //
  // 8. send our results back to main: hint C[startRow]
  // 
  count = chunk_size * N;
  int dest = 0;

  MPI_Send(N[start_row], count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  // 9. delete memory for A, B, C
  Delete2dMatrix(A);
  Delete2dMatrix(B);
  Delete2dMatrix(C);
}
