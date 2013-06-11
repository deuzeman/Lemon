#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int rank = 0;
  int size = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("Initialized node %d of %d.\n", rank + 1, size);

  int data[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int blocklens[8] = {1, 1, 1, 1, 1, 1, 1, 1};
//  int indices[8] = {7, 6, 5, 4, 3, 2, 1, 0};
  int indices[8] = {0, 1, 1, 1, 2, 2, 2, 3};

  MPI_Datatype reverse;

  MPI_Type_indexed(8, blocklens, indices, MPI_INT, &reverse);
  MPI_Type_commit(&reverse);
  MPI_Status status;

  if (rank == 0)
  {
    for (int idx = 0; idx < 8; ++idx)
      data[idx] = idx + 1;
    MPI_Send(data, 1, reverse, 1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Recv(data, 1, reverse, 0, 0, MPI_COMM_WORLD, &status);

  for (int ctr = 0; ctr < size; ++ctr)
  {
    if (ctr == rank)
    {
      printf("On node %d:\n  [", rank + 1);
      for (int idx = 0; idx < 7; ++idx)
        printf("%d, ", data[idx]);
      printf("%d]\n", data[7]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Type_free(&reverse);
  MPI_Finalize();
  return 0;
}
