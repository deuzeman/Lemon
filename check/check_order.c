/*****************************************************************************
 * LEMON v1.1                                                                *
 *                                                                           *
 * This file is part of the LEMON implementation of the SCIDAC LIME format.  *
 *                                                                           *
 * It is based directly upon the original c-lime implementation,             *
 * as maintained by C. deTar for the USQCD Collaboration,                    *
 * and inherits its license model and parts of its original code.            *
 *                                                                           *
 * LEMON is free software: you can redistribute it and/or modify             *
 * it under the terms of the GNU General Public License as published by      *
 * the Free Software Foundation, either version 3 of the License, or         *
 * (at your option) any later version.                                       *
 *                                                                           *
 * LEMON is distributed in the hope that it will be useful,                  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 * GNU General Public License for more details. You should have received     *
 * a copy of the GNU General Public License along with LEMON. If not,        *
 * see <http://www.gnu.org/licenses/>.                                       *
 *                                                                           *
 * LEMON was written for the European Twisted Mass Collaboration.            *
 * For support requests or bug reports, please contact                       *
 *    A. Deuzeman (deuzeman@itp.unibe.ch)                                    *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>

#include "../include/lemon.h"

typedef int site_t[4];

int main(int argc, char **argv)
{
  int write = 0;
  int mpi_dims[] = {0, 0, 0, 0};
  
  if (argc == 1)
    write = 1;

  if (argc == 5)
  {
    write = atoi(argv[1]);
    mpi_dims[1] = atoi(argv[2]);
    mpi_dims[2] = atoi(argv[3]);
    mpi_dims[3] = atoi(argv[4]);
  }
   
  int mpi_size;
  int rank;
  int mpi_coords[4];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  

  int periods[] = {1, 1, 1, 1};
  int const mapping[] = {0, 3, 2, 1};
  
 
  MPI_Comm cartesian;
  MPI_Dims_create(mpi_size, 4, mpi_dims);

  MPI_Cart_create(MPI_COMM_WORLD, 4, mpi_dims, periods, 1, &cartesian);
  
  MPI_Comm_rank(cartesian, &rank);
  MPI_Cart_coords(cartesian, rank, 4, mpi_coords);

  /* We work on a fixed 32^4 lattice volume. */
  int lat_dims[] = {32, 32, 32, 32};
  int lat_vol = 32 * 32 * 32 * 32;
 
  int loc_dims[4];
  for (int idx = 0; idx < 4; ++idx)
    loc_dims[idx] = lat_dims[idx] / mpi_dims[idx];
  int loc_vol = lat_vol / mpi_size;

  site_t *data = (site_t*)malloc(sizeof(site_t) * loc_vol);;
  int ctr = 0;
  for (int t = loc_dims[0] * mpi_coords[0]; t < loc_dims[0] * mpi_coords[0] + loc_dims[0]; ++t)
    for (int z = loc_dims[1] * mpi_coords[1]; z < loc_dims[1] * mpi_coords[1] + loc_dims[1]; ++z)
      for (int y = loc_dims[2] * mpi_coords[2]; y < loc_dims[2] * mpi_coords[2] + loc_dims[2]; ++y)
        for (int x = loc_dims[3] * mpi_coords[3]; x < loc_dims[3] * mpi_coords[3] + loc_dims[3]; ++x)
        {
          data[ctr][0] = t;
          data[ctr][1] = z;
          data[ctr][2] = y;
          data[ctr][3] = x;
          ++ctr;
        }

  MPI_File fp;  
  char const *type;
  
  if (!rank)
  {
    printf("Write = %d\n", write);
    printf("Dims  = %d, %d, %d, %d\n", mpi_dims[0], mpi_dims[1], mpi_dims[2], mpi_dims[3]);
  }
  if (write)
  {
    LemonWriter *w;
    LemonRecordHeader *h;
   
    int ME_flag=1, MB_flag=1, status=0;

    /* Start of code - writing */
    /* Note that the following is the only way to truncate the file with MPI */
    MPI_File_open(cartesian, "parallel_32x4.test", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);
    MPI_File_set_size(fp, 0);
    w = lemonCreateWriter(&fp, cartesian);

    h = lemonCreateHeader(MB_flag, ME_flag, "parallel-test", lat_vol * sizeof(site_t));
    status = lemonWriteRecordHeader(h, w);
    lemonDestroyHeader(h);

    lemonWriteLatticeParallelMapped(w, data, sizeof(site_t), lat_dims, mapping);

    lemonWriterCloseRecord(w);
    lemonDestroyWriter(w);
    MPI_File_close(&fp);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  /* Reading */
  
  LemonReader *r;

  MPI_File_open(cartesian, "parallel_32x4.test", MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
  r = lemonCreateReader(&fp, cartesian);

  if (lemonReaderNextRecord(r))
    fprintf(stderr, "Node %d reports: next record failed.\n", rank);

  type = lemonReaderType(r);
  if (strncmp(type, "parallel-test", 13))
    fprintf(stderr, "Node %d reports: wrong type read.\n", rank);

  lemonReadLatticeParallelMapped(r, data, sizeof(site_t), lat_dims, mapping);

  int failure = 0;
  ctr = 0;
  for (int t = loc_dims[0] * mpi_coords[0]; t < loc_dims[0] * (mpi_coords[0] + 1); ++t)
    for (int z = loc_dims[1] * mpi_coords[1]; z < loc_dims[1] * (mpi_coords[1] + 1); ++z)
      for (int y = loc_dims[2] * mpi_coords[2]; y < loc_dims[2] * (mpi_coords[2] + 1); ++y)
        for (int x = loc_dims[3] * mpi_coords[3]; x < loc_dims[3] * (mpi_coords[3] + 1); ++x)
        {
          failure += (data[ctr][0] != t);
          failure += (data[ctr][1] != z);
          failure += (data[ctr][2] != y);
          failure += (data[ctr][3] != x);
          ++ctr;
        }
  
  if (failure)
    fprintf(stderr, "Node %d reports: %d failures.\n", rank, failure);
  else
    fprintf(stderr, "Node %d reports data okay.\n", rank);

  lemonDestroyReader(r);
  MPI_File_close(&fp);
  MPI_Finalize();

  free(data);

  return(0);
}


