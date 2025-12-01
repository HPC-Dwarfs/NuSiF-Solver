/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include "profiler.h"
#include "comm.h"

double T[NUMREGIONS];
size_t C[NUMREGIONS];
static FILE *Filehandle;
static double CommVolume = 0.0;
static double Updates    = 0.0;
static double Dataset    = 0.0;

void initProfiler(CommType *c)
{
#ifdef PROFILING
  int tmpSize       = 0;
  double totalSize  = 0.0;

  int imaxLocal     = c->imaxLocal;
  int jmaxLocal     = c->jmaxLocal;
  int kmaxLocal     = c->kmaxLocal;
  double datapoints = (kmaxLocal + 2) * (jmaxLocal + 2) * (imaxLocal + 2);
  Updates           = kmaxLocal * jmaxLocal * imaxLocal;
  Dataset           = datapoints * 16.0 * 1.E-6;

  for (size_t i = 0; i < NUMREGIONS; i++) {
    T[i] = 0.0;
    C[i] = 0;
  }

#ifdef _MPI
  for (int i = 0; i < NDIRS; i++) {
    MPI_Type_size(c->sbufferTypes[i], &tmpSize);
    totalSize += tmpSize;
    MPI_Type_size(c->rbufferTypes[i], &tmpSize);
    totalSize += tmpSize;
  }

  CommVolume = totalSize * 1.E-6;
  char filename[MAX_STR_LENGTH];
  sprintf(filename, "profile-%d.txt", c->rank);
  Filehandle = fopen(filename, "we");
#endif /* ifdef _MPI */
#endif
}

void printProfile(CommType *c, int it)
{
#ifdef PROFILING
  double tmin[NUMREGIONS];
  double tmax[NUMREGIONS];
  double tavg[NUMREGIONS];

  commReduce(T, tmax, NUMREGIONS, MAX);
  commReduce(T, tmin, NUMREGIONS, MIN);
  commReduce(T, tavg, NUMREGIONS, SUM);

  if (commIsMaster(c)) {
    double mlups = (it * Updates) * 1.E-6;

    for (int i = 0; i < NUMREGIONS; i++) {
      tavg[i] /= c->size;
    }

    double walltime = tavg[SOLVER];

    printf("Solver Time: %.2f s Dataset size: %.2f MB Performance: %.2f MLUPS/s "
           "Performance total: %.2f MLUPS/s\n",
        walltime,
        Dataset,
        mlups / walltime,
        (mlups * c->size) / walltime);

    // printf(
    //     "min %11.2f  max %11.2f avg %11.2f\n", tmin[SOLVER], tmax[SOLVER], tavg[SOLVER]);
  }

#ifdef _MPI
  double dataset  = CommVolume * C[COMM];
  double walltime = T[COMM];
  fprintf(Filehandle,
      "Calls %lu Communication Time: %.2f s Data volume: %.2f MB Bandwidth: %.2f "
      "MB/s\n",
      C[COMM],
      walltime,
      dataset,
      dataset / walltime);
  fprintf(Filehandle,
      "Communication Time/Call: %.6f s Data volume/Call: %.3f"
      "MB\n",
      walltime / C[COMM],
      dataset / C[COMM]);
#endif /* ifdef _MPI */

  C[SOLVER] = 0;
  C[COMM]   = 0;
  T[SOLVER] = 0;
  T[COMM]   = 0;
#endif
}

void finalizeProfiler()
{
#ifdef PROFILING
  fclose(Filehandle);
#endif
}
