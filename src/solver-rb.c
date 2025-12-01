/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "comm.h"
#include "mpi_proto.h"
#include "parameter.h"
#include "solver.h"
#include "timing.h"
#include "util.h"

#define MAX_STR_LENGTH 30
typedef enum { SOLVER = 0, COMM, NUMREGIONS } RegionsType;

#ifdef VERBOSE
#define PROFILE(tag, call)                                                               \
  const double ts = getTimeStamp();                                                      \
  call;                                                                                  \
  T[tag] += (getTimeStamp() - ts);                                                       \
  C[tag]++;

static double CommVolume = 0.0;
static double T[NUMREGIONS];
static size_t C[NUMREGIONS];
static FILE *Filehandle;
#else
#define PROFILE(call) call;
#endif

static void initProfile(CommType *c)
{
  int tmpSize      = 0;
  double totalSize = 0.0;

  for (size_t i = 0; i < NUMREGIONS; i++) {
    T[i] = 0.0;
    C[i] = 0;
  }

  for (int i = 0; i < NDIRS; i++) {
    MPI_Type_size(c->sbufferTypes[i], &tmpSize);
    totalSize += tmpSize;
    MPI_Type_size(c->rbufferTypes[i], &tmpSize);
    totalSize += tmpSize;
  }

  CommVolume = totalSize * 1.E-6;
  char filename[MAX_STR_LENGTH];
  sprintf(filename, "profile-%d.txt", c->rank);
  Filehandle = fopen(filename, "w");
}

static void printProfile(Solver *s, int it)
{
  CommType *c = s->comm;
  double tmin[NUMREGIONS];
  double tmax[NUMREGIONS];
  double tavg[NUMREGIONS];

  MPI_Reduce(T, tmin, NUMREGIONS, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(T, tmax, NUMREGIONS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(T, tavg, NUMREGIONS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (commIsMaster(c)) {
    int imaxLocal       = s->comm->imaxLocal;
    int jmaxLocal       = s->comm->jmaxLocal;
    int kmaxLocal       = s->comm->kmaxLocal;
    double datapoints   = (kmaxLocal + 2) * (jmaxLocal + 2) * (imaxLocal + 2);
    double updatepoints = kmaxLocal * jmaxLocal * imaxLocal;
    double dataset      = datapoints * 16.0 * 1.E-6;
    double mlups        = (it * updatepoints) * 1.E-6;

    for (int i = 0; i < NUMREGIONS; i++) {
      tavg[i] /= c->size;
    }

    double walltime = tavg[SOLVER];

    printf("Solver Time: %.2f s Dataset size: %.2f MB Performance: %.2f MLUPS/s\n",
        walltime,
        dataset,
        mlups / walltime);

    printf(
        "min %11.2f  max %11.2f avg %11.2f\n", tmin[SOLVER], tmax[SOLVER], tavg[SOLVER]);
  }

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

  C[SOLVER] = 0;
  C[COMM]   = 0;
  T[SOLVER] = 0;
  T[COMM]   = 0;
}

void initSolver(Solver *s, Discretization *d, Parameter *p)
{
  s->eps     = p->eps;
  s->omega   = p->omg;
  s->itermax = p->itermax;
  s->grid    = &d->grid;
  s->comm    = &d->comm;
  s->problem = p->name;

  initProfile(s->comm);
}

void finalizeSolver()
{
  fclose(Filehandle);
}

double solve(Solver *s, double *p, const double *rhs)
{
  int imaxLocal = s->comm->imaxLocal;
  int jmaxLocal = s->comm->jmaxLocal;
  int kmaxLocal = s->comm->kmaxLocal;

  int imax      = s->grid->imax;
  int jmax      = s->grid->jmax;
  int kmax      = s->grid->kmax;

  double eps    = s->eps;
  int itermax   = s->itermax;
  double dx2    = s->grid->dx * s->grid->dx;
  double dy2    = s->grid->dy * s->grid->dy;
  double dz2    = s->grid->dz * s->grid->dz;
  double idx2   = 1.0 / dx2;
  double idy2   = 1.0 / dy2;
  double idz2   = 1.0 / dz2;

  double factor =
      s->omega * 0.5 * (dx2 * dy2 * dz2) / (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);
  double epssq = eps * eps;
  int it       = 0;
  double res   = 1.0;
  int pass, ksw, jsw, isw;

  double timeStart = getTimeStamp();
  while ((res >= epssq) && (it < itermax)) {
    ksw = 1;

    for (pass = 0; pass < 2; pass++) {
      jsw = ksw;
      PROFILE(COMM, commExchange(s->comm, p));

      for (int k = 1; k < kmaxLocal + 1; k++) {
        isw = jsw;
        for (int j = 1; j < jmaxLocal + 1; j++) {
          for (int i = isw; i < imaxLocal + 1; i += 2) {

            double r = RHS(i, j, k) -
                       ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                           (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                           (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2);

            P(i, j, k) -= (factor * r);
            res += (r * r);
          }
          isw = 3 - isw;
        }
        jsw = 3 - jsw;
      }
      ksw = 3 - ksw;
    }
#ifdef _MPI
    if (commIsBoundary(s->comm, FRONT)) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, j, 0) = P(i, j, 1);
        }
      }
    }

    if (commIsBoundary(s->comm, BACK)) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, j, kmaxLocal + 1) = P(i, j, kmaxLocal);
        }
      }
    }

    if (commIsBoundary(s->comm, BOTTOM)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, 0, k) = P(i, 1, k);
        }
      }
    }

    if (commIsBoundary(s->comm, TOP)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, jmaxLocal + 1, k) = P(i, jmaxLocal, k);
        }
      }
    }

    if (commIsBoundary(s->comm, LEFT)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          P(0, j, k) = P(1, j, k);
        }
      }
    }

    if (commIsBoundary(s->comm, RIGHT)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          P(imaxLocal + 1, j, k) = P(imaxLocal, j, k);
        }
      }
    }
#else
    for (int j = 1; j < jmax + 1; j++) {
      for (int i = 1; i < imax + 1; i++) {
        P(i, j, 0)        = P(i, j, 1);
        P(i, j, kmax + 1) = P(i, j, kmax);
      }
    }

    for (int k = 1; k < kmax + 1; k++) {
      for (int i = 1; i < imax + 1; i++) {
        P(i, 0, k)        = P(i, 1, k);
        P(i, jmax + 1, k) = P(i, jmax, k);
      }
    }

    for (int k = 1; k < kmax + 1; k++) {
      for (int j = 1; j < jmax + 1; j++) {
        P(0, j, k)        = P(1, j, k);
        P(imax + 1, j, k) = P(imax, j, k);
      }
    }
#endif
    commReduction(&res, SUM);
    res = res / (double)(imax * jmax * kmax);
#ifdef DEBUG
    if (commIsMaster(&s->comm)) {
      printf("%d Residuum: %e\n", it, res);
    }
#endif

    PROFILE(COMM, commExchange(s->comm, p));
    it++;
  }
  T[SOLVER] = getTimeStamp() - timeStart;

#ifdef VERBOSE
  if (commIsMaster(s->comm)) {
    printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
  }

  printProfile(s, it);
#endif

  return res;
}
