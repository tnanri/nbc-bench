#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef FUJITSU
#include <mpi-ext.h>
#endif
#include "omp.h"

#ifndef N 
#define N 500
#endif

#ifndef M
#define M 100000
#endif

#ifndef T
#define T 100
#endif

void start_alltoall(double *, double *, int, MPI_Request *);
void start_allreduce(double *, double *, int, MPI_Request *);
void start_bcast(double *, double *, int, MPI_Request *);
void finish_alltoall(MPI_Request *);
void finish_allreduce(MPI_Request *);
void finish_bcast(MPI_Request *);
void do_alltoall(double *, double *, int);
void do_allreduce(double *, double *, int);
void do_bcast(double *, double *, int);

void (*start_funcs[])(double *, double *, int, MPI_Request *) = {
    start_alltoall, 
    start_allreduce,
    start_bcast
};

void (*finish_funcs[])(MPI_Request *) = {
    finish_alltoall, 
    finish_allreduce,
    finish_bcast
};

void (*do_funcs[])(double *, double *, int) = {
    do_alltoall, 
    do_allreduce,
    do_bcast
};


void start_alltoall(double *sbuf, double *rbuf, int size, MPI_Request *reqs)
{
    MPI_Ialltoall(sbuf, size, MPI_DOUBLE, rbuf, size, MPI_DOUBLE, MPI_COMM_WORLD, &(reqs[0]));
}

void start_allreduce(double *sbuf, double *rbuf, int size, MPI_Request *reqs)
{
    MPI_Iallreduce(sbuf, rbuf, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &(reqs[0]));
}

void start_bcast(double *sbuf, double *rbuf, int size, MPI_Request *reqs)
{
    MPI_Ibcast(sbuf, size, MPI_DOUBLE, 0, MPI_COMM_WORLD, &(reqs[0]));
}

void finish_alltoall(MPI_Request *reqs)
{
    MPI_Status st;
    MPI_Wait(&(reqs[0]), &st);
}

void finish_allreduce(MPI_Request *reqs)
{
    MPI_Status st;
    MPI_Wait(&(reqs[0]), &st);
}

void finish_bcast(MPI_Request *reqs)
{
    MPI_Status st;
    MPI_Wait(&(reqs[0]), &st);
}

void do_alltoall(double *sbuf, double *rbuf, int size)
{
    MPI_Alltoall(sbuf, size, MPI_DOUBLE, rbuf, size, MPI_DOUBLE, MPI_COMM_WORLD);
}

void do_allreduce(double *sbuf, double *rbuf, int size)
{
    MPI_Allreduce(sbuf, rbuf, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void do_bcast(double *sbuf, double *rbuf, int size)
{
    MPI_Bcast(sbuf, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void do_comp_static(double *a, double *b, double *c, int n, int nn)
{
  int i, j, k;

#pragma omp parallel for private(j, k) schedule (static)
    for (i = 0; i < n; i++)
          for (k = 0; k < n; k++)
              for (j = 0; j < nn; j++)
                  c[i*nn+j] += a[i*n+k] * b[k*nn+j];

}

void do_comp_dynamic(double *a, double *b, double *c, int n, int nn)
{
    int i, j, k;

#pragma omp parallel for private(j, k) schedule (dynamic)
    for (i = 0; i < n; i++)
          for (k = 0; k < n; k++)
              for (j = 0; j < nn; j++)
                  c[i*nn+j] += a[i*n+k] * b[k*nn+j];

}

int main(int argc, char *argv[])
{
  int myrank, procs, thr;
  MPI_Request reqs[2];
  int i, j, r, n, m, t, f, p, q, d, w, nn;
  double *sbuf, *rbuf;
  double *a, *b, *c;
  double st, et, tcomm, tcomp, tall, rovlp;
  // char nm[128];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  // MPI_Get_processor_name(nm, &i);

#pragma omp parallel
  {
    thr = omp_get_num_threads();
  }

  // fprintf(stderr, "%d / %d %s\n", myrank, procs, nm);

  // setup parameters
  n = N; // size of array
  w = 1; // number of procs per node
  m = M; // size of message
  t = T; // repeat time
  f = 0; // communitation type (0: alltoall, 1: allreduce, 2: bcast)
  p = 0; // use progress pragma (0: no, 1: yes)
  q = 0; // use non-blocking (0: no, 1: yes)
  d = 0; // use OpenMP dynamic scheduling (0: no, 1: yes)

  i = 1;
  while (i < argc){
    if (strcmp(argv[i], "-n")==0){
      n = atoi(argv[i+1]);
      i+=2;
    }else if (strcmp(argv[i], "-w")==0){
      w = atoi(argv[i+1]);
      i+=2;
    }else if (strcmp(argv[i], "-m")==0){
      m = atoi(argv[i+1]);
      i+=2;
    }else if (strcmp(argv[i], "-t")==0){
      t = atoi(argv[i+1]);
      i+=2;
    }else if (strcmp(argv[i], "-f")==0){
      f = atoi(argv[i+1]);
      i+=2;
    }else if (strcmp(argv[i], "-p")==0){
      p = atoi(argv[i+1]);
      i+=2;
    }else if (strcmp(argv[i], "-q")==0){
      q = atoi(argv[i+1]);
      i+=2;
    }else if (strcmp(argv[i], "-d")==0){
      d = atoi(argv[i+1]);
      i+=2;
    }
  }

  //  fprintf(stderr, "n %d w %d m %d t %d f %d p %d q %d d %d  \n",
  //	 n, w, m, t, f, p, q, d);

  // nn: number of columns to be calculated by this process
  nn = n / w;
  if (nn * w < n) {
    fprintf(stderr, "Error: n is not dividable by the number of processes per node\n");
  }

  // prepare arrays
  sbuf = (double *)malloc(m*procs*sizeof(double));
  rbuf = (double *)malloc(m*procs*sizeof(double));
  a = (double *)malloc(n*n*sizeof(double));
  b = (double *)malloc(n*nn*sizeof(double));
  c = (double *)malloc(n*nn*sizeof(double));

  for (i = 0; i < m; i++)
      sbuf[i] = myrank * m + i;
  for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
          a[i * n + j] = (myrank - procs / 2.0) * 0.0001 + j * 0.0001 + i * n * 0.0001;
  for (i = 0; i < n; i++)
      for (j = 0; j < nn; j++)
          b[i * nn + j] = (myrank - procs / 2.0) * 0.0001;
  for (i = 0; i < n; i++)
      for (j = 0; j < nn; j++)
          c[i * nn + j] = 0.0;

  // warm-up
  MPI_Barrier(MPI_COMM_WORLD);

  for (i = 0; i < 10; i++) {
      if (q == 0) {
          do_funcs[f](sbuf, rbuf, m);
      } else {
          start_funcs[f](sbuf, rbuf, m, reqs);
          finish_funcs[f](reqs);
      }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  // communication only
  st = MPI_Wtime();
  for (r = 0; r < t; r++) {
      if (q == 0) {
          do_funcs[f](sbuf, rbuf, m);
      } else {
          start_funcs[f](sbuf, rbuf, m, reqs);
          finish_funcs[f](reqs);
      }
  }
  et = MPI_Wtime();
  tcomm = (et - st) / t;

  // computation only
  MPI_Barrier(MPI_COMM_WORLD);
  st = MPI_Wtime();
  for (r = 0; r < t; r++) {
      if (d == 0) 
	do_comp_static(a, b, c, n, nn);
      else
	do_comp_dynamic(a, b, c, n, nn);
  }
  et = MPI_Wtime();
  tcomp = (et - st) / t;

  // do both
  MPI_Barrier(MPI_COMM_WORLD);
  st = MPI_Wtime();
  for (r = 0; r < t; r++){
      if (q == 0) {
          do_funcs[f](sbuf, rbuf, m);
          if (d == 0) 
              do_comp_static(a, b, c, n, nn);
          else
              do_comp_dynamic(a, b, c, n, nn);
      } else {
          start_funcs[f](sbuf, rbuf, m, reqs);
#ifdef FUJITSU
          if (p == 1)
              FJMPI_Progress_start();
#endif
          if (d == 0) 
              do_comp_static(a, b, c, n, nn);
          else
              do_comp_dynamic(a, b, c, n, nn);
#ifdef FUJITSU
          if (p == 1)
              FJMPI_Progress_stop();
#endif
          finish_funcs[f](reqs);
      }
  }
  et = MPI_Wtime();
  tall = (et - st) / t;

  rovlp = ((tcomm - (tall - tcomp)) / tcomm) * 100.0;

  if (myrank==0)
      printf(" procs %d thr %d m %d n %d w %d nn %d d %d p %d q %d f %d comm %e comp %e all %e overlap %e %%\n",
             procs, thr, m, n, w, nn, d, p, q, f, tcomm, tcomp, tall, rovlp);

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
}

