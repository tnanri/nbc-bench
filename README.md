# nbc-bench
Benchmark program for measuring effect of communication overlap by non-blocking collectives.

Usage:
 - Modify makefile and type make.
 - Run:
    mpirun -np #procs ./ovlpbench [options]
    Available options:
      -n # : size of array for matrix multiplication
      -w # : number of procs per node
      -m # : size of message
      -t # : repeat time
      -f # : communitation type (0: alltoall, 1: allreduce, 2: bcast)
      -p # : use progress pragma (0: no, 1: yes) (Fujitsu MPI only)
      -q # : use non-blocking (0: no, 1: yes)
      -d # : use OpenMP dynamic scheduling (0: no, 1: yes)

Output: 
    procs : number of procs
    thr : number of threads per proc
    m : size of message
    n : size of array for matrix multiplication
    w : number of procs per node
    nn : number of columns of array per proc in matrix multiplication
    d : thread scheduling (0: static, 1: dynamic)
    p : usage of progress pragma
    q : usage of non-blocking collective
    f : function number (0: alltoall, 1: allreduce, 2: bcast)
    comm : time for communication only
    comp : time for computation only
    all : time with computation and communication
    overlap : ratio of overlap (= (comm - (all - comp)) / comm)

Note: 
 - n (size of array for matrix multiplication) should be chosen carefully
   according to the message size, to see sufficient effects of overlap.
