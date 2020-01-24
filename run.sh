#!/bin/sh
#PJM -L "rscgrp=fx-middle"
#PJM -L "node=48"
#PJM --mpi "proc=48"
#PJM --mpi "rank-map-bynode"
#PJM -L "elapse=05:00:00"
#PJM -j
#PJM -S

prog=./ovlpbench

proc=48
thr=32
T=20
M=131072
mode=3

export OMP_NUM_THREADS=${thr}
log=log.${PJM_JOBID}
/bin/rm -f ${log} ${log}.blocking ${log}.nbc-noassist ${log}.nbc-assist

for c in "0 1536" "0 3072" "1 384" "1 768" "2 192" "2 384"
do 
  f=`echo $c | awk '{print $1}'`    
  n=`echo $c | awk '{print $2}'`    
    
  echo "blocking"
  mpiexec -ofout ${log} -n ${proc} ${prog} -n ${n} -m ${M} -t ${T} -f ${f} -w 1
  cat ${log} >> ${log}.blocking
  /bin/rm -f ${log}

  echo "nbc noassist"
  mpiexec -ofout ${log} -n ${proc} ${prog} -n ${n} -m ${M} -t ${T} -f ${f} -w 1 -q 1 
  cat ${log} >> ${log}.nbc-noassist
  /bin/rm -f ${log}

  echo "nbc assist"
  mpiexec -ofout ${log} -n ${proc} --mca opal_progress_thread_mode ${mode} ${prog} -n ${n} -m ${M} -t ${T} -f ${f} -w 1 -q 1 
  cat ${log} >> ${log}.nbc-assist
  /bin/rm -f ${log}
done


