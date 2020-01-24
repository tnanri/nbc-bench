MPI = mpifccpx -Kfast
DEFS = -DN=200 -DM=1000000 -DT=10 -DFUJITSU

all: ovlpbench

ovlpbench: ovlpbench.c
	$(MPI) -Kopenmp $(DEFS)  ovlpbench.c -o ovlpbench

clean: 
	/bin/rm -f ovlpbench
