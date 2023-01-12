NCPUS=8
MEM=11800mb
TIME="71:58:02"

#!/bin/sh

#PBS -l walltime=${TIME}
#PBS -l select=1:ncpus=${NCPUS}:mem=${MEM}

$HOME/julia-1.6.4/bin/julia $HOME/Tethys.jl/src/hpc.jl


mkdir $WORK/$PBS_JOBID
cp * $WORK/$PBS_JOBID



