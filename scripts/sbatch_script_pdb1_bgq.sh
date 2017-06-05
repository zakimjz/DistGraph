#!/bin/bash
#SBATCH --job-name=MPI_GSPAN_BGQ
#SBATCH --partition=small
##SBATCH --ntasks=256

#this job requests x nodes
#SBATCH --nodes=4

# this job requests exclusive access to the nodes it is given
# this mean it will be the only job running on the node
#SBATCH --exclusive

# only request 1 MPI task per node
#SBATCH --ntasks-per-node=1

# and request 16 cpus per task for OpenMP threads
#SBATCH --cpus-per-task=8


#SBATCH -t 24:00:00
#SBATCH -D /gpfs/u/home/HPDM/HPDMtlkd/scratch/single-graph-miner/src/distributed
#SBATCH --mail-type=ALL
#SBATCH --mail-user=talukn@rpi.edu
#SBATCH -o joboutput.%J
#SBATCH -e error.%J.log

export OMP_WAIT_POLICY=ACTIVE
export BG_SMP_FAST_WAKEUP=YES
export OMP_PROC_BIND=TRUE

. ../../scripts/localvars.sh


export OMP_NUM_THREADS=8 #$SLURM_CPUS_PER_TASK
echo "num threads = $OMP_NUM_THREADS"

SUPPORT_LIST="20 19 18 17 16 15 14 13 12 11 10"
#SUPPORT_LIST="22"
FILETYPE_SEQ="-txt"
FILETYPE_PAR="-txt"
FILETYPE_PAR="-ladjp"

SUPPORT_LIST="200 190 180 170 160 150 140 130 120 110 100"
SUPPORT_LIST="200 180"
INPUT_GRAPH_SEQ="pdb1.stxt"
INPUT_GRAPH_PAR="pdb1_part_4"
NUM_PARTITIONS=4

RUNID=1
NUMPROC=$SLURM_NPROCS
NUMTHRDS=$OMP_NUM_THREADS #$SLURM_CPUS_PER_TASK
for MPI_EXECFILE in distgraph
do

   for support in ${SUPPORT_LIST}
   do
      PAR_LOGFILE="` get_par_log_filename $support ${INPUT_GRAPH_PAR} ${NUMPROC} ${NUMTHRDS} ${RUNID} `-${MPI_EXECFILE}"
      SEQ_LOGFILE="` get_seq_log_filename $support ${INPUT_GRAPH_SEQ} `"
       
      echo "${PAR_LOGFILE} ${SEQ_LOGFILE}"
      #srun --ntasks 1 ./gspan_seq ${FILETYPE_SEQ} ../../testdata/${INPUT_GRAPH_SEQ} ${support} > ${SEQ_LOGFILE}
      srun   --runjob-opts="--mapping TEDCBA --envs BG_COREDUMPDISABLED=0 BG_COREDUMPONEXIT=1 " ${MPI_EXECFILE} ${FILETYPE_PAR} ../../testdata/${INPUT_GRAPH_PAR} ${support} ${NUM_PARTITIONS} >  ${PAR_LOGFILE}

   done

done
