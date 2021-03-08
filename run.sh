#!/bin/bash

#SBATCH --mail-type=end
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --tasks-per-node=1
#SBATCH --time=18:00:00
#SBATCH --job-name=v2_2_DAD_AP

#Your model directory here
BASE_DIR=/scratch/users/qjin4@jhu.edu/v2.2_DAD
BUILD_DIR=$BASE_DIR
EXEC=$BASE_DIR/stoch3d

#Place params.txt in WORK_DIR before running2; output will also be saved to WORK_DIR
WORK_DIR=$BASE_DIR/output_AP
#If you want to save 3D VTK files, create a "stochout" subdirectory in WORK_DIR
#mkdir $WORK_DIR/stochout
cd $WORK_DIR

# load the modules
# . /etc/profile.modules
# module load openmpi/openmpi-1.6.5-gnu
# echo Modules loaded

#Do not recompile if using Intel because it is slow (see Notes.txt)
#cd $BUILD_DIR
#rm main.o
#make clean
# make

#echo pbs nodefile:
#cat $PBS_NODEFILE
#NPROCS=`wc -l < $PBS_NODEFILE`
#cat $PBS_NODEFILE > host.list

#Try this line if you get errors about limited virtual memory
#ulimit -l unlimited

echo Running job...
#export OMP_NUM_THREADS=$NPROCS
export OMP_NUM_THREADS=24
$EXEC
#mpirun --mca btl openib,self -np 1 -hostfile host.list $EXEC
# mpirun -np 1 -hostfile host.list $EXEC

#Try ^openib flag if you get errors about openfabrics:
#mpirun --mca btl ^openib -x OMP_NUM_THREADS=4 -np 1 -hostfile host.list $EXEC
