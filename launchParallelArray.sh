#!/bin/bash -l

#SBATCH -J R0.1
#SBATCH -D /home/mani/abalaji/3DCatalystModeling/
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -o %x_%A_%a.out
#SBATCH --mail-type=NONE
#SBATCH --mail-user=abalaji@stanford.edu
#SBATCH --array 1-13

ml apps/matlab

##########################
### Beginning of Execution
##########################

echo
echo The master node of this job is `hostname`
echo The working directory is `pwd`
echo The slurm task ID is $SLURM_ARRAY_TASK_ID
echo
echo -------------------------
echo Execution Starting Time:
echo `date`
echo -------------------------
echo

FOLDER1="/fastscratch/abalaji/newBase/rxn0.1/input.json"
PARAMETERNAME="targetCurrent"
PARAMETERVALUES="[-0.5E3 -0.6E3 -0.75E3 -0.8E3 -1.0E3 -1.25E3 -1.5E3 -1.75E3 -2.0E3 -2.25E3 -2.5E3 -2.75E3 -3E3]"
matlab -batch "runParallelJob('${FOLDER1}', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, '${PARAMETERNAME}', '${PARAMETERVALUES}')"

echo
echo -------------------------
echo Execution Ending Time:
echo `date`
echo -------------------------
echo

