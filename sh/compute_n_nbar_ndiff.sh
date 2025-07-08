#!/bin/bash -l
#SBATCH --job-name=n_sc1bd6
#SBATCH --output=results/outs_jacobian/%x_%A_%a.out
#SBATCH --error=results/outs_jacobian/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00
##SBATCH --array=0-10 

# STRINGS=(
# 'sc1b_d6' 'sc1b_d7'
# 'sc2b_d4' 'sc2b_d6'
# 'sc3b_d4' 'sc3b_d5' 'sc3b_d7')

# scenarios = [
# 'sc1bd6','sc1bd7',
# 'sc2bd4','sc2bd6',
# 'sc3bd4','sc3bd5','sc3bd7']

# Pick the string based on SLURM_ARRAY_TASK_ID
# S=${STRINGS[$SLURM_ARRAY_TASK_ID]}

py_file=/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/cocoa_photoz/codes/compute_n_nbar_ndiff_sc1bd6.py

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Slurm submit DIR is $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR
module purge > /dev/null 2>&1

# source /home/souzadio/miniconda3/etc/profile.d/conda.sh
# conda activate cocoapy310_env
# source start_cocoa.sh

module load slurm

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

# $CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} python $py_file --sc "$S"
$CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} python $py_file