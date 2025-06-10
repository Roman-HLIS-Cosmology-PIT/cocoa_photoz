#!/bin/bash -l
#SBATCH --job-name=JACOBIAN-LSST
#SBATCH --output=./cocoa_photoz/outs_jacobian/%x_%A_%a.out
#SBATCH --error=./cocoa_photoz/outs_jacobian/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00

yaml_file=./cocoa_photoz/yamls/lsst_y1/NZ_EVALUATE${SLURM_ARRAY_TASK_ID}.yaml

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

$CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run $yaml_file -f