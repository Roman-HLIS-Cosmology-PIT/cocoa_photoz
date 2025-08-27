#!/bin/bash
#SBATCH -J FISHER
#SBATCH -o jacobian_outs/stencil_5pt_log/outs/%x_%A_%a.out
#SBATCH -e jacobian_outs/stencil_5pt_log/outs/%x_%A_%a.err
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH -t 10-00:00:00
##SBATCH --mem=0

path=cocoa_photoz/yamls/roman_pca/jacobian_RR
yaml=$path/NZ_EVALUATE${SLURM_ARRAY_TASK_ID}.yaml

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Slurm submit DIR is $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR
module purge > /dev/null 2>&1
module load slurm

source /cm/shared/apps/miniforge/etc/profile.d/conda.sh
conda activate cocoa
source start_cocoa.sh

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

$CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run $yaml -f