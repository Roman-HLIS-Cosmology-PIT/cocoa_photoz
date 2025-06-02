#!/bin/bash -l
#SBATCH --job-name=JACOBIAN-YAML
#SBATCH --output=/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/outs_jacobian/%x_%A_%a.out
#SBATCH --error=/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/outs_jacobian/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=souzadio.jpl.nasa.gov

yaml_file=/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/NZ_EVALUATE1.yaml

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Slurm submit DIR is $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR
module purge > /dev/null 2>&1

source /home/souzadio/miniconda3/etc/profile.d/conda.sh
module load slurm
conda activate cocoapy310_env
source start_cocoa.sh

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

$CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run $yaml_file -f