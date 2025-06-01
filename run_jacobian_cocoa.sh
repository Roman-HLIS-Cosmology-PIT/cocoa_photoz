#!/bin/bash -l

#SBATCH --job-name=JACOBIAN-COCOA
#SBATCH --output=/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/outs_jacobian/%x_%A_%a.out
#SBATCH --error=/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/outs_jacobian/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00

###fisher_file=/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PZ/cocoa_photoz/Cocoa/projects/cocoa_des_y3_nzpca/04_fisher.py
###fisher_file=/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PZ/cocoa_photoz/Cocoa/projects/cocoa_des_y3_nzpca/05_fisher_pile3.py
fisher_file=/gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/jacobian.py

cd $SLURM_SUBMIT_DIR

module purge > /dev/null 2>&1

source /home/souzadio/miniconda3/etc/profile.d/conda.sh
module load slurm
conda activate /home/souzadio/miniconda3/envs/cocoapy310_env/
source /gpfs/scratch/pit-roman-hlis/Diogo/cocoapy310/Cocoa/start_cocoa.sh

export SLURM_TRES_PER_TASK=cpu=1
export OMP_NUM_THREADS=1
export I_MPI_SPIN_COUNT=1
export I_MPI_SHM_EAGER_THRESHOLD=4096
export SLURM_CPU_BIND=none
export SLURM_WHOLE=1

$CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} python $fisher_file
# $CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --report-bindings --mca vader,btl tcp,self --bind-to core --map-by numa:pe=${OMP_NUM_THREADS} python $fisher_file
# $CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} python $fisher_file