#!/bin/bash -l

#SBATCH --job-name=JACOBIAN-COCOA
#SBATCH --output=/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PZ/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/outs_jacobian/%x_%A_%a.out
#SBATCH --error=/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PZ/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/outs_jacobian/%x_%A_%a.err
#SBATCH --partition=extended-96core
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=7
#SBATCH --time=7-00:00:00

###fisher_file=/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PZ/cocoa_photoz/Cocoa/projects/cocoa_des_y3_nzpca/04_fisher.py
###fisher_file=/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PZ/cocoa_photoz/Cocoa/projects/cocoa_des_y3_nzpca/05_fisher_pile3.py
fisher_file=/gpfs/projects/MirandaGroup/Diogo/ROMAN-NZ-PROJECT/PZ/cocoapy310/Cocoa/projects/lsst_y1/cocoa_photoz/jacobian.py

cd $SLURM_SUBMIT_DIR
module purge > /dev/null 2>&1
source /gpfs/projects/MirandaGroup/Diogo/miniconda3/etc/profile.d/conda.sh
module load slurmb
module load texlive
conda activate cocoapy310_env
source start_cocoa
export SLURM_TRES_PER_TASK=cpu=1
export OMP_NUM_THREADS=1
export I_MPI_SPIN_COUNT=1
export I_MPI_SHM_EAGER_THRESHOLD=4096
export SLURM_CPU_BIND=none
export SLURM_WHOLE=1

# mpirun -n 40 python $fisher_file
$CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --report-bindings --mca vader,btl tcp,self --bind-to core --map-by numa:pe=${OMP_NUM_THREADS} python $fisher_file