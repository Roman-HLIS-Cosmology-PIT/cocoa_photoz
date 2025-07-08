#!/bin/bash -l
#SBATCH --job-name=PCA_same_gauss
#SBATCH --output=cocoa_photoz/results/outs_jacobian/%x_%A_%a.out
#SBATCH --error=cocoa_photoz/results/outs_jacobian/%x_%A_%a.err
###############SBATCH --ntasks-per-node=4
###############SBATCH --cpus-per-task=7
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00

# yaml_file=cocoa_photoz/yamls/roman_sc1bd4_g/SHIFT_MODEL_MCMC${SLURM_ARRAY_TASK_ID}.yaml
# yaml_file=cocoa_photoz/yamls/roman_sc1bd4_test_dv/PCA_MODEL_printdv_3x2_MCMC${SLURM_ARRAY_TASK_ID}.yaml
# yaml_file=cocoa_photoz/yamls/roman_sc1bd4_g/ROMAN_REAL_MCMC0_print_datavector_file_example1_3x2pt.modelvector.yaml
yaml_file=cocoa_photoz/yamls/roman_sc1bd4_g_18decimal_same_params_than_cosmocov_gauss_alphas_prior/PCA_MODEL_MCMC${SLURM_ARRAY_TASK_ID}.yaml
# yaml_file=cocoa_photoz/yamls/roman_real/ROMAN_REAL_MCMC${SLURM_ARRAY_TASK_ID}.yaml

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
# $CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} python $yaml_file -f