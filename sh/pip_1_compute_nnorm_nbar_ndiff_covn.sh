#!/bin/bash -l
#SBATCH --job-name=pip_1_sc1bd4
#SBATCH --output=./cocoa_photoz/results/outs/%x_%A_%a.out
#SBATCH --error=./cocoa_photoz/results/outs/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=7-00:00:00

roman_scenario=""
py_file=/gpfs/scratch/pit-roman-hlis/Diogo/cocoa/Cocoa/cocoa_photoz/codes/pip_1_compute_nnorm_nbar_ndiff_covn.py

# sc1bd4 - done  (jpl)
# sc1bd5 - done  (edge) / tranferred to (jpl)  
# sc1bd6 - done  (jpl)
# sc1bd7 - done  (jpl)
# sc2bd4 - done  (edge) / tranferred to (jpl)  
# sc2bd5 - done  (jpl)
# sc2bd6 - done  (edge) / tranferred to (jpl)    
# sc2bd7 - done  (jpl)   
# sc3bd4 - done  (edge) / tranferred to (jpl)      
# sc3bd5 - done  (jpl)  
# sc3bd7 - done  (edge) / tranferred to (jpl)       
# sc1bd4_dz001 - mismatch of information  (waiting Boyan's confirmation) 

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Slurm submit DIR is $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR
module purge > /dev/null 2>&1

conda activate cocoa
source start_cocoa.sh

module load slurm

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi

$CONDA_PREFIX/bin/mpirun -n ${SLURM_NTASKS} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} python $py_file --scenario "$roman_scenario"