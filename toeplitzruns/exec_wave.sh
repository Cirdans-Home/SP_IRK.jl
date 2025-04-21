#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=@THREAD@
#SBATCH --job-name=RK-@THREAD@-@STAGES@
#SBATCH -o /scratch/durastante/rungekutta_tests/SP_IRK.jl/toeplitzruns/logfiles/outerr/RK-@THREAD@-@STAGES@.txt
#SBATCH -e /scratch/durastante/rungekutta_tests/SP_IRK.jl/toeplitzruns/logfiles/outerr/RK-@THREAD@-@STAGES@-err.txt
#SBATCH --partition=gpu
#SBATCH --time=2-00:00:00 

module purge
module load julia/gpu-compiled
module list

method=@METHOD@
stages=@STAGES@

# export OPENBLAS_NUM_THREADS=128
export OMP_NUM_THREADS=@THREAD@
cd /scratch/durastante/rungekutta_tests/SP_IRK.jl 

srun julia --project --threads=@THREAD@  experiments/nonlinearwave.jl ${method} ${stages} >> toeplitzruns/logfiles/rk_${method}_t@THREAD@_s${stages}.txt 2>&1
