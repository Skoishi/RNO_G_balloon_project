#!/bin/bash
#SBATCH --job-name=rnog_raytracing
#SBATCH --account=PAS2608
#SBATCH --time=02:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=outputs/rnog_raytracing_%j.out
#SBATCH --error=outputs/rnog_raytracing_%j.err

# Load modules
module load gcc/12.3.0
module load python/3.12
module load cmake/3.25.2
module load root/6.28.06

# Activate virtual environment
source ~/my_rnog_env/bin/activate

# Create output directory
cd ~/RNO_G/July_balloon_Code/finished_code/theory
mkdir -p Ice_model

# Run Python script
python Ice_model1.4_balloon.py

deactivate