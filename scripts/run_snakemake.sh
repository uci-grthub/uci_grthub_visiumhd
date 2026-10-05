#!/bin/bash
#SBATCH --job-name=smk_controller
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/smk_controller_%j.out
#SBATCH --error=logs/smk_controller_%j.out

# Run the Snakemake controller as its own job, so it outlives the terminal that
# submitted it. The controller only submits and polls; each rule is its own sbatch job
# through the slurm profile. Arguments are passed straight to snakemake:
#   sbatch scripts/run_snakemake.sh output/squidpy/square_016um/clustered.h5ad
#   sbatch scripts/run_snakemake.sh squidpy_all --until squidpy_combine

set -euo pipefail

# conda stays on PATH for rules with a conda: directive; snakemake itself is the
# pixi global install in $HOME, which starts in seconds instead of minutes on /dfs9
source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh

"$HOME/.pixi/bin/snakemake" --profile slurm --rerun-incomplete "$@"
