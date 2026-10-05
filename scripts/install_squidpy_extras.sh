#!/bin/bash
#SBATCH --job-name=squidpy_extras
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/squidpy_extras_%j.out
#SBATCH --error=logs/squidpy_extras_%j.err

# Add the packages squidpy_cluster and squidpy_ccc need to the squidpy_dfs env:
# harmonypy, leidenalg and python-igraph from conda-forge, liana from PyPI.
# The core scverse stack is pinned to what is installed so the solver cannot move it.
# A snapshot of the env is written first; roll back with
#   conda install -n squidpy_dfs --revision <n>   (see `conda list -n squidpy_dfs --revisions`)

set -euo pipefail

source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate squidpy_dfs
export PYTHONNOUSERSITE=1

mkdir -p logs
conda list --explicit > "logs/squidpy_dfs_before_extras_${SLURM_JOB_ID:-local}.txt"

mamba install -y -n squidpy_dfs -c conda-forge \
    harmonypy leidenalg python-igraph \
    "numpy=2.3.5" "pandas=2.3.3" "scipy=1.16.3" \
    "anndata=0.12.0" "scanpy=1.11.5" "squidpy=1.8.1"

# the omnipath upgrade fixes the np.NAN crash under NumPy 2
pip install -U omnipath liana

python - <<'EOF'
import importlib.metadata as md
for name in ("numpy", "pandas", "scanpy", "squidpy", "anndata",
             "harmonypy", "leidenalg", "igraph", "liana", "omnipath"):
    print(f"{name:10s} {md.version(name)}")
import harmonypy, leidenalg, igraph, liana  # noqa: F401
EOF

# warm the CellChatDB cache; compute nodes may have no outbound network later
python -c "import liana as li; print('cellchatdb', li.resource.select_resource('cellchatdb').shape)"
