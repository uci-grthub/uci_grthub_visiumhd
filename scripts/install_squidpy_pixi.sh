#!/bin/bash
#SBATCH --job-name=squidpy_pixi
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/squidpy_pixi_%j.out
#SBATCH --error=logs/squidpy_pixi_%j.err

# Install the locked squidpy pixi environment (envs/squidpy/pixi.lock), check that the
# packages the rules import load, and warm the CellChatDB cache for squidpy_ccc.

set -euo pipefail

MANIFEST=envs/squidpy/pixi.toml
export PYTHONNOUSERSITE=1

pixi install --locked --manifest-path "$MANIFEST"

pixi run --manifest-path "$MANIFEST" python - <<'EOF'
import importlib.metadata as md
import sys
for name in ("numpy", "pandas", "anndata", "scanpy", "squidpy", "liana", "omnipath",
             "harmonypy", "leidenalg", "igraph", "pyarrow", "shapely", "tifffile"):
    print(f"{name:10s} {md.version(name)}")
import harmonypy, leidenalg, igraph, liana, squidpy, shapely, tifffile  # noqa: F401
import numpy
# fails if ~/.local site-packages leaked in
assert numpy.__file__.startswith(sys.prefix), numpy.__file__
print("imports ok, numpy from", numpy.__file__)
EOF

# compute nodes may have no outbound network later
pixi run --manifest-path "$MANIFEST" python -c \
    "import liana as li; print('cellchatdb', li.resource.select_resource('cellchatdb').shape)"
