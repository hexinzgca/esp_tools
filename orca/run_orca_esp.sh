#!/usr/bin/env bash
# Run an ORCA input in the container, then extract ESP on a 10 x natoms grid.
#
# Usage: ./run_orca_esp.sh wb97m-v_def2-tzvpd.inp [water.xyz]
# The input must use `* xyzfile ... water.xyz` (or the geometry passed as $2).
set -euo pipefail

INP=${1:-wb97m-v_def2-tzvpd.inp}
XYZ=${2:-water.xyz}
SIF=${ORCA_SIF:-/mnt/e/apptainer_sifs/build/orca-6.1.sif}
BASE=${INP%.inp}

if [[ ! -f "$INP" ]]; then
  echo "missing ORCA input: $INP" >&2
  exit 1
fi
if [[ ! -f "$XYZ" ]]; then
  echo "missing XYZ geometry: $XYZ" >&2
  exit 1
fi

# ORCA requires a full pathname when %pal nprocs > 1.  Oversubscription keeps
# small workstations from failing when the input requests more ranks than
# physical cores; lower %pal in the input if you want fewer processes.
apptainer exec --bind "$PWD:/work" "$SIF" \
  bash -c 'set -e; cd /work; export OMP_NUM_THREADS=1; export OMPI_MCA_rmaps_base_oversubscribe=1; ORCA=$(command -v orca); exec "$ORCA" "$1"' _ "$INP" \
  > "${BASE}.out" 2>&1

if [[ ! -f "${BASE}.gbw" ]]; then
  echo "ORCA did not produce ${BASE}.gbw; see ${BASE}.out" >&2
  exit 1
fi

PYTHON=${PYTHON:-python3}
"$PYTHON" "$(dirname "$0")/extract_esp.py" \
  --xyz "$XYZ" \
  --gbw "${BASE}.gbw" \
  --orca-container "$SIF" \
  --output "${BASE}_esp.json"

echo "wrote ${BASE}_esp.json"
