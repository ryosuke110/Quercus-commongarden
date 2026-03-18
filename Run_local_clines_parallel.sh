#!/usr/bin/env bash
# Run local cline estimation in parallel batches
# Author: Ryosuke Ito

set -euo pipefail

### Configuration ###
MAX_JOBS=8          # maximum concurrent jobs
TOTAL_BATCHES=22    # total SNP batches
SCRIPT="Local_cline_fit.R"

running=0

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log "Starting ${TOTAL_BATCHES} batches (max ${MAX_JOBS} concurrent)"

for ((batch=1; batch<=TOTAL_BATCHES; batch++))
do
  log "Launching batch ${batch}"

  Rscript --vanilla "${SCRIPT}" "${batch}" &

  ((running++))

  if (( running >= MAX_JOBS )); then
    wait -n
    ((running--))
  fi
done

wait

log "All batches finished"