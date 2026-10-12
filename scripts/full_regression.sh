#!/usr/bin/env bash
set -euo pipefail

echo "[1/4] Reset output directory"
rm -rf "test/output"
mkdir -p "test/output"

echo "[2/4] Build docker image"
docker build --platform linux/amd64 -t coinstacteam/vbm_pre .

echo "[3/4] Run simulator"
../coinstac/packages/coinstac-simulator/bin/coinstac-simulator -c 1 --loglevel silly

echo "[4/4] Compare outputs"
python3 compare_outputs.py \
  "test/output.zip" \
  "test/output/local0/simulatorRun/vbm_outputs.zip" \
  --strip-prefix "output/local0/simulatorRun/vbm_outputs" \
  --ignore "output/local0/simulatorRun/vbm_outputs.zip" \
  --ignore "output/local0/simulatorRun/wc1Re.png" \
  --tolerate "*.mat"

echo "Full regression completed successfully."
