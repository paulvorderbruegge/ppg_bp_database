#!/usr/bin/env bash
# Download all datasets (or only the ones given) to data/raw/<dataset>. Each download is verified with the
# checksums published by the source. Rerunning is safe: complete files are skipped or downloaded again.
#     scripts/download_all.sh                  # everything (~395 GB download, PulseDB alone 388 GB + 636 GB unpacked; see README for space)
#     scripts/download_all.sh ppgbp pttppg     # only these
#     PHYSIONET_USER=<username> scripts/download_all.sh mimic3   # MIMIC-III clinical data (credentialed, asks for the
#                                              # password, so not part of "everything"; see download_mimic3.sh)
# Needs: bash, curl, unzip, md5sum/sha1sum/sha256sum, and 7-Zip for PulseDB.
set -uo pipefail
cd "$(dirname "$0")"

datasets=("$@")
[ ${#datasets[@]} -eq 0 ] && datasets=(ppgbp pttppg bpgraphene mimicbp pwdb pulsedb)

failed=()
for d in "${datasets[@]}"; do
    echo "=== $d"
    [ -x "download_$d.sh" ] || { echo "unknown dataset: $d"; failed+=("$d"); continue; }
    "./download_$d.sh" || failed+=("$d")
done

if [ ${#failed[@]} -gt 0 ]; then
    echo "FAILED: ${failed[*]} (rerun: scripts/download_all.sh ${failed[*]})"
    exit 1
fi
echo "All downloads OK. Next: uv run python scripts/structure.py"
