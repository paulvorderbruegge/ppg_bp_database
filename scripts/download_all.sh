#!/usr/bin/env bash
# Download all datasets (or only the ones given) to data/raw/<dataset>. Each download is verified with the
# checksums published by the source. Rerunning is safe: complete files are skipped or downloaded again.
#     scripts/download_all.sh                  # everything (~45 GB, PulseDB alone ~40 GB; see README for space)
#     scripts/download_all.sh ppgbp pttppg     # only these
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
echo "All downloads OK. Next: uv run python scripts/convert.py"
