#!/usr/bin/env bash
# Download MIMIC-BP v2.2 (Sanches et al. 2024, Harvard Dataverse doi:10.7910/DVN/DBM1NF) to data/raw/mimicbp,
# verify (MD5), unpack.
# By default only what the converter uses (ppg, abp, labels, subject splits; ~520 MB).
# `scripts/download_mimicbp.sh --all` also downloads ecg.zip and resp.zip (+ ~400 MB).
set -euo pipefail
cd "$(dirname "$0")/.."  # project root, so it works from any folder
mkdir -p data/raw/mimicbp && cd data/raw/mimicbp

# Dataverse file ID, MD5, name (from the Dataverse API)
FILES="
7574732 dbb4ed1c3c39808dac4732f6ba2b153b ppg.zip
7574731 65135891d34753f05be54d320b38fda5 abp.zip
7574727 d83b4bd797f48250e6c22f72c50d2d78 labels.zip
7574628 781e86f459fc3bdd7bbbec8bec3d41f4 train_subjects.txt
7574624 55591cd0466bb33854306ab0c917be49 val_subjects.txt
7574623 9f301714c829373fef5e106564148135 test_subjects.txt
"
if [ "${1:-}" = "--all" ]; then
    FILES="$FILES
7574838 640060d10d35fcb91d8f1f20c4148985 ecg.zip
7574839 1203fa59c78fefbdc3c7f83151f0324b resp.zip
"
fi

while read -r id md5 name; do
    [ -z "$id" ] && continue
    if [ -f "$name" ] && echo "$md5  $name" | md5sum -c --status -; then
        echo "$name: already downloaded"
        continue
    fi
    curl -fL --retry 5 -# -o "$name" "https://dataverse.harvard.edu/api/access/datafile/$id"
    echo "$md5  $name" | md5sum -c -
done <<< "$FILES"

# Unpack every zip into its own folder (the converter finds the .npy files anywhere below data/raw/mimicbp)
for z in *.zip; do
    unzip -q -o "$z" -d "${z%.zip}"
done
# Archives and unpacked files together need twice the space. Once everything is unpacked: rm *.zip
